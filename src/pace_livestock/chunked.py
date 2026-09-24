"""Chromosome-chunked execution with bounded memory.

Every PACE quantity used by the default model is computed within one chromosome:
candidates are cis, a gene's denominator only contains elements on its chromosome and
activity and contact are resolved per element. Scoring whole chromosomes separately
and concatenating the results therefore reproduces a single genome-wide run, while
peak memory follows the largest chunk instead of the genome.

Settings that estimate one parameter from all chromosomes (functional eta labels, a
frozen eta calibrator, eQTL weak models and per-pair kappa fitting) are rejected.
"""

from __future__ import annotations

import copy
import csv
import gzip
import re
import shutil
from collections import defaultdict
from pathlib import Path

from .errors import PaceError
from .io.tables import read_table, write_table
from .provenance import digest, file_hash, output_directory, write_json

# Tables read once and passed unchanged to every chunk.
GLOBAL_TABLES = ("samples", "sources", "evidence")
# Output tables concatenated across chunks. Chunks hold disjoint chromosomes, so their
# elements, genes and promoters (and every key built from them) are disjoint too.
MERGED_TABLES = (
    "scores.tsv.gz",
    "gene_summary.tsv",
    "region_scores.tsv",
    "resolved_activity.tsv",
    "resolved_contacts.tsv",
    "promoter_weights.tsv",
    "multiomics_features.tsv.gz",
)
IDENTIFIED_TABLES = {"evidence.tsv": "evidence_id", "sources.tsv": "source_id"}
DEFAULT_CHUNK_PAIRS = 1_500_000


def natural_key(name: str):
    return [int(part) if part.isdigit() else part for part in re.split(r"(\d+)", name)]


def check_chunkable(cfg: dict) -> None:
    allocation = cfg["allocation"]
    if allocation["labels_path"] or allocation["calibrator_path"]:
        raise PaceError("Chromosome chunks cannot share one eta fit; run without chunking")
    if allocation.get("weak_model_path"):
        raise PaceError("eQTL weak models are applied genome-wide; run without chunking")
    if cfg["contact"]["reliability"] == "per_pair" and cfg["contact"]["kappa"] == "auto":
        raise PaceError("Per-pair shrinkage fits kappa genome-wide; set a fixed kappa to chunk")
    if cfg["inputs"].get("labels"):
        raise PaceError("Label tables are genome-wide; run without chunking")


def _stream(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    handle = opener(path, "rt", encoding="utf-8", newline="")
    return handle, csv.reader(handle, delimiter="\t")


def chromosome_maps(cfg: dict):
    """element -> chrom, promoter -> chrom and gene -> chrom from the catalog tables."""
    element, promoter, gene = {}, {}, {}
    for row in read_table(cfg["inputs"]["units"], required=["element_id", "chrom"]):
        element[row["element_id"]] = row["chrom"]
    for row in read_table(cfg["inputs"]["promoters"], required=["gene_id", "promoter_id", "chrom"]):
        promoter[row["promoter_id"]] = row["chrom"]
        if gene.setdefault(row["gene_id"], row["chrom"]) != row["chrom"]:
            raise PaceError("Gene promoters on multiple chromosomes are unsupported")
    return element, promoter, gene


def row_chromosome(name, header, maps):
    element, promoter, gene = maps
    index = {field: i for i, field in enumerate(header)}

    def lookup(mapping, key, what):
        if key not in mapping:
            raise PaceError(f"{name}: unknown {what} {key}")
        return mapping[key]

    if name in ("units", "promoters", "methylation"):
        i = index["chrom"]
        return lambda row: row[i]
    if name == "expression":
        i = index["gene_id"]
        return lambda row: lookup(gene, row[i], "gene_id")
    if name == "features":
        t, e = index["entity_type"], index["entity_id"]
        kinds = {"element": element, "promoter": promoter, "gene": gene}

        def feature(row):
            if row[t] == "edge":
                return lookup(element, row[e].split("|", 1)[0], "element_id")
            return lookup(kinds.get(row[t], {}), row[e], row[t])

        return feature
    i = index["element_id"]
    return lambda row: lookup(element, row[i], "element_id")


def plan_chunks(cfg: dict, maps, max_pairs: int) -> list[list[str]]:
    """Pack whole chromosomes into chunks of at most max_pairs candidate pairs."""
    counts = defaultdict(int)
    handle, reader = _stream(cfg["inputs"]["candidates"])
    with handle:
        header = next(reader)
        where = row_chromosome("candidates", header, maps)
        for row in reader:
            counts[where(row)] += 1
    chunks, current, size = [], [], 0
    for chrom in sorted(counts, key=natural_key):
        if current and size + counts[chrom] > max_pairs:
            chunks.append(current)
            current, size = [], 0
        current.append(chrom)
        size += counts[chrom]
    if current:
        chunks.append(current)
    return chunks


def split_inputs(cfg: dict, maps, chunks, work: Path) -> list[dict]:
    """Write each chromosome-specific input table once, split by chunk."""
    chunk_of = {chrom: i for i, group in enumerate(chunks) for chrom in group}
    configs = []
    for i in range(len(chunks)):
        (work / f"chunk_{i:04d}").mkdir(parents=True)
        configs.append(copy.deepcopy(cfg))
    for name, path in cfg["inputs"].items():
        if not path or name in GLOBAL_TABLES:
            continue
        handle, reader = _stream(path)
        outputs = {}
        try:
            with handle:
                header = next(reader, None)
                if header is None:
                    raise PaceError(f"inputs.{name}: empty file")
                where = row_chromosome(name, header, maps)
                for row in reader:
                    chunk = chunk_of.get(where(row))
                    if chunk is None:
                        # Chromosomes without any candidate pair contribute nothing.
                        continue
                    if chunk not in outputs:
                        target = work / f"chunk_{chunk:04d}" / f"{name}.tsv"
                        stream = open(target, "w", encoding="utf-8", newline="")
                        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
                        writer.writerow(header)
                        outputs[chunk] = (stream, writer)
                    outputs[chunk][1].writerow(row)
        finally:
            for stream, _ in outputs.values():
                stream.close()
        for i, chunk_cfg in enumerate(configs):
            target = work / f"chunk_{i:04d}" / f"{name}.tsv"
            if not target.exists():
                with open(target, "w", encoding="utf-8", newline="") as stream:
                    csv.writer(stream, delimiter="\t", lineterminator="\n").writerow(header)
            chunk_cfg["inputs"][name] = str(target)
    return configs


def merge_qc(values: list):
    """Sum counts, merge mappings and keep single values; conflicting values become lists."""
    first = values[0]
    if all(isinstance(v, dict) for v in values):
        keys = list(dict.fromkeys(k for v in values for k in v))
        return {k: merge_qc([v[k] for v in values if k in v]) for k in keys}
    if all(type(v) is int for v in values):
        return sum(values)
    distinct = []
    for value in values:
        if value not in distinct:
            distinct.append(value)
    return first if len(distinct) == 1 else distinct


def concatenate(sources: list[Path], target: Path) -> int:
    """Concatenate TSV files with a union header; returns the number of rows."""
    fields = []
    for path in sources:
        handle, reader = _stream(path)
        with handle:
            for field in next(reader, []):
                if field not in fields:
                    fields.append(field)
    headers = []
    for path in sources:
        handle, reader = _stream(path)
        with handle:
            headers.append(next(reader, []))
    if all(h == fields for h in headers):
        # Same columns everywhere: copy the data lines without parsing them.
        return _copy_lines(sources, target)
    n = 0
    opener = gzip.open if str(target).endswith(".gz") else open
    kwargs = {"compresslevel": 6} if str(target).endswith(".gz") else {}
    with opener(target, "wt", encoding="utf-8", newline="", **kwargs) as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        for path in sources:
            handle, reader = _stream(path)
            with handle:
                header = next(reader, [])
                position = [header.index(f) if f in header else None for f in fields]
                for row in reader:
                    writer.writerow(["" if i is None else row[i] for i in position])
                    n += 1
    return n


def _copy_lines(sources: list[Path], target: Path) -> int:
    n = 0
    zipped = str(target).endswith(".gz")
    out = gzip.open(target, "wb", compresslevel=6) if zipped else open(target, "wb")
    with out:
        for k, path in enumerate(sources):
            opener = gzip.open if str(path).endswith(".gz") else open
            with opener(path, "rb") as handle:
                header = handle.readline()
                if k == 0:
                    out.write(header)
                for line in handle:
                    out.write(line)
                    n += 1
    return n


def merge_identified(sources: list[Path], target: Path, identifier: str) -> None:
    """Union of evidence/source rows; repeated IDs must agree apart from their checksum.

    Rows are streamed in chunk order; only IDs occurring in several chunks (run-level
    sources, models and the core-formula record) are held back and written last.
    """
    headers = []
    for path in sources:
        handle, reader = _stream(path)
        with handle:
            headers.append(next(reader, []))
    fields = list(dict.fromkeys(f for h in headers for f in h))
    counts = defaultdict(int)
    for path, header in zip(sources, headers, strict=True):
        handle, reader = _stream(path)
        with handle:
            next(reader, None)
            column = header.index(identifier)
            for row in reader:
                counts[row[column]] += 1
    shared = {}
    with open(target, "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        for path, header in zip(sources, headers, strict=True):
            handle, reader = _stream(path)
            with handle:
                next(reader, None)
                for values in reader:
                    row = dict(zip(header, values, strict=True))
                    key = row[identifier]
                    if counts[key] == 1:
                        writer.writerow([row.get(f, "") for f in fields])
                    elif key not in shared:
                        shared[key] = {**row, "_checksums": [row.get("checksum")]}
                    else:
                        kept = shared[key]
                        for field, value in row.items():
                            if field != "checksum" and kept.get(field) != value:
                                raise PaceError(
                                    f"{target.name}: {key} differs between chunks ({field})"
                                )
                        kept["_checksums"].append(row.get("checksum"))
        for key in sorted(shared):
            row = shared[key]
            checksums = row.pop("_checksums")
            if len(set(checksums)) > 1:
                row["checksum"] = digest(checksums)
            writer.writerow([row.get(f, "") for f in fields])


def run_by_chromosome(
    cfg: dict,
    out,
    *,
    max_pairs: int = DEFAULT_CHUNK_PAIRS,
    dest: Path | None = None,
    relative_to=None,
    log=None,
    threads: int = 1,
):
    """Score chunk by chunk (in parallel with threads > 1) and publish one result folder."""
    if dest is None:
        with output_directory(out) as staging:
            return run_by_chromosome(
                cfg,
                out,
                max_pairs=max_pairs,
                dest=staging,
                relative_to=staging,
                log=log,
                threads=threads,
            )
    from .schemas import infer_metadata

    check_chunkable(cfg)
    global_checks(cfg)
    maps = chromosome_maps(cfg)
    chunks = plan_chunks(cfg, maps, max_pairs)
    if not chunks:
        raise PaceError("inputs.candidates: a nonempty table is required")
    work = dest / ".chunks"
    cfg = copy.deepcopy(cfg)
    inferred = global_metadata(cfg, infer_metadata)
    for name, rows in inferred.items():
        path = dest / f"inferred_{name}.tsv"
        write_table(path, rows)
        cfg["inputs"][name] = str(path)
    configs = split_inputs(cfg, maps, chunks, work)
    jobs = [
        (group, chunk_cfg, work / f"chunk_{i:04d}" / "result")
        for i, (group, chunk_cfg) in enumerate(zip(chunks, configs, strict=True))
    ]
    results = [None] * len(jobs)

    def report(i):
        if log:
            group = chunks[i]
            shown = ", ".join(group[:3]) + (f" +{len(group) - 3}" if len(group) > 3 else "")
            log(f"    chunk {i + 1}/{len(chunks)} done: {shown}")

    if threads > 1 and len(jobs) > 1:
        import concurrent.futures
        import multiprocessing

        context = multiprocessing.get_context("spawn")
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=min(threads, len(jobs)), mp_context=context
        ) as pool:
            futures = {pool.submit(score_chunk, *job): i for i, job in enumerate(jobs)}
            for future in concurrent.futures.as_completed(futures):
                i = futures[future]
                results[i] = future.result()
                report(i)
    else:
        for i, job in enumerate(jobs):
            results[i] = score_chunk(*job)
            report(i)
    parts = [work / f"chunk_{i:04d}" / "result" for i in range(len(chunks))]
    for name in MERGED_TABLES:
        concatenate([p / name for p in parts], dest / name)
    for name, identifier in IDENTIFIED_TABLES.items():
        merge_identified([p / name for p in parts], dest / name, identifier)
    for name in ("eta_calibration.json", "ml_feature_contract.json"):
        shutil.copyfile(parts[0] / name, dest / name)
    shutil.rmtree(work)
    qc = merge_qc([r["qc"] for r in results])
    qc["inferred_metadata"] = sorted(inferred)
    chunk_table = [
        {"chromosomes": r["chromosomes"], "n_candidates": r["n_candidates"]} for r in results
    ]
    qc["execution"] = {"mode": "by_chromosome", "chunks": chunk_table}
    manifest = merge_manifest(cfg, results)
    write_json(dest / "qc_report.json", qc)
    write_json(dest / "run_manifest.json", manifest)
    from .pipeline import write_config_and_report

    write_config_and_report(dest, cfg, qc, qc["n_candidates"], relative_to=relative_to)
    return {"qc": qc, "manifest": manifest, "eta_calibration": results[0]["eta"]}


def score_chunk(chromosomes, chunk_cfg, result_dir):
    """Score one chunk and write its result folder; runs in a worker process."""
    from .pipeline import compute, write_results

    result = compute(chunk_cfg, chunk=True)
    result_dir.mkdir()
    write_results(result_dir, chunk_cfg, result)
    return {
        "chromosomes": chromosomes,
        "qc": result["qc"],
        "manifest": result["manifest"],
        "eta": result["eta_calibration"],
        "n_candidates": len(result["scores"]),
    }


def global_metadata(cfg: dict, infer) -> dict:
    """Infer samples/sources once from distinct identifiers across all input tables."""
    needed = {
        "observed_activity": ("sample_id", "assay", "source_id", "normalization_id"),
        "observed_contacts": ("sample_id", "source_id", "normalization_id"),
        "expression": ("sample_id",),
        "methylation": ("sample_id", "assay"),
        "region_membership": ("source_id",),
        "evidence": ("source_id",),
    }
    if cfg["inputs"]["samples"] and cfg["inputs"]["sources"]:
        return {}
    tables = defaultdict(list)
    for name, fields in needed.items():
        path = cfg["inputs"].get(name)
        if not path:
            continue
        distinct = set()
        handle, reader = _stream(path)
        with handle:
            header = next(reader, [])
            index = [header.index(f) if f in header else None for f in fields]
            for row in reader:
                distinct.add(tuple(None if i is None else row[i] or None for i in index))
        tables[name] = [
            {f: (None if v == "NA" else v) for f, v in zip(fields, values, strict=True)}
            for values in sorted(distinct, key=lambda t: tuple(str(v) for v in t))
        ]
    for name in ("samples", "sources"):
        if cfg["inputs"][name]:
            tables[name] = read_table(cfg["inputs"][name])
    return infer(tables, cfg)


def merge_manifest(cfg: dict, results: list[dict]) -> dict:
    manifest = copy.deepcopy(max(results, key=lambda r: r["n_candidates"])["manifest"])
    chunk_ids = [
        {"chromosomes": r["chromosomes"], **r["manifest"]["universe_ids"]} for r in results
    ]
    merged_ids = {
        key: digest([(c["chromosomes"], c[key]) for c in chunk_ids])
        for key in results[0]["manifest"]["universe_ids"]
    }
    manifest["universe_ids"] = merged_ids
    manifest["comparison_contract"].update(merged_ids)
    manifest["chunk_universe_ids"] = chunk_ids
    manifest["config_hash"] = digest(cfg)
    manifest["input_hashes"] = {
        name: file_hash(path) for name, path in cfg["inputs"].items() if path
    }
    differences = {}
    reference = manifest["comparison_contract"]
    for r in results:
        for key, value in r["manifest"]["comparison_contract"].items():
            if key not in merged_ids and value != reference.get(key):
                differences.setdefault(key, []).append(r["chromosomes"])
    manifest["execution"] = {
        "mode": "by_chromosome",
        "chunks": [r["chromosomes"] for r in results],
        "contract_differences_by_chunk": differences,
        "note": "Scores equal a genome-wide run: candidates, denominators and resolution are cis.",
    }
    return manifest


def global_checks(cfg: dict) -> None:
    """Checks that concern the whole genome, run once before chunking."""
    from .run_options import _column_values

    panel = set(cfg["activity"]["panel"])
    available = set()
    for name in ("observed_activity", "resolved_activity"):
        path = cfg["inputs"].get(name)
        if path:
            available |= _column_values(path, "assay", limit=10**6)
    if not panel <= available:
        raise PaceError(
            "Measured mode requires tables for every declared activity assay; "
            "choose an explicit single-layer panel when appropriate"
        )
    if cfg["contact"]["mode"] == "observed" and not (
        cfg["inputs"].get("observed_contacts")
        or cfg["inputs"].get("resolved_contacts")
        or (cfg["contact"]["allow_prior_fallback"] and cfg["contact"]["prior_path"])
    ):
        raise PaceError(
            "Observed contact mode needs contact observations or an explicitly allowed prior fallback"
        )
