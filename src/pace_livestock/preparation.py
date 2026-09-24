"""Small reproducible preparation steps for experimental datasets."""

import math
from collections import defaultdict
from pathlib import Path

from .config import DEFAULTS, load_yaml
from .errors import PaceError
from .io.tables import integer, number, read_table, unique, write_table
from .provenance import file_hash, output_directory, read_json, write_json
from .schemas import SCHEMAS

MERGE_KEYS = {
    "observed_activity": ("element_id", "sample_id", "assay"),
    "observed_contacts": ("element_id", "promoter_id", "sample_id"),
    "expression": ("gene_id", "sample_id"),
    "methylation": ("chrom", "dyad_start0", "sample_id", "assay"),
}


def merge_tables(args):
    rows = []
    for path in args.inputs:
        rows.extend(read_table(path, required=SCHEMAS[args.table].split()))
    unique(rows, MERGE_KEYS[args.table], args.table)
    with output_directory(args.out) as dest:
        write_table(
            dest / f"{args.table}.tsv", rows, fields=None if rows else SCHEMAS[args.table].split()
        )
        write_json(
            dest / "preparation_report.json",
            {
                "operation": "merge_tables",
                "rows": len(rows),
                "inputs": {str(Path(p).resolve()): file_hash(p) for p in args.inputs},
            },
        )
    return {"output": args.out, "rows": len(rows)}


def prepare_pairs(args):
    rows = build_pairs(args.catalog_dir)
    with output_directory(args.out) as dest:
        write_table(dest / "pairs.tsv", rows)
    return {"output": args.out, "pairs": len(rows)}


def build_pairs(catalog_dir) -> list[dict]:
    """Unique element-TSS pairs that need a contact value, sorted by key."""
    root = Path(catalog_dir)
    units = read_table(root / "units.tsv", required=SCHEMAS["units"].split())
    promoters = read_table(root / "promoters.tsv", required=SCHEMAS["promoters"].split())
    candidates = read_table(root / "candidates.tsv", required=SCHEMAS["candidates"].split())
    unique(units, ("element_id",), "units")
    unique(promoters, ("gene_id", "promoter_id"), "promoters")
    unique(candidates, ("element_id", "gene_id"), "candidates")
    elements = {r["element_id"]: r for r in units}
    genes = defaultdict(list)
    for p in promoters:
        genes[p["gene_id"]].append(p)
    rows = {}
    for edge in candidates:
        e, gene = edge["element_id"], edge["gene_id"]
        if e not in elements or gene not in genes:
            raise PaceError("Candidate references an unknown unit or gene")
        for p in genes[gene]:
            unit = elements[e]
            if p["chrom"] != unit["chrom"]:
                raise PaceError("Only cis contact pairs are supported")
            key = e, p["promoter_id"]
            row = {
                "element_id": e,
                "promoter_id": p["promoter_id"],
                "chrom": unit["chrom"],
                "anchor0": integer(unit["anchor0"], "anchor0"),
                "tss0": integer(p["tss0"], "tss0"),
            }
            if key in rows and rows[key] != row:
                raise PaceError("One physical promoter has inconsistent coordinates")
            rows[key] = row
    return [rows[k] for k in sorted(rows)]


def normalize_activity(args):
    rows = read_table(args.counts, required=["element_id", "sample_id", "assay", "count"])
    libs = read_table(args.library_sizes, required=["sample_id", "library_size"])
    units = read_table(args.units, required=["element_id", "start", "end"])
    unique(rows, MERGE_KEYS["observed_activity"], "counts")
    unique(libs, ("sample_id",), "library_sizes")
    unique(units, ("element_id",), "units")
    sizes = {r["sample_id"]: number(r["library_size"], "library_size", minimum=1) for r in libs}
    widths = {
        r["element_id"]: integer(r["end"], "end") - integer(r["start"], "start") for r in units
    }
    if any(v <= 0 for v in widths.values()):
        raise PaceError("Activity windows must have positive widths")
    if not args.window_id and len(set(widths.values())) != 1:
        raise PaceError("Variable-width regions require an explicit --window-id")
    window = args.window_id or f"grid:{next(iter(widths.values()))}:mean"
    output = []
    for r in rows:
        if r["sample_id"] not in sizes or r["element_id"] not in widths:
            raise PaceError("Count has an unknown sample or scoring unit")
        count = number(r["count"], "count", minimum=0, missing=True)
        total = sizes[r["sample_id"]]
        if count > total:
            raise PaceError("Window count exceeds its declared filtered library size")
        output.append(
            {
                "element_id": r["element_id"],
                "sample_id": r["sample_id"],
                "assay": r["assay"],
                "signal": count * 1e6 / total / widths[r["element_id"]],
                "measurement_status": "observed" if math.isfinite(count) else "unmeasured",
                "callable_fraction": 1.0 if math.isfinite(count) else 0.0,
                "unit": "CPM_per_bp",
                "normalization_id": args.normalization_id,
                "window_id": window,
            }
        )
    with output_directory(args.out) as dest:
        write_table(dest / "observed_activity.tsv", output)
        write_json(
            dest / "normalization.json",
            {
                "method": "count_per_million_filtered_fragments_per_bp",
                "library_sizes": libs,
                "input_hashes": {
                    name: file_hash(getattr(args, name))
                    for name in ("counts", "library_sizes", "units")
                },
                "note": "Library sizes must count all filtered fragments, not only fragments in candidate regions.",
            },
        )
    return {"output": args.out, "rows": len(output)}


def promoter_weights(args):
    promoters = read_table(args.promoters, required=SCHEMAS["promoters"].split())
    signals = read_table(args.signals, required=["promoter_id", "signal"])
    unique(signals, ("promoter_id",), "promoter signals")
    unique(promoters, ("gene_id", "promoter_id"), "promoters")
    lookup = {r["promoter_id"]: number(r["signal"], "promoter signal", minimum=0) for r in signals}
    genes = defaultdict(list)
    for p in promoters:
        if p["promoter_id"] not in lookup:
            raise PaceError(
                "Every planned TSS requires a measured promoter signal before weights are frozen"
            )
        genes[p["gene_id"]].append(p)
    source = f"{args.assay}:{args.normalization_id}:{file_hash(args.signals)}"
    for gene, group in genes.items():
        total = math.fsum(lookup[p["promoter_id"]] for p in group)
        if total == 0 and args.zero_policy == "error":
            raise PaceError(
                f"All promoter signals are zero for {gene}; choose an explicit --zero-policy equal if justified"
            )
        for p in group:
            p.update(
                pi=lookup[p["promoter_id"]] / total if total else 1 / len(group),
                pi_source=source if total else "all_zero_equal:" + source,
            )
    with output_directory(args.out) as dest:
        write_table(dest / "promoters.tsv", promoters)
        write_json(
            dest / "preparation_report.json",
            {
                "assay": args.assay,
                "normalization_id": args.normalization_id,
                "signal_sha256": file_hash(args.signals),
                "zero_policy": args.zero_policy,
            },
        )
    return {"output": args.out, "genes": len(genes)}


def init_project(args):
    inputs, catalog = {}, {}
    if args.catalog_dir:
        root = Path(args.catalog_dir).resolve()
        if not root.is_dir():
            raise PaceError("--catalog-dir must exist")
        if (root / "run_catalog_config.yaml").is_file():
            catalog = load_yaml(root / "run_catalog_config.yaml")["catalog"]
            if catalog.get("chrom_sizes_path"):
                catalog["chrom_sizes_path"] = str(root / catalog["chrom_sizes_path"])
        for name in DEFAULTS["inputs"]:
            options = [p for p in (root / f"{name}.tsv", root / f"{name}.tsv.gz") if p.is_file()]
            if len(options) > 1:
                raise PaceError(f"Ambiguous catalog table: {name}")
            if options:
                inputs[name] = str(options[0])
    contact = {"mode": "observed"}
    if args.contact_prior:
        p = Path(args.contact_prior).resolve()
        prior = read_json(p / "manifest.json" if p.is_dir() else p)
        contact = {
            "prior_path": str(p),
            "mode": "observed" if inputs.get("observed_contacts") else "prior_only",
            "allow_prior_fallback": True,
            **{
                k: prior.get(k)
                for k in ("scale", "resolution", "normalization_id", "balancing", "window_id")
            },
        }
    if args.prior_preset:
        if args.contact_prior:
            raise PaceError("Choose a fitted prior OR an explicit human reference baseline")
        contact = {"mode": "prior_only", "prior_preset": args.prior_preset}
        inputs.pop("observed_contacts", None)
        inputs.pop("resolved_contacts", None)
    cfg = {
        "context": {"species": args.species, "assembly": args.assembly, "context_id": args.tissue},
        "target_level": args.target_level,
        "activity": {"panel": args.panel},
        "inputs": inputs,
        "contact": contact,
    }
    if catalog:
        cfg["catalog"] = catalog
    with output_directory(args.out) as dest:
        required = ["units", "promoters", "candidates", "observed_activity"]
        if args.target_level == "population_mean":
            # Several animals: the donor of every library must be declared.
            required += ["samples", "sources"]
        if (
            contact["mode"] == "observed"
            and "observed_contacts" not in inputs
            and not args.contact_prior
        ):
            required.append("observed_contacts")
        (dest / "data").mkdir()
        missing = []
        for name in required:
            if name not in inputs:
                inputs[name] = f"data/{name}.tsv"
                write_table(dest / inputs[name], [], fields=SCHEMAS[name].split())
                missing.append(name)
        script = run_script(cfg)
        (dest / "run.sh").write_text(script, encoding="utf-8")
        (dest / "README.md").write_text(
            "# PACE project\n\n"
            "This template contains no invented experimental values. Fill the listed input tables first.\n\n"
            + "Tables to prepare: "
            + (", ".join(missing) or "none; verify the existing tables")
            + ".\n\n"
            "Use `pace catalog`, `pace activity`, `pace contacts`, `pace merge` and `pace fit-prior` as needed.\n\n"
            "Then, from this directory, check and score with `bash run.sh`, which runs:\n\n```bash\n"
            + script.split("\n", 2)[2]
            + "```\n\n"
            "Library-size normalization applies to raw counts only. A human contact prior is an explicitly transferred, unvalidated baseline.\n",
            encoding="utf-8",
        )
    return {
        "output": args.out,
        "script": str(Path(args.out) / "run.sh"),
        "tables_to_prepare": missing,
    }


INPUT_FLAGS = {
    "observed_activity": "--activity",
    "observed_contacts": "--contacts",
    "resolved_activity": "--resolved-activity",
    "resolved_contacts": "--resolved-contacts",
    "support_bounds": "--support-bounds",
}


def run_options(cfg: dict) -> list[str]:
    """Command-line options equivalent to a project's settings."""
    import shlex

    context = cfg["context"]
    words = [
        "--species",
        context["species"],
        "--assembly",
        context["assembly"],
        "--tissue",
        context["context_id"],
        "--target-level",
        cfg["target_level"],
        "--panel",
        *cfg["activity"]["panel"],
    ]
    for name, path in cfg["inputs"].items():
        if name == "region_membership":
            continue
        words += [INPUT_FLAGS.get(name, "--" + name.replace("_", "-")), path]
    catalog = cfg.get("catalog", {})
    for key, flag in (
        ("width_bp", "--unit-width"),
        ("offset_bp", "--grid-offset"),
        ("candidate_radius_bp", "--candidate-radius"),
        ("chrom_sizes_path", "--chrom-sizes"),
    ):
        if catalog.get(key) is not None:
            words += [flag, str(catalog[key])]
    if catalog.get("include_promoter_units") is False:
        words.append("--no-include-promoters")
    contact = cfg["contact"]
    words += ["--contact-mode", contact["mode"]]
    if contact.get("prior_path"):
        words += ["--contact-prior", contact["prior_path"]]
    if contact.get("allow_prior_fallback"):
        words.append("--allow-prior-fallback")
    if contact.get("prior_preset"):
        words += ["--prior-preset", contact["prior_preset"]]
    return [shlex.quote(str(w)) for w in words]


def run_script(cfg: dict) -> str:
    options = " \\\n  ".join(" ".join(pair) for pair in _pairs(run_options(cfg)))
    return (
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        f"pace validate \\\n  {options}\n"
        f"pace run \\\n  {options} \\\n  -o results\n"
    )


def _pairs(words):
    """Group an option with its values for readable line breaks."""
    group = []
    for word in words:
        if word.startswith("--") and group:
            yield group
            group = []
        group.append(word)
    if group:
        yield group
