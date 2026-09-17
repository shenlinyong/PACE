"""Reference FASTA and normalized VCF/BCF readers (no variant calling or liftover)."""

from __future__ import annotations

import gzip
import re
from bisect import bisect_left
from collections import defaultdict
from pathlib import Path

from ..errors import PaceError


class Reference:
    """Read regular uncompressed FASTA by byte offsets without loading the genome."""

    def __init__(self, path):
        self.path = Path(path)
        self.index = {}
        self.handle = open(path, "rb")
        chrom, offset, length, bases, width, short = None, 0, 0, 0, 0, False
        try:
            while True:
                position = self.handle.tell()
                line = self.handle.readline()
                if not line or line.startswith(b">"):
                    if chrom is not None:
                        if not length:
                            raise PaceError(f"Empty FASTA sequence: {chrom}")
                        self.index[chrom] = (offset, length, bases, width)
                    if not line:
                        break
                    chrom = line[1:].split()[0].decode()
                    if chrom in self.index:
                        raise PaceError(f"Duplicate FASTA chromosome: {chrom}")
                    offset, length, bases, width, short = self.handle.tell(), 0, 0, 0, False
                else:
                    sequence = line.rstrip(b"\r\n")
                    if chrom is None or not sequence or set(sequence.upper()) - set(b"ACGTN"):
                        raise PaceError("Reference must be uncompressed FASTA containing A/C/G/T/N")
                    if not bases:
                        bases, width, offset = len(sequence), len(line), position
                    elif (
                        short
                        or len(sequence) > bases
                        or (len(sequence) == bases and len(line) != width and line.endswith(b"\n"))
                    ):
                        raise PaceError("FASTA lines must have regular width except the final line")
                    short = len(sequence) < bases or not line.endswith(b"\n")
                    length += len(sequence)
            if not self.index:
                raise PaceError("Empty FASTA")
        except Exception:
            self.handle.close()
            raise

    @property
    def sizes(self):
        return {k: v[1] for k, v in self.index.items()}

    def fetch(self, chrom, start, end):
        if chrom not in self.index or start < 0 or end < start or end > self.index[chrom][1]:
            raise PaceError(f"FASTA interval outside reference: {chrom}:{start}-{end}")
        offset, _, bases, width = self.index[chrom]
        result = []
        while start < end:
            n = min(end - start, bases - start % bases)
            self.handle.seek(offset + (start // bases) * width + start % bases)
            result.append(self.handle.read(n).decode("ascii").upper())
            start += n
        return "".join(result)

    def close(self):
        self.handle.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def read_variants(path, *, sample_id=None):
    if str(path).endswith(".bcf"):
        try:
            import pysam
        except ImportError as exc:
            raise PaceError("BCF requires pace-livestock[io]") from exc
        with pysam.VariantFile(path) as source:
            names = list(source.header.samples)
            sample = _choose_sample(names, sample_id)
            return [
                {
                    "chrom": r.chrom,
                    "pos0": r.pos - 1,
                    "ref": r.ref.upper(),
                    "alts": tuple(r.alts or ()),
                    "gt": tuple(r.samples[sample].get("GT") or (None,)),
                    "phased": bool(r.samples[sample].phased),
                    "phase_set": r.samples[sample].get("PS"),
                    "end": r.stop,
                    "filter": ";".join(r.filter.keys()) or ".",
                }
                for r in source
            ]
    opener = gzip.open if str(path).endswith(".gz") else open
    rows, sample_index = [], None
    with opener(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                fields = line.strip().split("\t")
                sample = _choose_sample(fields[9:], sample_id)
                sample_index = fields.index(sample)
                continue
            if line.startswith("#") or not line.strip():
                continue
            if sample_index is None:
                raise PaceError("VCF needs a #CHROM header and selected sample")
            fields = line.strip().split("\t")
            if len(fields) <= sample_index:
                raise PaceError("Malformed VCF sample record")
            fmt = dict(zip(fields[8].split(":"), fields[sample_index].split(":")))
            gt_text = fmt.get("GT", ".")
            separator = "|" if "|" in gt_text else "/"
            try:
                gt = tuple(None if a == "." else int(a) for a in gt_text.split(separator))
                info = dict(x.split("=", 1) for x in fields[7].split(";") if "=" in x)
                pos = int(fields[1]) - 1
                end = int(info.get("END", pos + len(fields[3])))
            except ValueError as exc:
                raise PaceError("Invalid VCF coordinate, GT or END") from exc
            rows.append(
                {
                    "chrom": fields[0],
                    "pos0": pos,
                    "ref": fields[3].upper(),
                    "alts": tuple(fields[4].split(",")),
                    "gt": gt,
                    "phased": separator == "|",
                    "phase_set": fmt.get("PS"),
                    "end": end,
                    "filter": fields[6],
                }
            )
    if sample_index is None:
        raise PaceError("VCF needs a #CHROM header and selected sample")
    return rows


def _choose_sample(names, sample):
    if len(set(names)) != len(names):
        raise PaceError("VCF sample names must be unique")
    if sample is None and len(names) == 1:
        return names[0]
    if sample not in names:
        raise PaceError("Specify genome.sample_id matching exactly one VCF sample")
    return sample


def is_structural(variant, *, max_indel=50):
    """Classify alleles carried by the selected sample, not every ALT in the record."""
    return any(
        a.startswith("<")
        or any(c in a for c in "[]*.")
        or abs(len(a) - len(variant["ref"])) > max_indel
        for a in selected_alts(variant)
    )


def selected_alts(variant):
    """Return only called non-reference alleles; missing GT remains unknown.

    Callers must handle ``None`` in GT separately. It must never be converted to
    reference or used as an ALT index.
    """
    gt = variant.get("gt", ())
    if not gt:
        raise PaceError("VCF record requires a selected sample genotype")
    if any(a is not None and (a < 0 or a > len(variant["alts"])) for a in gt):
        raise PaceError("VCF genotype allele index outside ALT list")
    return tuple(variant["alts"][a - 1] for a in sorted({a for a in gt if a}))


def relationship_affected(variant):
    """Whether a reported call prevents using unchanged reference distances.

    Unknown/filtered calls cannot establish that the relationship is intact.
    SNVs and unselected alternate alleles do not change genomic distance.
    """
    alts = selected_alts(variant)
    return (
        any(a is None for a in variant["gt"])
        or variant.get("filter", "PASS") not in ("PASS", ".")
        or is_structural(variant)
        or any(len(a) != len(variant["ref"]) for a in alts)
    )


def breakend_endpoints(variant):
    """Return remote zero-based BND endpoints relevant to the selected genotype.

    With a missing genotype, possible remote endpoints are retained solely to
    mark uncertainty; this does not assert that the sample carries the BND.
    """
    selected = selected_alts(variant)
    alleles = variant["alts"] if any(a is None for a in variant["gt"]) else selected
    endpoints = set()
    for allele in alleles:
        if not any(c in allele for c in "[]"):
            continue
        match = re.fullmatch(
            r"(?:[ACGTNacgtn]+)?([\[\]])([^\[\]]+):([1-9][0-9]*)\1(?:[ACGTNacgtn]+)?",
            allele,
        )
        if match is None:
            raise PaceError(f"Malformed VCF breakend ALT: {allele}")
        endpoints.add((match.group(2), int(match.group(3)) - 1))
    return sorted(endpoints)


class VariantIndex:
    """Index selected calls and both BND endpoints, including spanning SVs."""

    def __init__(self, variants, *, include_remote_breakends=True):
        self.rows, self.starts, self.max_ends = {}, {}, {}
        by_chrom = defaultdict(list)
        for row in variants:
            # END may describe an ALT absent from the sample (e.g. C,<DEL> with
            # GT=1/1). It must not make that sample's SNV span the unused deletion.
            end = (
                row["end"]
                if is_structural(row) or any(a is None for a in row["gt"])
                else row["pos0"] + len(row["ref"])
            )
            by_chrom[row["chrom"]].append({**row, "end": end})
            for chrom, pos0 in breakend_endpoints(row) if include_remote_breakends else []:
                by_chrom[chrom].append(
                    {
                        **row,
                        "chrom": chrom,
                        "pos0": pos0,
                        "end": pos0 + 1,
                        "remote_breakend": True,
                        "source_chrom": row["chrom"],
                        "source_pos0": row["pos0"],
                    }
                )
        for chrom, rows in by_chrom.items():
            rows.sort(key=lambda x: (x["pos0"], x["end"]))
            self.rows[chrom], self.starts[chrom] = rows, [x["pos0"] for x in rows]
            current, ends = 0, []
            for row in rows:
                current = max(current, row["end"])
                ends.append(current)
            self.max_ends[chrom] = ends

    def query(self, chrom, start, end):
        rows = self.rows.get(chrom, [])
        hi = bisect_left(self.starts.get(chrom, []), end)
        lo = bisect_left(self.max_ends.get(chrom, []), start + 1)
        return [r for r in rows[lo:hi] if r["end"] > start]
