"""BED0/GTF adapters with explicit chromosome aliases and physical TSS deduplication."""

import gzip
import re
from collections import defaultdict

from ..errors import PaceError
from .tables import integer


def text_open(path):
    return (
        gzip.open(path, "rt", encoding="utf-8")
        if str(path).endswith(".gz")
        else open(path, encoding="utf-8")
    )


def read_bed(path, *, source_id: str, aliases=None):
    rows = []
    with text_open(path) as handle:
        for n, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 3:
                raise PaceError(f"BED line {n} needs at least three tab-delimited columns")
            chrom = (aliases or {}).get(fields[0], fields[0])
            start, end = integer(fields[1], "BED start"), integer(fields[2], "BED end")
            if start >= end:
                raise PaceError(f"BED line {n}: end must exceed start")
            rows.append(
                {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "region_id": fields[3] if len(fields) > 3 else f"{chrom}:{start}-{end}",
                    "source_id": source_id,
                }
            )
    return rows


def read_gtf(path, *, aliases=None):
    coords, transcripts = set(), []
    with text_open(path) as handle:
        for n, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9:
                raise PaceError(f"GTF line {n}: expected nine fields")
            if fields[2] != "transcript":
                continue
            start, end = (
                integer(fields[3], "GTF start", minimum=1),
                integer(fields[4], "GTF end", minimum=1),
            )
            if end < start or fields[6] not in ("+", "-"):
                raise PaceError(f"GTF line {n}: invalid interval/strand")
            attr = dict(re.findall(r'(\S+)\s+"([^"]+)"', fields[8]))
            if not attr.get("gene_id") or not attr.get("transcript_id"):
                raise PaceError(
                    "GTF transcript requires gene_id and transcript_id; versions are retained"
                )
            chrom = (aliases or {}).get(fields[0], fields[0])
            tss = start - 1 if fields[6] == "+" else end - 1
            coords.add((attr["gene_id"], chrom, tss, fields[6]))
            transcripts.append(
                {
                    "transcript_id": attr["transcript_id"],
                    "gene_id": attr["gene_id"],
                    "promoter_id": f"{chrom}:{tss}:{fields[6]}",
                }
            )
    counts = defaultdict(int)
    for gene, *_ in coords:
        counts[gene] += 1
    promoters = [
        {
            "gene_id": gene,
            "promoter_id": f"{chrom}:{tss}:{strand}",
            "chrom": chrom,
            "tss0": tss,
            "strand": strand,
            "pi": 1 / counts[gene],
            "pi_source": "equal_physical_tss",
        }
        for gene, chrom, tss, strand in sorted(coords)
    ]
    if not promoters:
        raise PaceError(
            "GTF has no transcript records; provide a transcript annotation or promoters.tsv"
        )
    return promoters, transcripts
