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
                    **({"strand": fields[5]} if len(fields) > 5 else {}),
                }
            )
    return rows


BIOTYPE_KEYS = ("gene_biotype", "gene_type", "transcript_biotype", "transcript_type")


def read_chrom_sizes(path, *, aliases=None) -> dict[str, int]:
    """Chromosome lengths from a chrom.sizes file, a FASTA .fai index or a bigWig header.

    Text files may have a `chrom length` header line or none (UCSC/.fai style).
    """
    aliases = aliases or {}
    if str(path).lower().endswith((".bw", ".bigwig")):
        try:
            import pyBigWig
        except ImportError as exc:
            raise PaceError("Reading sizes from a bigWig requires the io extra") from exc
        with pyBigWig.open(str(path)) as bw:
            return {aliases.get(k, k): int(v) for k, v in bw.chroms().items()}
    sizes = {}
    with text_open(path) as handle:
        for n, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                fields = line.split()
            if n == 1 and fields[:2] == ["chrom", "length"]:
                continue
            if len(fields) < 2:
                raise PaceError(f"{path}:{n}: expected chromosome name and length")
            chrom = aliases.get(fields[0], fields[0])
            if chrom in sizes:
                raise PaceError(f"{path}:{n}: duplicate chromosome {chrom}")
            sizes[chrom] = integer(fields[1], f"{path}:{n} chromosome length", minimum=1)
    if not sizes:
        raise PaceError(f"{path}: no chromosome lengths found")
    return sizes


def read_gtf(path, *, aliases=None, gene_types=None):
    """Distinct physical TSSs from GTF transcript lines.

    gene_types keeps only transcripts whose gene_biotype/gene_type/transcript_biotype/
    transcript_type attribute is listed (for example protein_coding or mRNA).
    """
    wanted = set(gene_types or ())
    annotated = False
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
            if wanted:
                types = {attr[k] for k in BIOTYPE_KEYS if k in attr}
                annotated |= bool(types)
                if not types & wanted:
                    continue
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
    if wanted and not annotated:
        raise PaceError(
            "--gene-types was given but no transcript line has a gene_biotype, gene_type, "
            "transcript_biotype or transcript_type attribute"
        )
    if not promoters:
        raise PaceError(
            "GTF has no transcript records; provide a transcript annotation or promoters.tsv"
        )
    return promoters, transcripts
