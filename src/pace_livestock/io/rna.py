"""RNA annotation; gene TPM is never interpreted as measured transcript usage."""

import math
from collections import defaultdict

from ..errors import PaceError
from .tables import number, unique


def transcript_tpm_to_gene(rows, mapping):
    unique(mapping, ("transcript_id",), "transcript mapping")
    lookup = {r["transcript_id"]: r for r in mapping}
    unique(rows, ("transcript_id", "sample_id"), "transcript expression")
    values = defaultdict(float)
    incomplete = set()
    for row in rows:
        if row["transcript_id"] not in lookup:
            raise PaceError("Transcript TPM has no gene/physical TSS mapping")
        gene = lookup[row["transcript_id"]]["gene_id"]
        key = gene, row["sample_id"]
        status = row.get("status", "observed")
        if status not in {
            "observed",
            "unmeasured",
            "low_coverage",
            "unmappable",
            "invalid",
            "not_applicable",
        }:
            raise PaceError("Invalid transcript expression status")
        values.setdefault(key, 0.0)
        if status != "observed":
            incomplete.add(key)
            continue
        values[key] += number(row["tpm"], "transcript TPM", minimum=0)
    return [
        {
            "gene_id": gene,
            "sample_id": sample,
            "tpm": math.nan if (gene, sample) in incomplete else value,
            "status": "unmeasured" if (gene, sample) in incomplete else "observed",
        }
        for (gene, sample), value in sorted(values.items())
    ]
