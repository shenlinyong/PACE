"""RNA annotation; gene TPM is never interpreted as measured transcript usage."""

from collections import defaultdict

from ..errors import PaceError
from .tables import number, unique


def transcript_tpm_to_gene(rows, mapping):
    unique(mapping, ("transcript_id",), "transcript mapping")
    lookup = {r["transcript_id"]: r for r in mapping}
    unique(rows, ("transcript_id", "sample_id"), "transcript expression")
    values = defaultdict(float)
    for row in rows:
        if row["transcript_id"] not in lookup:
            raise PaceError("Transcript TPM has no gene/physical TSS mapping")
        gene = lookup[row["transcript_id"]]["gene_id"]
        values[gene, row["sample_id"]] += number(row["tpm"], "transcript TPM", minimum=0)
    return [
        {"gene_id": gene, "sample_id": sample, "tpm": value, "status": "observed"}
        for (gene, sample), value in sorted(values.items())
    ]
