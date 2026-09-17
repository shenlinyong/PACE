"""Reproducible sparse-candidate/kernel benchmark; no biological performance claim."""

import json
import resource
import time

from pace_livestock.catalog import candidate_edges
from pace_livestock.core import score
from pace_livestock.provenance import environment


def main():
    units = [
        dict(element_id=f"E{i}", chrom="chrSynthetic", anchor0=i * 1000 + 249) for i in range(20000)
    ]
    promoters = [
        dict(gene_id=f"G{i}", chrom="chrSynthetic", tss0=i * 20000 + 10000) for i in range(1000)
    ]
    start = time.perf_counter()
    edges = candidate_edges(units, promoters, radius=50000)
    candidate_seconds = time.perf_counter() - start
    for row in edges:
        row.update(A_used=1.0, Cbar=2.0)
    start = time.perf_counter()
    scored, summary = score(edges, eta=1)
    score_seconds = time.perf_counter() - start
    print(
        json.dumps(
            dict(
                n_units=len(units),
                n_genes=len(promoters),
                n_edges=len(scored),
                n_gene_summaries=len(summary),
                candidate_seconds=candidate_seconds,
                score_seconds=score_seconds,
                peak_rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024,
                measurement_scope="Linux process RSS; sparse candidates and numerical kernel only",
                environment=environment(),
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
