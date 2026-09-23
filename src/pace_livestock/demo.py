"""Deterministic offline fixture with measured activity and contact tables."""

from pathlib import Path

import yaml

from .io.tables import write_table
from .provenance import output_directory
from .schemas import SCHEMAS


def create_example(destination, regime="measured"):
    if regime != "measured":
        raise ValueError("PACE demonstrations require measured activity")
    dest = Path(destination)
    dest.mkdir(parents=True, exist_ok=True)
    context = {"species": "synthetic", "assembly": "toy_assembly", "context_id": "toy_tissue"}
    units = [
        {
            "element_id": f"E{i + 1}",
            "chrom": "chrToy",
            "start": 5000 + i * 500,
            "end": 5500 + i * 500,
            "anchor0": 5249 + i * 500,
            "element_roles": "enhancer",
            "canonical_catalog_id": "toy_grid",
        }
        for i in range(3)
    ]
    promoters = [
        {
            "gene_id": gene,
            "promoter_id": p,
            "chrom": "chrToy",
            "tss0": tss,
            "strand": "+",
            "pi": 1,
            "pi_source": "synthetic_equal",
        }
        for gene, p, tss in [("G1", "P1", 10000), ("G2", "P2", 14000)]
    ]
    candidates = [
        {
            "element_id": u["element_id"],
            "gene_id": p["gene_id"],
            "candidate_universe_id": "toy_candidates",
        }
        for u in units
        for p in promoters
    ]
    write_table(dest / "units.tsv", units)
    write_table(dest / "promoters.tsv", promoters)
    write_table(dest / "candidates.tsv", candidates)
    source = dict(
        source_id="toy_source",
        path_or_accession="generated:PACE-synthetic-fixture",
        source_type="synthetic",
        assembly="toy_assembly",
        processing_method="deterministic_fixture",
        normalization_id="toy_mean",
        checksum="synthetic_recipe",
    )
    write_table(dest / "sources.tsv", [source])
    write_table(dest / "evidence.tsv", [], fields=SCHEMAS["evidence"].split())
    cfg = {
        "schema_version": "pace-1",
        "run_id": f"toy_{regime}",
        "regime": regime,
        "execution_profile": "demonstration",
        "estimand": "bulk_proxy",
        "target_level": "individual",
        "context": context,
        "inputs": {
            k: f"{k}.tsv" for k in ("units", "promoters", "candidates", "sources", "evidence")
        },
        "catalog": {"include_promoter_units": False},
        "activity": {"panel": ["ATAC", "H3K27ac"]},
        "contact": {"scale": "toy_contact", "mode": "observed"},
        "allocation": {"eta": "auto"},
        "seed": 17,
    }
    samples = [
        {
            "sample_id": f"S_{assay}",
            "donor_id": "toy_animal",
            "assay": assay,
            "biological_replicate": "1",
            "technical_replicate": "1",
            **context,
            "source_id": "toy_source",
        }
        for assay in ("ATAC", "H3K27ac", "Hi-C")
    ]
    signals = [
        {
            "element_id": u["element_id"],
            "sample_id": f"S_{assay}",
            "assay": assay,
            "signal": value,
            "measurement_status": "observed",
            "callable_fraction": 1,
            "unit": "toy_signal",
            "normalization_id": "toy_mean",
            "window_id": "grid:500:mean",
        }
        for u, value in zip(units, [4, 2, 1], strict=True)
        for assay in ("ATAC", "H3K27ac")
    ]
    contacts = [
        {
            "element_id": u["element_id"],
            "promoter_id": p["promoter_id"],
            "sample_id": "S_Hi-C",
            "contact_value": value,
            "measurement_status": "observed",
            "bin_pair_id": f"{u['element_id']}:{p['promoter_id']}",
            "scale": "toy_contact",
            "resolution": 500,
            "source_id": "toy_source",
        }
        for p, cs in zip(promoters, [[3, 2, 2], [1, 2, 6]], strict=True)
        for u, value in zip(units, cs, strict=True)
    ]
    for name, rows in (
        ("samples", samples),
        ("observed_activity", signals),
        ("observed_contacts", contacts),
    ):
        write_table(dest / f"{name}.tsv", rows)
        cfg["inputs"][name] = f"{name}.tsv"
    (dest / "config.yaml").write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")
    return dest / "config.yaml"


def demo(regime, out):
    from .pipeline import run

    final = Path(out).resolve()
    with output_directory(final) as dest:
        config = create_example(dest / "inputs", regime)
        result = run(config, dest / "results")
        # The demo transaction moves the persistent inputs and outputs together.
        for file in (dest / "results").iterdir():
            if file.suffix in (".json", ".yaml", ".md"):
                file.write_text(
                    file.read_text(encoding="utf-8").replace(str(dest), str(final)),
                    encoding="utf-8",
                )
    return {"regime": regime, "n_candidates": len(result["scores"]), "output": str(final)}
