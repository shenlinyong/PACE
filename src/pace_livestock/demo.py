"""Small, deterministic, genuinely offline fixtures for all three public modes."""

from pathlib import Path

import numpy as np
import yaml

from .io.tables import write_table
from .provenance import file_hash, output_directory, write_json
from .schemas import SCHEMAS


def create_example(destination, regime):
    if regime not in ("measured", "hybrid", "genome_only"):
        raise ValueError("Unknown demonstration regime")
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
        "allocation": {"eta": 0},
        "seed": 17,
    }
    if regime != "genome_only":
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
        if regime == "hybrid":
            signals[-1]["signal"], signals[-1]["measurement_status"] = None, "unmeasured"
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
    if regime != "measured":
        model_dir = dest / "models" / "sequence"
        model_dir.mkdir(parents=True)
        np.savez(
            model_dir / "weights.npz",
            coefficients=np.array([[4.0, 1.0, 2.0, 1.0], [1.0, 3.0, 1.0, 2.0]]),
            bias=np.array([0.2, 0.1]),
        )
        meta = {
            "species": "synthetic",
            "assembly": "toy_assembly",
            "context_id": "toy_tissue",
            "is_synthetic": True,
            "target_level": "individual",
            "training_sources": ["deterministic_fixture_recipe"],
            "calibration_sources": [],
            "test_sources": [],
            "validation": {},
        }
        manifest = {
            **meta,
            "model_id": "synthetic_fixed_sequence",
            "kind": "sequence",
            "adapter": "fixed_linear",
            "assays": ["ATAC", "H3K27ac"],
            "input_length": 8192,
            "output_window": 500,
            "output_type": "normalized_signal",
            "signal_unit": "toy_signal",
            "normalization_id": "toy_mean",
            "weights_file": "weights.npz",
            "weights_sha256": file_hash(model_dir / "weights.npz"),
        }
        write_json(model_dir / "manifest.json", manifest)
        reference = list("ACGT" * 5000)
        reference[5500:6000] = "A" * 500
        reference[6000:6500] = "C" * 500
        seq = "".join(reference)
        (dest / "genome.fa").write_text(
            ">chrToy\n" + "\n".join(seq[i : i + 80] for i in range(0, len(seq), 80)) + "\n",
            encoding="utf-8",
        )
        (dest / "sample.vcf").write_text(
            "##fileformat=VCFv4.2\n##contig=<ID=chrToy,length=20000>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttoy_animal\n"
            "chrToy\t5601\t.\tA\tC\t60\tPASS\t.\tGT:PS\t0|1:1\n",
            encoding="utf-8",
        )
        (dest / "callable.bed").write_text("chrToy\t0\t20000\n", encoding="utf-8")
        write_table(dest / "ploidy.tsv", [{"chrom": "chrToy", "ploidy": 2}])
        cfg["sequence"] = {"model_path": "models/sequence"}
        cfg["genome"] = {
            "individual_id": "toy_animal",
            "reference_path": "genome.fa",
            "variant_path": "sample.vcf",
            "callability_path": "callable.bed",
            "ploidy_path": "ploidy.tsv",
            "sample_id": "toy_animal",
        }
        if regime == "genome_only":
            prior_dir = dest / "models" / "contact"
            prior_dir.mkdir()
            write_json(
                prior_dir / "manifest.json",
                {
                    **meta,
                    "kind": "contact_prior",
                    "model_id": "synthetic_contact",
                    "a": 1.0,
                    "gamma": 1.0,
                    "d_min": 1000.0,
                    "d_ref": 1000.0,
                    "scale": "toy_contact",
                    "resolution": 500,
                },
            )
            cfg["contact"].update(
                {"mode": "prior_only", "prior_path": "models/contact", "allow_prior_fallback": True}
            )
        else:
            fusion_dir = dest / "models" / "fusion"
            fusion_dir.mkdir()
            write_json(
                fusion_dir / "manifest.json",
                {
                    **meta,
                    "kind": "fusion",
                    "model_id": "synthetic_fusion",
                    "calibration_target": "individual_state",
                    "signal_unit": "toy_signal",
                    "normalization_id": "toy_mean",
                    "output_window": 500,
                    "strata": {
                        "default": {
                            assay: {"identifiable": True, "weight": 0.5, "scale": 1.0, "n": 20}
                            for assay in ("ATAC", "H3K27ac")
                        }
                    },
                },
            )
            cfg["fusion"] = {"calibrator_path": "models/fusion"}
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
