"""Regenerate explicitly synthetic canonical examples; never fetch data or real weights."""

from pathlib import Path
from tempfile import TemporaryDirectory

import yaml

from pace_livestock.config import load_config
from pace_livestock.demo import create_example
from pace_livestock.io.tables import read_table, write_table
from pace_livestock.pipeline import compute
from pace_livestock.provenance import clean


def classifier_contract():
    """Match the published H3K4me1 example to the actual measured pipeline."""
    with TemporaryDirectory() as directory:
        cfg = load_config(create_example(Path(directory) / "inputs"))
        samples = read_table(cfg["inputs"]["samples"])
        samples.append({**samples[0], "sample_id": "S_H3K4me1", "assay": "H3K4me1"})
        write_table(cfg["inputs"]["samples"], samples)
        observations = read_table(cfg["inputs"]["observed_activity"])
        for element, value in (("E1", 0), ("E2", 5), ("E3", 10)):
            observations.append(
                {
                    **observations[0],
                    "element_id": element,
                    "sample_id": "S_H3K4me1",
                    "assay": "H3K4me1",
                    "signal": value,
                }
            )
        write_table(cfg["inputs"]["observed_activity"], observations)
        return clean(compute(cfg)["manifest"]["ml_feature_contract"])


def main(destination=None):
    root = Path(destination) if destination else Path(__file__).resolve().parents[1] / "examples"
    create_example(root / "measured")
    contact_examples(root / "contact")
    advanced = root / "training"
    advanced.mkdir(parents=True, exist_ok=True)

    def yaml_file(name, content):
        (advanced / name).write_text(yaml.safe_dump(content, sort_keys=False), encoding="utf-8")

    meta = dict(
        species="synthetic",
        assembly="toy_assembly",
        context_id="toy_tissue",
        target_level="individual",
        is_synthetic=True,
    )
    contact = [
        dict(
            bin_pair_id=f"{d}_{i}", distance_bp=d, contact_value=v, split="train", region_id="train"
        )
        for d, c in [(1000, 8), (2000, 4), (4000, 2)]
        for i, v in enumerate([c, 0])
    ]
    contact.append(
        dict(
            bin_pair_id="test",
            distance_bp=3000,
            contact_value=4 / 3,
            split="test",
            region_id="test",
        )
    )
    write_table(advanced / "contacts.tsv", contact)
    yaml_file(
        "contact.yaml",
        dict(
            **meta,
            model_id="synthetic_fitted_contact",
            data="contacts.tsv",
            scale="toy_contact",
            resolution=500,
            bin_edges=[500, 1500, 3000, 5000],
            d_ref=1000,
            d_min=500,
        ),
    )
    learning = []
    for split, n in [("train", 6), ("calibration", 2), ("test", 2)]:
        for group in range(n):
            for label in (0, 1):
                identifier = f"{split}_{group}_{label}"
                learning.append(
                    dict(
                        element_id=identifier,
                        gene_id="G",
                        assayed_region_id=identifier,
                        mapping_count=1,
                        group_id=f"{split}_{group}",
                        split=split,
                        label_status="enhancing_positive" if label else "powered_negative",
                        effect_direction="down" if label else "none",
                        A_used=1,
                        Cbar=1,
                        distance_bp=100,
                        pace_score=0.5,
                        regime="measured",
                        activity_sources="observed",
                        contact_sources="observed",
                        H3K4me1=10 * label,
                    )
                )
    write_table(advanced / "learning.tsv", learning)
    yaml_file(
        "learning.yaml",
        dict(
            model_id="synthetic_classifier",
            data="learning.tsv",
            is_synthetic=True,
            context={k: meta[k] for k in ("species", "assembly", "context_id")},
            extra_features=["H3K4me1"],
            penalties=[[0.01, 0.01], [0.02, 0.01]],
            folds=3,
            seed=17,
            calibrate=True,
            feature_contract=classifier_contract(),
        ),
    )
    analysis = root / "analysis"
    analysis.mkdir(exist_ok=True)
    labels = [
        dict(
            label_id=f"L{i}",
            assayed_region_id=e,
            gene_id=g,
            context_id="toy_tissue",
            perturbation_type="synthetic_inhibition",
            effect_direction="down" if i % 2 else "none",
            effect_size=-1 if i % 2 else 0,
            label_status="enhancing_positive" if i % 2 else "powered_negative",
            assay_id="toy_assay",
            group_id=e,
            source_id="toy_source",
        )
        for i, (e, g) in enumerate([(e, g) for e in ("E1", "E2", "E3") for g in ("G1", "G2")])
    ]
    write_table(analysis / "labels.tsv", labels)
    write_table(
        analysis / "membership.tsv",
        [
            dict(
                region_id=e,
                element_id=e,
                source_id="toy_source",
                membership_rule="synthetic_one_to_one",
            )
            for e in ("E1", "E2", "E3")
        ],
    )
    (analysis / "benchmark.yaml").write_text(
        yaml.safe_dump(
            dict(
                run_config="../measured/config.yaml",
                labels="labels.tsv",
                region_membership="membership.tsv",
                stratify=["gene_id"],
            )
        ),
        encoding="utf-8",
    )
    (analysis / "comparison.yaml").write_text(
        yaml.safe_dump(dict(left="../../results/measured", right="../../results/measured_repeat")),
        encoding="utf-8",
    )
    write_table(
        analysis / "replicates.tsv",
        [
            dict(run_path=f"../../results/{p}", donor_id="toy_animal", replicate_type="technical")
            for p in ("measured", "measured_repeat")
        ],
    )
    (analysis / "stability.yaml").write_text("replicates: replicates.tsv\n", encoding="utf-8")


def contact_examples(root):
    """Small, synthetic four-chromosome inputs for the complete contact workflow."""
    create_example(root)
    cfg = yaml.safe_load((root / "config.yaml").read_text())
    originals = {
        k: read_table(root / f"{k}.tsv")
        for k in ("units", "promoters", "candidates", "observed_activity", "observed_contacts")
    }
    batches = {k: [] for k in originals}
    variants, motifs = [], []
    for chrom in ("c1", "c2", "c3", "c4"):
        for row, start in zip(originals["units"], [5000, 20000, 30000], strict=True):
            batches["units"].append(
                {
                    **row,
                    "element_id": chrom + row["element_id"],
                    "chrom": chrom,
                    "start": start,
                    "end": start + 500,
                    "anchor0": start + 249,
                }
            )
        for row, tss in zip(originals["promoters"], [10000, 40000], strict=True):
            batches["promoters"].append(
                {
                    **row,
                    "gene_id": chrom + row["gene_id"],
                    "promoter_id": chrom + row["promoter_id"],
                    "chrom": chrom,
                    "tss0": tss,
                }
            )
        batches["candidates"] += [
            {**r, "element_id": chrom + r["element_id"], "gene_id": chrom + r["gene_id"]}
            for r in originals["candidates"]
        ]
        batches["observed_activity"] += [
            {
                **r,
                "element_id": chrom + r["element_id"],
                "signal": {"E1": 1.0, "E2": 2.8, "E3": 1.0}[r["element_id"]],
            }
            for r in originals["observed_activity"]
        ]
        for row in originals["observed_contacts"]:
            batches["observed_contacts"].append(
                {
                    **row,
                    "element_id": chrom + row["element_id"],
                    "promoter_id": chrom + row["promoter_id"],
                    "bin_pair_id": chrom + row["bin_pair_id"],
                    "raw_count": int(row["contact_value"]),
                    "count_to_contact": 1.0,
                    "chrom": chrom,
                    "anchor0": {"E1": 5249, "E2": 20249, "E3": 30249}[row["element_id"]],
                    "tss0": {"P1": 10000, "P2": 40000}[row["promoter_id"]],
                    "normalization_id": "toy_mean",
                    "balancing": "unbalanced",
                    "window_id": "bin_pair",
                }
            )
        variants += [
            dict(
                variant_id=chrom + gene,
                chrom=chrom,
                pos0=pos,
                gene_id=chrom + gene,
                pip=0.9,
                signal_id="signal",
            )
            for gene, pos in [("G1", 5000), ("G2", 30000)]
        ]
        motifs += [
            dict(chrom=chrom, start=x, end=x + 10, strand=strand, strength=0.8)
            for x, strand in [(24900, "-"), (25100, "+")]
        ]
    for key, rows in batches.items():
        write_table(root / f"{key}.tsv", rows)
    write_table(root / "eqtl.tsv", variants)
    write_table(root / "motifs.tsv", motifs)

    def config(name, content):
        (root / name).write_text(yaml.safe_dump(content, sort_keys=False), encoding="utf-8")

    config("boundaries.yaml", dict(motifs="motifs.tsv"))
    prior = dict(
        **cfg["context"],
        target_level=cfg["target_level"],
        is_synthetic=True,
        model_id="synthetic_boundary_prior",
        scale="toy_contact",
        resolution=500,
        normalization_id="toy_mean",
        balancing="unbalanced",
        window_id="bin_pair",
        boundaries="boundary_asset/boundaries.tsv",
        d_ref=1000,
        d_min=1000,
    )
    config("prior.yaml", dict(**prior, a=1, gamma=1, beta=0, kappa=2))
    config(
        "fit-hic.yaml",
        dict(
            **prior,
            data="observed_contacts.tsv",
            gamma_grid=[0.5, 1.0, 1.5],
            beta_grid=[0.0, 0.5, 1.0],
            test_chromosomes=["c4"],
        ),
    )
    cfg["contact"].update(
        mode="shrinkage",
        reliability="per_pair",
        prior_path="prior_asset",
        normalization_id="toy_mean",
        balancing="unbalanced",
        window_id="bin_pair",
    )
    config("fuse.yaml", cfg)
    cfg["contact"].update(mode="prior_only", reliability=None)
    cfg["inputs"]["observed_contacts"] = None
    config("config.yaml", cfg)
    config(
        "fit-labels.yaml",
        dict(
            run_config="config.yaml",
            labels="eqtl.tsv",
            **cfg["context"],
            gamma_grid=[1.0],
            beta_grid=[0.0, 0.5],
            eta_grid=[0.0, 0.5, 1.0],
            aggregation="independent_signals",
            test_chromosomes=["c4"],
        ),
    )
    cfg["allocation"]["weak_model_path"] = "weak_fit/weak_model.json"
    config("calibrated.yaml", cfg)


if __name__ == "__main__":
    main()
