"""Regenerate explicitly synthetic canonical examples; never fetch data or real weights."""

from pathlib import Path

import yaml

from pace_livestock.demo import create_example
from pace_livestock.io.tables import write_table


def main():
    root = Path(__file__).resolve().parents[1] / "examples"
    for regime in ("measured", "hybrid", "genome_only"):
        create_example(root / regime, regime)
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
    fusion = [
        dict(
            assay=assay,
            quality_stratum="default",
            observed=9,
            predicted=3,
            target=40**0.5 - 1,
            split="calibration" if i < 20 else "test",
            group_id=f"{assay}_{i}",
            input_donor=f"D{i}",
            target_donor=f"D{i}",
            input_measurement_id=f"{assay}_low_{i}",
            target_measurement_id=f"{assay}_independent_{i}",
        )
        for assay in ("ATAC", "H3K27ac")
        for i in range(24)
    ]
    write_table(advanced / "fusion.tsv", fusion)
    yaml_file(
        "fusion.yaml",
        dict(
            **meta,
            model_id="synthetic_fitted_fusion",
            data="fusion.tsv",
            calibration_target="individual_state",
            signal_unit="toy_signal",
            normalization_id="toy_mean",
            output_window=500,
            scales={"ATAC": 1, "H3K27ac": 1},
            minimum_samples=20,
            measurement_design="synthetic_independent_same_donor_measurements",
        ),
    )
    sequences = []
    for split, n in [("train", 4), ("validation", 2), ("test", 2)]:
        for i in range(n):
            sequences.append(
                dict(
                    sequence_id=f"{split}_{i}",
                    sequence=("ACGT" if i % 2 else "AAAA") * 2048,
                    split=split,
                    group_id=f"{split}_{i}",
                    chrom=f"{split}_{i}",
                    start=5000,
                    end=5500,
                    ATAC=i + 1,
                    H3K27ac=None if i == 0 else i + 2,
                )
            )
    write_table(advanced / "sequence.tsv", sequences)
    yaml_file(
        "sequence.yaml",
        dict(
            **meta,
            model_id="synthetic_cnn",
            data="sequence.tsv",
            assays=["ATAC", "H3K27ac"],
            signal_unit="toy_signal",
            normalization_id="toy_mean",
            input_length=8192,
            output_window=500,
            seed=17,
            epochs=2,
            batch_size=2,
            channels=[8, 8, 8],
            threads=1,
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
    (analysis / "variants.yaml").write_text(
        yaml.safe_dump(
            dict(run_config="../genome_only/config.yaml", variants="../genome_only/sample.vcf")
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


if __name__ == "__main__":
    main()
