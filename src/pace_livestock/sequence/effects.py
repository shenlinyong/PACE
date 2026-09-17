"""Single-variant REF/ALT scenarios; these are not reconstructed individual genomes."""

import math

from ..config import load_config, operation_config
from ..errors import PaceError
from ..evidence.assets import load_asset
from ..io.tables import write_table
from ..io.variants import Reference, VariantIndex, read_variants
from ..provenance import output_directory, write_json
from ..schemas import load_tables
from .genome import build_window
from .model import predict_windows


def variant_effects_command(path, out):
    cfg = operation_config(
        path,
        allowed={"run_config", "variants", "sample_id"},
        required=["run_config", "variants"],
        paths=["run_config", "variants"],
    )
    run_cfg = load_config(cfg["run_config"])
    tables = load_tables(run_cfg)
    model = load_asset(run_cfg["sequence"]["model_path"], run_cfg, kind="sequence")
    if not run_cfg["genome"]["reference_path"]:
        raise PaceError("Variant scenarios require a reference FASTA")
    variants = read_variants(cfg["variants"], sample_id=cfg.get("sample_id"))
    # Scenario alleles are substitutions at the physical VCF record, not at a
    # synthetic remote endpoint used to mask individual contact relationships.
    index = VariantIndex(variants, include_remote_breakends=False)
    output = []
    with Reference(run_cfg["genome"]["reference_path"]) as reference:
        for unit in tables["units"]:
            flank = (model["input_length"] - (unit["end"] - unit["start"])) // 2
            local = index.query(unit["chrom"], unit["start"] - flank, unit["end"] + flank)
            for variant in local:
                for allele, alt in enumerate(variant["alts"], 1):
                    ref_window = build_window(
                        reference,
                        unit,
                        [],
                        input_length=model["input_length"],
                        ploidy=1,
                        callable_fraction=1,
                        reference_only=True,
                    )
                    alt_window = build_window(
                        reference,
                        unit,
                        [{**variant, "gt": (allele,), "phased": True}],
                        input_length=model["input_length"],
                        ploidy=1,
                        callable_fraction=1,
                        reference_only=True,
                        context_margin=abs(len(alt) - len(variant["ref"]))
                        if not alt.startswith("<")
                        else 0,
                    )
                    ref_predictions = predict_windows(
                        [ref_window], model, max_n_fraction=run_cfg["sequence"]["max_n_fraction"]
                    )
                    alt_predictions = predict_windows(
                        [alt_window], model, max_n_fraction=run_cfg["sequence"]["max_n_fraction"]
                    )
                    for a, b in zip(ref_predictions, alt_predictions, strict=True):
                        output.append(
                            {
                                "element_id": unit["element_id"],
                                "chrom": variant["chrom"],
                                "pos0": variant["pos0"],
                                "ref": variant["ref"],
                                "alt": alt,
                                "assay": a["assay"],
                                "reference_signal": a["predicted_value"],
                                "alternate_signal": b["predicted_value"],
                                "delta_signal": b["predicted_value"] - a["predicted_value"],
                                "full_delta_pace": math.nan,
                                "scenario": "single_variant_on_reference_context",
                                "status": b["status"],
                                "reason": b["reason"]
                                if b["status"] != "resolved"
                                else "sequence_only_scenario",
                            }
                        )
    with output_directory(out) as dest:
        write_table(
            dest / "variant_effects.tsv",
            output,
            fields=None
            if output
            else ["element_id", "chrom", "pos0", "ref", "alt", "assay", "delta_signal", "reason"],
        )
        write_json(
            dest / "report.json",
            {
                "n_variant_assay_scenarios": len(output),
                "individual_reconstruction": False,
                "full_delta_pace": "not_claimed; use compatible complete run comparisons for full deltas",
            },
        )
    return output
