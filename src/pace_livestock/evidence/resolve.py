"""Resolve real observations, quantitative predictions and explicit contact priors."""

from __future__ import annotations

import math
from collections import defaultdict

from ..core import bulk_mean
from ..errors import PaceError
from ..io.variants import is_structural
from ..provenance import digest
from .contact import distance_prior, shrink
from .fusion import resolve_signal


def aggregate_observations(rows, samples, value_column, *, minimum_callable=0):
    """Equal technical means within biological replicate, then donor means.

    This is a declared equal-weight policy for already normalized signals. Raw count
    pooling must happen upstream. No donor is counted once per technical replicate.
    """
    by_donor = defaultdict(lambda: defaultdict(list))
    used_samples = []
    for row in rows:
        if (
            row["measurement_status"] != "observed"
            or row.get("callable_fraction", 1) < minimum_callable
        ):
            continue
        sample = samples[row["sample_id"]]
        value = row[value_column]
        if not math.isfinite(value):
            continue
        by_donor[sample["donor_id"]][sample["biological_replicate"]].append(value)
        used_samples.append(row["sample_id"])
    donor_means = []
    for reps in by_donor.values():
        donor_means.append(float(bulk_mean([float(bulk_mean(v)) for v in reps.values()])))
    return (float(bulk_mean(donor_means)) if donor_means else math.nan, sorted(set(used_samples)))


def resolve_activity(t, cfg, *, predictions, sequence_asset=None, fusion_asset=None):
    samples = {r["sample_id"]: r for r in t["samples"]}
    if cfg["target_level"] == "individual" and len({r["donor_id"] for r in t["samples"]}) > 1:
        raise PaceError("target_level=individual cannot combine observations from multiple donors")
    if cfg["target_level"] == "individual" and cfg["genome"]["variant_path"] and t["samples"]:
        if {r["donor_id"] for r in t["samples"]} != {cfg["genome"]["individual_id"]}:
            raise PaceError(
                "Hybrid observations and individual genotype must identify the same donor"
            )
    obs, pred, imported = defaultdict(list), {}, {}
    for row in t["observed_activity"]:
        obs[row["element_id"], row["assay"]].append(row)
    for row in predictions:
        key = row["element_id"], row["assay"]
        if key in pred:
            raise PaceError("Duplicate prediction for a unit/assay; average copies before import")
        pred[key] = row
    for row in t["resolved_activity"]:
        imported[row["element_id"], row["assay"]] = row
    if fusion_asset:
        wanted = "individual_state" if cfg["target_level"] == "individual" else "population_mean"
        if fusion_asset.get("calibration_target") != wanted:
            raise PaceError("Fusion calibration target does not match the run estimand")
    result = []
    for unit in sorted(t["units"], key=lambda r: r["element_id"]):
        for assay in cfg["activity"]["panel"]:
            key = unit["element_id"], assay
            expected_window = f"grid:{cfg['catalog']['width_bp']}:mean"
            observations = obs.get(key, [])
            p = pred.get(key)
            contract = {(r["unit"], r["normalization_id"], r["window_id"]) for r in observations}
            if len(contract) > 1:
                raise PaceError(
                    f"Incompatible observed activity units/normalization/window for {key}"
                )
            if cfg["catalog"]["profile"] == "canonical_grid" and any(
                r["window_id"] != expected_window for r in observations
            ):
                raise PaceError(f"Observed activity window mismatch for {key}")
            if p:
                if sequence_asset is None or p["model_id"] != sequence_asset["model_id"]:
                    raise PaceError("Predictions require a matching sequence model manifest")
                wanted = (
                    sequence_asset["signal_unit"],
                    sequence_asset["normalization_id"],
                    expected_window,
                )
                if (p["unit"], p["normalization_id"], p["window_id"]) != wanted or (
                    contract and next(iter(contract)) != wanted
                ):
                    raise PaceError(
                        f"Sequence and observed units/normalization/window mismatch for {key}"
                    )
                contract.add(wanted)
            observed, used = aggregate_observations(
                observations,
                samples,
                "signal",
                minimum_callable=cfg["activity"]["minimum_callable_fraction"],
            )
            predicted = p["predicted_value"] if p and p["status"] == "resolved" else math.nan
            calibration = None
            if fusion_asset:
                calibration = (
                    fusion_asset.get("strata", {})
                    .get(cfg["fusion"]["quality_stratum"], {})
                    .get(assay)
                )
                if calibration is None:
                    raise PaceError(f"Fusion asset lacks assay/quality stratum for {assay}")
                if contract and fusion_asset.get("signal_unit") != next(iter(contract))[0]:
                    raise PaceError("Fusion signal unit mismatch")
                if contract and fusion_asset.get("normalization_id") != next(iter(contract))[1]:
                    raise PaceError("Fusion normalization mismatch")
                if fusion_asset.get("output_window") != cfg["catalog"]["width_bp"]:
                    raise PaceError("Fusion target window mismatch")
            value, source, w = resolve_signal(
                observed, predicted, regime=cfg["regime"], calibration=calibration
            )
            signal_unit, normalization, window_id = (
                next(iter(contract))
                if contract
                else (
                    sequence_asset["signal_unit"] if sequence_asset else None,
                    sequence_asset["normalization_id"] if sequence_asset else None,
                    expected_window,
                )
            )
            evidence_type = "aggregate" if source == "observed" and len(used) > 1 else source
            evidence_id = digest(
                {
                    "key": key,
                    "samples": used,
                    "source": source,
                    "model": sequence_asset.get("model_id") if sequence_asset else None,
                }
            )
            row = {
                "element_id": key[0],
                "assay": assay,
                "observed_value": observed,
                "predicted_value": predicted,
                "resolved_value": value,
                "evidence_id": evidence_id,
                "evidence_type": evidence_type,
                "observation_sample_id": used[0] if len(used) == 1 and w > 0 else None,
                "parent_evidence_ids": ";".join(f"observation:{key[0]}:{s}:{assay}" for s in used)
                if w > 0
                else None,
                "model_id": sequence_asset["model_id"] if sequence_asset and w < 1 else None,
                "calibrator_id": fusion_asset["model_id"]
                if fusion_asset and source == "fused"
                else None,
                "fusion_weight": w,
                "resolution_status": "resolved" if math.isfinite(value) else "unresolved",
                "reason": "resolved"
                if math.isfinite(value)
                else p.get("reason", "missing_assay")
                if p
                else "missing_assay",
                "unit": signal_unit,
                "normalization_id": normalization,
                "window_id": window_id,
                "aggregation_method": "equal_technical_then_biological_then_donor_mean"
                if len(used) > 1
                else "single_source",
                "n_observation_samples": len(used),
                "structural_status": p.get("structural_status", "not_assessed")
                if p
                else "not_assessed",
            }
            if key in imported:
                row = validate_imported_activity(
                    imported[key], cfg, sequence_asset, fusion_asset, expected_window
                )
            result.append(row)
    # Freeze normalization across units, not merely across sources at a single unit.
    for assay in cfg["activity"]["panel"]:
        contracts = {
            (r["unit"], r.get("normalization_id"), r["window_id"])
            for r in result
            if r["assay"] == assay and r["resolution_status"] == "resolved"
        }
        if len(contracts) > 1:
            raise PaceError(f"Activity scale/normalization differs across units for {assay}")
    return result


def validate_imported_activity(row, cfg, sequence_asset, fusion_asset, expected_window):
    row = dict(row)
    source = row["evidence_type"]
    if cfg["regime"] == "measured" and source not in ("observed", "aggregate"):
        raise PaceError("measured mode cannot import predicted/fused activity")
    if cfg["regime"] == "genome_only" and source != "sequence_prediction":
        raise PaceError("genome_only activity imports must be sequence_prediction")
    if cfg["catalog"]["profile"] == "canonical_grid" and row["window_id"] != expected_window:
        raise PaceError("Imported activity window differs from canonical target")
    if source in ("sequence_prediction", "fused"):
        if (
            not sequence_asset
            or row["model_id"] != sequence_asset["model_id"]
            or row["unit"] != sequence_asset["signal_unit"]
        ):
            raise PaceError("Imported prediction requires a matching quantitative model")
        if row.get("normalization_id") != sequence_asset["normalization_id"]:
            raise PaceError("Imported prediction normalization_id mismatch")
    if source == "fused" and (not fusion_asset or row["calibrator_id"] != fusion_asset["model_id"]):
        raise PaceError("Imported fused activity requires its matching calibrator")
    if row["resolution_status"] != "resolved":
        row["resolved_value"] = math.nan
    row.setdefault("structural_status", "not_assessed")
    return row


def resolve_contacts(t, cfg, *, prior_asset=None, variants=None):
    samples = {r["sample_id"]: r for r in t["samples"]}
    units = {r["element_id"]: r for r in t["units"]}
    promoters = defaultdict(list)
    for p in t["promoters"]:
        promoters[p["gene_id"]].append(p)
    observed, imported, shared = defaultdict(list), {}, {}
    for row in t["observed_contacts"]:
        observed[row["element_id"], row["promoter_id"]].append(row)
        key = row["sample_id"], row["bin_pair_id"], row["resolution"]
        value = row["contact_value"] if row["measurement_status"] == "observed" else None
        if key in shared and shared[key] != value:
            raise PaceError("Rows sharing one sample/bin pair disagree on contact measurement")
        shared[key] = value
    for row in t["resolved_contacts"]:
        imported[row["element_id"], row["promoter_id"]] = row
    if prior_asset and prior_asset.get("scale") != cfg["contact"]["scale"]:
        raise PaceError("Contact prior scale does not match observed/run scale")
    mode = cfg["contact"]["mode"]
    if mode in ("prior_only", "shrinkage") and not prior_asset:
        raise PaceError(f"contact.mode={mode} requires a contact prior asset")
    r = cfg["contact"]["reliability"]
    if mode == "shrinkage" and (r is None or not cfg["contact"]["reliability_source"]):
        raise PaceError("Shrinkage needs explicit reliability and reliability_source")
    output, done = [], set()
    for edge in t["candidates"]:
        unit = units[edge["element_id"]]
        for p in promoters[edge["gene_id"]]:
            key = unit["element_id"], p["promoter_id"]
            if key in done:
                continue
            done.add(key)
            distance = abs(unit["anchor0"] - p["tss0"])
            obs, used = aggregate_observations(observed[key], samples, "contact_value")
            prior = distance_prior(distance, prior_asset) if prior_asset else math.nan
            same_bin = any(
                unit["anchor0"] // row["resolution"] == p["tss0"] // row["resolution"]
                for row in observed[key]
            )
            near = same_bin or distance < cfg["contact"]["near_diagonal_bp"]
            reason = "resolved"
            if near:
                value, source, weight = (
                    (prior, "contact_prior", 0.0)
                    if cfg["contact"]["near_diagonal_policy"] == "prior_or_unresolved"
                    else (math.nan, "observed", math.nan)
                )
                reason = (
                    "near_diagonal_prior" if math.isfinite(value) else "near_diagonal_unresolved"
                )
            elif mode == "prior_only":
                value, source, weight = prior, "contact_prior", 0.0
            elif mode == "shrinkage":
                value, source, weight = shrink(obs, prior, r), "fused", r
            elif math.isfinite(obs):
                value, source, weight = obs, "observed", 1.0
            elif cfg["contact"]["allow_prior_fallback"]:
                value, source, weight = prior, "contact_prior", 0.0
            else:
                value, source, weight = math.nan, "observed", math.nan
            structural = "not_assessed"
            if variants:
                between = variants.query(
                    unit["chrom"], min(unit["start"], p["tss0"]), max(unit["end"], p["tss0"] + 1)
                )
                affected = [
                    v
                    for v in between
                    if any(a is None or a != 0 for a in v["gt"])
                    and (is_structural(v) or any(len(a) != len(v["ref"]) for a in v["alts"]))
                ]
                if affected:
                    value, reason, structural = (
                        math.nan,
                        "altered_distance_or_tss_unsupported",
                        "reported_variant_affects_relationship",
                    )
            row = {
                "element_id": key[0],
                "promoter_id": key[1],
                "observed_value": obs,
                "prior_value": prior,
                "resolved_value": value,
                "evidence_id": digest([key, source, used]),
                "evidence_type": "aggregate" if source == "observed" and len(used) > 1 else source,
                "observation_sample_id": used[0] if len(used) == 1 and weight > 0 else None,
                "parent_evidence_ids": ";".join(f"contact:{key[0]}:{key[1]}:{s}" for s in used)
                if weight > 0
                else None,
                "prior_id": prior_asset["model_id"] if prior_asset and weight < 1 else None,
                "reliability": weight,
                "reliability_source": cfg["contact"]["reliability_source"],
                "resolved_mode": "prior_only"
                if weight == 0
                else "observed"
                if weight == 1
                else mode,
                "bin_pair_id": ";".join(sorted({r["bin_pair_id"] for r in observed[key]})) or None,
                "resolution_status": "resolved" if math.isfinite(value) else "unresolved",
                "reason": reason
                if math.isfinite(value) or reason != "resolved"
                else "missing_contact",
                "scale": cfg["contact"]["scale"],
                "distance_bp": distance,
                "structural_status": structural,
            }
            if key in imported:
                imported_row = dict(imported[key])
                if imported_row["scale"] != cfg["contact"]["scale"]:
                    raise PaceError("Imported contact scale mismatch")
                if imported_row["evidence_type"] in ("contact_prior", "fused") and (
                    not prior_asset or imported_row["prior_id"] != prior_asset["model_id"]
                ):
                    raise PaceError("Imported contacts need a matching prior asset")
                if (
                    cfg["regime"] == "genome_only"
                    and imported_row["evidence_type"] != "contact_prior"
                ):
                    raise PaceError("genome_only contact imports must use the declared prior")
                imported_row["distance_bp"] = distance
                imported_row["structural_status"] = structural
                if imported_row["resolution_status"] != "resolved" or structural != "not_assessed":
                    imported_row["resolved_value"] = math.nan
                row = imported_row
            output.append(row)
    return sorted(output, key=lambda r: (r["element_id"], r["promoter_id"]))
