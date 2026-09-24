"""Resolve measured activity and explicitly declared contact evidence."""

from __future__ import annotations

import json
import math
from collections import defaultdict

from ..boundary_prior import (
    contact_prior,
    fit_kappa,
    load_boundaries,
    posterior_contact,
    unique_count_pairs,
)
from ..core import bulk_mean
from ..errors import PaceError
from ..io.tables import integer, number
from ..provenance import digest
from .contact import distance_prior, shrink


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


def resolve_activity(t, cfg):
    """Aggregate compatible measurements; missing required assays stay unavailable."""
    samples = {r["sample_id"]: r for r in t["samples"]}
    if cfg["target_level"] == "individual" and len({r["donor_id"] for r in t["samples"]}) > 1:
        raise PaceError("target_level=individual cannot combine observations from multiple donors")
    observed, imported = defaultdict(list), {}
    for row in t["observed_activity"]:
        observed[row["element_id"], row["assay"]].append(row)
    for row in t["resolved_activity"]:
        imported[row["element_id"], row["assay"]] = row
    result = []
    for unit in sorted(t["units"], key=lambda r: r["element_id"]):
        for assay in cfg["activity"]["panel"]:
            key = unit["element_id"], assay
            expected_window = f"grid:{cfg['catalog']['width_bp']}:mean"
            observations = observed[key]
            contract = {(r["unit"], r["normalization_id"], r["window_id"]) for r in observations}
            if len(contract) > 1:
                raise PaceError(
                    f"Incompatible observed activity units/normalization/window for {key}"
                )
            if cfg["catalog"]["profile"] == "canonical_grid" and any(
                r["window_id"] != expected_window for r in observations
            ):
                raise PaceError(f"Observed activity window mismatch for {key}")
            value, used = aggregate_observations(
                observations,
                samples,
                "signal",
                minimum_callable=cfg["activity"]["minimum_callable_fraction"],
            )
            signal_unit, normalization, window = (
                next(iter(contract)) if contract else (None, None, expected_window)
            )
            source = "aggregate" if len(used) > 1 else "observed"
            row = {
                "element_id": key[0],
                "assay": assay,
                "observed_value": value,
                "resolved_value": value,
                "evidence_id": digest({"key": key, "samples": used, "source": source}),
                "evidence_type": source,
                "observation_sample_id": used[0] if len(used) == 1 else None,
                "parent_evidence_ids": ";".join(
                    f"observation:{key[0]}:{sample}:{assay}" for sample in used
                )
                or None,
                "resolution_status": "resolved" if math.isfinite(value) else "unresolved",
                "reason": "resolved" if math.isfinite(value) else "missing_assay",
                "unit": signal_unit,
                "normalization_id": normalization,
                "window_id": window,
                "aggregation_method": "equal_technical_then_biological_then_donor_mean"
                if len(used) > 1
                else "single_source",
                "n_observation_samples": len(used),
                "structural_status": "not_assessed",
            }
            if key in imported:
                row = validate_imported_activity(imported[key], cfg, expected_window)
            result.append(row)
    for assay in cfg["activity"]["panel"]:
        contracts = {
            (r["unit"], r.get("normalization_id"), r["window_id"])
            for r in result
            if r["assay"] == assay and r["resolution_status"] == "resolved"
        }
        if len(contracts) > 1:
            raise PaceError(f"Activity scale/normalization differs across units for {assay}")
    return result


def validate_imported_activity(row, cfg, expected_window):
    row = dict(row)
    if row["evidence_type"] not in ("observed", "aggregate"):
        raise PaceError("PACE requires observed or aggregated measured activity")
    if not all(row.get(k) for k in ("unit", "normalization_id", "window_id")):
        raise PaceError("Imported activity requires unit, normalization_id and window_id")
    if cfg["catalog"]["profile"] == "canonical_grid" and row["window_id"] != expected_window:
        raise PaceError("Imported activity window differs from canonical target")
    if row["resolution_status"] != "resolved":
        row["resolved_value"] = math.nan
    row.setdefault("structural_status", "not_assessed")
    return row


def contact_measurement_contract(t, cfg, *, prior_asset=None):
    """Freeze the measurement definition before pooling or comparing contacts.

    A scale name alone does not make different Hi-C resolutions interchangeable.
    Optional normalization/balancing/window identifiers remain unspecified for
    legacy tables, but an unspecified value never matches an explicitly different
    measurement definition. Config values can declare the contract for every row.
    """
    fields = ("resolution", "scale", "normalization_id", "balancing", "window_id")
    declared = {field: cfg["contact"].get(field) for field in fields}
    sources = {row["source_id"]: row for row in t.get("sources", [])}

    def normalize(row, *, source, fallback=None):
        result = {}
        for field in fields:
            value = row.get(field)
            if value is None:
                value = declared[field]
            if value is None and fallback is not None:
                value = fallback[field]
            if field == "resolution" and value is not None:
                value = integer(value, f"{source}.resolution", minimum=1)
            if declared[field] is not None and value != declared[field]:
                raise PaceError(f"Contact {field} differs from the run measurement contract")
            result[field] = value
        return result

    observed = []
    for row in t["observed_contacts"]:
        source_normalization = sources.get(row.get("source_id"), {}).get("normalization_id")
        row_normalization = row.get("normalization_id") or declared["normalization_id"]
        if (
            source_normalization is not None
            and row_normalization is not None
            and source_normalization != row_normalization
        ):
            raise PaceError("Contact normalization_id differs from its declared source")
        observed.append(
            normalize(
                {**row, "normalization_id": row_normalization or source_normalization},
                source="observed contact",
            )
        )
    contracts = {tuple(row[field] for field in fields) for row in observed}
    if len(contracts) > 1:
        raise PaceError(
            "Incompatible observed contact resolution/scale/normalization/balancing/window; "
            "harmonize contacts before aggregation"
        )
    observed_contract = observed[0] if observed else None
    prior_contract = normalize(prior_asset, source="contact prior") if prior_asset else None
    if prior_contract and prior_asset.get("resolution") is None:
        raise PaceError("Contact prior must declare its fitted resolution")
    if observed_contract and prior_contract and observed_contract != prior_contract:
        raise PaceError(
            "Contact prior resolution/scale/normalization/balancing/window does not match "
            "observed contacts"
        )
    baseline = observed_contract or prior_contract
    imported_contracts = []
    for row in t["resolved_contacts"]:
        imported = normalize(row, source="imported contact", fallback=baseline)
        if imported["resolution"] is None:
            raise PaceError("Imported contact requires a declared resolution")
        if baseline is not None and imported != baseline:
            raise PaceError("Imported contact measurement contract differs from run contacts")
        imported_contracts.append(imported)
    if len({tuple(row[field] for field in fields) for row in imported_contracts}) > 1:
        raise PaceError("Imported contacts have incompatible measurement contracts")
    if baseline is None and imported_contracts:
        baseline = imported_contracts[0]
    return {"contract_version": 1, **(baseline or normalize({}, source="contact"))}


def resolve_contacts(t, cfg, *, prior_asset=None):
    measurement = contact_measurement_contract(t, cfg, prior_asset=prior_asset)
    samples = {r["sample_id"]: r for r in t["samples"]}
    units = {r["element_id"]: r for r in t["units"]}
    promoters = defaultdict(list)
    for p in t["promoters"]:
        promoters[p["gene_id"]].append(p)
    promoter_lookup = {p["promoter_id"]: p for p in t["promoters"]}
    same_tss_bin = defaultdict(list)
    observed, imported, shared = defaultdict(list), {}, {}
    for row in t["observed_contacts"]:
        observed[row["element_id"], row["promoter_id"]].append(row)
        p = promoter_lookup[row["promoter_id"]]
        same_tss_bin[row["element_id"], p["chrom"], p["tss0"] // row["resolution"]].append(row)
        key = row["sample_id"], row["bin_pair_id"], row["resolution"]
        value = row["contact_value"] if row["measurement_status"] == "observed" else None
        if cfg["contact"]["reliability"] == "per_pair" and value is not None:
            value = (value, row.get("raw_count"), row.get("count_to_contact"))
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
    pseudocount_policy = cfg["contact"]["pseudocount"]
    if pseudocount_policy == "powerlaw" and (not prior_asset or mode != "observed"):
        raise PaceError(
            "Power-law pseudocount requires observed mode and a compatible fitted prior"
        )
    r = cfg["contact"]["reliability"]
    if mode == "shrinkage" and (r is None or not cfg["contact"]["reliability_source"]):
        raise PaceError("Shrinkage needs explicit reliability and reliability_source")
    per_pair = r == "per_pair"
    boundaries = load_boundaries(prior_asset) if prior_asset else None
    kappa, kappa_source = None, None
    if per_pair:
        if not t["observed_contacts"]:
            raise PaceError(
                "per_pair shrinkage requires raw observed_contacts, including when reusing resolved contacts"
            )
        kappa = cfg["contact"]["kappa"]
        if kappa == "auto":
            kappa = prior_asset.get("kappa")
            kappa_source = "prior_asset"
        else:
            kappa_source = "configuration"
        if kappa is None:
            raw = [
                {
                    **row,
                    "chrom": units[row["element_id"]]["chrom"],
                    "anchor0": units[row["element_id"]]["anchor0"],
                    "tss0": promoter_lookup[row["promoter_id"]]["tss0"],
                }
                for row in t["observed_contacts"]
            ]
            fitting = unique_count_pairs(
                raw, resolution=measurement["resolution"], boundaries=boundaries
            )
            expected = [
                distance_prior(x["distance_bp"], prior_asset)
                * math.exp(-prior_asset.get("beta", 0) * x["boundary_strength"])
                / x["count_to_contact"]
                for x in fitting
            ]
            kappa = fit_kappa([x["raw_count"] for x in fitting], expected)["kappa"]
            kappa_source = "run_unique_bin_pairs"
        kappa = number(kappa, "kappa", minimum=0)
        if kappa == 0:
            raise PaceError("kappa must be positive")
    output, done = [], set()
    for edge in t["candidates"]:
        unit = units[edge["element_id"]]
        for p in promoters[edge["gene_id"]]:
            key = unit["element_id"], p["promoter_id"]
            if key in done:
                continue
            done.add(key)
            distance = abs(unit["anchor0"] - p["tss0"])
            observations = observed[key]
            # Reuse the same measured bin without counting duplicate TSS queries as replicates.
            if not observations and measurement["resolution"]:
                candidates = same_tss_bin[
                    key[0], p["chrom"], p["tss0"] // measurement["resolution"]
                ]
                by_sample = {}
                for candidate in candidates:
                    sid = candidate["sample_id"]
                    if sid in by_sample and (
                        candidate["contact_value"] != by_sample[sid]["contact_value"]
                        and candidate["measurement_status"]
                        == by_sample[sid]["measurement_status"]
                        == "observed"
                    ):
                        raise PaceError("Same-bin TSS contact measurements disagree")
                    if sid not in by_sample or candidate["measurement_status"] == "observed":
                        by_sample[sid] = candidate
                observations = list(by_sample.values())
            obs, used = aggregate_observations(observations, samples, "contact_value")
            resolution = measurement["resolution"]
            left, right = unit["anchor0"], p["tss0"]
            if per_pair:
                left, right = (
                    x // resolution * resolution + resolution // 2 for x in (left, right)
                )
            prior = (
                contact_prior(unit["chrom"], left, right, prior_asset, boundaries)
                if prior_asset
                else math.nan
            )
            posteriors = []
            same_bin = (
                resolution is not None and unit["anchor0"] // resolution == p["tss0"] // resolution
            )
            near = same_bin or distance < cfg["contact"]["near_diagonal_bp"]
            reason = "resolved"
            near_method = None
            if near:
                near_policy = cfg["contact"]["near_diagonal_policy"]
                value, source, weight = (
                    (prior, "contact_prior", 0.0)
                    if near_policy != "unresolved"
                    else (math.nan, "observed", math.nan)
                )
                if near_policy == "prior_or_neighbor" and not math.isfinite(value) and same_bin:
                    corrected = [
                        {
                            **r,
                            "contact_value": r["near_diagonal_value"],
                            "measurement_status": "observed",
                        }
                        for r in observations
                        if r.get("near_diagonal_method") == "neighbor_max"
                        and r.get("near_diagonal_value") is not None
                        and math.isfinite(r["near_diagonal_value"])
                    ]
                    value, used = aggregate_observations(corrected, samples, "contact_value")
                    source, weight = "aggregate", 1.0
                    near_method = "neighbor_max" if math.isfinite(value) else None
                reason = (
                    "near_diagonal_neighbor_max"
                    if near_method
                    else "near_diagonal_prior"
                    if math.isfinite(value)
                    else "near_diagonal_unresolved"
                )
            elif mode == "prior_only":
                value, source, weight = prior, "contact_prior", 0.0
            elif per_pair:
                for observation in observations:
                    if observation["measurement_status"] != "observed":
                        continue
                    posterior = posterior_contact(
                        observation.get("raw_count"),
                        observation.get("count_to_contact"),
                        prior,
                        kappa,
                    )
                    if not math.isclose(
                        observation["contact_value"],
                        float(observation["raw_count"]) * float(observation["count_to_contact"]),
                        rel_tol=1e-8,
                        abs_tol=1e-12,
                    ):
                        raise PaceError("Raw counts/factor do not reproduce the measured contact")
                    posteriors.append(
                        {
                            **posterior,
                            "sample_id": observation["sample_id"],
                            "raw_count": float(observation["raw_count"]),
                            "count_to_contact": float(observation["count_to_contact"]),
                            "measurement_status": "observed",
                        }
                    )
                if posteriors:
                    value, used = aggregate_observations(posteriors, samples, "resolved_value")
                    weight, _ = aggregate_observations(posteriors, samples, "reliability")
                    source = "fused"
                else:
                    value, source, weight = prior, "contact_prior", 0.0
                    reason = "unavailable_bin_prior"
            elif mode == "shrinkage":
                value, source, weight = shrink(obs, prior, r), "fused", r
            elif math.isfinite(obs):
                value, source, weight = obs, "observed", 1.0
            elif cfg["contact"]["allow_prior_fallback"]:
                value, source, weight = prior, "contact_prior", 0.0
            else:
                value, source, weight = math.nan, "observed", math.nan
            structural = "not_assessed"
            row = {
                "element_id": key[0],
                "promoter_id": key[1],
                "observed_value": obs,
                "prior_value": prior,
                "prior_coordinate_policy": "bin_centers"
                if per_pair
                else prior_asset.get(
                    "fitting_coordinate_policy",
                    prior_asset.get("prior_coordinate_policy", "genomic_anchors"),
                )
                if prior_asset
                else None,
                "resolved_value": value,
                "evidence_id": digest([key, source, used]),
                "evidence_type": "aggregate" if source == "observed" and len(used) > 1 else source,
                "observation_sample_id": used[0] if len(used) == 1 and weight > 0 else None,
                "parent_evidence_ids": ";".join(
                    sorted(
                        {
                            f"contact:{r['element_id']}:{r['promoter_id']}:{r['sample_id']}"
                            for r in observations
                            if r["sample_id"] in used
                        }
                    )
                )
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
                "bin_pair_id": ";".join(sorted({r["bin_pair_id"] for r in observations})) or None,
                "resolution_status": "resolved" if math.isfinite(value) else "unresolved",
                "reason": reason
                if math.isfinite(value) or reason != "resolved"
                else "missing_contact",
                "scale": cfg["contact"]["scale"],
                **{k: v for k, v in measurement.items() if k != "contract_version"},
                "distance_bp": distance,
                "near_diagonal_method": near_method,
                "shared_tss_bin": any(r["promoter_id"] != key[1] for r in observations),
                "pseudocount_value": 0.0,
                "structural_status": structural,
            }
            if per_pair:
                row.update(
                    kappa=kappa,
                    kappa_source=kappa_source,
                    prior_coordinate_policy="bin_centers",
                    posterior_samples=json.dumps(posteriors, sort_keys=True),
                )
            if key in imported and per_pair:
                old = imported[key]
                for field in ("resolved_value", "reliability", "kappa", "prior_value"):
                    if not math.isclose(
                        number(old.get(field), field), row[field], rel_tol=1e-10, abs_tol=1e-12
                    ):
                        raise PaceError(
                            "Imported per-pair posterior differs from raw-data recomputation"
                        )
                if (
                    old.get("prior_id") != row["prior_id"]
                    or old.get("posterior_samples") != row["posterior_samples"]
                ):
                    raise PaceError(
                        "Imported per-pair posterior has different prior/sample provenance"
                    )
            if key in imported and not per_pair:
                imported_row = dict(imported[key])
                if imported_row["scale"] != cfg["contact"]["scale"]:
                    raise PaceError("Imported contact scale mismatch")
                imported_source = imported_row["evidence_type"]
                allowed_sources = (
                    {"contact_prior"}
                    if mode == "prior_only"
                    else {"fused"}
                    if mode == "shrinkage"
                    else {"observed", "aggregate"}
                )
                if (near and cfg["contact"]["near_diagonal_policy"] != "unresolved") or (
                    mode == "observed" and cfg["contact"]["allow_prior_fallback"]
                ):
                    allowed_sources.add("contact_prior")
                if mode == "observed" and pseudocount_policy != "none" and prior_asset:
                    allowed_sources.add("regularized")
                if imported_source not in allowed_sources:
                    raise PaceError(
                        "Imported contact evidence conflicts with contact.mode/fallback policy"
                    )
                if imported_row["resolution_status"] == "resolved":
                    imported_weight = number(
                        imported_row["reliability"],
                        "imported contact reliability",
                        minimum=0,
                        maximum=1,
                    )
                    expected_weight = (
                        0.0
                        if imported_source == "contact_prior"
                        else 1.0
                        if imported_source in ("observed", "aggregate", "regularized")
                        else r
                    )
                    if imported_weight != expected_weight:
                        raise PaceError(
                            "Imported contact reliability conflicts with the declared evidence policy"
                        )
                    expected_mode = (
                        "regularized"
                        if imported_source == "regularized"
                        else "prior_only"
                        if imported_weight == 0
                        else "observed"
                        if imported_weight == 1
                        else "shrinkage"
                    )
                    if imported_row["resolved_mode"] != expected_mode:
                        raise PaceError(
                            "Imported contact resolved_mode conflicts with its reliability"
                        )
                    imported_row["reliability"] = imported_weight
                if imported_row["evidence_type"] in ("contact_prior", "fused", "regularized") and (
                    not prior_asset or imported_row["prior_id"] != prior_asset["model_id"]
                ):
                    raise PaceError("Imported contacts need a matching prior asset")
                imported_row["distance_bp"] = distance
                imported_row.update(
                    {k: v for k, v in measurement.items() if k != "contract_version"}
                )
                imported_row["structural_status"] = structural
                if imported_row["resolution_status"] != "resolved" or structural != "not_assessed":
                    imported_row["resolved_value"] = math.nan
                    imported_row["resolution_status"] = "unresolved"
                    if structural != "not_assessed":
                        imported_row["reason"] = reason
                if near and (
                    cfg["contact"]["near_diagonal_policy"] == "unresolved"
                    or not (
                        imported_row["evidence_type"] == "contact_prior"
                        or (
                            cfg["contact"]["near_diagonal_policy"] == "prior_or_neighbor"
                            and imported_row["evidence_type"] == "aggregate"
                            and imported_row.get("near_diagonal_method") == "neighbor_max"
                        )
                    )
                ):
                    imported_row["resolved_value"] = math.nan
                    imported_row["resolution_status"] = "unresolved"
                    imported_row["reason"] = "near_diagonal_unresolved"
                row = imported_row
            # Add a distance-dependent pseudocount once, only on the fitted measurement scale.
            if (
                mode == "observed"
                and pseudocount_policy != "none"
                and prior_asset
                and not near
                and row["resolution_status"] == "resolved"
                and row["evidence_type"] in ("observed", "aggregate", "regularized")
            ):
                pc = cfg["contact"]["pseudocount_strength"] * min(
                    prior, distance_prior(cfg["contact"]["pseudocount_distance_bp"], prior_asset)
                )
                if row["evidence_type"] == "regularized":
                    raw = number(row.get("observed_value"), "regularized observed_value", minimum=0)
                    if not math.isclose(
                        number(row.get("pseudocount_value"), "pseudocount_value", minimum=0),
                        pc,
                        rel_tol=1e-10,
                    ) or not math.isclose(row["resolved_value"], raw + pc, rel_tol=1e-10):
                        raise PaceError(
                            "Imported pseudocount differs from the declared prior or raw contact"
                        )
                elif pc > 0:
                    raw = row["resolved_value"]
                    row.update(
                        observed_value=raw,
                        resolved_value=raw + pc,
                        prior_id=prior_asset["model_id"],
                        prior_value=prior,
                        evidence_type="regularized",
                        resolved_mode="regularized",
                        reason="powerlaw_pseudocount",
                    )
                row["pseudocount_value"] = pc
            row["evidence_id"] = digest({"key": key, "resolution": row})
            output.append(row)
    return sorted(output, key=lambda r: (r["element_id"], r["promoter_id"]))
