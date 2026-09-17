"""Conservative fixed-target haplotype construction for SNVs and flanking short indels."""

from __future__ import annotations

from bisect import bisect_right
from collections import defaultdict

from ..errors import PaceError
from ..io.bed_gtf import read_bed
from ..io.tables import integer, read_table, unique
from ..io.variants import Reference, VariantIndex, is_structural, read_variants


class Callability:
    def __init__(self, intervals):
        by_chrom = defaultdict(list)
        for row in intervals:
            by_chrom[row["chrom"]].append((row["start"], row["end"]))
        self.intervals, self.ends = {}, {}
        for chrom, rows in by_chrom.items():
            merged = []
            for start, end in sorted(rows):
                if merged and start <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
                else:
                    merged.append((start, end))
            self.intervals[chrom], self.ends[chrom] = merged, [r[1] for r in merged]

    def fraction(self, chrom, start, end):
        rows = self.intervals.get(chrom, [])
        i = bisect_right(self.ends.get(chrom, []), start)
        overlap = 0
        for left, right in rows[i:]:
            if left >= end:
                break
            overlap += max(0, min(end, right) - max(start, left))
        return overlap / (end - start)


def build_window(
    reference,
    unit,
    variants,
    *,
    input_length,
    ploidy,
    callable_fraction,
    unrecorded_policy="require_callable",
    reference_only=False,
    context_margin=0,
):
    target_width = unit["end"] - unit["start"]
    if input_length < target_width or (input_length - target_width) % 2:
        raise PaceError("Input and target lengths need a symmetric integer flank")
    flank = (input_length - target_width) // 2
    left, right = unit["start"] - flank - context_margin, unit["end"] + flank + context_margin
    result = {
        "element_id": unit["element_id"],
        "output_target_id": unit["element_id"],
        "status": "resolved",
        "reason": "reference_context_prediction" if reference_only else "individual_sequence",
        "callable_fraction": callable_fraction,
        "assumed_reference_fraction": 0.0,
        "structural_status": "not_assessed",
        "sequences": [],
        "additional_context_bp_each_side": context_margin,
    }

    def unresolved(reason):
        return {**result, "status": "unresolved", "reason": reason, "sequences": []}

    if left < 0 or right > reference.sizes.get(unit["chrom"], 0):
        return unresolved("input_window_out_of_bounds")
    if not reference_only and callable_fraction < 1:
        if unrecorded_policy == "require_callable":
            return unresolved("uncallable_input_window")
        result["assumed_reference_fraction"] = 1 - callable_fraction
    local = sorted(variants, key=lambda r: r["pos0"])
    phase_sets, prev_end = set(), -1
    for v in local:
        pos, end = v["pos0"], v["pos0"] + len(v["ref"])
        if len(v["gt"]) != ploidy:
            return unresolved("ploidy_mismatch")
        if any(a is None for a in v["gt"]):
            return unresolved("missing_genotype")
        if any(a < 0 or a > len(v["alts"]) for a in v["gt"]):
            raise PaceError("VCF genotype allele index outside ALT list")
        if is_structural(v) and any(v["gt"]):
            result["structural_status"] = "reported_sv_unsupported"
            return unresolved("reported_sv_unsupported")
        if pos < left or end > right:
            return unresolved("variant_crosses_input_boundary")
        if reference.fetch(unit["chrom"], pos, end) != v["ref"]:
            raise PaceError(f"VCF REF mismatch at {unit['chrom']}:{pos + 1}")
        if v["filter"] not in ("PASS", "."):
            return unresolved("filtered_variant")
        if len(set(v["gt"])) > 1:
            if not v["phased"]:
                return unresolved("unphased_heterozygote")
            # Missing PS is a single local implicit block; mixed known/unknown blocks are rejected.
            phase_sets.add(v["phase_set"] if v["phase_set"] is not None else "implicit_local_block")
        if pos < prev_end:
            return unresolved("overlapping_variants")
        prev_end = end
        for allele in v["gt"]:
            alt = v["ref"] if allele == 0 else v["alts"][allele - 1]
            if set(alt.upper()) - set("ACGTN"):
                return unresolved("unsupported_allele")
            if len(alt) != len(v["ref"]) and pos < unit["end"] and end > unit["start"]:
                return unresolved("indel_changes_output_target")
    if len(phase_sets) > 1:
        return unresolved("unlinked_phase_blocks")
    base = reference.fetch(unit["chrom"], left, right)
    for h in range(ploidy):
        # Positive added context permits fixed-length recentering after flanking insertions.
        delta_before = sum(
            (len(v["alts"][v["gt"][h] - 1]) - len(v["ref"]))
            for v in local
            if v["gt"][h] and v["pos0"] + len(v["ref"]) <= unit["start"]
        )
        delta_total = sum(
            (len(v["alts"][v["gt"][h] - 1]) - len(v["ref"])) for v in local if v["gt"][h]
        )
        crop_start = context_margin + delta_before
        if crop_start < 0 or context_margin + delta_total < delta_before:
            # Deletions need extra sequence outside the validated/callable input window.
            # Keep the output unavailable instead of padding with reference assumptions.
            return unresolved("indel_requires_additional_context")
        sequence = base
        for v in reversed(local):
            allele = v["gt"][h]
            alt = v["ref"] if allele == 0 else v["alts"][allele - 1]
            i = v["pos0"] - left
            sequence = sequence[:i] + alt + sequence[i + len(v["ref"]) :]
        sequence = sequence[crop_start : crop_start + input_length]
        if len(sequence) != input_length:
            return unresolved("indel_window_length_mismatch")
        result["sequences"].append(sequence)
    result["mapping_status"] = (
        "fixed_target_flanking_indel"
        if any(len(a) != len(v["ref"]) for v in local for a in v["alts"])
        else "identity_or_snv"
    )
    return result


def prepare_windows(cfg: dict, units: list[dict], *, input_length: int):
    genome = cfg["genome"]
    if not genome["reference_path"]:
        raise PaceError("Sequence inference requires genome.reference_path")
    variants = (
        read_variants(genome["variant_path"], sample_id=genome["sample_id"])
        if genome["variant_path"]
        else []
    )
    index = VariantIndex(variants)
    callable_regions = (
        read_bed(genome["callability_path"], source_id="callability")
        if genome["callability_path"]
        else []
    )
    callable_index = Callability(callable_regions)
    reference_only = genome["variant_path"] is None
    if not reference_only and not genome["individual_id"]:
        raise PaceError("Individual variants require genome.individual_id")
    ploidies = {}
    if genome["ploidy_path"]:
        rows = read_table(genome["ploidy_path"], required=["chrom", "ploidy"])
        unique(rows, ("chrom",), "ploidy")
        ploidies = {r["chrom"]: integer(r["ploidy"], "ploidy", minimum=1) for r in rows}
    output = []
    with Reference(genome["reference_path"]) as reference:
        for v in variants:
            if (
                v["chrom"] not in reference.sizes
                or v["pos0"] < 0
                or v["end"] > reference.sizes[v["chrom"]]
            ):
                raise PaceError("VCF coordinates or chromosome names do not match reference")
        for unit in units:
            ploidy = 1 if reference_only else ploidies.get(unit["chrom"])
            if ploidy not in (1, 2):
                raise PaceError(f"Explicit haploid/diploid ploidy required for {unit['chrom']}")
            flank = (input_length - (unit["end"] - unit["start"])) // 2
            left, right = unit["start"] - flank, unit["end"] + flank
            margin = 0
            # Include newly exposed variant context until the required indel margin stabilizes.
            # A finite chromosome bounds expansion; unsupported SVs never drive this expansion.
            while True:
                local = index.query(unit["chrom"], left - margin, right + margin)
                needed = sum(
                    max(abs(len(a) - len(v["ref"])) for a in v["alts"])
                    for v in local
                    if not is_structural(v) and any(v["gt"])
                )
                if needed <= margin:
                    break
                margin = needed
                if left - margin < 0 or right + margin > reference.sizes.get(unit["chrom"], 0):
                    break
            fraction = (
                1.0
                if reference_only
                else callable_index.fraction(unit["chrom"], left - margin, right + margin)
            )
            item = build_window(
                reference,
                unit,
                index.query(unit["chrom"], left - margin, right + margin),
                input_length=input_length,
                ploidy=ploidy,
                callable_fraction=fraction,
                unrecorded_policy=genome["unrecorded_site_policy"],
                reference_only=reference_only,
                context_margin=margin,
            )
            if item["structural_status"] == "not_assessed" and genome["sv_assessed"]:
                item["structural_status"] = "no_reported_sv_in_window"
            output.append(item)
    return output, index
