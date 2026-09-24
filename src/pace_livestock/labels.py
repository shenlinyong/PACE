"""One rule for functional perturbation labels, shared by benchmarks, training and eta."""

from .errors import PaceError

DIRECTIONS = ("down", "up", "none")


def classify_label(row) -> tuple[int | None, str]:
    """Return (label, reason): 1 for an enhancing positive, 0 for a powered negative.

    An enhancing positive requires effect_direction=down (the knockdown lowered the
    gene) and a powered negative requires effect_direction=none. Upregulation may be
    repressive or indirect and is excluded. Directions outside down/up/none are data
    errors (for example typos) and are rejected, never silently relabelled.
    """
    direction, status = row.get("effect_direction"), row.get("label_status")
    if direction not in DIRECTIONS:
        raise PaceError(
            f"effect_direction must be one of {', '.join(DIRECTIONS)}; received {direction!r} "
            f"(label {row.get('label_id') or row.get('assayed_region_id') or ''})"
        )
    if direction == "up":
        return None, "potential_repressive_or_complex"
    if status == "enhancing_positive":
        return (1, "used") if direction == "down" else (None, "positive_requires_downregulation")
    if status == "powered_negative":
        return (0, "used") if direction == "none" else (None, "negative_requires_no_effect")
    return None, "unusable_label_status"
