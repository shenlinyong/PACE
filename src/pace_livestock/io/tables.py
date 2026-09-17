"""Lossless TSV ingestion: explicit NA, no silent type coercion."""

from __future__ import annotations

import csv
import gzip
import math
import re
from pathlib import Path

from ..errors import PaceError


def read_table(path: str | Path, *, required=()) -> list[dict]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)):
            raise PaceError(f"{path}: duplicate column name")
        missing = set(required) - set(fields)
        if missing:
            raise PaceError(f"{path}: missing columns {sorted(missing)}")
        rows = []
        for n, row in enumerate(reader, 2):
            if None in row or any(v is None for v in row.values()):
                raise PaceError(f"{path}:{n}: inconsistent number of TSV fields")
            rows.append({k: None if v in ("NA", "") else v for k, v in row.items()})
    return rows


def write_table(path: str | Path, rows: list[dict], *, fields=None) -> None:
    fields = fields or list(dict.fromkeys(k for row in rows for k in row))
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    k: "NA" if v is None or (isinstance(v, float) and not math.isfinite(v)) else v
                    for k, v in row.items()
                }
            )


def number(value, name: str, *, missing=False, minimum=None, maximum=None) -> float:
    if value is None and missing:
        return math.nan
    if isinstance(value, bool):
        raise PaceError(f"{name}: booleans are not numbers")
    try:
        result = float(value)
    except (ValueError, TypeError) as exc:
        raise PaceError(f"{name}: expected a number, received {value!r}") from exc
    if (
        not math.isfinite(result)
        or (minimum is not None and result < minimum)
        or (maximum is not None and result > maximum)
    ):
        raise PaceError(f"{name}: invalid value {value!r}; finite range [{minimum}, {maximum}]")
    return result


def integer(value, name: str, *, minimum=0) -> int:
    if isinstance(value, bool) or not re.fullmatch(r"[+-]?\d+", str(value)):
        raise PaceError(f"{name}: expected an integer, received {value!r}")
    result = int(value)
    if result < minimum:
        raise PaceError(f"{name}: must be >= {minimum}, received {result}")
    return result


def unique(rows: list[dict], keys: tuple[str, ...], name: str) -> None:
    seen = set()
    for row in rows:
        key = tuple(row.get(k) for k in keys)
        if None in key or key in seen:
            raise PaceError(f"{name}: missing or duplicate key {key}")
        seen.add(key)
