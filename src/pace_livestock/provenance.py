"""Deterministic hashes, strict JSON, atomic output transactions and runtime provenance."""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import math
import os
import platform
import shutil
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator

from .errors import PaceError


def clean(value):
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, Path):
        return str(value)
    if hasattr(value, "item"):
        return clean(value.item())
    return value


def digest(value) -> str:
    return hashlib.sha256(
        json.dumps(clean(value), sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()
    ).hexdigest()


def file_hash(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path: str | Path, value) -> None:
    Path(path).write_text(
        json.dumps(clean(value), indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False)
        + "\n",
        encoding="utf-8",
    )


def read_json(path: str | Path) -> dict:
    def reject_constant(value):
        raise PaceError(f"{path}: non-standard JSON number {value}")

    return json.loads(Path(path).read_text(encoding="utf-8"), parse_constant=reject_constant)


@contextmanager
def output_directory(path: str | Path) -> Iterator[Path]:
    """Commit all outputs together; refuse an existing path, including an empty directory."""
    final = Path(path).resolve()
    if final.exists():
        raise PaceError(f"Output already exists: {final}; choose a new --out directory")
    final.parent.mkdir(parents=True, exist_ok=True)
    temp = Path(tempfile.mkdtemp(prefix=f".{final.name}-", dir=final.parent))
    try:
        yield temp
        if final.exists():
            raise PaceError(f"Output appeared during computation: {final}")
        os.rename(temp, final)
    finally:
        if temp.exists():
            shutil.rmtree(temp)


def environment() -> dict:
    names = [
        "pace-livestock",
        "numpy",
        "PyYAML",
        "torch",
        "scikit-learn",
        "cooler",
        "pyBigWig",
        "pysam",
    ]
    packages = {}
    for name in names:
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            pass
    return {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "packages": packages,
        "threads": {
            k: os.environ.get(k)
            for k in ["OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"]
        },
    }
