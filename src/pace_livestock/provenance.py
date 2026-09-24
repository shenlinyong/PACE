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
import uuid
from contextlib import contextmanager
from contextvars import ContextVar
from pathlib import Path
from typing import Iterator

from .errors import PaceError

_FORCE_OUTPUT = ContextVar("pace_force_output", default=False)


@contextmanager
def output_policy(*, force=False):
    token = _FORCE_OUTPUT.set(force)
    try:
        yield
    finally:
        _FORCE_OUTPUT.reset(token)


def clean(value):
    kind = type(value)
    # Fast path for the scalar types that make up almost every table row.
    if value is None or kind is str or kind is int or kind is bool:
        return value
    if kind is float:
        return value if math.isfinite(value) else None
    if isinstance(value, dict):
        out = {}
        for k, v in value.items():
            vt = type(v)
            if v is None or vt is str or vt is int or vt is bool:
                out[k if type(k) is str else str(k)] = v
            elif vt is float:
                out[k if type(k) is str else str(k)] = v if math.isfinite(v) else None
            else:
                out[k if type(k) is str else str(k)] = clean(v)
        return out
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


def digest_rows(rows) -> str:
    """Order-sensitive SHA-256 of many records without building one large JSON string."""
    h = hashlib.sha256()
    for row in rows:
        h.update(
            json.dumps(
                clean(row), sort_keys=True, separators=(",", ":"), ensure_ascii=False
            ).encode()
        )
        h.update(b"\n")
    return h.hexdigest()


def file_hash(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def software_hash() -> str:
    """Hash the installed Python implementation, independently of its installation path."""
    root = Path(__file__).parent
    return digest(
        [(str(path.relative_to(root)), file_hash(path)) for path in sorted(root.rglob("*.py"))]
    )


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
    """Publish complete results; explicit force preserves a previous PACE result as a backup."""
    if Path(path).is_symlink():
        raise PaceError("Output must not be a symbolic link")
    final = Path(path).resolve()
    previous = None
    if final.exists():
        if not _FORCE_OUTPUT.get():
            raise PaceError(
                f"Output already exists: {final}; choose a new --out directory or --force"
            )
        marker = final / "pace_output.json"
        legacy = final / "run_manifest.json"
        recognized = marker.is_file() and read_json(marker).get("format") == "pace-output-1"
        recognized |= legacy.is_file() and read_json(legacy).get("schema_version") == "pace-1"
        if not recognized:
            raise PaceError("--force only replaces an identifiable PACE output directory")
        previous = final.stat()
    final.parent.mkdir(parents=True, exist_ok=True)
    temp = Path(tempfile.mkdtemp(prefix=f".{final.name}-", dir=final.parent))
    try:
        yield temp
        write_json(temp / "pace_output.json", {"format": "pace-output-1"})
        if previous is not None:
            current = final.stat()
            if (current.st_dev, current.st_ino) != (previous.st_dev, previous.st_ino):
                raise PaceError("Output directory changed during computation")
            backup = final.with_name(final.name + ".backup-" + uuid.uuid4().hex[:8])
            os.rename(final, backup)
            try:
                os.rename(temp, final)
            except OSError:
                os.rename(backup, final)
                raise
        elif final.exists():
            raise PaceError(f"Output appeared during computation: {final}")
        else:
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
