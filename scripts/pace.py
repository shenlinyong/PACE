#!/usr/bin/env python3
"""Compatibility launcher for the single installed PACE command."""

import sys

try:
    from pace_livestock.cli import main
except ModuleNotFoundError as exc:
    if exc.name != "pace_livestock":
        raise
    sys.exit("Install PACE with Python >=3.11: python -m pip install .; then run pace --help")

if __name__ == "__main__":
    raise SystemExit(main())
