#!/usr/bin/env bash
# Compatibility name for the current isolated PACE installer.
set -euo pipefail
pace_root="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
exec "$pace_root/install.sh" "$@"
