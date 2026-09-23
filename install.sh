#!/usr/bin/env bash
# Install an isolated PACE environment and executable in a user-selected prefix.
set -euo pipefail

pace_prefix="${HOME}/.local"
pace_python="python3"
pace_extras=""
while (($#)); do
    case "$1" in
        --prefix|--python|--extras)
            if (($# < 2)); then
                printf 'Missing value for %s\n' "$1" >&2
                exit 2
            fi
            case "$1" in
                --prefix) pace_prefix="$2" ;;
                --python) pace_python="$2" ;;
                --extras) pace_extras="$2" ;;
            esac
            shift 2 ;;
        -h|--help)
            printf '%s\n' 'Usage: ./install.sh [--prefix ~/.local] [--python python3] [--extras io,ml]'
            printf '%s\n' 'Requires Python >=3.11. Installs into PREFIX/share/pace/venv and PREFIX/bin/PACE.'
            exit 0 ;;
        *) printf 'Unknown option: %s\n' "$1" >&2; exit 2 ;;
    esac
done

if [[ ! "$pace_extras" =~ ^(io|ml|dev)(,(io|ml|dev))*$ && -n "$pace_extras" ]]; then
    printf '%s\n' 'Extras must be a comma-separated subset of io,ml,dev.' >&2
    exit 2
fi
"$pace_python" -c 'import sys; sys.exit("PACE requires Python >=3.11; use --python /path/to/python") if sys.version_info < (3,11) else None'
pace_root="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p -- "$pace_prefix"
pace_prefix="$(CDPATH= cd -- "$pace_prefix" && pwd)"
pace_env="$pace_prefix/share/pace/venv"
pace_bin="$pace_prefix/bin/PACE"
if [[ -e "$pace_bin" || -L "$pace_bin" ]]; then
    if [[ ! -L "$pace_bin" || "$(readlink -- "$pace_bin")" != "$pace_env/bin/PACE" ]]; then
        printf 'Refusing to replace an unmanaged executable: %s\n' "$pace_bin" >&2
        exit 2
    fi
fi
if [[ -e "$pace_env" && ! -f "$pace_env/.pace-managed" ]]; then
    printf 'Refusing to modify an unmanaged environment: %s\n' "$pace_env" >&2
    exit 2
fi
if [[ ! -e "$pace_env" ]]; then
    "$pace_python" -m venv "$pace_env"
    touch "$pace_env/.pace-managed"
fi
pace_requirement="$pace_root"
if [[ -n "$pace_extras" ]]; then
    pace_requirement="$pace_root[$pace_extras]"
fi
"$pace_env/bin/python" -m pip install --upgrade "$pace_requirement"
"$pace_env/bin/PACE" --version
mkdir -p -- "$pace_prefix/bin"
ln -sfn -- "$pace_env/bin/PACE" "$pace_bin"
printf 'Installed: %s\n' "$pace_bin"
printf 'If needed, add to PATH: export PATH="%s/bin:$PATH"\n' "$pace_prefix"
