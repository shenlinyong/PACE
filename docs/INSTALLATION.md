# Installation

PACE requires Python 3.11 or newer. Install from this repository; no PyPI publication
is implied. A base installation runs canonical table scoring and the offline
synthetic examples without a GPU.

```bash
git clone https://github.com/shenlinyong/PACE.git
cd PACE
python -m venv .venv
source .venv/bin/activate
python -m pip install .
PACE --version
PACE run --help
```

Alternatively use the isolated prefix installer:

```bash
./install.sh --python python3 --prefix "$HOME/.local"
export PATH="$HOME/.local/bin:$PATH"
PACE --version
```

Choose a Python >=3.11 interpreter with `--python /path/to/python` when the system
Python is older. The installer creates `PREFIX/share/pace/venv` and
`PREFIX/bin/PACE`, refusing to replace unrelated environments or executables.
`setup.sh` delegates to the same installer. For an existing Conda environment,
install the package with its Python; `environment.yml` creates a minimal current
PACE environment when invoked from the repository root.

| Extra | Installation in an activated environment | Purpose |
|---|---|---|
| IO | `python -m pip install '.[io]'` | bigWig, cool/mcool, BCF and associated adapters |
| Sequence | `python -m pip install '.[sequence]'` | Quantitative CNN training/inference |
| ML ecosystem | `python -m pip install '.[ml]'` | Optional scientific ML dependencies |
| Development | `python -m pip install -e '.[dev,io,sequence,ml]'` | Full current tests and build tools |

The prefix installer accepts `--extras io,sequence,ml`. Select a PyTorch CPU/GPU
build suitable for your environment. Installation may need network access for
packages; scoring does not download weights or invoke remote prediction services.
FASTA/TSV input and the deterministic synthetic demos use only base dependencies.

Bedtools, samtools, peak callers and read-alignment pipelines can prepare upstream
data; they are not required to execute the installed PACE scoring kernel. No
Snakemake workflow or automatic raw-read analysis is distributed in the current
branch. See [input preparation](input_preparation.md) and [quick start](QUICKSTART.md).
