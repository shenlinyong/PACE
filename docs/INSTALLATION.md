# Download and installation

The supported installation is from this GitHub repository. Commands below target
Linux servers and macOS; on Windows use a Linux environment such as WSL2 or Docker.
Python 3.11 or newer is required; the supplied Conda environment uses Python 3.12.
Installing dependencies requires access to the configured package repositories.
Bundled examples subsequently run offline, without downloading data or weights.

After installation, start with the mode that matches your data:

```bash
PACE measured --help
PACE hybrid --help
PACE genome --help
```

These commands show a short copyable example before the complete parameter reference. The examples use filenames such as `animal_001_liver_ATAC_H3K27ac.tsv` to make the expected assay and sample identity visible; replace them with your own files. PACE accepts one, two, three or more biological replicates when each sample is registered in `samples.tsv` and the same `sample_id` is used in the measurement tables.

## 1. Download the source

With Git:

```bash
git clone https://github.com/shenlinyong/PACE.git
cd PACE
git rev-parse HEAD
```

Record the printed commit for reproducibility. `main` changes as development
continues. To repeat an analysis, check out the exact commit recorded in that
analysis with `git checkout YOUR_RECORDED_COMMIT` before installing. The placeholder
is a commit you recorded, not a release identifier supplied by this manual.

On a server with wget but without Git:

```bash
wget -O PACE-main.tar.gz https://github.com/shenlinyong/PACE/archive/refs/heads/main.tar.gz
tar -xzf PACE-main.tar.gz
cd PACE-main
```

This is GitHub's source archive, not an installer or a model weight download.
For a fixed revision, GitHub also serves
`https://github.com/shenlinyong/PACE/archive/<FULL_COMMIT>.tar.gz`; substitute a
verified full commit. There is no claimed published PyPI package, Docker registry
image or livestock model download in this documentation.

## 2. Option A: Conda

Run from the downloaded repository root:

```bash
conda env create -f environment.yml
conda activate pace
PACE --version
PACE run --help
```

The environment installs PACE with genomic IO and ML dependencies. Conda is an
external prerequisite; obtain it through your institution's supported installation.
If you prefer an existing environment:

```bash
conda create -n pace python=3.12 pip -c conda-forge
conda activate pace
python -m pip install '.[io,ml]'
```

Real CNN training/inference additionally needs the sequence extra. For CPU work:

```bash
python -m pip install 'torch>=2.6,<3' --index-url https://download.pytorch.org/whl/cpu
python -m pip install '.[sequence]'
```

For GPU use, install a PyTorch build appropriate for the server's supported CUDA
configuration before the sequence extra. A GPU is not required for PACE's scoring
kernel or synthetic example models. GPU drivers are managed outside PACE.

## 3. Option B: Python virtual environment

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install '.[io,ml]'
PACE --version
```

A minimal `python -m pip install .` is sufficient for prepared TSV scoring and the
three synthetic demos. Use the table below to install only the required additions.
Always use the activated environment's `python -m pip`.

| Installation | Purpose |
|---|---|
| `python -m pip install .` | Core scoring, TSV/FASTA/VCF, synthetic examples |
| `python -m pip install '.[io]'` | bigWig, cool/mcool, BCF readers |
| `python -m pip install '.[sequence]'` | Quantitative sequence CNN training/inference |
| `python -m pip install '.[ml]'` | Optional scientific ML ecosystem |
| `python -m pip install -e '.[dev,io,ml,sequence]'` | Editable development, complete regression suite |

The independent prefix installer is another option:

```bash
./install.sh --python python3 --prefix "$HOME/.local" --extras io,ml
export PATH="$HOME/.local/bin:$PATH"
PACE --version
```

It creates its own environment below the prefix and refuses to replace unrelated
executables. It is unnecessary if the Conda or virtual-environment installation
already works. `setup.sh` delegates to the same installer.

## 4. Option C: Docker

Docker must already be installed and usable by your account. Build the image from
the source you downloaded; `pace:local` is a local image name:

```bash
docker build -t pace:local .
docker run --rm pace:local --version
```

The default image includes `io` and `ml`. To include CNN dependencies:

```bash
docker build --build-arg PACE_EXTRAS=io,ml,sequence -t pace:sequence .
```

The image entrypoint is `PACE`. Mount the project directory as `/work`; using the
host UID/GID makes created files belong to the calling Linux account:

```bash
docker run --rm --user "$(id -u):$(id -g)" \
  -v "$PWD":/work -w /work pace:local \
  measured --config examples/measured/config.yaml --out results/docker_measured

docker run --rm --user "$(id -u):$(id -g)" \
  -v "$PWD":/work -w /work pace:local \
  hybrid --config examples/hybrid/config.yaml --out results/docker_hybrid

docker run --rm --user "$(id -u):$(id -g)" \
  -v "$PWD":/work -w /work pace:local \
  genome --config examples/genome_only/config.yaml --out results/docker_genome
```

All YAML inputs must refer to paths visible **inside the container**. Relative
paths work unchanged when the whole project is mounted. A host path such as
`/home/user/data` is not visible unless separately mounted. The container does not
run read alignment or fetch a species model. Docker disk needs depend mainly on
whether PyTorch and its dependencies are included.

## 5. Check the installation

From the repository root in the activated environment:

```bash
PACE measured --config examples/measured/config.yaml --out results/check_measured
PACE hybrid --config examples/hybrid/config.yaml --out results/check_hybrid
PACE genome --config examples/genome_only/config.yaml --out results/check_genome
```

Use new output paths on each attempt: successful outputs are not overwritten.
Expect `scores.tsv.gz`, `gene_summary.tsv`, `qc_report.json`, `run_manifest.json`
and `resolved_config.yaml`. These fixtures test installation and workflow, not
biological accuracy. Continue with the [tutorial](TUTORIAL.md) or the
[complete Chinese manual](USER_GUIDE.zh-CN.md).

For developers:

```bash
python -m pytest -q
python scripts/check_public_docs.py
```

Dependencies and supported adapters are listed in [pyproject.toml](../pyproject.toml).
Raw FASTQ analysis, peak calling, variant calling, phasing and general genome
liftover are upstream tasks. No Snakemake raw-read workflow is distributed in the
current implementation. See [troubleshooting](TROUBLESHOOTING.md) for setup errors.
