# Download and installation

PACE requires Python >=3.11. Run these commands on Linux/macOS; WSL2 provides a Linux environment on Windows. A GPU is unnecessary.

## Download

```bash
git clone https://github.com/shenlinyong/PACE.git
cd PACE
git rev-parse HEAD
```

Record the commit for the analysis. Without git:

```bash
wget -O PACE-main.tar.gz https://github.com/shenlinyong/PACE/archive/refs/heads/main.tar.gz
tar -xzf PACE-main.tar.gz
cd PACE-main
```

## Conda

Install Conda/Miniforge first, then from the repository root:

```bash
conda env create -f environment.yml
conda activate pace
pace --version
pace demo --out results/install_check
```

The environment contains Python, the installed package and common experimental IO dependencies. Updating the checkout also requires reinstalling the package with the active environment's Python.

## Python virtual environment

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install '.[io,ml]'
pace --help
```

| Installation | Purpose |
|---|---|
| `python -m pip install .` | Core scoring, standard tables and the offline example |
| `python -m pip install '.[io]'` | bigWig and cool/mcool readers |
| `python -m pip install '.[ml]'` | Optional scientific ML dependencies |
| `python -m pip install -e '.[dev,io,ml]'` | Editable development and the complete supported test suite |

An isolated prefix installation is also available:

```bash
./install.sh --prefix "$HOME/.local" --python python3 --extras io,ml
export PATH="$HOME/.local/bin:$PATH"
pace --version
```

## Docker

The repository supplies a Dockerfile; build the image locally:

```bash
docker build -t pace:local .
docker run --rm pace:local --version
docker run --rm --user "$(id -u):$(id -g)" -v "$PWD":/work -w /work pace:local \
  run --config examples/measured/config.yaml --out results/docker_measured
```

The image entry point is already `pace`. The mounted directory must contain the input files and configuration. Paths inside the container use `/work`, not inaccessible host paths. Outputs are written into the host's mounted directory. Use a new output directory on each run.

## Verify and update

```bash
pace validate --config examples/measured/config.yaml
pace run --config examples/measured/config.yaml --out results/measured_check
```

After an authorized repository update, reinstall with `python -m pip install '.[io,ml]'` in the chosen environment, or rebuild the Docker image. Do not mix results from different versions without checking their measurement and comparison compatibility checks.

[Tutorial](TUTORIAL.md) · [Chinese manual](USER_GUIDE.zh-CN.md) · [Troubleshooting](TROUBLESHOOTING.md)
