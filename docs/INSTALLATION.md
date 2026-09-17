# Installation

PACE (Prediction of Activity-based regulatory Connections for Enhancers) runs locally or on a Linux compute server. The recommended setup uses Conda and the supplied environment files.

## Choose an environment

| Your starting point | Environment | Includes |
| --- | --- | --- |
| Quantified tables, called peaks or processed signal tracks | `pace`, from `environment.yml` | Scoring, genomic interval tools, direct commands, QC and tests |
| Aligned reads requiring automated peak calling | `pace-workflow`, from `workflow/envs/pace-env.yml` | Everything above plus MACS2 and Snakemake |
| Binary `.hic` or `.cool` contact files | Either environment plus the applicable optional reader | `hic-straw` or `cooler` |

Choose one of the first two environments for your analysis; the workflow environment already contains the direct-workflow dependencies. BEDPE contact tables need no extra reader.

## Before installation

Install **Bash, Git and Conda**. Conda supplies Python; a separate system Python installation is unnecessary. [Miniforge](https://github.com/conda-forge/miniforge#install) provides Conda with conda-forge configured. Follow its instructions for your operating system and CPU architecture, then reopen the shell or initialize Conda for that shell.

On a managed compute server, use the administrator-provided Conda installation if available. To install Conda yourself on Linux x86-64, the following follows the [Miniforge installation instructions](https://github.com/conda-forge/miniforge#install). It requires `curl` and uses an interactive installer, which lets you choose the installation directory:

```bash
curl -fL -o Miniforge3-Linux-x86_64.sh \
  https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
```

Accept shell initialization and open a new Bash terminal. Skip this step when `conda --version` already works. Git is needed to clone the repository; an extracted source ZIP can be used instead by entering its top-level directory. PACE's Conda environments do not require administrator privileges.

The documented validation target is **Linux x86-64**. On Windows, use a Linux environment such as WSL2. macOS installations depend on availability of the chosen bioinformatics packages and are not part of the recorded Linux validation. A GPU is not required. Genome-wide memory and runtime depend on the number of candidate–TSS pairs; pilot one chromosome before submitting a full-genome job. The small bundled examples do not establish a genome-wide resource requirement.

```bash
git --version
conda --version
```

## Recommended environment

```bash
git clone https://github.com/shenlinyong/PACE.git
cd PACE
CONDA_CHANNEL_PRIORITY=strict conda env create -f environment.yml
conda activate pace
```

The recipe uses **conda-forge before bioconda**, excludes default channels, and requests strict priority for the installation command. This follows [Bioconda's channel guidance](https://bioconda.github.io/#with-conda) without editing a user's global Conda configuration.

| Package | Requested version | Used for |
| --- | --- | --- |
| Python | 3.11 | Running the command-line programs |
| NumPy | ≥1.26, <3 | Numerical scoring |
| pandas | ≥2.2, <3 | Genomic and sample tables |
| PyYAML | ≥6, <7 | YAML configuration |
| SciPy | ≥1.11, <2 | Supporting numerical and analysis utilities |
| matplotlib-base | ≥3.8, <4 | QC plots without a desktop GUI |
| pyBigWig | ≥0.3.22, <0.4 | Exact regional summaries of bigWig tracks |
| bedtools | ≥2.31, <3 | Candidate interval operations and read counting |
| samtools | ≥1.19, <2 | BAM indexing and reference preparation |
| pytest | ≥8, <10 | Regression tests |
| pip | ≥24 | Optional Python packages |

All entries are installed by `conda env create`; users do not need to install them individually. The numerical-table interface itself only uses Python, NumPy and pandas. The recommended environment also supports direct genomic-file workflows and all shipped tests.

Verify the active environment:

```bash
python --version
python -c "import numpy, pandas, yaml, scipy, matplotlib, pyBigWig; print('Imports OK')"
bedtools --version
samtools --version
python scripts/pace.py --help
python -m pytest tests -q
python scripts/smoke_test.py --output-dir results/smoke
```

`Imports OK` confirms the Python modules are available; bedtools and samtools should each print a version. The regression suite should pass **48 tests**, and the small-file check should end with `PASS`. These checks exercise installed software; they do not benchmark a species or tissue.

For a batch job without interactive activation, use Conda's environment runner:

```bash
conda run -n pace python scripts/pace.py \
  --pairs example_quantified/candidates.tsv \
  --activity-config example_quantified/activity.json \
  --output results/installation_check/predictions.tsv
```

This should write six gene-level rows. If your environment has a different name or was created with `--prefix`, replace `-n pace` with the matching name or `-p /path/to/env`.

Continue with the [quick start](QUICKSTART.md).

## Optional peak-calling and Snakemake environment

Create a separate environment if you want Snakemake to run MACS2 and the downstream workflow on **already aligned reads**:

```bash
CONDA_CHANNEL_PRIORITY=strict conda env create -f workflow/envs/pace-env.yml
conda activate pace-workflow
snakemake --version
macs2 --version
```

This includes the packages above plus **MACS2 2.x**, **Snakemake 7.32.4**, **PuLP 2.7**, and **setuptools <81**. The pins preserve the workflow's command-line and package-resource API compatibility. MACS2 is unnecessary when candidate peaks have already been called; Snakemake is unnecessary for the direct commands.

Run the workflow inside this environment as shown in the [tutorial](TUTORIAL.md#optional-snakemake-workflow). FASTQ alignment and raw Hi-C processing are upstream tasks. No aligner, Hi-C matrix builder or RNA quantifier is installed implicitly.

## Optional contact readers

BEDPE inputs need no additional reader. To read binary contact matrices, install the applicable package in the environment you will run:

```bash
# For .hic files:
conda install --override-channels -c conda-forge -c bioconda \
  --strict-channel-priority hic-straw

# For .cool files:
conda install --override-channels -c conda-forge -c bioconda \
  --strict-channel-priority cooler
```

Check imports with `python -c 'import hicstraw'` or `python -c 'import cooler'`. Binary-reader installation does not supply the expected contact curve or local QC. Measured-contact scoring still requires the metadata described in [Inputs](INPUTS.md#contact-metadata). The BEDPE path is included in the small-file smoke check; this does not validate every real `.hic` or `.cool` dataset.

Optional preprocessing tools such as [deepTools](https://deeptools.readthedocs.io/en/develop/content/installation.html) can generate normalized bigWig tracks. They are not required for the bundled examples and their normalization must be chosen for the experiment.

## Installation helper and existing environments

`bash setup.sh` creates `pace` from `environment.yml`. `bash setup.sh --workflow` creates `pace-workflow`; `bash setup.sh --prefix /path/to/env` creates a named-location environment. The helper uses Conda and does not install into system Python.

For an environment you deliberately want to update:

```bash
CONDA_CHANNEL_PRIORITY=strict conda env update -n pace -f environment.yml
```

For an existing Python environment, the direct-workflow Python dependencies can be installed with `python -m pip install -r requirements.txt`. External tools still need separate installation. `requirements-core.txt` is for prepared-table scoring only; `requirements-test.txt` adds test dependencies. Conda is the supported complete installation route described here.

## Save the exact software environment

After activating the environment used for the analysis:

```bash
mkdir -p results/provenance
git rev-parse HEAD > results/provenance/pace_commit.txt
conda env export --no-builds > results/provenance/environment.yml
conda list --explicit > results/provenance/conda-explicit.txt
```

The repository recipe specifies compatible ranges. The explicit export records the resolved packages for reproducing the same platform. Preserve the configuration, sample sheet, reference assembly/accessions and unfiltered predictions with these files. See [Troubleshooting](TROUBLESHOOTING.md) if installation or execution fails.

## Recreate the tested Linux builds

The checked Linux x86-64 package builds and MD5 hashes are also provided as explicit Conda lock files. Choose this route when you need those exact builds:

```bash
conda create -n pace --file environment-linux-64.lock.txt
# Optional workflow environment:
conda create -n pace-workflow --file workflow/envs/pace-linux-64.lock.txt
```

Use a new environment name if it already exists. These files are platform-specific; the YAML recipes remain the installation entry point for a fresh dependency resolution.
