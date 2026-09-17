# PACE

[![Canonical software tests](https://github.com/shenlinyong/PACE/actions/workflows/ci.yml/badge.svg)](https://github.com/shenlinyong/PACE/actions/workflows/ci.yml)

**Auditable activity–contact support for livestock enhancer–gene research.**

The installable Python package implements the September 2026 canonical-grid contract with
three evidence regimes: **measured**, **hybrid** and **genome_only**. It keeps quantitative
signals, evidence provenance, missingness and normalization backgrounds explicit. Author and
maintainer: **申林用 (Linyong Shen), Northwest A&F University**.

[New software manual](docs/software.md) · [中文说明](README.zh-CN.md) ·
[Input preparation](docs/input_preparation.md) · [Training](docs/training.md) ·
[Comparisons](docs/comparison.md) · [Validation](docs/validation.md)

## Installable canonical interface

Python 3.11+; base execution needs NumPy and PyYAML, with no GPU or runtime download.

```bash
git clone https://github.com/shenlinyong/PACE.git
cd PACE
python -m pip install .
pace-livestock demo --regime measured --out results/demo_measured
pace-livestock demo --regime hybrid --out results/demo_hybrid
pace-livestock demo --regime genome_only --out results/demo_genome
```

For real genomic binary files install `'.[io]'`; for quantitative CNN training/inference install
`'.[sequence]'`. [Configuration and command reference](docs/software.md#commands).

The score is `A × Cbar × B^eta`, normalized within the actually scoreable candidates of each
gene; default **eta=0**. The panel stays fixed, true zeros stay zero, required missing contacts
stay NA, and no arbitrary residual mass is added. `compare` recomputes a common denominator
from raw support and distinguishes complete from conditional changes. RNA and auxiliary
marks remain annotations unless explicitly used by an independently fitted classifier.

All bundled examples and weights are **synthetic software fixtures**. Real species/tissue
weights and independent biological validation are **not supplied**. Software test success does
not establish individual-effect accuracy. The score is a relative support share, not a causal
probability or expression fold change. See [scope and limitations](docs/limitations.md).

```bash
pace-livestock validate --config examples/measured/config.yaml
pace-livestock run --config examples/measured/config.yaml --out results/measured
pace-livestock capabilities --config examples/genome_only/config.yaml
```

Runs export all candidate scores, per-gene denominators, resolved evidence, nine-omics features,
QC, input/model hashes and actual configuration. Existing output directories are never overwritten.
The model, [data dictionary](docs/data_dictionary.md), and [parameter reference](docs/parameters.md)
describe the scientific contract. [CITATION.cff](CITATION.cff) records software authorship; cite the
exact commit used. No publication DOI or formal software release is implied.

## Existing region-based interface

**Everything below documents the preserved `scripts/pace.py` and `workflow/` interface.** Its
historical defaults (including eta=1 and missing-assay handling) differ from the canonical package
above. Its commands and existing analyses remain available for reproduction. Do not apply its
defaults or thresholds to canonical scores, or compare the two without a new compatible analysis.

### PACE (Prediction of Activity-based regulatory Connections for Enhancers)

**Predict enhancer–gene links from chromatin activity, promoter contact and gene annotation.**

[Installation](docs/INSTALLATION.md) · [Quick start](docs/QUICKSTART.md) · [Tutorial](docs/TUTORIAL.md) · [ABC comparison](docs/ABC_COMPARISON.md) · [Parameters](docs/PARAMETERS.md) · [中文手册](docs/README_zh.md)

PACE is a general framework for ranking candidate enhancer–gene links across species using matching genome references and annotations. It combines available activity measurements, adjusts contact evidence by its reliability, and integrates distinct transcription start sites (TSSs). Each prediction includes a relative score and a separate assessment of input evidence.

These features make PACE particularly well suited to livestock datasets, where assay coverage is often uneven, tissue-matched Hi-C is limited, and promoter annotation is incomplete. The same features apply to other species with comparable data constraints. Pig, cattle and chicken are application examples; the scoring model contains no livestock-specific species restriction.

## PACE and the original ABC model

The comparison below refers to [Fulco et al. (2019)](https://doi.org/10.1038/s41588-019-0538-0) and the [NG2019 implementation](https://github.com/EngreitzLab/ABC-Enhancer-Gene-Prediction-20250314-archive/tree/NG2019).

ABC provides an activity–contact baseline when those inputs can be estimated reliably. PACE extends that formulation to make uneven coverage, variable contact reliability and incomplete annotation explicit in scoring and reporting.

| Component | Original ABC formulation | PACE modification | Why it helps with livestock data |
| --- | --- | --- | --- |
| **Contact** | Activity is multiplied by the chosen contact estimate; averaged Hi-C and distance alternatives were already considered | Adjust observed/expected contact toward a distance prior according to local reliability | Low-quality or unavailable matched Hi-C has less influence; source and fallback remain visible |
| **Activity** | The original two-assay profile uses the geometric mean of accessibility and H3K27ac | Combine effective measurements with declared scales and weights; distinguish `NA` from measured `0` | A missing assay does not erase an available measurement or become evidence of inactivity |
| **Candidate targets** | Normalize activity × contact independently for each gene | First allocate support across an enhancer's candidate genes through $B(E,G)$ | Reduces raw support for broadly shared enhancers, aiming to improve specificity in gene-dense regions |
| **Promoters** | Contact is defined relative to the selected gene promoter | Deduplicate transcript TSSs and average contacts with gene-level TSS weights | Integrates alternative promoters without inflating support through repeated annotation records |
| **Evidence** | The relative ABC score summarizes activity–contact support | Report score separately from $Q$, component qualities and reason codes | A candidate can remain available while incomplete evidence is labelled `provisional` |
| **Unscored candidates** | The normalization is defined over candidate activity–contact products | Use finite support in $\mathcal E^{\mathrm{obs}}(G)$, retain unscored rows and allow independent residual support $U(G)$ | An uncomputable contribution remains missing; omitted regulatory support is explicitly acknowledged |

For the two models, the normalization rules are

$$
\mathrm{ABC}(E,G)=\frac{A_{\mathrm{ABC}}(E)C_{\mathrm{ABC}}(E,G)}{\sum_{e\in\mathcal E(G)}A_{\mathrm{ABC}}(e)C_{\mathrm{ABC}}(e,G)},
$$

$$
\mathit{PACE}(E,G)=\frac{A(E)C(E,G)B(E,G)^\eta}{\sum_{e\in\mathcal E^{\mathrm{obs}}(G)}A(e)C(e,G)B(e,G)^\eta+U(G)}.
$$

In PACE, $C(E,G)$ is already reliability-adjusted and averaged across distinct TSSs; it is the adjusted contact $\overline C(E,G)$ in a direct ABC comparison. The default allocation exponent is $\eta=1$. Unknown $U(G)$ is computationally zero and remains flagged; PACE does not estimate missing regulatory support automatically.

Target allocation and multi-TSS integration build on [generalized ABC / STARE](https://doi.org/10.1093/bioinformatics/btad062). The [comparison guide](docs/ABC_COMPARISON.md) explains all six changes, numerical examples and attribution. Suitability for incomplete data is a design property; improved accuracy requires a matched biological benchmark.

## Installation

Use **Linux, Bash, Git and Conda**. [Miniforge](https://github.com/conda-forge/miniforge#install) supplies Conda. The environment file installs Python and the required packages together.

```bash
git clone https://github.com/shenlinyong/PACE.git
cd PACE
CONDA_CHANNEL_PRIORITY=strict conda env create -f environment.yml
conda activate pace
python scripts/pace.py --help
```

The `pace` environment includes Python 3.11, NumPy, pandas, PyYAML, SciPy, matplotlib, pyBigWig, bedtools, samtools, pytest and pip. It supports numerical tables, called peaks and processed genomic signals. A GPU and training labels are unnecessary.

| Additional task | Additional packages | Installation |
| --- | --- | --- |
| Automated peak calling from aligned reads | MACS2, Snakemake, compatible PuLP and setuptools | [Separate `pace-workflow` Conda environment](docs/INSTALLATION.md#optional-peak-calling-and-snakemake-environment) |
| Read `.hic` contact files | `hic-straw` | [Optional contact readers](docs/INSTALLATION.md#optional-contact-readers) |
| Read `.cool` contact files | `cooler` | [Optional contact readers](docs/INSTALLATION.md#optional-contact-readers) |

[Package versions and roles, Linux lock files, installation checks and environment export](docs/INSTALLATION.md)

## Quick start

Run from the repository root with the `pace` environment active:

```bash
python scripts/pace.py \
  --pairs example_quantified/candidates.tsv \
  --activity-config example_quantified/activity.json \
  --output results/quickstart/predictions.tsv
```

Expected console output:

```text
Wrote 6 gene-level edges; all scores are uncalibrated.
```

The nine enhancer–TSS input rows become six enhancer–gene rows. Two retain missing scores because their activity was not measured. Inspect the results:

```bash
python - <<'PYCODE'
import pandas as pd
p = pd.read_csv('results/quickstart/predictions.tsv', sep='\t')
print(p[['TargetGeneEnsemblID', 'start', 'PACE.Score',
         'contact_state', 'evidence_status']].to_string(index=False))
PYCODE
```

The [quick-start guide](docs/QUICKSTART.md) includes expected scores and filtering commands. The [worked examples](docs/WORKED_EXAMPLES.md) show what changes when an assay or contact measurement is missing.

To run the bundled genomic-read example:

```bash
bash example/run_example_direct.sh
```

This constructs candidates, quantifies signals, scores pairs and writes filtered predictions and QC to `example/results/Example_Sample/`. The unfiltered table contains **12,000 enhancer–gene pairs**. Both examples use synthetic inputs to demonstrate software behavior.

## Run your own data

| Starting data | Required inputs | Use |
| --- | --- | --- |
| Quantified measurements | Enhancer coordinates, stable gene IDs, TSSs, positive contact priors, activity or named assay signals | [`scripts/pace.py`](docs/INPUTS.md#quantified-candidate-table) |
| Called peaks and processed signals | BED/narrowPeak, chromosome sizes, matched GTF/TSS catalogue, accessibility signal | [Step-by-step tutorial](docs/TUTORIAL.md#step-by-step-genomic-file-example) |
| Aligned accessibility reads | BAM/tagAlign, reference files, sample sheet and YAML configuration | [Snakemake tutorial](docs/TUTORIAL.md#optional-snakemake-workflow) |

Use one assembly throughout and BED0 coordinates for genomic tables. `prepare_tss.py` converts GTF coordinates. Supply `NA` for an unavailable measurement and `0` for measured zero. Measured Hi-C is optional; qualified contact use requires compatible expected contacts, reliability and source metadata. FASTQ alignment and raw Hi-C processing are upstream steps.

For a new species, replace the reference files and review the candidate catalogue, signal scales and contact prior. [The tutorial](docs/TUTORIAL.md#replace-the-example-with-your-species-and-tissue) explains these choices, including livestock-specific preparation considerations.

## Interpret the output

Keep the **unfiltered prediction table** as the analysis record. Apply selection thresholds after scoring so that filtering does not redefine the normalization background.

| Field | Interpretation |
| --- | --- |
| `PACE.Score` | Relative support within the supplied gene background; `NA` means unscorable |
| `contact_gene`, `target_share`, `raw_support` | Adjusted gene contact, enhancer-side allocation and unnormalized support |
| `TargetGeneTSSs`, `n_tss` | Distinct TSSs included in the gene-level result |
| `contact_state` | Prior-only, matched, surrogate or contact evidence requiring further QC |
| `evidence_status`, `evidence_reasons` | Whether input evidence is sufficient, provisional or insufficient, and why |
| `score_scope`, `unscored_candidates` | Residual-support status and supplied candidates lacking finite support |

A score is neither a regulatory probability nor a calibrated false discovery rate. The example cutoff **0.02** needs independent calibration for the intended species and tissue. `provisional` candidates may still be useful for follow-up; `sufficient_input_evidence` describes the inputs, not functional validation. [Full output reference](docs/IO_FORMATS.md)

## Documentation and verification

| Guide | Contents |
| --- | --- |
| [ABC comparison](docs/ABC_COMPARISON.md) | Six modifications, their rationale and relevance to livestock datasets |
| [Parameters](docs/PARAMETERS.md) | Every setting, default, actual control point and expected effect |
| [Installation](docs/INSTALLATION.md) | Conda setup, package roles, optional dependencies and exact Linux builds |
| [Quick start](docs/QUICKSTART.md) / [Tutorial](docs/TUTORIAL.md) | Runnable examples, expected outputs and real-data preparation |
| [Worked examples](docs/WORKED_EXAMPLES.md) | Missing assays, missing versus zero contact, and allocation ablation |
| [Inputs](docs/INPUTS.md) / [Outputs](docs/IO_FORMATS.md) | File formats, missing values, scores and evidence states |
| [Equations](docs/FORMULA.md) / [Notation](docs/NOTATION.md) | Mathematical definitions and their software fields |
| [Troubleshooting](docs/TROUBLESHOOTING.md) | Installation errors, empty outputs and unexpected scores |
| [中文手册](docs/README_zh.md) | 中文模型比较、参数说明、安装与使用教程 |

```bash
python -m pytest tests -q
python scripts/smoke_test.py --output-dir results/smoke
```

These checks cover the scoring kernel and small-file execution. [Validation](VALIDATION.md) distinguishes software verification from biological evidence.

## Citation and support

Author and maintainer: **shenlinyong, 申林用 (Linyong Shen), Northwest A&F University**.

Cite the PACE repository and the exact commit used; a publication DOI is not assigned in this repository. Cite the underlying [ABC](https://doi.org/10.1038/s41588-019-0538-0) and [generalized ABC](https://doi.org/10.1093/bioinformatics/btad062) methods when discussing those components.

Report problems through [GitHub Issues](https://github.com/shenlinyong/PACE/issues), including the command, environment and a small reproducible input. PACE uses the [MIT licence](LICENSE); third-party attribution appears in [AUTHORS.md](AUTHORS.md).
