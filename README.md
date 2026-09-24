# PACE

**Enhancer–gene prediction for livestock and poultry.**

PACE ranks candidate regulatory elements for each gene in a tissue using measured activity and contact evidence. It is built on the Activity-by-Contact (ABC) model and adapted to the data that livestock labs actually have: incomplete epigenomic panels, sparse or low-resolution Hi-C, rough transcript annotation, and only one or a few animals per tissue.

The input formats are species-independent and require a matching reference, annotation and experimental activity data. Pig, cattle, sheep, goat and chicken are the main use cases.

```bash
git clone https://github.com/shenlinyong/PACE.git && cd PACE
conda env create -f environment.yml && conda activate pace
pace demo --out results/demo
```

---

## Contents

1. [How it works](#how-it-works)
2. [Is PACE right for your data?](#is-pace-right-for-your-data)
3. [Installation](#installation)
4. [Test run](#test-run)
5. [Running on your own data](#running-on-your-own-data)
   - [Step 0. Collect your files](#step-0-collect-your-files)
   - [Step 1. Build the candidate catalog](#step-1-build-the-candidate-catalog)
   - [Step 2. Extract activity signal](#step-2-extract-activity-signal)
   - [Step 3. Extract Hi-C contacts](#step-3-extract-hi-c-contacts)
   - [Step 4. Write the sample and source tables](#step-4-write-the-sample-and-source-tables)
   - [Step 5. Score](#step-5-score)
6. [No Hi-C? Low-resolution Hi-C?](#no-hi-c-low-resolution-hi-c)
7. [Reading the results](#reading-the-results)
8. [Comparing animals, tissues or treatments](#comparing-animals-tissues-or-treatments)
9. [RNA-seq, methylation and other data](#rna-seq-methylation-and-other-data)
10. [Configuration reference](#configuration-reference)
11. [Input table reference](#input-table-reference)
12. [The formula in full](#the-formula-in-full)
13. [How PACE differs from ABC](#how-pace-differs-from-abc)
14. [Troubleshooting](#troubleshooting)
15. [Citation](#citation)

---

## How it works

For every candidate element–gene pair within 5 Mb, PACE computes:

```math
\boxed{
\mathrm{PACE}(E,G)=
\frac{A_\star(E)\,\overline C(E,G)}
{\displaystyle\sum_{e\in\mathcal E(G)} A_\star(e)\,\overline C(e,G)}
}
```

- **Activity** is how open and active the element is, taken from ATAC-seq, DNase-seq and/or H3K27ac ChIP-seq. When you have two assays, PACE uses their geometric mean.
- **Contact** is how often the element touches the gene's promoter in 3D, taken from Hi-C. If a gene has several distinct transcription start sites (TSSs), PACE measures contact to each one and combines them.

For a complete candidate set with positive total support, the scores for one gene add up to 1. A score of 0.3 means 30% of the modeled activity-weighted contact in that set, not 30% of gene expression. Incomplete genes have NA primary scores by default; available-subset scores are separate.

Two consequences worth knowing:

- Multiplying a whole assay by a constant (for example, a different sequencing depth) cancels out within a gene. So ATAC and H3K27ac do not need to be on the same scale.
- A score is a share, not an absolute strength. An element's score can go down simply because a neighbouring element got stronger.

## Is PACE right for your data?

**Good fit**

- You have ATAC-seq, DNase-seq or H3K27ac ChIP-seq (or CUT&Tag) from a livestock tissue.
- You want to connect GWAS hits, eQTLs or selection signals to target genes.
- You have one or a few animals, not dozens.
- Your Hi-C is 10–25 kb resolution, shallow, from another tissue, or missing entirely.

**Not a good fit**

- You only have FASTQ files. PACE starts from processed files (peaks, bigWig, cool/mcool). Do alignment, peak calling and Hi-C processing first.
- You have no activity data at all. PACE does not predict enhancer activity from DNA sequence.
- You need a causal probability. PACE ranks candidates. Functional validation is still needed to prove regulation.

## Installation

Requirements: Linux or macOS, Python 3.11 or newer. No GPU.

**Conda (recommended)**

```bash
git clone https://github.com/shenlinyong/PACE.git
cd PACE
conda env create -f environment.yml
conda activate pace
pace --version
```

**pip**

```bash
git clone https://github.com/shenlinyong/PACE.git
cd PACE
python3 -m venv .venv
source .venv/bin/activate
pip install '.[io,ml]'
```

The extras are optional. `io` adds pyBigWig, cooler and pandas, which you need to read bigWig and Hi-C files. `ml` adds scikit-learn and scipy for the optional classifier. `pip install .` alone is enough for the demo.

**Docker**

```bash
docker build -t pace:local .
docker run --rm --user "$(id -u):$(id -g)" -v "$PWD":/work -w /work pace:local \
  demo --out results/docker_demo
```

Any path you pass to the container must be inside the mounted directory.

Use `pace --help` to list commands and `pace <command> --help` for their full options.
`pace-livestock` is an alias of `pace`.

**No git?**

```bash
wget -O PACE-main.tar.gz https://github.com/shenlinyong/PACE/archive/refs/heads/main.tar.gz
tar -xzf PACE-main.tar.gz && cd PACE-main
```

For published analyses, record the exact version you used:

```bash
git rev-parse HEAD
```

## Test run

```bash
pace demo --out results/demo
```

If `results/demo/results/scores.tsv.gz` appears, the installation works. The demo data are synthetic. They test the software, not biology.

To see what a full configuration looks like and run it:

```bash
pace validate --config examples/measured/config.yaml
pace run      --config examples/measured/config.yaml --out results/example
```

`validate` checks inputs and computes a scoring preview without writing a result directory. By default `--out` must be new. `--force` replaces a recognized PACE result only after retaining a sibling backup.

## Running on your own data

The whole pipeline looks like this:

| Step | Input | Output |
|---|---|---|
| 1. Catalog | Peaks, GTF and chromosome sizes | Units, promoters and candidates |
| 2. Activity | ATAC/DNase/H3K27ac bigWigs | observed_activity.tsv |
| 3. Contact | Hi-C matrix and query pairs | observed_contacts.tsv |
| 4. Metadata | Library identities and data sources | samples.tsv and sources.tsv |
| 5. Scoring | Prepared inputs and settings | scores.tsv.gz and QC reports |

A tidy project layout helps. The examples below assume:

| Path within my_project | Contents |
|---|---|
| data/ | Processed BED, GTF, bigWig and cool/mcool files |
| prepared/ | Input tables |
| results/ | PACE outputs |
| run.yaml | Analysis settings |

All commands are run from `my_project/`.

### Step 0. Collect your files

Every file must use the **same genome assembly and the same chromosome names**. For example, mixing `chr1` and `1` prevents matching coordinates.

| File                                          | Required?                   | How to get it                                                |
| --------------------------------------------- | --------------------------- | ------------------------------------------------------------ |
| Peaks or candidate regions (BED)              | Yes                         | MACS2/MACS3 on ATAC or DNase; merge peaks from all samples of the tissue |
| Gene annotation (GTF) with `transcript` lines | Yes                         | Ensembl or NCBI for your assembly                            |
| Chromosome sizes                              | Yes                         | `cut -f1,2 genome.fa.fai > data/chrom_sizes.tsv`, then add a header line `chrom	length` |
| ATAC-seq or DNase-seq bigWig                  | At least one activity assay | deepTools `bamCoverage --normalizeUsing CPM`                 |
| H3K27ac bigWig                                | At least one activity assay | same                                                         |
| Hi-C `.mcool`                                 | Recommended                 | HiC-Pro, distiller or Juicer (then `hic2cool`)               |

Supported activity panels: ATAC alone, DNase alone, H3K27ac alone, ATAC + H3K27ac, or DNase + H3K27ac. Evaluate the chosen panel on independent data; adding a low-quality assay does not guarantee better rankings.

### Step 1. Build the candidate catalog

Save as `prepare_catalog.yaml`:

```yaml
kind: catalog
bed: data/peaks.bed               # candidate regions
gtf: data/annotation.gtf          # gene annotation
chrom_sizes: data/chrom_sizes.tsv
source_id: pig_liver_peaks        # any short name; used to trace results back to this file
width: 500                        # non-overlapping grid width in bp; not ABC peak-centered windows
offset: 0
radius: 5000000                   # search for candidates up to 5 Mb from each TSS
include_promoters: true           # retain promoter units in the normalization background
```

Run:

```bash
pace prepare --config prepare_catalog.yaml --out prepared/catalog
```

You get:

- `units.tsv`: the candidate elements
- `promoters.tsv`: one row per distinct physical TSS (transcripts sharing a TSS are counted once)
- `candidates.tsv`: which elements are candidates for which genes
- `region_membership.tsv`: which input region each element came from

The grid includes only cells overlapping input regions or required promoters, not the whole genome. Broad or noisy peaks still enlarge the denominator.

`radius` is a search window, not a biological claim. 5 Mb matches ABC. If you change it, check that your top predictions are stable.

### Step 2. Extract activity signal

Run once per bigWig. Save as `prepare_atac.yaml`:

```yaml
kind: bigwig
track: data/pig1_ATAC.bw
units: prepared/catalog/units.tsv
sample_id: pig1_ATAC              # must match samples.tsv in Step 4
assay: ATAC                       # ATAC, DNase or H3K27ac
unit: normalized_signal
normalization_id: cpm             # a label for how the bigWig was normalized; keep it identical across samples
window_id: grid:500:mean
missing_is_measured_zero: false
minimum_callable_fraction: 0.8    # skip elements where less than 80% of bases have data
```

```bash
pace prepare --config prepare_atac.yaml --out prepared/pig1_atac
```

Copy the file, change `track`, `sample_id`, `assay` and the output folder, and run it again for H3K27ac.

About `missing_is_measured_zero`: bigWig files often leave out regions with no reads. That can mean "measured, zero signal" or "not measured". If you are not sure how your bigWig was made, leave it `false`.

Now merge all activity tables into one:

```bash
pace merge-tables --table observed_activity \
  --inputs prepared/pig1_atac/observed_activity.tsv \
           prepared/pig1_h3k27ac/observed_activity.tsv \
  --out prepared/activity
```

### Step 3. Extract Hi-C contacts

First list every element–promoter pair that needs a contact value:

```bash
pace prepare-pairs --catalog-dir prepared/catalog --out prepared/pairs
```

Then save as `prepare_hic.yaml`:

```yaml
kind: cooler
contact: data/pig1.mcool::/resolutions/5000   # pick the finest resolution your data support
pairs: prepared/pairs/pairs.tsv
resolution: 5000
balanced: true                    # true if the file has ICE/KR balancing weights
missing_pixels_are_zero: true
sample_id: pig1_HiC
source_id: pig1_HiC_source
scale: balanced_contact           # a label; must match contact.scale in run.yaml
normalization_id: ice             # a label; must match contact.normalization_id in run.yaml
```

```bash
pace prepare --config prepare_hic.yaml --out prepared/pig1_hic
```

Things that will give wrong results:

- Observed/expected (O/E) matrices. They have the distance effect removed, which is exactly what PACE needs.
- p-values or loop calls instead of contact frequencies.
- Averaging contacts from different resolutions or normalizations.

For `.hic` input, first convert to cool/mcool using a compatible converter and verify available resolutions and balancing weights; a filename extension does not determine the stored resolution.

### Step 4. Write the sample and source tables

`prepared/samples.tsv` has one row per sequencing library. Columns are tab-separated:

```
sample_id	donor_id	assay	biological_replicate	technical_replicate	species	assembly	context_id	source_id
pig1_ATAC	pig1	ATAC	1	1	pig	Sscrofa11.1	liver	pig1_ATAC_source
pig1_K27ac	pig1	H3K27ac	1	1	pig	Sscrofa11.1	liver	pig1_K27ac_source
pig1_HiC	pig1	HiC	1	1	pig	Sscrofa11.1	liver	pig1_HiC_source
```

- `donor_id` is the animal. Technical replicates of the same library share a `biological_replicate` number.
- `species`, `assembly` and `context_id` must be spelled exactly as in `run.yaml`.

`prepared/sources.tsv` records where each file came from, so the analysis can be traced later:

```
source_id	path_or_accession	source_type	assembly	processing_method	normalization_id	checksum
pig_liver_peaks	data/peaks.bed	bed	Sscrofa11.1	MACS3 callpeak -q 0.01	NA	NA
pig1_ATAC_source	data/pig1_ATAC.bw	bigwig	Sscrofa11.1	bowtie2 + bamCoverage	cpm	NA
pig1_K27ac_source	data/pig1_K27ac.bw	bigwig	Sscrofa11.1	bowtie2 + bamCoverage	cpm	NA
pig1_HiC_source	data/pig1.mcool	mcool	Sscrofa11.1	HiC-Pro + cooler balance	ice	NA
```

If you have public data, put the accession (for example a GEO or SRA ID) in `path_or_accession`. `pace validate` will tell you if any value is not accepted.

### Step 5. Score

Save as `run.yaml`. Lines marked `# CHANGE` are the ones you must edit. Check measurement settings against your inputs. Omitted schema, input mode and scoring-target fields use built-in defaults.

```yaml
run_id: pig1_liver                        # CHANGE: a name for this analysis
execution_profile: research
target_level: individual                  # one animal: individual. Several animals averaged: population_mean

context:
  species: pig                            # CHANGE
  assembly: Sscrofa11.1                   # CHANGE
  context_id: liver                       # CHANGE: tissue or cell type

inputs:
  units: prepared/catalog/units.tsv
  region_membership: prepared/catalog/region_membership.tsv
  promoters: prepared/catalog/promoters.tsv
  candidates: prepared/catalog/candidates.tsv
  samples: prepared/samples.tsv
  sources: prepared/sources.tsv
  observed_activity: prepared/activity/observed_activity.tsv
  observed_contacts: prepared/pig1_hic/observed_contacts.tsv   # CHANGE, or remove if you have no Hi-C

catalog:
  profile: canonical_grid
  width_bp: 500
  offset_bp: 0
  include_promoter_units: true
  chrom_sizes_path: prepared/catalog/chrom_sizes.tsv
  candidate_radius_bp: 5000000

activity:
  panel: [ATAC, H3K27ac]                  # CHANGE: e.g. [ATAC], [H3K27ac] or [DNase, H3K27ac]
  minimum_callable_fraction: 0.8

contact:
  mode: observed                          # observed, prior_only or shrinkage (see next section)
  scale: balanced_contact                 # CHANGE: same label as Step 3
  resolution: 5000                        # CHANGE: same as Step 3
  normalization_id: ice                   # CHANGE: same label as Step 3
  balancing: balanced
  window_id: bin_pair
  near_diagonal_policy: prior_or_neighbor
  pseudocount: auto

scoring:
  partial_policy: withhold

promoters:
  weights: equal                          # use "provided" only if promoters.tsv has trusted pi values

allocation:
  eta: auto

multiomics:
  mode: annotate

seed: 17
```

Paths inside a YAML file are relative to **the YAML file's folder**, not to where you type the command. If you put `run.yaml` in a subfolder, add `../` to the paths.

Check, then run:

```bash
pace validate --config run.yaml
pace run      --config run.yaml --out results/pig1_liver
```

`pace measured` does exactly the same thing as `pace run`.

**Several animals?** Score each animal separately first (`target_level: individual`) and check that the top predictions agree between animals. Then, if you want a tissue-level map, run once more with all animals in `samples.tsv` and `target_level: population_mean`. PACE averages technical replicates within each biological replicate, then biological replicates within each animal, then gives every animal equal weight. Technical replicates never count as extra animals.

## No Hi-C? Low-resolution Hi-C?

Hi-C is the data type livestock projects most often lack. PACE has three contact modes:

| `contact.mode` | Uses                                                     | When to use                               |
| -------------- | -------------------------------------------------------- | ----------------------------------------- |
| `observed` | Measured Hi-C, optional matching-prior regularization | A contact matrix with suitable coverage and resolution |
| `shrinkage` | Observed contact and a compatible distance prior | An independently chosen observation weight and matching measurement scales |
| `prior_only`   | Distance-decay prior only                                | No usable Hi-C                            |

The distance-decay prior has the form

```
contact_prior(d) = a × ( max(d, d_min) / d_ref ) ^ (−γ)
```

where `d` is the element–TSS distance. Point `contact.prior_path` at a prior fitted to your species. Use `pace fit-prior --cooler ... --species ... --assembly ... --tissue ... --out ...` to fit a prior including callable zero pixels. The lower distance defaults to the matrix resolution. A different tissue requires `contact.allow_cross_context_prior: true`; source and target contexts are retained and the transfer is labeled unvalidated. The species, assembly, target level and measurement definitions must match. Tissue transfer requires evaluation; it is not assumed equivalent to target-tissue Hi-C.

Without any usable Hi-C, the explicit research baseline `prior_preset: abc_human` with `mode: prior_only` uses the human ABC reference gamma=1.024238616787792 on a relative scale. Omit contact tables for this baseline. It is not a fitted livestock parameter. See [the model](docs/ADVANCED.md#fitting-and-transferring-a-prior).

For `shrinkage`, set `contact.reliability` to the weight given to observed contacts (1 = observed only, 0 = prior only) and `contact.reliability_source` to a short note saying how you chose it.

Results from `prior_only` and `shrinkage` runs are labelled as such in the output. Report the mode you used in your methods.

Two settings matter most with low-resolution Hi-C:

- `near_diagonal_policy: prior_or_neighbor` (default). Elements in the same Hi-C bin as the TSS would otherwise get the diagonal value, which mostly reflects bin self-contact. PACE uses the prior if one is given, or else the recorded maximum of valid neighboring contacts. If neither is available, contact remains NA.
- `pseudocount: auto` (default). Adds a distance-based term to finite off-diagonal measurements only when a compatible fitted prior is available. Without one, observations are unchanged. Missing values are not converted to measured zeros.

## Reading the results

The output folder contains:

| File                         | What's in it                                                 |
| ---------------------------- | ------------------------------------------------------------ |
| `scores.tsv.gz`              | **Main result.** One row per element–gene pair               |
| `gene_summary.tsv`           | One row per gene: number of candidates, how many could be scored, coverage, denominator |
| `qc_report.json`             | Overall quality summary. Read this first                     |
| `report.md`                  | Short human-readable summary of the run                      |
| `resolved_activity.tsv` | Measured values and optional assay offsets, kept separate |
| `promoter_weights.tsv` | Original/effective TSS weights, retained weight and filtering reasons |
| `resolved_contacts.tsv`      | The exact contact value used for every element–TSS pair, and where it came from |
| `resolved_config.yaml`       | Your configuration with every default filled in              |
| `run_manifest.json`          | Input file hashes, software version and environment. Useful for the methods section |
| `eta_calibration.json`       | Whether the optional allocation term was used (normally not) |
| `multiomics_features.tsv.gz` | Extra annotations, if you supplied RNA-seq, methylation etc. |

Key columns in `scores.tsv.gz`:

| Column                   | Meaning                                                      |
| ------------------------ | ------------------------------------------------------------ |
| `element_id`, `gene_id`  | The pair                                                     |
| `pace_score`             | The relative score, 0–1; withheld for incomplete genes by default                 |
| `A_used`                 | Activity of the element                                      |
| `Cbar`                   | Contact between the element and the gene's promoters         |
| `support`, `log_support` | Activity × contact before dividing by the gene total         |
| `scoreable`              | Whether the pair has calculable support; the full gene denominator may still be incomplete                            |
| `normalization_status`   | `complete`, `partial`, `zero_support` or `empty` (see below) |
| `reason`                 | Why a value is NA                                            |

Look at the top of the file:

```bash
zcat results/pig1_liver/scores.tsv.gz | head -5 | column -t
```

Or in Python:

```python
import pandas as pd
s = pd.read_csv("results/pig1_liver/scores.tsv.gz", sep="\t")
top = s[s.normalization_status == "complete"].dropna(subset=["pace_score"]).sort_values("pace_score", ascending=False)
print(top[["element_id", "gene_id", "pace_score", "A_used", "Cbar"]].head(20))
```

**How to interpret the results**

- **Check coverage before scores.** `normalization_status` tells you whether every planned candidate for the gene could be scored (`complete`) or only some (`partial`). With `partial_policy: withhold` (default), partial genes get no score, because dividing by an incomplete total would inflate the remaining candidates.
- **NA is not zero.** NA means the value could not be computed. Zero means the assay measured zero signal; it does not prove the absence of biological activity or contact. Never replace NA with 0.
- **A score of 1 is not automatically a strong enhancer.** If only one candidate has any support, it gets 1 by definition.
- **Choose your own threshold.** Published ABC cutoffs depend on the assay panel, candidate definition and human perturbation benchmark. Calibrate a livestock threshold on independent functional labels where available; eQTL colocalisation provides complementary association evidence, not interchangeable causal labels. Filter after scoring, never before, so the gene totals stay the same.
- **A score is not a probability.** It ranks candidates for follow-up.

## Comparing animals, tissues or treatments

Do not subtract two score tables directly. Because scores are shares, a change in one element shifts every other element of the same gene. Use `pace compare`, which rescores both runs over the candidates they have in common.

Save as `compare.yaml`:

```yaml
left: results/pig1_liver
right: results/pig2_liver
minimum_common_units: 2
allow_eta_difference: false
allow_evidence_difference: false
```

```bash
pace compare --config compare.yaml --out results/pig2_vs_pig1
```

Differences are right minus left. Both runs must use the same context, candidate catalog, effective promoter definition and compatible measurements (including resolution and normalization). Cross-tissue runs are not accepted as full comparable effects by this command. When interpreting a difference, look at activity, contact, `support` and the gene total together, not only at `pace_score`.

## RNA-seq, methylation and other data

The main score uses only activity (ATAC/DNase/H3K27ac) and contact (Hi-C). Other data can be added as annotations, so they appear next to each prediction without changing the score:

| Data                                      | Input                                                        | What PACE does with it                               |
| ----------------------------------------- | ------------------------------------------------------------ | ---------------------------------------------------- |
| RNA-seq                                   | `inputs.expression` (`gene_id`, `sample_id`, `tpm`, `status`) | Gene expression annotation                           |
| DNA methylation (WGBS, RRBS)              | `inputs.methylation` (CpG counts)                            | Element and promoter methylation, with coverage kept |
| H3K4me1, H3K4me3, H3K27me3, H3K9me3, CTCF | `inputs.features`                                            | Chromatin state annotations                          |

Why RNA-seq is not multiplied into the score: a gene's expression would multiply both the top and the bottom of the fraction, so it cancels. Use it to filter out genes that are not expressed in your tissue.

Uncovered CpGs are reported as missing, not as unmethylated.

**Optional classifier.** If you have functional labels (for example CRISPRi results), you can train a separate classifier with `multiomics.mode: ml`. It writes `pace_ml_score` and, if calibration passes, `pace_ml_probability`. These columns are kept separate from `pace_score` and are never averaged with it.

## Configuration reference

Only the settings you are likely to change are listed. `resolved_config.yaml` in every output folder shows the full set with defaults filled in.

**General**

| Setting                                       | Default      | What it does                                                 |
| --------------------------------------------- | ------------ | ------------------------------------------------------------ |
| `run_id`                                      | `pace`       | Name of the analysis                                         |
| `target_level`                                | `individual` | `individual` needs exactly one animal. `population_mean` averages animals with equal weight |
| `context.species`, `.assembly`, `.context_id` | none         | Must match every input table exactly                         |
| `seed`                                        | `17`         | Random seed for anything stochastic                          |

**Catalog**

| Setting                          | Default          | What it does                                                 |
| -------------------------------- | ---------------- | ------------------------------------------------------------ |
| `catalog.profile`                | `canonical_grid` | `canonical_grid` uses fixed-width elements. `provided_regions` uses your own non-overlapping regions as-is |
| `catalog.width_bp`               | `500`            | Element width                                                |
| `catalog.candidate_radius_bp`    | `5000000`        | How far from a TSS to look for candidates                    |
| `catalog.include_promoter_units` | `true`           | Whether promoters of other genes compete as candidates       |

**Activity**

| Setting                              | Default            | What it does                                                 |
| ------------------------------------ | ------------------ | ------------------------------------------------------------ |
| `activity.panel`                     | `[ATAC, H3K27ac]`  | Which assays to use. Fixed for the whole run                 |
| `activity.combine`                   | `geometric_equal`  | Geometric mean of the assays in the panel                    |
| `activity.minimum_callable_fraction` | `0.0`              | Elements with less bigWig coverage than this become NA. 0.8 is a reasonable starting point |
| `activity.replicate_aggregation`     | `equal_donor_mean` | Technical → biological → animal averaging                    |

If an element is missing one assay of a two-assay panel, its activity is NA. PACE never quietly switches to a single assay for some elements, because that would mix two different scales in one gene's total.

**Contact**

| Setting                              | Default             | What it does                                                 |
| ------------------------------------ | ------------------- | ------------------------------------------------------------ |
| `contact.mode`                       | `observed`          | `observed`, `shrinkage` or `prior_only`                      |
| `contact.resolution`                 | none                | Hi-C bin size. Must match Step 3                             |
| `contact.scale` | `depth_normalized_contact` | Must match the actual measurement scale in Step 3 |
| `contact.normalization_id` | none | Must match the processing used in Step 3 |
| `contact.prior_path`                 | none                | Distance-decay prior, needed for `shrinkage` and `prior_only` |
| `contact.reliability`                | none                | Weight on observed contacts in `shrinkage` mode, 0–1         |
| `contact.near_diagonal_policy`       | `prior_or_neighbor` | How to handle elements in the same bin as the TSS            |
| `contact.pseudocount`                | `auto`              | Distance-based pseudocount for sparse matrices               |
| `contact.pseudocount_distance_bp`    | `5000`              | Distance the pseudocount is anchored at                      |
| `contact.allow_prior_fallback`       | `false`             | If `true`, missing observed contacts fall back to the prior  |

**Scoring and promoters**

| Setting                  | Default    | What it does                                                 |
| ------------------------ | ---------- | ------------------------------------------------------------ |
| `scoring.partial_policy` | `withhold` | Genes where not every candidate could be scored get no score |
| `promoters.weights`      | `provided` | `provided` uses the `pi` column of `promoters.tsv`. `equal` weights every distinct TSS of a gene equally. Use `equal` unless you have reliable TSS usage data (for example from CAGE) |

**Experimental allocation**

| Setting          | Default | What it does                                                 |
| ---------------- | ------- | ------------------------------------------------------------ |
| `allocation.eta` | `auto`  | Exponent on an optional cross-gene contact-share term; a biological competition interpretation has not been established. `auto` stays at 0 unless you supply functional labels that support a higher value. Leave it alone unless you have CRISPR data |

## Input table reference

All tables are tab-separated with a header row, UTF-8, optionally gzipped. Do not rename an `.xlsx` file to `.tsv`. Coordinates are 0-based, half-open, like BED. Write `NA` for missing values and `0` only for measured zeros.

Steps 1–3 generate most of these for you. You normally only write `samples.tsv` and `sources.tsv` by hand.

| Table                    | Columns                                                      |
| ------------------------ | ------------------------------------------------------------ |
| `units`                  | `element_id chrom start end anchor0 element_roles canonical_catalog_id` |
| `region_membership`      | `region_id element_id source_id membership_rule`             |
| `promoters`              | `gene_id promoter_id chrom tss0 strand pi pi_source`         |
| `candidates`             | `element_id gene_id candidate_universe_id`                   |
| `samples`                | `sample_id donor_id assay biological_replicate technical_replicate species assembly context_id source_id` |
| `sources`                | `source_id path_or_accession source_type assembly processing_method normalization_id checksum` |
| `observed_activity`      | `element_id sample_id assay signal measurement_status callable_fraction unit normalization_id window_id` |
| `observed_contacts`      | `element_id promoter_id sample_id contact_value measurement_status bin_pair_id scale resolution source_id` |
| `expression` (optional)  | `gene_id sample_id tpm status`                               |
| `methylation` (optional) | `chrom dyad_start0 methylated_count total_count sample_id assay` |
| `features` (optional)    | `entity_type entity_id feature_name value evidence_id status` |
| `labels` (optional)      | `label_id assayed_region_id gene_id context_id perturbation_type effect_direction effect_size label_status assay_id group_id source_id` |

`measurement_status` is one of `observed`, `unmeasured`, `low_coverage`, `unmappable`, `invalid`, `not_applicable`.

## The formula in full

The default model is the main equation above, with the complete planned candidate
set in its denominator. For a fixed assay panel and distinct TSSs:

```math
A_\star(E)=\left[\prod_{m\in\mathcal M}x_m(E)\right]^{1/|\mathcal M|},\qquad
\overline C(E,G)=\sum_{t\in\mathcal T(G)}\pi(t\mid G)\widetilde C(E,t).
```

ATAC=4 and H3K27ac=9 give activity 6. Ctilde includes the configured contact policy:
observations, near-diagonal correction, optional additive regularization, prior-only
contact or explicit shrinkage. A missing required TSS gives NA by default.

Optional `activity.pseudocounts` adds per-assay offsets after aggregation, without
filling NA. `promoters.minimum_weight` and `promoters.missing_policy: drop_missing`
select one shared TSS set per gene and renormalize its weights. At least
`promoters.minimum_retained_weight` (default 0.9) of the original weight must remain;
otherwise the gene stays unscoreable. Filtered scores describe only the selected
TSS definition. Both options are off by default; see [worked settings](docs/ADVANCED.md#sparse-activity-and-alternative-tsss).

The experimental allocation extension multiplies support by B(E,G)^eta, where B is
the element's contact share across its candidate genes. `eta=auto` uses zero without
eligible independent functional evidence. This extension is separate from the
main formula; see [equations and calibration](docs/ADVANCED.md#experimental-target-allocation).

## How PACE differs from ABC

The activity–contact normalization follows ABC. The implementations also depend
on candidate selection, contact processing and promoter definitions:

| Component | PACE implementation | Comparison requirement |
|---|---|---|
| Activity | Fixed single/two-assay panel; measured zero and NA are distinct | Match assays, windows and normalization |
| Candidates | Unique grid cells overlapping peaks and required promoters | Align candidates before comparing rankings or thresholds |
| TSSs | Weighted distinct TSSs, optional gene-wide filtering | Keep annotations and TSS definitions fixed |
| Replicates | Technical, biological and donor means in that order | Match the biological quantity being estimated |
| Sparse contact | Recorded correction, prior fitting and explicit transfer | Match contact scale and evaluate transfer independently |
| Missing denominator | Primary score withheld; conditional score and bounds separate | Report coverage alongside predictive accuracy |
| Additional omics | Annotations or a separately evaluated classifier | Evaluate each model output independently |

With identical activity, contact, one TSS and a complete candidate set, the default
normalization agrees with the ABC-style formula. That algebraic agreement does not
imply identical rankings from different preprocessing pipelines. Software tests do
not establish superiority over ABC. See [comparison methods](docs/ADVANCED.md#comparisons-and-benchmarks).

## Troubleshooting

**`pace: command not found`**
The environment is not active. Run `conda activate pace`, or `source .venv/bin/activate` if you used pip.

**File not found, but the file exists**
Paths in YAML files are relative to the YAML file's own folder. Paths on the command line are relative to your current folder. In Docker, the file must be inside the mounted folder.

**`Output directory already exists`**
Choose a new `--out` folder, or use `--force` to preserve the existing PACE result as a sibling backup before replacement.

**Chromosome name mismatch**
Your BED, GTF, bigWig and mcool disagree on names (`chr1` vs `1`, `chrMT` vs `MT`). Rename in one place so all files match.

**Most genes are `partial` or NA**
Usually Hi-C is missing for some TSSs. Look at the `reason` column in `scores.tsv.gz`. Check near-diagonal corrections and available contact coverage. Consider an appropriate resolution, an explicitly allowed prior fallback, or the documented gene-wide TSS filter. Equal weights do not repair missing contacts. Do not replace NA with 0.

**Everything is NA for one assay**
Check that `sample_id` and `assay` in `observed_activity.tsv` match `samples.tsv`, and that `activity.panel` uses the same assay names.

**`scale` or `normalization_id` mismatch**
The labels in `run.yaml` must match those used in Step 3. Changing a label does not convert data. If two Hi-C files really are on different scales, reprocess them the same way.

**`eta_calibration.json` says eta is 0**
That is normal. It stays 0 unless you supplied suitable functional labels.

**Only one animal. Is that OK?**
Yes. Use `target_level: individual`. You just can't say anything about variation between animals.

**My species isn't pig, cattle or chicken**
The readers accept species-independent formats. Use matching measured activity, annotation and reference data; a compatible contact prior is needed when contact is unavailable. Performance must be assessed for your species and tissue.

Still stuck? Open an [issue](https://github.com/shenlinyong/PACE/issues) and include the output of `pace --version`, your `run.yaml`, and `qc_report.json`.

## Citation

PACE is not yet published. If you use it, please cite this repository with the commit you used (`git rev-parse HEAD`), and the ABC papers it builds on:

> Fulco CP, Nasser J, Jones TR, et al. Activity-by-contact model of enhancer–promoter regulation from thousands of CRISPR perturbations. *Nature Genetics* 51, 1664–1669 (2019).
>
> Nasser J, Bergman DT, Fulco CP, et al. Genome-wide enhancer maps link risk variants to disease genes. *Nature* 593, 238–243 (2021).

## License and contact

MIT License. Written and maintained by Linyong Shen, Northwest A&F University. Questions and bug reports: [GitHub Issues](https://github.com/shenlinyong/PACE/issues).

[Advanced use](docs/ADVANCED.md): multiomics processing, classifier training, eta calibration, worked calculations, validation and complete defaults.
