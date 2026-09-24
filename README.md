# PACE

**Enhancer–gene prediction for livestock and poultry.**

PACE ranks candidate regulatory elements for each gene in a tissue using measured activity and contact evidence. It is built on the Activity-by-Contact (ABC) model and adapted to the data that livestock labs actually have: incomplete epigenomic panels, sparse or low-resolution Hi-C, rough transcript annotation, and only one or a few animals per tissue.

The input formats are species-independent and require a matching reference, annotation and experimental activity data. Pig, cattle, sheep, goat and chicken are the main use cases.

Every command is driven by ordinary command-line options, like `bedtools`, `samtools` or `macs3`. No configuration file is needed.

There are three ways to run it, depending on the Hi-C you have: Hi-C of the sample (A), Hi-C of the same species (B), or none (C). See [Which case am I?](#which-case-am-i-three-data-situations). The example below is case A.

```bash
git clone https://github.com/shenlinyong/PACE.git && cd PACE
conda env create -f environment.yml && conda activate pace

pace predict -b peaks.bed -g genes.gtf --atac atac.bw --h3k27ac h3k27ac.bw \
  --hic pig1.mcool --hic-resolution 10000 \
  --species pig --assembly Sscrofa11.1 --tissue liver -t 4 -o pig1_liver
```

---

## Contents

1. [How it works](#how-it-works)
2. [Is PACE right for your data?](#is-pace-right-for-your-data)
3. [Installation](#installation)
4. [Test run](#test-run)
5. [Quick start: one command](#quick-start-one-command)
6. [Step by step](#step-by-step)
7. [Several animals](#several-animals)
8. [No Hi-C? Low-resolution Hi-C?](#no-hi-c-low-resolution-hi-c)
9. [Whole genomes: memory and threads](#whole-genomes-memory-and-threads)
10. [Reading the results](#reading-the-results)
11. [Comparing animals, tissues or treatments](#comparing-animals-tissues-or-treatments)
12. [RNA-seq, methylation and other data](#rna-seq-methylation-and-other-data)
13. [Command reference](#command-reference)
14. [Input table reference](#input-table-reference)
15. [The formula in full](#the-formula-in-full)
16. [How PACE differs from ABC](#how-pace-differs-from-abc)
17. [Troubleshooting](#troubleshooting)
18. [Citation](#citation)

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
  demo -o results/docker_demo
```

Any path you pass to the container must be inside the mounted directory.

`pace` alone lists every command; `pace COMMAND --help` shows its options with defaults and an example.
`pace-livestock` is an alias of `pace`.

**No git?**

```bash
wget -O PACE-main.tar.gz https://github.com/shenlinyong/PACE/archive/refs/heads/main.tar.gz
tar -xzf PACE-main.tar.gz && cd PACE-main
```

For published analyses, record the exact version you used:

```bash
pace --version
git rev-parse HEAD
```

## Test run

```bash
pace demo -o results/demo
```

If `results/demo/results/scores.tsv.gz` appears, the installation works. The demo data are synthetic. They test the software, not biology.

The example folders contain small scripts that use only command-line options:

```bash
cd examples/measured && bash run.sh      # score prepared tables
cd ../contact && bash run.sh             # boundary priors, sparse Hi-C and eQTL tuning
```

`pace validate` checks inputs and computes a scoring preview without writing a result directory. By default `-o/--out` must be a new folder. `--force` replaces a recognized PACE result only after keeping the old one as a sibling backup.

## Quick start: one command

`pace predict` runs the whole pipeline: it builds candidate elements from your peaks, quantifies each bigWig, extracts Hi-C contacts, and scores every element–gene pair.

Every file must use the **same genome assembly and the same chromosome names** (`chr1` and `1` do not match).

| Input | Option | How to get it |
|---|---|---|
| Peaks or candidate regions (BED) | `-b/--peaks` | MACS2/MACS3 on ATAC or DNase; several files (for example one per animal) are united |
| Gene annotation (GTF) with `transcript` lines | `-g/--gtf` | Ensembl or NCBI for your assembly |
| ATAC-seq, DNase-seq and/or H3K27ac bigWig | `--atac`, `--dnase`, `--h3k27ac` | deepTools `bamCoverage --normalizeUsing CPM` |
| Hi-C `.cool` or `.mcool` (recommended) | `--hic` and `--hic-resolution` | HiC-Pro, distiller or Juicer (then `hic2cool`); balance with `cooler balance` |
| Chromosome sizes (optional) | `-c/--chrom-sizes` | `chrom.sizes`, `genome.fa.fai` or any bigWig; default: the first bigWig |
| Species, assembly, tissue | `--species --assembly --tissue` | Free text, recorded in every output |

### Which case am I? Three data situations

Activity data (peaks + at least one ATAC, DNase or H3K27ac bigWig) is always required. What differs between users is **contact** (Hi-C) data. Pick the one row that matches you:

| Case | Your Hi-C situation | Option | Contact used for scoring | Parameters to set by hand |
|---|---|---|---|---|
| **A** | Hi-C from the same animal/tissue you score | `--hic FILE` | Measured contacts; a power law fitted on the same map fills masked bins | None |
| **B** | No Hi-C for this sample, but **any** Hi-C of the same species and genome assembly (another tissue, another animal, a public dataset) | `--prior DIR` from `pace fit-prior` | Distance power law with γ fitted on **your species' own** Hi-C | None (γ is fitted) |
| **C** | No Hi-C for this species at all | `--abc-prior` | Human ABC power law, γ = 1.024 | None (γ is built in) |

Prefer A over B over C. The contact prior is `contact(d) = a × (max(d, d_min)/d_ref)^(−γ)`: contact falls as a power of the element–TSS distance `d`. The amplitude `a` cancels in PACE scores, so only γ (how fast contact decays) changes the ranking. See [No Hi-C? Low-resolution Hi-C?](#no-hi-c-low-resolution-hi-c) for details.

**Case A – Hi-C of this sample**

```bash
pace predict \
  -b data/liver_peaks.bed -g data/Sus_scrofa.Sscrofa11.1.gtf.gz \
  --atac data/pig1_ATAC.bw --h3k27ac data/pig1_H3K27ac.bw \
  --hic data/pig1.mcool --hic-resolution 10000 \
  --species pig --assembly Sscrofa11.1 --tissue liver \
  -t 4 -o results/pig1_liver
```

`--hic-resolution` picks the bin size inside an `.mcool` file (a single-resolution `.cool` needs none). PACE also fits the distance power law of this same map, as ABC does, and uses it for a small distance-based pseudocount and for element–TSS pairs in masked (unbalanceable) Hi-C bins. Each such pair is labelled `contact_prior` in `resolved_contacts.tsv`. Add `--strict-contacts` to leave those pairs NA instead.

**Case B – Hi-C of the same species, but not of this sample**

Step 1, once per species/assembly: fit γ from the Hi-C you have. Name the tissue the Hi-C came from.

```bash
pace fit-prior --cooler public/pig_liver.mcool::/resolutions/10000 \
  --species pig --assembly Sscrofa11.1 --tissue liver -o priors/pig_liver
```

Step 2: score your samples with that prior.

```bash
pace predict -b data/liver_peaks.bed -g data/genes.gtf --atac data/pig2_ATAC.bw \
  --prior priors/pig_liver \
  --species pig --assembly Sscrofa11.1 --tissue liver -o results/pig2_liver
```

The fitted γ is stored in `priors/pig_liver/manifest.json`; `distance_bins.tsv` holds the observed and fitted decay for a plot. Species and assembly must match exactly. A prior from **another tissue** (for example liver Hi-C used for muscle) is allowed but must be requested explicitly; build the tables step by step and add `--allow-cross-context-prior` to `pace run` (see [Transfer between tissues](docs/ADVANCED.md#transfer-between-tissues)). The output then records the transfer as unvalidated.

**Case C – no Hi-C for this species**

```bash
pace predict -b data/liver_peaks.bed -g data/genes.gtf --atac data/pig1_ATAC.bw \
  --abc-prior \
  --species pig --assembly Sscrofa11.1 --tissue liver -o results/pig1_liver
```

This uses the human ABC reference γ = 1.024238616787792. It was validated with CRISPR perturbations in human cells, not in livestock, so the outputs are labelled `unvalidated_for_target_context`. State this in your methods, and switch to case B as soon as any Hi-C of your species becomes available.

Several bigWigs for one assay (`--atac a.bw b.bw`) are averaged as replicates of the same animal. Supported activity panels: ATAC, DNase or H3K27ac alone, ATAC + H3K27ac, or DNase + H3K27ac. Evaluate the chosen panel on independent data; adding a low-quality assay does not guarantee better rankings.

The output folder holds the results (`scores.tsv.gz` and the files described in [Reading the results](#reading-the-results)) and every intermediate table under `prepared/`. The folder is self-contained: `resolved_config.yaml` records every setting with paths relative to the folder, so it can be moved or archived with a publication.

Useful options (`pace predict --help` lists all):

| Option | Default | Meaning |
|---|---|---|
| `-w/--width` | 500 | Element width on a fixed genome grid (bp) |
| `-r/--radius` | 5000000 | Pair each gene with elements up to this distance from its TSS |
| `--gene-types` | all | Keep only these biotypes, for example `protein_coding` (Ensembl) or `mRNA` (NCBI) |
| `--skip-unlisted-chroms` | off | Ignore peaks and genes on contigs missing from the chromosome sizes |
| `--min-callable-fraction` | 0 | Elements with less bigWig coverage become NA; 0.8 is a common choice |
| `--missing-as-zero` | off | Treat bases absent from the bigWig as measured zero signal |
| `--partial-policy` | withhold | Genes with unscored candidates: `withhold` or `conditional` |
| `-t/--threads` | 1 | Chromosome chunks scored in parallel |

`radius` is a search window, not a biological claim. 5 Mb matches ABC. If you change it, check that your top predictions are stable.

## Step by step

`predict` is a shortcut for the commands below. Run them yourself to reuse a catalog across animals, to combine replicates, or to inspect each step.

| Step | Command | Output |
|---|---|---|
| 1. Candidates | `pace catalog` | `units.tsv`, `promoters.tsv`, `candidates.tsv`, `region_membership.tsv` |
| 2. Activity | `pace activity` (once per bigWig), then `pace merge` | `observed_activity.tsv` |
| 3. Contact | `pace contacts` (Hi-C) and/or `pace fit-prior` | `observed_contacts.tsv`, prior folder |
| 4. Scores | `pace run` | `scores.tsv.gz` and QC reports |

```bash
# 1. candidate elements and distinct TSSs
pace catalog -b data/peaks.bed -g data/genes.gtf -c data/genome.fa.fai -o prepared/catalog

# 2. activity per bigWig, then one table
pace activity -i data/pig1_ATAC.bw    -a ATAC    -d prepared/catalog -o prepared/pig1_atac
pace activity -i data/pig1_H3K27ac.bw -a H3K27ac -d prepared/catalog -o prepared/pig1_k27ac
pace merge -t observed_activity \
  -i prepared/pig1_atac/observed_activity.tsv prepared/pig1_k27ac/observed_activity.tsv \
  -o prepared/activity

# 3. Hi-C contacts for every element-TSS pair, plus the distance prior of the same map
pace contacts -i data/pig1.mcool -r 10000 -d prepared/catalog -o prepared/pig1_hic
pace fit-prior --cooler data/pig1.mcool::/resolutions/10000 \
  --scale balanced_contact --normalization-id cooler_weight \
  --species pig --assembly Sscrofa11.1 --tissue liver -o prepared/pig1_prior

# 4. score
pace run -d prepared/catalog \
  --activity prepared/activity/observed_activity.tsv \
  --contacts prepared/pig1_hic/observed_contacts.tsv \
  --contact-prior prepared/pig1_prior --allow-prior-fallback \
  --species pig --assembly Sscrofa11.1 --tissue liver \
  --by-chromosome -t 4 -o results/pig1_liver
```

Notes:

- `pace catalog` builds a non-overlapping grid of `--width` bp cells overlapping your peaks and all TSSs, not ABC peak-centred windows. `promoters.tsv` has one row per distinct physical TSS (transcripts sharing a TSS are counted once). Broad or noisy peaks enlarge the denominator.
- `pace activity`: bigWig files often leave out regions with no reads. That can mean "measured, zero signal" or "not measured". If you are not sure how your bigWig was made, do not use `--missing-as-zero`. The sample name defaults to the file name; `--normalization-id` is a label for how the bigWig was normalized and must be the same for files you want to compare.
- `pace contacts` needs contact frequencies (balanced by default, or raw with `--no-balance`). Observed/expected matrices, p-values or loop calls give wrong results, because they remove the distance effect PACE needs. Do not mix resolutions or normalizations. For `.hic` input, first convert to cool/mcool with a compatible converter and check the stored resolutions and balancing weights.
- `pace run` without `--samples` treats every sample as a replicate of one animal and writes the tables it inferred (`inferred_samples.tsv`, `inferred_sources.tsv`) to the result folder. The contact scale label is read from the contact table.

## Several animals

Score each animal separately first and check that the top predictions agree between animals. Then, for a tissue-level map, run once with all animals and a sample table that says which animal each library came from:

```
sample_id	donor_id	assay	biological_replicate	technical_replicate	species	assembly	context_id	source_id
pig1_ATAC	pig1	ATAC	1	1	pig	Sscrofa11.1	liver	pig1_ATAC
pig2_ATAC	pig2	ATAC	1	1	pig	Sscrofa11.1	liver	pig2_ATAC
pig1_HiC	pig1	HiC	1	1	pig	Sscrofa11.1	liver	pig1_HiC
```

```bash
pace run -d prepared/catalog --activity prepared/all_activity.tsv \
  --contacts prepared/all_contacts.tsv --samples prepared/samples.tsv \
  --target-level population_mean \
  --species pig --assembly Sscrofa11.1 --tissue liver -o results/liver_population
```

`species`, `assembly` and `context_id` in the table must be spelled exactly as on the command line. PACE averages technical replicates within each biological replicate, then biological replicates within each animal, then gives every animal equal weight. Technical replicates never count as extra animals. A `sources.tsv` table recording where each file came from (GEO/SRA accession, processing, checksum) is optional; `pace validate` reports any value that is not accepted.

## No Hi-C? Low-resolution Hi-C?

Hi-C is the data type livestock projects most often lack. PACE has three contact modes (`pace run --contact-mode`):

| Mode | Uses | When to use |
| --- | --- | --- |
| `observed` | Measured Hi-C, optional matching-prior regularization | A contact matrix with suitable coverage and resolution |
| `shrinkage` | Observed contact and a compatible distance prior | Sparse maps: a fixed observation weight, or per-pair Gamma–Poisson shrinkage (`pace fuse`) |
| `prior_only` | Distance-decay prior only | No usable Hi-C |

The distance-decay prior has the form

```
contact_prior(d) = a × ( max(d, d_min) / d_ref ) ^ (−γ)
```

where `d` is the element–TSS distance. Fit it on Hi-C from your species with `pace fit-prior`, which counts callable zero pixels:

```bash
pace fit-prior --cooler pig_liver.mcool::/resolutions/10000 \
  --species pig --assembly Sscrofa11.1 --tissue liver -o pig_prior
pace predict ... --prior pig_prior        # or: pace run ... --contact-prior pig_prior --contact-mode prior_only
```

The lower distance defaults to the matrix resolution. A prior from a different tissue requires `--allow-cross-context-prior`; source and target contexts are retained and the transfer is labelled unvalidated. The species, assembly, target level and measurement definitions must match. Tissue transfer requires evaluation; it is not assumed equivalent to target-tissue Hi-C.

Without any usable Hi-C, `--abc-prior` (`pace run --prior-preset abc_human --contact-mode prior_only`) uses the human ABC reference gamma=1.024238616787792 on a relative scale. It is not a fitted livestock parameter. See [the model](docs/ADVANCED.md#fitting-and-transferring-a-prior).

For `shrinkage`, either give a fixed weight with `--contact-reliability 0.5 --reliability-source <how it was chosen>`, or use `pace fuse`, which estimates a weight per pair from raw counts and the prior (`pace contacts` exports the counts and conversion factors). See [boundary priors and weak calibration](docs/ADVANCED.md#boundary-priors-and-weak-calibration).

Results from `prior_only` and `shrinkage` runs are labelled as such in the output. Report the mode you used in your methods.

Two settings matter most with low-resolution Hi-C:

- Near-diagonal pairs (default `prior_or_neighbor`). Elements in the same Hi-C bin as the TSS would otherwise get the diagonal value, which mostly reflects bin self-contact. PACE uses the prior if one is given, or else the recorded maximum of valid neighbouring contacts. If neither is available, contact remains NA.
- `--pseudocount auto` (default). Adds a distance-based term to finite off-diagonal measurements only when a compatible fitted prior is available. Without one, observations are unchanged. Missing values are not converted to measured zeros.

**Boundary-aware priors and eQTL tuning.** Optional commands add a CTCF-boundary attenuation to the prior (`pace boundaries`, `pace prior`, `pace fit-hic`), per-pair Gamma–Poisson shrinkage (`pace fuse`) and chromosome-held-out tuning of gamma, beta and eta on eQTL fine-mapping (`pace fit-labels`, then `pace run --weak-model`). `examples/contact/run.sh` runs all of them on synthetic tables. These components have software tests, not an established livestock accuracy gain. PIP mass is association evidence, not a CRISPR response or a calibrated link probability. Equations, input formats and limitations are in [Advanced use](docs/ADVANCED.md#boundary-priors-and-weak-calibration).

## Whole genomes: memory and threads

Everything in the default model is computed within a chromosome, so PACE can score whole chromosomes separately and concatenate them without changing any score. `pace predict` does this by default; add `--by-chromosome` to `pace run`. `-t/--threads N` scores N chromosome chunks in parallel.

Measured on synthetic data at livestock density (60 peaks and 8 genes per Mb, 5 Mb radius): about 1.05 million element–gene pairs per 120 Mb took 2.5 minutes with 4 threads and 1.6 GB of memory per thread. A whole livestock genome has roughly 20 million pairs: expect tens of minutes with 4–8 threads, and memory set by the largest chromosome (roughly 6 GB per million pairs of that chromosome). Lower `--radius` reduces both.

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

Do not subtract two score tables directly. Because scores are shares, a change in one element shifts every other element of the same gene. Use `pace compare`, which rescores both runs over the candidates they have in common:

```bash
pace compare --left results/pig1_liver --right results/pig2_liver -o results/pig2_vs_pig1
```

Differences are right minus left. Both runs must use the same context, candidate catalog, effective promoter definition, compatible measurements (including resolution and normalization) and the same execution mode (both chunked or both unchunked). Cross-tissue runs are not accepted as full comparable effects by this command. When interpreting a difference, look at activity, contact, `support` and the gene total together, not only at `pace_score`.

## RNA-seq, methylation and other data

The main score uses only activity (ATAC/DNase/H3K27ac) and contact (Hi-C). Other data can be added as annotations, so they appear next to each prediction without changing the score:

| Data | Prepare with | Add to `pace run` with | What PACE does with it |
| --- | --- | --- | --- |
| RNA-seq | `pace expression --tpm transcripts.tsv -d prepared/catalog` | `--expression` | Gene expression annotation |
| DNA methylation (WGBS, RRBS) | `pace methylation --counts cpg.tsv -d prepared/catalog` | `--methylation` | Element and promoter methylation, with coverage kept |
| H3K4me1, H3K4me3, H3K27me3, H3K9me3, CTCF peaks | `pace features -b ctcf.bed -d prepared/catalog` | `--features` | Chromatin state annotations |

Why RNA-seq is not multiplied into the score: a gene's expression would multiply both the top and the bottom of the fraction, so it cancels. Use it to filter out genes that are not expressed in your tissue.

Uncovered CpGs are reported as missing, not as unmethylated.

**Optional classifier.** If you have functional labels (for example CRISPRi results), you can train a separate classifier with `pace train` and apply it with `pace run --ml-model` or `pace predict-ml`. It writes `pace_ml_score` and, if calibration passes, `pace_ml_probability`. These columns are kept separate from `pace_score` and are never averaged with it.

## Command reference

`pace COMMAND --help` shows every option with its default.

| Command | Purpose | Key options |
|---|---|---|
| `predict` | Peaks + GTF + bigWig (+ Hi-C) to scores in one step | `-b -g --atac/--dnase/--h3k27ac --hic/--prior/--abc-prior -o` |
| `catalog` | Candidate elements, TSSs and pairs | `-b -g -c -w -r --gene-types -o` |
| `activity` | Signal per element from a bigWig | `-i -a -d -s -o` |
| `contacts` | Element–TSS contacts from cool/mcool | `-i -r -d --no-balance -o` |
| `merge` | Combine tables with duplicate checks | `-t -i -o` |
| `run` / `validate` | Score prepared tables / check them | `-d --activity --contacts --contact-prior --species --assembly --tissue -o` |
| `fit-prior` | Distance prior from a cool/mcool map | `--cooler --species --assembly --tissue -o` |
| `boundaries`, `prior`, `fit-hic` | Boundary-aware priors | see `--help` |
| `fuse` | Per-pair shrinkage of sparse contacts | same as `run`, plus `--contact-prior` |
| `fit-labels` | eQTL tuning of gamma, beta, eta | `--run --eqtl --test-chromosomes -o` |
| `compare`, `stability`, `benchmark` | Downstream evaluation | `--left --right`, `--replicates`, `--run --labels --membership` |
| `fit-eta`, `train`, `predict-ml` | Functional-label models | see `--help` |
| `features`, `methylation`, `expression` | Annotation tables | see `--help` |
| `pairs`, `normalize-activity`, `promoter-weights`, `init`, `demo` | Helpers | see `--help` |

Main `pace run` settings:

| Option | Default | What it does |
|---|---|---|
| `--target-level` | `individual` | `individual` needs exactly one animal. `population_mean` averages animals with equal weight |
| `--panel` | assays present | Assays used for activity: `ATAC`, `DNase`, `H3K27ac`, or accessibility + `H3K27ac` |
| `--contact-mode` | `observed` | `observed`, `shrinkage` or `prior_only` |
| `--contact-prior` | none | Distance-decay prior, needed for `shrinkage` and `prior_only` |
| `--allow-prior-fallback` | off | Missing observed contacts use the prior |
| `--pseudocount` | `auto` | Distance-based pseudocount for sparse matrices |
| `--partial-policy` | `withhold` | Genes where not every candidate could be scored get no score |
| `--eta` | `auto` | Experimental allocation exponent; `auto` stays 0 without a calibrator |
| `--by-chromosome`, `-t` | off, 1 | Chunked scoring and parallel chunks |

Every setting and its default is listed in [Advanced use](docs/ADVANCED.md#complete-configuration-defaults); `resolved_config.yaml` in every result folder records the complete set used. Without `--panel`, the panel is the set of ATAC, DNase and H3K27ac assays present in the activity table.

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
Paths on the command line are relative to your current folder. In Docker, the file must be inside the mounted folder.

**`Output already exists`**
Choose a new `-o` folder, or use `--force` to keep the existing PACE result as a sibling backup before replacement.

**Chromosome name mismatch / "not in the chromosome sizes"**
Your BED, GTF, bigWig and mcool disagree on names (`chr1` vs `1`, `chrMT` vs `MT`). Rename in one place so all files match, or give a renaming table with `--aliases`. If the extra names are unplaced contigs you do not need, add `--skip-unlisted-chroms`.

**Most genes are `partial` or NA**
Usually Hi-C is missing for some TSSs or elements (masked bins). Look at the `reason` column in `scores.tsv.gz`. `pace predict --hic` uses the distance prior of the same map for masked bins; with `pace run`, give `--contact-prior` and `--allow-prior-fallback`. Equal weights do not repair missing contacts. Do not replace NA with 0.

**"requires tables for every declared activity assay"**
`--panel` names an assay that is not in the activity table. Check the `assay` column, or leave `--panel` out to use the assays present.

**`scale` or `normalization_id` mismatch**
All contact tables and the prior must describe the same measurement. Changing a label does not convert data. If two Hi-C files really are on different scales, reprocess them the same way.

**`eta_calibration.json` says eta is 0**
That is normal. It stays 0 without a calibrator, and fitted models can also select 0. Supply functional labels with `fit-eta`, or use `fit-labels` and `pace run --weak-model` for separately marked eQTL weak calibration.

**Only one animal. Is that OK?**
Yes. That is the default (`--target-level individual`). You just can't say anything about variation between animals.

**Out of memory**
Use chunked scoring (default in `predict`, `--by-chromosome` in `run`), fewer `--threads`, or a smaller `--radius`.

**My species isn't pig, cattle or chicken**
The readers accept species-independent formats. Use matching measured activity, annotation and reference data; a compatible contact prior is needed when contact is unavailable. Performance must be assessed for your species and tissue.

Still stuck? Open an [issue](https://github.com/shenlinyong/PACE/issues) and include the output of `pace --version`, the command you ran, and `qc_report.json`.

## Citation

PACE is not yet published. If you use it, please cite this repository with the commit you used (`git rev-parse HEAD`), and the ABC papers it builds on:

> Fulco CP, Nasser J, Jones TR, et al. Activity-by-contact model of enhancer–promoter regulation from thousands of CRISPR perturbations. *Nature Genetics* 51, 1664–1669 (2019).
>
> Nasser J, Bergman DT, Fulco CP, et al. Genome-wide enhancer maps link risk variants to disease genes. *Nature* 593, 238–243 (2021).

## License and contact

MIT License. Written and maintained by Linyong Shen, Northwest A&F University. Questions and bug reports: [GitHub Issues](https://github.com/shenlinyong/PACE/issues).

[Advanced use](docs/ADVANCED.md): multiomics processing, classifier training, eta calibration, worked calculations, validation and complete defaults.
