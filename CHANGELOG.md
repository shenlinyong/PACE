# Changelog

## 0.8.0

- Every command takes ordinary command-line options (bedtools/samtools style); no configuration file is needed. YAML input is still read through a hidden `--config` to reproduce older runs.
- New `pace predict`: peaks + GTF + bigWig (+ Hi-C, a fitted prior or the human ABC power law) to scores in one command. With Hi-C it fits the distance power law of the same map for the pseudocount and for masked bins, as ABC does (`--strict-contacts` disables the fallback).
- New `catalog`, `activity`, `contacts`, `merge`, `features`, `methylation` and `expression` commands replace `prepare --config`.
- Samples and sources tables are inferred for single-animal runs; the contact scale and activity panel are read from the input tables.
- Chromosome sizes may be a chrom.sizes file, a FASTA index or a bigWig; `--gene-types` filters biotypes; several BED files are united; `--skip-unlisted-chroms` drops contigs absent from the sizes.
- Chromosome-chunked scoring (`run --by-chromosome`, default in `predict`) with `-t/--threads` parallel chunks; results equal a single run.
- About 2–5× faster scoring and Hi-C extraction; bounded memory per chunk. The core-formula evidence record no longer lists every parent ID in one cell.
- `pace init` writes a `run.sh` with command-line options instead of `config.yaml`.

## 0.7.0

- Add CTCF motif scanning, candidate-boundary priors and Poisson-profile Hi-C fitting.
- Export raw cooler counts and balancing factors; resolve per-pair Gamma–Poisson posteriors with masked-bin fallback and sample provenance.
- Add separate eQTL weak calibration of gamma, beta and eta with chromosome-separated selection, final holdout and explicit zero-eta fallback.
- Connect `boundaries`, `fit-hic`, `prior`, `fuse` and `fit-labels` to the main pipeline; include an offline four-chromosome example.
- Keep functional calibration separate and preserve true fixed-eta benchmark ablations.


- Consolidate user documentation into README and `docs/ADVANCED.md`; preserve advanced methods, labels, examples and full defaults.
- Remove obsolete installation wrappers and dependency snapshots; install with pip, Conda or Docker.
- Check the README formula, advanced defaults and local documentation links in CI.

## 0.6.0

- Add optional assay-specific activity pseudocounts without filling missing measurements.
- Add gene-wide TSS filtering with a retained-weight threshold and original/effective weight output.
- Allow explicit same-species, same-assembly contact-prior transfer between tissues; retain source context and clear target-validation claims.
- Start cooler prior fits at the matrix resolution by default, including 10–25 kb maps.
- Present activity × contact normalization as the main formula; document allocation as experimental.
- Add the lowercase `pace-livestock` alias, remove a case-only duplicate documentation filename and update citation metadata.
- Report incomplete genes as partial even when all measured supports are zero.
- Remove review-draft documents and update the bilingual method and workflow documentation.

## 0.5.0

- Withhold incomplete primary scores by default; report conditional scores and support sensitivity intervals.
- Add matched-scale power-law contact pseudocounts, recorded neighbor correction and same-bin TSS query reuse.
- Fit contact priors directly from sparse coolers, counting valid zero opportunities; provide an explicit unvalidated human shape baseline.
- Add region summaries, promoter-signal weights, CPM density normalization, project initialization, pairs and merge commands.
- Use the lowercase pace executable; force output replacement preserves a prior-result backup.
- Add contact-power benchmark control and a controlled-degradation validation protocol.


## 0.4.0

- Restrict activity to measured ATAC, DNase or H3K27ac and valid experimental aggregates.
- Remove hybrid/genome-only modes, sequence/fusion training and prediction, individual reconstruction and variant-effect scenarios.
- Remove their examples, optional deep-learning dependencies and obsolete commands/configuration keys; fail explicitly on unsupported configurations.
- Preserve measured multiomics annotations, multiple promoters, optional allocation and declared contact priors.
- Apply same-bin contact checks to imported evidence and validate imported assay/sample identity.
- Bind auxiliary assay definitions and methylation platform names to classifier inputs; report held-out raw scores, calibrated probabilities and deployment coverage separately.
- Update installation, formulas, CLI, manuscript scope and supported test/CI workflows.

This is a software scope change. It does not establish higher biological predictive accuracy.

## 0.3.0

### Scientific input and scoring checks

- Check contact resolution, normalization, balancing and measurement windows before
  pooling observations, applying a prior, importing resolved contacts or comparing runs.
- Preserve finite log support when an ordinary floating-point support underflows.
- Keep legitimate chromosome-boundary exclusions consistent between catalog
  preparation and loading. Preserve catalog identity when files are relocated.
- Use only the alleles selected by a sample's genotype, and check both ends of a
  reported breakend. Do not treat an unused alternative allele as a structural change.
- Bind exported sequence predictions to their reference, sample, variants,
  callability, ploidy, target coordinates, model and processing policy. Revalidate
  imported evidence and retain invalid-window masks.

### Learning and multi-omics

- Require independent grouped validation before accepting an automatically fitted
  nonzero allocation exponent; otherwise retain zero. Old unvalidated allocation
  artifacts must be refitted.
- Keep repeated perturbations and edges together in actual cross-validation folds.
- Bind ML assets to feature definitions, candidate rules and evidence policies;
  return explicit out-of-scope results for incompatible inputs.
- Exclude failed RNA observations, aggregate repeated annotations by donor, and
  preserve complete feature-to-evidence provenance.
- Add promoter methylation summaries, reference CpG coverage, configurable coverage
  thresholds, and separate WGBS/RRBS features. Expose BED motif strand features.

### Installation and use

- Add a Docker build and a Conda environment with common genomic IO dependencies.
- Replace the quick-start and model documentation with executable examples,
  mode-specific configurations, input preparation instructions, expanded equations
  and a Chinese user manual.
- Extend continuous integration to the regression suite, wheel installation,
  Docker and Conda workflows.

The bundled data and weights remain synthetic software examples. This release
does not establish predictive accuracy for a livestock species or tissue. Existing
run outputs without the new measurement contracts should be rerun before comparison.

## 0.2.1 — 2026-09-17

- Fix GitHub formula rendering with supported upright names and fenced math blocks;
  guard against rejected macros and Markdown consuming LaTeX escapes.
- Consolidate the public repository around one current scoring implementation.
- Rewrite the homepage, formula, notation, comparison, installation, tutorial and
  output references; remove incompatible equations, parameters and unsupported
  biological validation statements from the current branch.
- Retire superseded region/ML workflows, their examples and their model-specific
  tests. Preserve their source and results through Git history, not parallel current APIs.
- Make compatibility launchers call the installed package; update environment and
  dependency instructions. Add public-documentation and launcher checks to CI.

## 0.2.0 — 2026-09-17

- Add continuous eta, automatic functional-label calibration, frozen reuse and
  training/test isolation with independent numerical tests.
- Add PACE/pace executables, three mode aliases, direct file arguments and a prefix
  installer. Export eta applicability and fitting provenance with each run.

## 0.1.0 — 2026-09-17

- Introduce the installable src package, shared measured/hybrid/genome pipeline,
  strict schemas, provenance, quantitative CNN, contact/fusion fitting, separate
  classifier, comparisons and functional benchmarks.
- Add synthetic examples, wheel distribution, numerical tests and hosted CI.

Earlier implementations remain available in Git history. These release notes describe software behavior, not demonstrated biological superiority.
