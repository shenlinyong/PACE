# Computational methods

PACE (Prediction of Activity-based regulatory Connections for Enhancers) predicts relative support for cis enhancer–gene links using a shared scoring kernel. The [equations](FORMULA.md) are the mathematical specification; [parameter rationale](PARAMETERS.md) identifies defaults and available controls.

## Candidate and promoter preparation

Use a fixed, same-assembly enhancer catalogue and complete available gene annotation before score filtering. The supplied narrowPeak adapter selects peaks by `signalValue`, extends summits and applies optional blacklist filtering. Its TSS-region option identifies overlaps; it does not add missing promoter candidates or merge overlapping intervals. If a study requires a union of enhancer and promoter windows, construct that catalogue upstream and use it consistently for all models.

`prepare_tss.py` reads gene/transcript GTF records, converts positions to BED0, keeps distinct TSSs per stable gene ID and labels gene-boundary fallback. It retains all provided biotypes. File predictors generate same-chromosome pairs whose interval-midpoint-to-TSS distance is strictly below the configured window. Numerical-table scoring instead uses the candidate pairs supplied by the analyst.

## Activity and contact

Activity combines nonnegative, appropriately scaled signals through the missing-aware shifted geometric mean. Unknown measurements and measured zero remain distinct. Assay priors and quality inputs have separate meanings. The primary configuration uses accessibility and available H3K27ac with equal priors; optional inhibitory attenuation is disabled.

Contact combines a positive distance prior with qualified observed/expected contact through supplied local reliability. Metadata must establish compatible units and source provenance. Missing prerequisites cause prior fallback rather than automatic acceptance of an unscaled observation. Distinct promoter contacts are averaged using gene-level TSS-use weights defined before windowing.

## Target allocation and gene normalization

Each enhancer's contact is allocated over its candidate genes. Raw support is activity × gene contact × target allocation raised to its configured exponent. The final score divides this support by finite gene-level support plus optional independently supplied residual mass. Unknown residual mass remains flagged. Unscorable candidates remain visible in the unfiltered result.

The activity–contact basis is from [Fulco et al. (2019)](https://doi.org/10.1038/s41588-019-0538-0); enhancer-centred allocation and the multiple-promoter principle draw on [Hecker et al. (2023)](https://doi.org/10.1093/bioinformatics/btad062). PACE's exact operations and defaults are not interchangeable with every ABC or STARE release. See [the comparison](ABC_COMPARISON.md).

## Evidence and implementation

Input evidence is reported separately from relative score. Activity, contact, promoter and catalogue quality determine an operational evidence state, with reason codes for unknown or low support. RNA provides stable-ID context and optional downstream eligibility filtering; it does not multiply structural scores.

The table CLI and file adapters call [workflow/scripts/pace_core.py](../workflow/scripts/pace_core.py). All gene-level candidates are written before optional output thresholding. Runtime dependencies and environment records are described in [Installation](INSTALLATION.md); [Validation](../VALIDATION.md) distinguishes numerical/software verification from biological evaluation.

## Reproducible biological use

Record source accessions, assembly/annotation versions, preprocessing, replicate handling, candidate/TSS construction, signal scales, contact QC, model settings and the exact commit. Calibration or comparison must preserve the declared data budget and include omitted predictions in coverage reporting. Software examples are synthetic, and a successful run does not by itself demonstrate functional enhancer regulation.
