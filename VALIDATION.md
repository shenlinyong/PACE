# Validation scope

Software verification and biological evaluation answer different questions.

## Reproducible software checks

On 17 September 2026, all 48 standalone tests passed with no skips. All eight smoke commands passed; the maximum absolute score error against the independent formula was 1.12 × 10⁻¹⁶.

Run `python -m pytest tests -q` from the repository root. The tests cover activity aggregation, gene normalization, contact reliability, multi-TSS behavior, missing measurements, command-line adapters and expression weighting after normalization in the archived sensitivity calculator. Read-count tests require bedtools; bigWig tests require pyBigWig.

Run `python scripts/smoke_test.py --output-dir results/smoke` for the eight-command small-file workflow. It creates three enhancers and two genes, checks both prediction interfaces against an independent formula calculation, and exercises measured-zero and missing BEDPE contacts. The expected output contains six unfiltered and four filtered links.

The bundled five-step read example also completed successfully, producing 12,000 unfiltered enhancer–gene predictions. Its generated outputs are not tracked as biological evidence.

The earlier workspace validation comprised 56 tests, including 12 analysis-specific tests outside this software repository. This repository ships the original 44 standalone software cases plus four legacy-expression regression cases. Its own passing test count must be obtained from the command above, rather than inferred from the larger workspace suite.

## Biological scope

The associated analysis used 10,356 human development records, 1,986 historically inspected external records and 24 livestock maps. Comparators included ABC, distance, native-input TargetFinder, released JEME catalogues and explicitly adapted GATv2EPI and EPIPDLF models. Inputs and training targets differed. Historically inspected external data do not constitute an untouched prospective validation set.

The livestock analysis used distance-prior scoring for 49,962,511 candidate pairs and retained 603,804 distal links at a descriptive 0.02 cutoff. Genetic-support comparisons, contact sensitivity and reporter-locus reanalysis have their own limitations. These data, large prediction maps and manuscript materials are not bundled as software test fixtures.

Scores and evidence indices are not calibrated probabilities. Current tests do not validate a complete FASTQ-to-prediction run, all Snakemake scheduling, or every real binary Hi-C/Cooler pathway. No data-repository DOI is implied by this software release.
