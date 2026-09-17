# Current software verification

Executed local checks on **2026-09-17** for software **0.2.1**. The tests use synthetic
fixtures, including actual small genomic file formats. They verify software behavior;
no real livestock weights, independent functional benchmark or individual-effect
accuracy is claimed.

## Executed checks

| Check | Result |
|---|---|
| Current test suite | **111 passed** |
| Numerical references | Hand-derived endpoint and fractional eta scores, known interior calibration optimum, activity, TSS, fusion, bulk-order and CpG calculations pass |
| Evidence and normalization | Missingness/zero distinction, fixed panel and B universe, partial/common denominators and overflow checks pass |
| Calibration | Context and split isolation, zero fallback, bounded fitting, frozen reuse, incompatible scope and invalid artifact rejection pass |
| Genomic IO | Small bigWig and sparse cooler fixtures, BED/GTF, CpG, RNA and variant/callability behavior pass |
| Sequence and ML | Actual CNN gradient/update/save/load and grouped classifier/preprocessing checks pass |
| Installed command | Three direct-file modes, file paths with spaces, YAML overrides and output protection pass |
| Compatibility launchers | scripts/pace.py invokes the same fractional-scoring implementation; setup.sh delegates to the current installer |
| Public documentation | Main displayed equations match FORMULA.md; documented defaults match configuration; relative links resolve; retired model terms and trees are rejected |
| Source style | Ruff lint/format, shell syntax and git patch checks pass |

The preceding mixed-interface suite had 153 tests. This revision retains its 105
current-model tests, removes 48 tests together with the retired implementation, and
adds six public-documentation/launcher checks. No current-model test was removed.

One cooler 0.10.4 warning concerns NumPy timedelta construction; it does not change
the tested outputs. The local environment uses Python 3.13.12, NumPy 2.5.3 and
PyYAML 6.0.3, with optional IO and sequence dependencies installed. The
[recorded dependency snapshot](../requirements-tested-python313.txt) identifies the
research environment; it is not a claim that every dependency is required for scoring.

## Verification map

| Contract | Implementation | Test modules |
|---|---|---|
| Primary formula and evidence | core/scoring.py, config.py, schemas.py, pipeline.py | canonical_core, canonical_pipeline |
| Continuous eta | learning/allocation.py | canonical_allocation |
| Genomic inputs and units | catalog/, io/, sequence/genome.py | canonical_adapters, canonical_genome |
| Quantitative model and supervised learning | sequence/, learning/model.py | canonical_sequence, canonical_learning |
| Preparation, comparison and evaluation | operations.py, evaluation/ | canonical_operations, canonical_pipeline |
| Public commands and documentation | cli.py, run_options.py, scripts/check_public_docs.py | canonical_cli, canonical_documentation |

Each numerical assertion names an independent mathematical expectation or invariant.
Small synthetic training runs test correct optimization and persistence rather than
an arbitrary accuracy target. [Worked examples](WORKED_EXAMPLES.md) and the
[review](../PACE_Review_and_Validation.md) give the reference calculations.

## Distribution and hosted checks

Source and wheel distributions use the same PACE package. Distribution acceptance
runs all three modes outside the checkout, verifies the installed source fingerprint,
and exercises eta fitting and frozen reuse. The
[GitHub workflow](https://github.com/shenlinyong/PACE/actions/workflows/ci.yml) checks
Python 3.11–3.13, minimum core dependencies, IO/sequence/ML extras, current regression,
public documentation and installed-wheel behavior. Its status is attached to each
source commit; a previous commit's result is not validation of a new revision.

## Biological validation still required

Cross-locus quantitative prediction, within-locus individual changes and enhancer–gene
link prediction are separate tasks. Each requires suitable independent evidence,
applicable sampling and explicit coverage accounting. QTL association is auxiliary
support, not functional ground truth. No historical dataset counts, unbundled study
results or synthetic checks establish current-model biological performance.
See [limitations](limitations.md).
