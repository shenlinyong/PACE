# Software verification record

This record describes executed software checks on **2026-09-17**. It is not a biological
validation study. All new examples, labels, variants and model weights used here are synthetic.
No true livestock quantitative weights, paired-animal effects or independent perturbation
benchmarks were supplied in this task.

## Executed local checks

| Check | Observed result |
|---|---|
| Complete repository test suite | **127 passed**, including retained region-interface tests |
| Canonical package statement coverage | **78.5%** (2,008 of 2,558 statements); a diagnostic, not a biological-quality score |
| Independent mathematical fixtures | Exact rational eta 0/1 scores, activity, multi-TSS, fusion, bulk-order and CpG examples passed |
| Missingness/denominators | True zero, missing assay/contact, fixed B universe, partial/common normalization and overflow passed |
| Real small-file IO | bigWig missing/negative signal, sparse cool queries, BED/GTF, stranded CpG and RNA mappings passed |
| Sequence handling | REF, GT, ploidy, phase sets, callability, multiallelic SNV, fixed-target indel and reported SV checks passed |
| Actual CNN | CPU gradients/loss finite, masked heads respected, parameters updated, safetensors loading and central-position response passed |
| Grouped classifier | Fold isolation, label exclusions, train-only preprocessing, actual auxiliary-feature use and calibration passed |
| Offline modes | measured, hybrid and genome_only ran through the same scoring kernel |
| Actual command workflows | Training → CNN predictions → genome scoring; prior/fusion fitting; train/predict ML; compare/stability/benchmark/variant scenarios succeeded |
| Distribution | sdist and wheel built successfully; clean wheel-only environment ran all three demos outside the source checkout |
| Offline guarantee in wheel test | NumPy/PyYAML-only environment; torch/pandas absent; socket connection calls disabled during all demos |
| Style/patch checks | Ruff lint/format and git diff --check passed |
| Mutation checks | In isolated copies, deleting B allocation or replacing geometric by arithmetic activity caused independent reference tests to fail |

One dependency warning was observed: cooler 0.10.4 uses a NumPy timedelta construction that
NumPy 2.5.3 deprecates. It did not change the tested outputs. The initial legacy-suite attempt
reported a missing matplotlib dependency; after installing that documented legacy dependency,
the complete suite passed. No legacy mathematical function was changed to make these tests pass.

## Requirements-to-evidence map

| Development contract | Main implementation | Test evidence |
|---|---|---|
| C/E: schemas, source identity, window/assembly consistency, numerical core | config.py, schemas.py, core/scoring.py | test_canonical_core.py, test_canonical_pipeline.py |
| C: grid, de-duplication, sparse candidates, format adapters | catalog/, io/ | test_canonical_adapters.py |
| D/F: three regimes, asset scope, same bulk order, safe quantitative models | evidence/, sequence/, pipeline.py | pipeline/genome/sequence tests |
| E: common denominator, full/conditional Delta, compositional interpretation | evaluation/compare.py | pipeline counterexamples and actual compare/stability commands |
| F: priors, independent-target fusion | evidence/contact.py, evidence/fusion.py, operations.py | test_canonical_operations.py |
| G: labels, groups, fixed transforms, elastic net, calibration | learning/model.py | test_canonical_learning.py |
| H: functional-label baselines, coverage and AP | evaluation/benchmark.py, metrics.py | operation and metric tests |
| I/J: packaging, offline use, preserved history, documentation | pyproject.toml, CI, demos, bilingual guides | local distribution test and GitHub workflow |

Each numerical reference test names the independent source or analytical construction. The
test constants are not regenerated from current software output. Tests deliberately avoid
requiring a small synthetic training run to reach an arbitrary prediction-accuracy threshold.

## Measured synthetic scale

`scripts/measure_canonical_performance.py` was run with 20,000 units and 1,000 genes. It emitted
99,880 sparse cis edges: candidate generation **0.615 s**, eta-1 scoring **1.228 s**, process peak
RSS **135.8 MiB**. This measures candidate generation plus the numerical kernel, including
interpreter/process memory; it does **not** measure genomic track IO, CNN inference, end-to-end
whole-genome analysis or biological accuracy. Single measurements are descriptive, not a
performance guarantee.

Environment: Python 3.13.12, Linux x86_64 / glibc 2.28; NumPy 2.5.3, PyYAML 6.0.3,
PyTorch 2.14.0+cpu, cooler 0.10.4, pyBigWig 0.3.26, pysam 0.24.1. Full local package versions are
recorded in [requirements-tested-python313.txt](../requirements-tested-python313.txt).
The core CI matrix covers Python 3.11–3.13 and separately tests minimum core dependencies;
IO, sequence, ML and legacy regression jobs are separate. Current hosted results are linked
from the [Actions workflow](https://github.com/shenlinyong/PACE/actions/workflows/ci.yml).

## Scientific validation still required

Cross-locus quantitative prediction, within-locus individual changes and enhancer–gene link
prediction remain **separate unvalidated tasks** until appropriate real held-out evidence is
provided. Population scope, functional-label direction/power, calibration independence,
candidate discovery and end-to-end missing-positive accounting need study-specific evaluation.
QTL association is auxiliary evidence. No accuracy improvement or general livestock readiness
is inferred from the tests in this record. See [limitations](limitations.md).
