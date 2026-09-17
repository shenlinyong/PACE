# PACE changes

## 0.2.1 — 2026-09-17

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

See [migration](docs/migration.md) for access to superseded source history. These
release notes describe software behavior, not demonstrated biological superiority.
