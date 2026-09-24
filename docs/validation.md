# Verification and release acceptance

Software verification and biological validation answer different questions. The test suite uses synthetic data and small real file formats to check implementation behavior. It supplies no measured accuracy claim for an animal species or tissue.

## Reproduce software checks

```bash
python -m pip install -e '.[dev,io,ml]'
python -m pytest tests -q
python scripts/check_public_docs.py
ruff check src tests scripts
ruff format --check src tests scripts
python -m build
```

Inspect skipped optional IO tests and the [CI run for the exact commit](https://github.com/shenlinyong/PACE/actions/workflows/ci.yml). A configured CI job alone does not mean it has passed.

## Required behavioral coverage

| Area | Invariant |
|---|---|
| Formula | Independent hand calculations; fixed assay and candidate sets; exact zeros and NA |
| Measured scope | Retired modes/options/configurations fail; unavailable activity is never imputed |
| Contact | Compatible resolutions and sources; same-bin imports use the same policy as raw inputs |
| Imports | Real sample identities and correct assay; measured activity only |
| Allocation | Zero fallback, grouped validation, bounded fit and test-set isolation |
| Comparison | Common denominators; underflow-safe reloads; partial results remain conditional |
| Multiomics | Experimental annotations, measurement definitions, CpG coverage and RNA status |
| Classifier | Grouped splits, inference scope, separate score and probability evaluation |
| Packaging | Supported Python versions, wheel outside checkout, measured demo, Conda and Docker |
| Documentation | Shared total equation, current defaults, working links and executable interfaces |

## Biological acceptance

Use independent functional perturbation labels suitable for the claim. Record candidate selection, sampling design, tissue and experiment. Do not treat input contact data or associations as independent functional truth. Separate fitting, selection and final evaluation; report coverage and unavailable tested positives. Comparisons across families, breeds or tissues require a split appropriate to that claim. See [limitations](limitations.md).
