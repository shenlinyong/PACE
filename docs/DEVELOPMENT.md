# Development and verification

The Python package contains pure scoring kernels, evidence resolution, file adapters, preparation commands and reporting. Genomic measurements remain distinct from priors and independently trained classifier outputs. Keep changes to evidence policies reflected in frozen calibration, classifier and comparison contracts.

```bash
python -m pip install -e '.[dev,io,ml]'
python -m pytest tests -q
python scripts/check_public_docs.py
ruff check src tests scripts
ruff format --check src tests scripts
python -m build
```

Independent numerical tests cover incomplete denominators, support bounds, zero counts, near-diagonal correction, same-bin TSS reuse, contact-prior scales and sparse valid-pair counting. File-based tests exercise actual bigWig/cooler adapters where their optional dependencies are installed. Core-only CI intentionally skips optional readers.

CI also installs wheels outside the checkout and runs the Conda environment and Docker image. Installation checks do not establish biological performance. [The validation protocol](ROBUSTNESS_VALIDATION.md) describes the required external analyses.
