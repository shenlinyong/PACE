# Contributing to PACE

Report bugs with the command, a minimal synthetic input, the software commit and the relevant
QC reason. Please do not attach private animal or human data, credentials, or large model files.

```bash
python -m pip install -e '.[dev,io,sequence,ml]'
python -m pytest tests -q
ruff check src/pace_livestock tests/test_canonical_*.py
ruff format --check src/pace_livestock tests/test_canonical_*.py
python -m build
```

Numerical changes need an independent analytical expectation or a documented invariant.
Missing values and biological zeros must remain distinct. Preprocessing may only be fitted
inside the training split; real validation claims must cite actual evidence. Keep IO separate
from the mathematical functions. Loading models must not execute pickle or arbitrary code.

New canonical features belong in `src/pace_livestock`. The existing scripts and workflow are
a separate region-based interface; changes to their behavior require explicit regression
evidence and documentation. Additions should preserve their API and license notices.

Use small focused commits and include tests and relevant documentation. The project follows
the numerical-correctness and reproducibility workflow of
[research-software-engineering](https://github.com/a-attia/scicomp-research-skills/tree/8435b16d91972c4f31b006de7bcacf1f5eb47e8e/skills/research-software-engineering),
with a Scientific Python src-layout and PEP 621 package metadata.
