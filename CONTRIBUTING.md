# Contributing to PACE

Report bugs with the command, a minimal synthetic input, the software commit and the relevant
QC reason. Please do not attach private animal or human data, credentials, or large model files.

```bash
python -m pip install -e '.[dev,io,ml]'
python -m pytest tests -q
ruff check src/pace_livestock tests/test_canonical_*.py
ruff format --check src/pace_livestock tests/test_canonical_*.py
python -m build
```

Numerical changes need an independent analytical expectation or a documented invariant.
Missing values and biological zeros must remain distinct. Preprocessing may only be fitted
inside the training split; real validation claims must cite actual evidence. Keep IO separate
from the mathematical functions. Loading models must not execute pickle or arbitrary code.

All scientific implementations belong in `src/pace_livestock`. Command aliases and compatibility
launchers must call this same package. Keep the public formula, examples and schemas consistent;
run `python scripts/check_public_docs.py` before publishing documentation changes. Superseded
implementations are accessible only through Git history. Preserve licence and attribution notices.

Use small focused commits and include tests and relevant documentation. The project follows
the numerical-correctness and reproducibility workflow of
[research-software-engineering](https://github.com/a-attia/scicomp-research-skills/tree/8435b16d91972c4f31b006de7bcacf1f5eb47e8e/skills/research-software-engineering),
with a Scientific Python src-layout and PEP 621 package metadata.
