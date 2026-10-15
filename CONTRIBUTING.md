# Contributing to PACE

Report bugs with the command, a minimal input, the software commit and the relevant QC reason. Do not attach private data or credentials.

## Development checks

```bash
python -m pip install -e '.[dev,io,ml]'
python -m pytest tests -q
ruff check src/pace_livestock tests
ruff format --check src/pace_livestock tests
```

Numerical changes need an independent hand calculation. Missing values and measured zeros must stay distinct. Keep file readers separate from mathematical functions, and keep the README formulas, defaults and schemas consistent with the code. User instructions belong in [README](README.md); methods belong in [Advanced use](docs/ADVANCED.md).
