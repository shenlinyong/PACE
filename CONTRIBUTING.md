# Contributing to PACE

Report bugs with the command, a minimal synthetic input, the software commit and
the relevant QC reason. Do not attach private data, credentials or large model files.

## Development checks

From the repository root:

```bash
python -m pip install -e '.[dev,io,ml]'
python -m pytest tests -q
python scripts/check_public_docs.py
ruff check src/pace_livestock tests scripts
ruff format --check src/pace_livestock tests scripts
python -m build
```

Numerical changes need an independent analytical expectation or a documented
invariant. Missing values and measured zeros must remain distinct. Fit preprocessing
only inside the training split; support biological validation claims with actual
evidence. Keep file readers separate from mathematical functions. Loading models
must not execute pickle or arbitrary code.

All scientific implementations belong in `src/pace_livestock`. The installed
`pace` and `pace-livestock` commands and `python -m pace_livestock` use the same
package. Keep formulas, examples and schemas consistent. User instructions belong
in [README](README.md); detailed methods belong in [Advanced use](docs/ADVANCED.md).
The documentation checker verifies the README equation, the full defaults in
ADVANCED and local links, including section anchors.

Use small, focused commits with relevant tests and documentation. Preserve license
and attribution notices. Earlier implementations remain available in Git history.

## Developer utilities

- `python scripts/generate_canonical_examples.py` regenerates the bundled synthetic
  examples. Review the resulting diff before committing it; the script does not
  download experiments or real trained weights.
- `python scripts/measure_canonical_performance.py` measures sparse candidate
  generation and the scoring kernel. Its JSON records timings, Linux process RSS
  and the environment. It does not measure biological accuracy or full-pipeline
  throughput.

## Release checks

Run the regression suite with optional IO dependencies installed, inspect skipped
tests, and complete the checks above. Build and install the wheel outside the
checkout, then run the measured demo and an explicitly chosen fractional eta.
Check the Conda, Docker and supported-Python CI jobs for the exact commit.

Review the diff for unintended files, private data, generated results and large
assets. Record the tested environment and known limitations in release notes.
Push without rewriting history and tag a release only after its checks pass.
A source commit, package version, GitHub release, container publication and DOI
are separate records; cite only artifacts that exist. Biological model assets
also need their validation scope and redistribution license.

The [validation protocol](docs/ADVANCED.md#validation) describes independent
analyses needed for biological claims. Software tests alone do not establish
livestock prediction accuracy.
