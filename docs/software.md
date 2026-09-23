# Software scope and workflow

The package has one measured-activity implementation in `src/pace_livestock`.
`PACE`, `pace` and `pace-livestock` call the same command interface.

1. [Install](INSTALLATION.md) using Conda, venv or Docker.
2. [Prepare](input_preparation.md) processed experimental tracks, contact data and a fixed catalog.
3. [Validate and run](TUTORIAL.md) a measured-data configuration.
4. Check missingness, normalization background, evidence provenance and [formula interpretation](FORMULA.md).
5. Use [independent functional labels](training.md) for optional parameter fitting or classification, and [benchmark](comparison.md) on held-out data.

The software reports all planned candidates, uses transactional output directories and saves input/software hashes. Invalid configurations fail explicitly. [Migration](migration.md) lists incompatible API changes; historical implementations are retained in Git history.

Tests and demonstrations establish software behavior on controlled fixtures. No real livestock accuracy or independently validated universal weights are claimed. See [limitations](limitations.md).
