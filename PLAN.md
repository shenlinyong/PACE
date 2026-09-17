# PACE development record

Updated 2026-09-17. Maintainer: 申林用 (Linyong Shen), Northwest A&F University.
Public repository: https://github.com/shenlinyong/PACE.

## Current model and software

The current branch provides one implementation in `src/pace_livestock`. All PACE
command aliases and the compatibility launcher call it. The mathematical source
of truth is [FORMULA.md](docs/FORMULA.md), with the detailed
[model](docs/model.md), [calibration contract](docs/eta_calibration.md) and
[independent calculations](PACE_Review_and_Validation.md).

Measured, hybrid and genome-only modes share the same bulk-proxy scoring order.
Automatic eta uses zero without suitable functional data, otherwise a bounded
training/calibration estimate; test labels do not fit eta. Main normalization uses
only actual scoreable support. Separate classifier outputs and QC retain their roles.

## Public documentation correction

The user identified incompatible formulas and outdated information still visible
in the public repository. Software 0.2.1 removes the retired code/configuration/
example paths from the current branch and rewrites the public documentation.
Git history preserves prior source and results. The
[content audit](notes/public_docs_audit.md) records preserved, replaced and removed
material; [migration](docs/migration.md) links the archival snapshot.

The local current-model suite passes 111 tests, including six new public-documentation
and launcher checks. Packaging and hosted CI checks accompany publication. Final executed results are recorded in
[validation](docs/validation.md). No real biological validation or livestock weights
are claimed from these synthetic software tests.

## Development practices

The work uses research-software-engineering, agent-resource-discipline and
human-facing-doc-authoring from scicomp-research-skills, with independent numerical
references, deterministic seeds, safe model artifacts and source hashing.
[Implementation decisions](notes/impl_canonical.md) and
[continuous-eta decisions](notes/impl_eta_cli.md) retain the numerical rationale.
