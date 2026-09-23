# PACE development record

Maintainer: Linyong Shen (申林用), Northwest A&F University.

The current package focuses on measured activity and declared contact evidence.
The source of truth is [the formula](docs/FORMULA.md). Entry points, configuration,
dependencies, examples and documentation use this scope. Removed APIs and historical
source are documented in [migration](docs/migration.md).

Verification covers independent numerical expectations, explicit failures for
unsupported inputs, measured-data round trips, contact integrity, grouped functional
validation, experimental annotations and package installation. Run the checks in
[validation](docs/validation.md) for the exact commit. Biological accuracy requires
appropriate independent data and is not inferred from passing software tests.
