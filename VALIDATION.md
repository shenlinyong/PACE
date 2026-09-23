# PACE verification and biological scope

The current [software verification record](docs/validation.md) reports executed
numerical, schema, training, CLI, installation and CI checks. Independent mathematical
expectations are exercised by the tests described in the
[development guide](docs/DEVELOPMENT.md).

Bundled example data and learned model weights are synthetic software fixtures.
The optional `abc_human` contact-shape preset is a published human reference,
explicitly unvalidated for a target livestock context. Real livestock learned
weights, independent individual-effect validation and functional-link accuracy
are not bundled. Historical analyses using different software definitions are
not evidence of performance for the current formula.

Run `python -m pytest tests -q` after installing development dependencies. The
[GitHub workflow](https://github.com/shenlinyong/PACE/actions/workflows/ci.yml) also
checks public documentation consistency and installed package behavior.
