# PACE verification and biological scope

The current [software verification record](docs/validation.md) reports executed
numerical, schema, training, CLI, installation and CI checks. Independent mathematical
expectations are recorded in [the review](PACE_Review_and_Validation.md).

All bundled examples and model weights are synthetic software fixtures. Real
livestock weights, independent individual-effect validation and functional-link
accuracy are not bundled. Historical analyses using different software definitions
are not evidence of performance for the current formula.

Run `python -m pytest tests -q` after installing development dependencies. The
[GitHub workflow](https://github.com/shenlinyong/PACE/actions/workflows/ci.yml) also
checks public documentation consistency and installed package behavior.
