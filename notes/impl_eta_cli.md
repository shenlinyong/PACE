# Continuous allocation and the PACE command

The September 2026 user revision replaces the initial binary allocation switch with
a bounded continuous exponent. The numerical default remains zero. In automatic
mode, supplied, applicable functional training/calibration labels can estimate the
exponent; held-out test labels never enter that estimate.

## Numerical contract

Support is `A * Cbar * B**eta`, evaluated as `log(A) + log(Cbar) + eta*log(B)`.
At eta zero the implementation skips B, including its missing-data requirements.
Positive eta retains the fixed candidate-gene universe and true-zero behavior.

The estimator minimizes a gene-balanced pairwise logistic ranking loss on positive
versus powered-negative enhancer perturbations within each gene. Log support
differences are affine in eta; the objective is convex on [0, 1], so an
endpoint/derivative-bisection solver gives a deterministic bounded estimate. This
is a ranking surrogate, not a probability model for the PACE composition. No
performance improvement or optimal biological exponent is assumed. Three
informative genes are the conservative engineering default for attempting a fit,
not a biological sample-size justification. Fewer genes, no identifiable allocation
contrast, or no usable pairs produce an explicitly recorded zero fallback.

Label context, direction, power/status, gene/element/group separation and nonzero
scoreable evidence are checked before fitting. Invalid schemas or leaking splits
raise errors; inapplicable evidence is reported and excluded. A saved calibrator
records fitting support, input hashes, scope and objective; reuse must match the
scientific contract. Final performance needs independent functional data.

## Command contract and acceptance

`PACE --mode measured|hybrid|genome --out DIR` accepts either a run YAML or direct
input/model options. `PACE measured`, `PACE hybrid`, and `PACE genome` are aliases.
CLI paths resolve from the working directory; paths inside YAML resolve from the
YAML directory. Overrides are validated by the same configuration code. The
existing `pace-livestock` entry point remains usable.

Acceptance covers independently derived fractional scores and an interior eta
optimum, endpoint/zero/missing cases, split leakage and artifact mismatch, all
three modes with flags alone, installation outside the checkout, documentation,
full regression tests and hosted CI. Synthetic fixtures verify computation only.

Continuous calibration was introduced in software 0.2.0. Its numerical and behavior
tests remain in the current suite, including fractional-score and interior-optimum fixtures.
Direct flag-only measured, hybrid and genome modes, prefix installation and
wheel execution outside the checkout are exercised. The [validation record](../docs/validation.md)
and GitHub workflow track distribution and hosted regression results.

The subsequent public-documentation consolidation replaces duplicated initial
contracts with maintained entry points. Current formula/parameter pages and both
READMEs describe the continuous rule. Retired implementations and their defaults
are available only through Git history; see [the public audit](public_docs_audit.md).
