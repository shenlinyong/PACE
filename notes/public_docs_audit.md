# Public model and documentation consolidation

The user requested complete removal of incorrect and superseded information from
the current repository on 2026-09-17. Adding scope notices above incompatible
formulas did not achieve that. This revision makes the current branch describe
and execute one scientific model. Earlier commits remain available in Git history.

## Content inventory and disposition

| Material | Previous location | Destination or reason for removal |
|---|---|---|
| Identity, author, repository, licence and citation | README, AUTHORS, CITATION | Preserve in current entry points and attribution files |
| Installation, commands, input/output guide links | README and uppercase guides | Rewrite around the installed PACE executable and current schemas |
| Main score and allocation exponent | README, FORMULA, NOTATION, METHODS, Chinese guide | One authoritative FORMULA page; all displayed equations use the scoreable support sum and continuous eta policy |
| Activity, contacts, TSS and missing-data interpretation | Formula, comparison, parameter guides | Replace shifted/missing-assay aggregation and implicit priors with the implemented fixed-panel and qualified-contact rules |
| ABC comparison and numeric examples | ABC_COMPARISON, WORKED_EXAMPLES | Preserve algebraic comparisons with explicit scope; use current independent analytical fixtures |
| Optional ML and additional omics | ML_INTEGRATION, old configs and scripts | Document the separate current classifier and annotation roles; remove obsolete weighting paths |
| Old install recipes, Snakemake and raw-read commands | setup, environment, workflow, example | Remove incompatible workflow; setup delegates to the current installer, preparation docs define processed-input boundaries |
| Wrong-model code and generated predictions | legacy, old scripts, example_quantified, workflow | Remove from the current branch; retain original files in Git history, without recomputing or relabelling old results |
| Old numerical regression tests | Non-canonical test modules | Retire with the removed implementation; retain every current model test and add public-contract checks |
| Duplicated development/model appendices | Root model/development documents | Replace with maintained entry points and implementation requirements; retain valid independent hand calculations in the review |
| Historical biological counts and performance claims | VALIDATION and old README | Remove from current validation because evidence for the current model is not bundled |
| Model scope, limitations, training and scientific checks | Current lowercase guides | Audit against current source and keep supported content; remove stale implementation-status statements |
| Development history and completed work | PLAN, CHANGELOG, implementation notes | Rewrite as current status plus concise dated change history; link the archived commit instead of reproducing superseded claims |

## Acceptance

Audit tracked Markdown, configuration, code, examples and entry points. Public
formula displays must agree with the numerical kernel; automatic eta uses zero
without applicable calibration and fits only eligible train/calibration labels.
Current guides must not offer old commands, output fields or validation claims.
Run the current model suite, documentation/link checks, installation/build checks
and hosted CI, and verify the published default-branch files after pushing.

Local result: 111 current tests pass, including six new documentation/launcher checks.
All 105 existing current-model tests were retained; 48 tests specific to the removed
implementation were retired with it. Public links, displayed equations and configuration
defaults pass the automated audit. Distribution and hosted results are tracked in
[the verification record](../docs/validation.md) and the linked CI workflow.
