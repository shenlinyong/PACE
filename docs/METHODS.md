# Implemented methods

PACE uses one [mathematical definition](FORMULA.md) for measured, hybrid and genome
inputs. Detailed biological assumptions and research extensions are in
[the model description](model.md); executable contracts are in [parameters](parameters.md).

## Shared calculation

1. Freeze nonoverlapping scoring units, distinct TSSs and candidate enhancer–gene edges.
2. Resolve each declared assay from compatible observations, quantitative predictions
   or an independently calibrated combination. Preserve missing values and true zeros.
3. Form the equal geometric activity of the fixed panel. For bulk sequence evidence,
   average each assay across copies before forming activity.
4. Resolve contact from the declared observed/prior policy and average distinct TSS
   contacts with frozen promoter weights. Retain distance background in contact.
5. Use eta zero by default without applicable calibration; otherwise fit eligible
   functional training/calibration labels or reuse their frozen artifact.
6. Multiply activity, contact and optional allocation, then normalize over the
   actual scoreable candidates of each gene. Export coverage, denominator and provenance.

## Inputs and fitted assets

Measured mode needs qualified activity and contact or an explicit compatible prior.
Hybrid mode additionally permits matching quantitative sequence evidence; actual
fusion needs a target-specific calibrator. Genome mode needs trained quantitative
sequence weights and contact prior, reference/individual sequence and a frozen
candidate catalog. No universal species weights or contact constants are supplied.

The sequence CNN, contact-prior fit, signal-fusion fit, eta calibration and separate
elastic-net classifier have distinct targets. Their data splits and validity checks
are documented in [training](training.md). Additional omics are annotations unless
explicitly used by that classifier. The primary score remains a support composition.

## Reproducibility and evaluation

Runs export input/model/software hashes, resolved configuration, evidence tables,
actual eta and normalization identities. Tests compare independent hand calculations,
missing-data invariants and installed command behavior. Independent biological
accuracy still requires appropriate functional and individual-level datasets.
See [validation](validation.md), [comparisons](comparison.md) and [limitations](limitations.md).
