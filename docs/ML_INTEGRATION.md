# Additional omics and supervised learning

RNA, auxiliary histone marks, CTCF and methylation are annotations by default.
They do not multiply the primary [PACE formula](FORMULA.md). With appropriate
functional labels, the separate elastic-net classifier can use specified features.
Its outputs are `pace_ml_score` and, after appropriate calibration, `pace_ml_probability`.

[Training](training.md) describes grouped tuning, training-only preprocessing and
fixed inference with `pace train` / `pace predict-ml`. This classifier is separate
from [continuous eta fitting](eta_calibration.md); neither converts the primary
support share into a universal causal probability.
