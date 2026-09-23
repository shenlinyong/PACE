# Implemented methods

PACE constructs a frozen set of non-overlapping regulatory units, deduplicated physical TSSs and cis candidate links. Qualified measured signals are aggregated by technical replicate, biological replicate and donor. A fixed ATAC/DNase/H3K27ac panel produces activity by an equal geometric mean. Missing required layers remain unavailable and measured zeros are preserved.

Contacts use a predeclared observed/prior policy on compatible scales and resolutions, then fixed promoter weights. Optional allocation uses the complete candidate-gene contact set and an exponent in [0,1]; automatic selection falls back to zero without eligible independent functional evidence. The primary gene-wise score is withheld when planned support is incomplete; conditional scores and support-assumption sensitivity bounds are reported separately. Matched priors can regularize sparse observed contacts, and same-bin contacts use a declared prior or measured-neighbor correction.

RNA, other histone marks, CTCF and methylation are annotations by default. Optional supervised outputs are evaluated separately from the primary formula. Comparisons recompute common denominators, and final functional tests remain separate from parameter selection.

See [equations](FORMULA.md), [input preparation](input_preparation.md), [training](training.md), [comparison](comparison.md), and [limitations](limitations.md). These implementation details do not assert superior biological accuracy.
