# PACE parameters

| Parameter | Default | Interpretation |
| --- | --- | --- |
| Activity mode | `missing_geometric` | Missing-aware shifted geometric mean |
| Activity priors | equal for enabled accessibility and H3K27ac | Declared design choices |
| Target allocation exponent | 1 | Enhancer-centric allocation |
| Inhibitory attenuation strength | 0 | No inhibitory attenuation in primary predictions |
| Residual support | computationally 0, marked unknown | No estimator is supplied |
| Evidence threshold | 0.5 | Operational input-quality cutoff, not FDR |
| Descriptive score cutoff | 0.02 | Requires independent calibration for decision use |
| Candidate cis window | 5 Mb | Strict midpoint-to-TSS distance bound |

Distance-prior constants and the exact equations are given in [METHODS.md](METHODS.md). Unknown quality remains unknown. RNA-expression options retained for compatibility do not multiply the current primary score. Historical expression weights, when comparing the earlier model, are applied after normalization as defined in [NOTATION.md](NOTATION.md).
