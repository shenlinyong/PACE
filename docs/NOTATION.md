# Mathematical notation and software fields

PACE (Prediction of Activity-based regulatory Connections for Enhancers) uses $\mathit{PACE}(E,G)$ for the mathematical score. Entity-dependent quantities use function notation.

| Symbol | Meaning | Software representation |
| --- | --- | --- |
| $E$, $G$, $t$, $i$ | Enhancer, gene, TSS, assay | Coordinates, stable target ID, `TargetGeneTSS`, named signal |
| $S(E,i)$ | Preprocessed assay signal | Input signal column |
| $a_i$ | Fixed positive assay scale | Activity JSON `scale` |
| $x(E,i)$ | Scaled signal $S(E,i)/a_i$ | Computed before aggregation |
| $m(E,i)$ | Measurement present/absent | Finite value versus `NA` |
| $w_i$, $q(E,i)$ | Assay prior and measurement quality | JSON `weight`, `quality_column` |
| $A(E)$ | Enhancer activity | `activity` |
| $\kappa$, $I(E)$ | Inhibitory strength and aggregate fraction | API `inhibition_strength`, output `repression` |
| $P(d)$ | Positive contact prior | `contact_prior` |
| $H(E,t)$, $D(d)$ | Observed and expected contact | `contact_observed`, `contact_expected` |
| $\lambda(E,t)$ | Contact reliability used in mixing | Input `contact_reliability`, computed `contact_weight` |
| $\pi(G,t)$ | TSS-use weight | `tss_weight` |
| $C_{\mathrm{adj}}(E,t)$ | Adjusted TSS-level contact | `contact` before aggregation |
| $C(E,G)$ | Gene-level contact | `contact_gene` |
| $B(E,G)$ | Enhancer-centred target share | `target_share` |
| $\eta$ | Allocation exponent | `competition_power` |
| $R(E,G)$ | Raw activity/contact/allocation support | `raw_support` |
| $\mathcal E^{\mathrm{obs}}(G)$ | Supplied candidates with finite support | Defines `observed_mass` |
| $U(G)$ | Independently supplied residual support | `unassigned_mass` |
| $\mathit{PACE}(E,G)$ | Relative gene-normalized score | `PACE.Score` |
| $Q_A,Q_C,Q_T,Q_{\mathrm{cat}}$ | Activity, contact, TSS and catalogue QC | `activity_quality`, `contact_quality`, `tss_quality`, `catalogue_quality` |
| $Q(E,G)$ | Minimum evidence component, unknown retained | `evidence_quality` |

Distance units are bp. Signal, observed-contact and residual-support units must be documented with the input profile. See [equations](FORMULA.md), [parameters](PARAMETERS.md) and [outputs](IO_FORMATS.md).
