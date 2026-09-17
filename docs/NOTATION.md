# Model notation

Entity-dependent quantities use function notation throughout the current manuscript, model documentation, responses and Fig. 1. The score is always $\mathit{PACE}(E,G)$; model versions are identified by their surrounding text.

| Symbol | Definition |
| --- | --- |
| $E$, $G$ | Focal enhancer and gene |
| $e$, $g$, $t$, $i$ | Summation indices for enhancers, genes, TSSs and assays |
| $A(E)$ | Enhancer activity |
| $S(E,i)$, $x(E,i)$ | Assay signal and scaled signal; $x(E,i)=S(E,i)/a_i$ |
| $m(E,i)$, $q(E,i)$, $w_i$, $a_i$ | Observation indicator, measurement quality, assay prior and positive assay scale |
| $\kappa$ | Overall attenuation strength; `inhibition_strength` in the core API; zero disables attenuation |
| $\mathrm{Scale}$, $\mathrm{Pseudocount}$, $\gamma$ | Distance-prior constants, retaining the previous formula names |
| $I(E)$ | Inhibitory score, retaining the previous model's symbol; the current optional extension averages observed, scaled context values |
| $d(E,t)$, $H(E,t)$, $\lambda(E,t)$ | Enhancer–TSS distance, observed contact and reliability |
| $P(d)$, $D(d)$ | Distance prior and compatible expected observed contact |
| $C_{\mathrm{adj}}(E,t)$ | Reliability-adjusted enhancer–TSS contact |
| $\mathcal T(G)$, $\pi(G,t)$ | Distinct TSS set and TSS-use weights for a gene |
| $C(E,G)$ | Enhancer–gene contact; the current model uses the weighted average over distinct TSSs |
| $\mathcal G(E)$, $B(E,G)$ | Candidate target-gene set and enhancer-centric target allocation |
| $R(E,G)$ | Raw pair support, $A(E)C(E,G)B(E,G)^\eta$ |
| $\mathcal E(G)$, $\mathcal E^{\mathrm{obs}}(G)$ | Candidate-enhancer set and subset with finite raw support |
| $U(G)$ | Nonnegative independently estimated residual support; unknown and computationally zero in the primary configuration |
| $Q(E,G)$ | Separate input-evidence index |
| $Q_A(E)$, $Q_C(E,G)$, $Q_T(G)$, $Q_{\mathrm{cat}}(G)$ | Activity, observed-contact, TSS and catalogue quality |
| $W_{\mathrm{expr}}(G)$ | Gene expression weight in the intended previous model only |

Equal symbols identify equal roles, not equal preprocessing. In particular, old PACE, current PACE and ABC can supply different activity/contact values. The current scoring kernel and numerical predictions are unchanged by this notation revision.

The intended previous score is

$$
\mathit{PACE}(E,G)=W_{\mathrm{expr}}(G)
\frac{A(E)C(E,G)}{\sum_{e\in\mathcal E(G)}A(e)C(e,G)}.
$$

For a positive denominator its within-gene sum is $W_{\mathrm{expr}}(G)$. A strictly positive expression weight preserves within-gene ranking. A zero weight ties all finite scores at zero.

The current score is

$$
\mathit{PACE}(E,G)=\frac{A(E)C(E,G)B(E,G)^\eta}
{\sum_{e\in\mathcal E^{\mathrm{obs}}(G)}A(e)C(e,G)B(e,G)^\eta+U(G)}.
$$

The earlier expression-weighted formula is included for model comparison. Reintroducing its weight into the current primary score would change cross-gene rankings and thresholds and would require new evaluation; it is not an algebraic rewrite of the existing model.

## Previous-to-current notation and implementation

| Earlier or inconsistent notation | Canonical notation | Code field or argument |
| --- | --- | --- |
| Earlier indexed or upright score notation | $\mathit{PACE}(E,G)$ | `PACE.Score`; identify the model version in prose |
| $A_e$ | $A(E)$ | `activity` |
| $C_{eg}$, $\bar C_{eg}$ | $C(E,G)$ | `contact_gene` |
| $R_e$, $R_{\mathrm{inh}}(E)$ | $I(E)$ | `repression` |
| $S_i$ (previous assay signal), $s(E,k)$ | $S(E,i)$ | input assay signal; enhancer and assay are explicit |
| $B_{eg}$, $R_{eg}$, $U_g$, $Q_{eg}$ | $B(E,G)$, $R(E,G)$, $U(G)$, $Q(E,G)$ | `target_share`, `raw_support`, `unassigned_mass`, `evidence_quality` |
| $\pi_{gt}$, $\lambda_{et}$ | $\pi(G,t)$, $\lambda(E,t)$ | TSS-use weights and observed-contact reliability |
| $s$, $d_0$ in the distance prior | $\mathrm{Scale}$, $\mathrm{Pseudocount}$ | distance scale and distance pseudocount |

The previous inhibitory formula used $\beta$ as the H3K27me3-specific weight.
It must not also denote the new overall attenuation strength. The latter is
$\kappa$, and the current activity factor is $\exp(-\kappa I(E))$.
At $\kappa=0$ this factor is exactly one. The previous activity used $1-I(E)$;
these are different model operations, even though $I(E)$ has the same inhibitory-score role.
The previous inhibitory formula also reused $\gamma$ for an H3K9me3 weight;
in the current model $\gamma$ denotes only the distance-decay exponent.
The current API uses named assay-weight dictionaries, so it does not need to
reuse either Greek letter for individual inhibitory-assay weights.

The previous $E_{5Mb}$ is the candidate set at a fixed 5-Mb window.
The current $\mathcal E(G)$ makes the gene dependence explicit and allows
the declared configurable window; $\mathcal E^{\mathrm{obs}}(G)$ is its
finite-support subset. These sets must not be conflated when candidates are missing.

Software identifiers, filenames and raw-data columns remain stable interface names;
the table maps them to the manuscript notation. PACE as a software/model name is
distinct from the mathematical score $\mathit{PACE}(E,G)$.
