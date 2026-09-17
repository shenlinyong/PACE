# PACE mathematical definition

This is the authoritative formula for the installed PACE software. All three
evidence modes use this definition. [Parameters](parameters.md), [input schemas](data_dictionary.md)
and [functional calibration](eta_calibration.md) specify how data enter it.

## Total score

```math
\boxed{
\mathrm{PACE}(E,G)=
\frac{A_\star(E)\,\overline C(E,G)\,[B(E,G)]^{\eta_{\mathrm{used}}}}
{\displaystyle\sum_{e\in\mathcal E^{\mathrm{score}}(G)}
 A_\star(e)\,\overline C(e,G)\,[B(e,G)]^{\eta_{\mathrm{used}}}}
}
```

Write the unnormalized support as

```math
S(E,G)=A_\star(E)\,\overline C(E,G)\,[B(E,G)]^{\eta_{\mathrm{used}}}.
```

Here E is a fixed, nonoverlapping scoring unit and G is a target gene.
The normalization set consists of the planned candidates whose support can be
computed under the declared evidence policy. The denominator is exactly their
support sum. It has no added residual contribution or numerical pseudocount.

In automatic mode,

```math
\eta_{\mathrm{used}}=
\begin{cases}
0, & \text{without applicable or sufficiently informative functional calibration},\\
\widehat\eta\in[0,1], & \text{after fitting eligible training/calibration evidence}.
\end{cases}
```

A matching frozen calibration supplies the same fitted value without refitting.
Explicit fixed numeric exponents in [0,1] are also supported. Final test labels
never enter parameter estimation. The [calibration objective and checks](eta_calibration.md)
are separate from the primary score's interpretation.

## Activity

```math
A_\star(E)=\left[\prod_{m\in\mathcal M}x_{\star,m}(E)\right]^{1/|\mathcal M|}.
```

The panel is fixed for a run: ATAC, DNase, H3K27ac, ATAC+H3K27ac or
DNase+H3K27ac. Each signal is nonnegative, quantitative and matched to the same
scoring window and background. ATAC and DNase are not counted twice as independent
accessibility layers. A required missing layer is unresolved; the program does
not recompute a smaller panel for that element. A measured zero remains zero.

Each assay is resolved from qualified measured, predicted or calibrated fused
evidence before calculating activity. Bulk sequence predictions are averaged
across copies **per assay before** the geometric mean. Optional observed/predicted
fusion operates in a calibrated log1p signal space; the activity formula itself
has no log shift. RNA and additional omics do not multiply primary support.

## Contact and promoters

```math
\widetilde C(E,t)=r(E,t)C_{\mathrm{obs}}(E,t)
+[1-r(E,t)]C_{\mathrm{prior}}(E,t),\qquad 0\le r(E,t)\le1,
```

```math
\overline C(E,G)=\sum_{t\in\mathcal T(G)}\pi(t\mid G)\widetilde C(E,t),
\qquad \pi(t\mid G)\ge0,\qquad\sum_t\pi(t\mid G)=1.
```

Distinct physical TSSs are deduplicated before weighting. Weights are fixed for
the gene and independent of E. Missing contact for a positive-weight TSS makes
the gene contact unresolved; weights are not renormalized over available TSSs.
At r=0 or r=1 only the used source is required. A prior must match the declared
species, assembly, tissue, contact scale and target; none is silently substituted
from an unrelated species. Near-diagonal contact follows the explicit run policy.

## Cross-gene allocation

```math
B(E,G)=\frac{\overline C(E,G)}
{\displaystyle\sum_{H\in\mathcal G(E)}\overline C(E,H)}.
```

The candidate-gene set is frozen before scoring. Positive eta requires all
necessary contacts in that set; missing genes are not dropped to improve B.
At eta=0 the implementation skips B entirely, avoiding an undefined zero power.
This is a target-allocation factor, not a physical conservation law or a boundary score.

## Zeros, missingness and interpretation

| Condition | Result |
|---|---|
| Required activity/contact is missing | NA support with a reason |
| Required evidence is known and support is zero | Zero support, retained as scoreable |
| Activity is known and all necessary element contacts are known and zero | Zero support even though B is undefined |
| Gene support sum is zero or has no scoreable candidates | NA PACE score |
| Only part of the planned gene background is scoreable | Conditional score with `normalization_status=partial` |

Positive support is evaluated as `log(A) + log(Cbar) + eta*log(B)` and normalized
in log space. Missingness is never repaired by adding an arbitrary small constant.
The reported score is a relative support share, not a gene-expression fraction,
causal probability or calibrated false discovery rate. Compare samples on common
scoreable backgrounds and inspect activity, support and gene totals alongside PACE.

Implementation: [scoring kernel](../src/pace_livestock/core/scoring.py),
[calibration](../src/pace_livestock/learning/allocation.py).
Independent references: [worked examples](WORKED_EXAMPLES.md) and
[mathematical review](../PACE_Review_and_Validation.md).
