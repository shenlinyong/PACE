# Comparing PACE with activity–contact baselines

This page describes mathematical comparisons to the implemented baselines; it does
not claim higher biological accuracy or a complete reproduction of another package.
The current PACE equation is defined in [FORMULA.md](FORMULA.md).

## What the eta-zero limit establishes

For the same activity, contact and scoreable candidates, setting eta to zero gives

```math
\mathrm{PACE}_{\eta=0}(E,G)=
\frac{A_\star(E)\overline C(E,G)}
{\displaystyle\sum_{e\in\mathcal E^{\mathrm{score}}(G)}A_\star(e)\overline C(e,G)}.
```

With one TSS, this is an activity–contact normalization rule. It does not imply that
candidate generation, measurement processing, contact priors or missingness are
identical to an external ABC implementation.

| Component | Current PACE behavior | Requirement for a fair comparison |
|---|---|---|
| Activity | Fixed-panel equal geometric mean | Match assays, units, windows and missingness |
| Contact | Declared observations/prior policy and weighted distinct TSSs | Match contact scale and evidence context |
| Allocation | B raised to a fixed or training/calibration-estimated eta in [0,1] | Freeze the exponent before final evaluation |
| Candidates | Explicit canonical catalog and frozen candidate edges | Align the candidate universe and disclose exclusions |
| Normalization | Sum of actually scoreable support within each gene | Recompute a common denominator rather than joining normalized outputs |
| Other omics | Annotation or a separate supervised classifier | Evaluate classifier and formula outputs separately |

## Built-in benchmark methods

`PACE benchmark` evaluates eta-zero and eta-one controls, negative distance and an
explicit `ABC_style_single_TSS` control. That control selects the physical TSS with
smallest coordinate, then promoter ID, before evaluation. It reuses resolved
contact evidence; it is not a claimed reproduction of an external software pipeline.
Configured calibrated PACE or fixed continuous eta is also reported. Calibration
and evaluation genes/elements/groups must not overlap.

External score tables require method version and configuration provenance. They
can be evaluated on their available coverage; a common-denominator comparison is
not claimed without the underlying support needed to reconstruct it. Use eligible
functional positives and powered negatives, report average precision and coverage,
and choose any decision threshold using training/calibration data. No generic
numeric cutoff establishes a validated regulatory link.

See [benchmark configuration](comparison.md), [worked calculations](WORKED_EXAMPLES.md)
and [eta calibration](eta_calibration.md).
