# Continuous allocation calibration

PACE accepts any finite `eta` in [0, 1]. The default configuration is `eta: auto`:
without applicable functional labels or a matching frozen artifact, the exponent
used numerically is **zero**. Explicit numeric eta values never trigger fitting.
The [primary formula](FORMULA.md) and the independent elastic-net classifier have
different roles; fitting eta does not train that classifier.

## Formula and objective

```math
\mathrm{PACE}(E,G)=\frac{A_\star(E)\overline C(E,G)B(E,G)^{\eta_{used}}}
{\sum_{e\in\mathcal E^{score}(G)}A_\star(e)\overline C(e,G)B(e,G)^{\eta_{used}}}.
```

In automatic mode, eta_used is zero unless the candidate estimate passes grouped
held-out selection and stability checks. The accepted value is bounded to [0,1]. At exactly zero, B and its missing-data requirements
are skipped. For every positive eta, the fixed candidate-gene universe is required
to define B. True-zero support remains zero; a zero normalization total is NA.

The estimator uses functional enhancer inhibition/deletion experiments with
direction-qualified positives and powered negatives **within each gene**. Define
`a = log(A*Cbar)` and `b = log(B)`. For a positive p and negative n for the same gene,
their log PACE ratio equals their log support ratio; the gene denominator cancels.
The declared ranking objective is

```math
\widehat\eta=\arg\min_{0\le\eta\le1}
\frac1{|\mathcal G_{fit}|}\sum_{G\in\mathcal G_{fit}}
\frac1{|P_G||N_G|}\sum_{p\in P_G,n\in N_G}
\log\left[1+\exp\{-[(a_p-a_n)+\eta(b_p-b_n)]\}\right].
```

Each informative gene has equal weight, and each positive-negative pair within a
gene has equal weight. The objective is convex. Endpoint derivative checks and
40 derivative-bisection iterations estimate the minimum on [0, 1]; this is not a
grid restricted to 0 and 1. Pair calculations use bounded memory blocks.

This is a ranking surrogate, **not** a Bernoulli likelihood of PACE or evidence
that its value is a causal probability. There is no assumed universal optimal eta
and no guarantee of improved held-out performance. Validate ranking and coverage
on independent functional experiments in the intended population and tissue.

## Required functional-label table

`--eta-labels` takes a tab-delimited file with these columns:

| Columns | Contract |
|---|---|
| `label_id` | Unique nonempty assay-label identity |
| `element_id`, `gene_id` | One frozen canonical element–gene edge; do not copy a multi-element perturbation onto individual tiles |
| `species`, `assembly`, `context_id`, `target_level` | Must match the scored scientific context |
| `perturbation_type` | `CRISPRi`, `deletion`, or `enhancer_inhibition`; `synthetic_inhibition` is demonstration-only |
| `label_status`, `effect_direction` | `enhancing_positive` with `down`, or `powered_negative` with `none` |
| `split` | train, calibration or test; train fits when present, calibration then confirms without refitting; calibration-only data require grouped validation |
| `group_id` | Predeclared independent locus/experimental group; no group may cross splits |
| `assay_id`, `source_id` | Nonempty experimental assay and source identifiers |

See the [synthetic table](../examples/training/eta_labels.tsv). This is a dedicated
calibration table, distinct from the region-level benchmark schema.
Convert a region label only when its canonical element mapping is unambiguous.
RNA abundance and association/QTL annotations do not constitute functional labels.
Negative status must reflect adequate assay power, not merely a nonsignificant test.

Gene IDs, element IDs and group IDs cannot cross train/calibration/test splits.
Duplicate labels, repeated usable edge labels, missing fields or leaking splits
raise errors. Context mismatch, unsupported perturbations, uncertain/repressive
labels and zero/unresolved support are excluded with recorded reasons. Test rows
are excluded before considering their functional outcomes. No test outcome affects
the fitted parameter. No pseudocount is added to make zero support fit a log model.

Every fitting gene needs at least one eligible positive, one eligible negative,
and an allocation contrast greater than numerical tolerance (1e-12 in log B).
Default `minimum_genes: 3` is an engineering guard, not a sample-size calculation.
Fewer informative genes cause eta=0 with a recorded fallback; each CV training
subset must also meet this requirement. The minimum may be predeclared as an integer
at least 2. These checks do not supply a biological sample-size calculation.

## Deployment validation after candidate fitting

The convex fit above proposes a candidate; it does not by itself authorize a
nonzero deployed exponent. Shared group_id, gene_id or element_id joins records
into connected independent components. At least `minimum_groups` (default 3)
components are required. Up to `validation_folds` (default 5) grouped folds refit
the continuous candidate using only each training fold.

The selection procedure tests 21 predeclared shrinkage multipliers from 0 to 1 in
steps of 0.05, multiplying each fold's separately fitted continuous candidate. It
uses held-out per-gene macro-average precision (AP), chooses the smallest multiplier
within one standard error of the best mean AP, then requires paired mean AP gain
over zero to exceed its standard error and positive gain in at least 80% of held-out
groups by default. This is conservative internal model selection, not an external
accuracy claim or a confidence interval for eta.

The candidate is then refitted on eligible fitting labels and multiplied by the
selected shrinkage. If explicit train rows exist, calibration rows are reserved
for a separate confirmation and cannot retune the candidate. If only calibration
rows are supplied, they still require the same grouped validation. Test outcomes
never enter fitting, selection or confirmation. Failed or uninformative validation
returns the zero default with reasons.

Artifacts distinguish `candidate_eta` from the actual `eta`, and record
`cross_validation`, `estimation_label_ids`, `confirmation_label_ids` and
`deployment_validated`. The current reusable schema is `pace-eta-2`; an old artifact
without deployment evidence must be recalibrated rather than relabelled. Explicit
manual exponents remain research ablations and are not automatically validated.

## Configuration and reuse

```yaml
allocation:
  eta: auto
  labels_path: functional_labels.tsv
  calibrator_path: null
  minimum_genes: 3
  minimum_groups: 3
  validation_folds: 5
  minimum_positive_fraction: 0.8
  missing_policy: fixed_gene_set
```

Alternatively set `calibrator_path` to a saved `eta_calibration.json` and remove
`labels_path`. The two sources are mutually exclusive. For predeclared ablation,
set `eta: 0`, `eta: 1` or another finite value in the interval and remove both paths.
CLI `--eta-labels`/`--eta-model` select automatic mode even when an older YAML
explicitly specified eta=0. Explicit `--eta NUMBER` selects manual mode.

Every run exports `eta_calibration.json`; the run manifest repeats its metadata and
stores the numeric exponent in `comparison_contract.eta`. Fitted artifacts record
the objective, fitted/boundary/fallback status, label and fitting-evidence hashes,
used genes/elements/groups, excluded-label reasons, baseline/fitted loss and scope.
Frozen reuse checks scope hashes, context, target level, regime, panel/scales,
catalog/candidate/promoter identities, contact policy and quantitative asset hashes.
Synthetic calibration cannot enter a research/validated run. An unchanged frozen
calibrator may be reused across compatible individuals, not across incompatible
tissues, candidate universes or measurement/model contracts.

For between-animal comparisons, fit once and reuse the same frozen artifact.
Refitting separately can confound genetic/evidence changes with parameter changes.
The benchmark command reports calibrated PACE in addition to its fixed baselines
when a calibration source is configured, and rejects evaluation labels that
overlap fitting genes, elements or groups. Calibrated in-sample scores are training
outputs; reporting them as independent performance is invalid.
