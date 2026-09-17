# Comparison with the original ABC model

PACE (Prediction of Activity-based regulatory Connections for Enhancers) retains the activity–contact interpretation of enhancer–gene scoring and changes how incomplete measurements, alternative promoters and competing targets enter that calculation.

## Reference definition

The comparator here is Fulco et al. (2019), [Activity-by-contact model of enhancer–promoter regulation from thousands of CRISPR perturbations](https://doi.org/10.1038/s41588-019-0538-0), with the [NG2019 code](https://github.com/EngreitzLab/ABC-Enhancer-Gene-Prediction-20250314-archive/tree/NG2019). The ABC maintainers [identify that branch separately from later releases](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction). Pin a commit when running a comparison.

The original score is

$$
\mathrm{ABC}(E,G)=\frac{A_{\mathrm{ABC}}(E)C_{\mathrm{ABC}}(E,G)}{\sum_{e\in\mathcal E(G)}A_{\mathrm{ABC}}(e)C_{\mathrm{ABC}}(e,G)}.
$$

The original profile combined accessibility and H3K27ac through a geometric mean. It already offered measured, averaged and power-law contact alternatives. Therefore, use without tissue-matched Hi-C is **not by itself a new PACE capability**. This manual compares defined computations rather than treating every later ABC feature as absent.

## Design changes and their consequences

| Change | What PACE computes | Why introduce it? | Relevance to livestock data | Boundary of the claim |
| --- | --- | --- | --- | --- |
| Missing-aware activity | Weighted mean in `log1p` space, then `expm1`; `NA` omitted and measured zero retained | Avoid discarding an available assay or silently replacing an absent assay by zero | Different tissues may have different assays or patchy coverage | It cannot recover an unmeasured enhancer; shifted activity can increase false positives |
| Explicit assay scales and priors | Divide each signal by a fixed positive scale; default equal priors | Keep numerical units and chosen assay importance explicit | Enables declared single-assay and two-assay input regimes | Equal weights are design choices, not livestock-trained optima |
| Reliability-aware contact | Distance prior multiplied by a mixture of 1 and observed/expected contact | Control how much qualified contact evidence modifies the prior | Allows sparse measured contact or a declared surrogate map to be used transparently | Expected contact and reliability must be supplied independently; no automatic tissue correction |
| Distinct promoter integration | Weighted average across deduplicated TSSs | Represent multiple promoters without making transcript count itself a multiplier | Alternative and incompletely annotated promoters matter when combining animal annotations | Missing promoters remain missing; uniform TSS use is an assumption |
| Enhancer-centred target allocation | `target_share = C(E,G) / sum_g C(E,g)` | Account for several candidate genes sharing one enhancer | Makes candidate target context part of the model | Missing gene annotation also changes the allocation denominator |
| Residual-support input | Add independently estimated `unassigned_mass` to gene support | Allow a declared estimate of support outside observed candidates | Makes incomplete candidate coverage explicit | Default is unknown, computationally zero; the software does not estimate missing regulatory mass |
| Separate evidence assessment | Minimum of four input-quality components plus reason codes | Avoid treating a large normalized score as evidence that all inputs are adequate | Helps prioritize follow-up when assays and annotation have uneven support | This index is not a confidence probability or an FDR estimate |
| RNA as context | Stable-ID expression annotation and optional downstream filtering | Preserve the structural score independently of target expression | Retains low-expression and missing-expression cases for inspection | The primary score does not learn cross-tissue expression covariance |

The target-allocation and multiple-TSS ideas are based on Hecker et al. (2023), [The adapted Activity-By-Contact model for enhancer–gene assignment and its application to single-cell data](https://doi.org/10.1093/bioinformatics/btad062), and its [STARE implementation](https://github.com/schulzlab/STARE). PACE's TSS contact is a weighted **average**, with weights assigned before the candidate window is applied. The attribution concerns the underlying ideas; the PACE formula and implementation should be evaluated as their own model.

## Settings retained as starting points

The 5-Mb cis window, approximately 500-bp summit-centred candidates and activity–contact gene normalization retain the ABC-style experimental setup. PACE does not claim that those values were optimized for pig, cattle or chicken. Likewise, the shipped distance exponent, distance offset and score cutoff are starting settings, not fitted species-specific constants. Their exact values, control points and reasons appear in [Parameters](PARAMETERS.md).

## What makes the workflow suitable for domestic-animal studies?

Suitability here means that the inputs and outputs accommodate common data constraints. Supply a species-matched genome and annotation, use the measured assays you actually have, retain all distinct annotated TSSs, and carry contact provenance and unknown quality into the result. No human gene-symbol list or pretrained human neural-network weights are required by the primary scoring kernel.

For pig and cattle, use genome sizes, chromosome naming and annotation versions appropriate to the chosen assembly. For chicken, do not reuse a mammalian effective genome size in MACS2, and keep the selected microchromosome/scaffold policy consistent across references and signals. No liftover or cross-species transfer is performed automatically. These are data-preparation requirements, not hard-coded species presets.

The associated project includes livestock maps, human perturbation comparisons and genetic-support analyses; [Validation](../VALIDATION.md) states their scope. The software tests demonstrate correct execution and specified numerical behavior. They do not establish that PACE universally outperforms ABC, or that an input regime works equally well in every tissue.

## A fair ABC comparison

Use the same assembly, accessible data budget, candidate background and independently determined eligibility rules. Run each declared native or adapted input profile transparently. Freeze thresholds on a separate validation set. Retain missing predictions when calculating positive coverage and recall; do not restrict evaluation to pairs both methods scored.

The table CLI accepts `--competition-power 0` as an allocation ablation. With one TSS, zero residual support and exactly the same activity/contact inputs, this recovers ABC-style score normalization. It **does not** recreate the full original ABC pipeline: activity preprocessing, candidate construction and contact processing can still differ.

For a two-method evaluation command and the required experimental-label schema, see [the tutorial](TUTORIAL.md#compare-predictions-and-annotate-genetic-support).
