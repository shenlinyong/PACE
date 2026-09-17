> **Interface scope:** This retained guide documents the legacy region-based scripts/workflow.
> For the installable canonical-grid package and its September 2026 defaults, see [software.md](software.md).

# PACE compared with the original ABC model

PACE is a general enhancer–gene scoring framework that is particularly well suited to the incomplete data conditions common in livestock studies. Its changes address measurement availability, contact reliability, promoter representation and candidate coverage. None depends on a livestock-specific gene list or species label.

This page compares the original ABC formulation with the implemented PACE model. Numerical examples illustrate the calculations; they are not biological performance estimates.

## Reference and equations

The comparator is [Fulco et al., Nature Genetics (2019)](https://doi.org/10.1038/s41588-019-0538-0), using the [NG2019 code](https://github.com/EngreitzLab/ABC-Enhancer-Gene-Prediction-20250314-archive/tree/NG2019). The [ABC documentation](https://abc-enhancer-gene-prediction.readthedocs.io/en/latest/) distinguishes this implementation from later releases.

$$
\mathrm{ABC}(E,G)=\frac{A_{\mathrm{ABC}}(E)C_{\mathrm{ABC}}(E,G)}{\sum_{e\in\mathcal E(G)}A_{\mathrm{ABC}}(e)C_{\mathrm{ABC}}(e,G)}.
$$

$$
\mathit{PACE}(E,G)=\frac{A(E)\overline C(E,G)B(E,G)^\eta}{\sum_{e\in\mathcal E^{\mathrm{obs}}(G)}A(e)\overline C(e,G)B(e,G)^\eta+U(G)}.
$$

Here $\overline C(E,G)$ denotes reliability-adjusted, TSS-averaged contact. It is **the same quantity** called $C(E,G)$ in the [scoring specification](FORMULA.md) and `contact_gene` in software output. The bar distinguishes it from the ABC contact estimate on this comparison page; it introduces no additional operation.

ABC offers a direct activity–contact formulation for reliably estimated inputs. PACE adds explicit rules for inputs with uneven availability or reliability and for incomplete promoter and candidate catalogues. This makes PACE a useful fit for many livestock datasets, while remaining applicable across species.

## Six changes and why they matter

| Original ABC treatment | PACE treatment | Reason for the change | Relevance to livestock data |
| --- | --- | --- | --- |
| Use the chosen contact estimate in activity × contact | Reliability-adjusted $\overline C(E,G)$ | Let qualified observations influence the score in proportion to independently assessed reliability | Matched Hi-C may be absent or vary in quality across tissues and conditions |
| Original activity profile combines accessibility and H3K27ac geometrically | Aggregate effective measurements in $A(E)$, with explicit `NA` and measured `0` | Keep available signal when another assay is missing, without declaring absence of activity | Multi-assay coverage often differs among species, breeds and tissues |
| Normalize separately for each gene | Add enhancer-centred allocation $B(E,G)$ | Reduce raw support when one enhancer distributes contact across several candidate genes | Provides a specificity mechanism for ambiguous links in gene-dense regions |
| Use contact to the selected promoter | Deduplicate TSSs and average contacts using $\pi(G,t)$ | Represent alternative promoters without a transcript-count multiplier | Annotation sources may contain repeated transcripts or different promoter alternatives |
| Report a relative score | Report score alongside $Q$ and evidence state | Separate rank from the adequacy of its supporting measurements | Incomplete evidence can be retained as `provisional` for further study |
| Normalize candidate activity–contact products | Define $\mathcal E^{\mathrm{obs}}(G)$, retain unscored rows and accept independent $U(G)$ | Preserve the distinction between zero support, uncomputable support and support outside the supplied catalogue | Missing assays and incomplete candidate discovery remain visible |

### 1. Contact reliability: how much should a measurement alter the prior?

Original ABC considered tissue-specific Hi-C, average Hi-C and distance-based contact alternatives. PACE's additional mechanism is explicit reliability adjustment at the enhancer–TSS level:

$$
C_{\mathrm{adj}}(E,t)=P(d)\left[1-\lambda(E,t)+\lambda(E,t)\frac{H(E,t)}{D(d)}\right].
$$

| Quantity | Role | How it enters the implementation |
| --- | --- | --- |
| $P(d)>0$ | Contact expected from genomic distance before using pair-specific evidence | `contact_prior` |
| $H(E,t)\geq0$ | Observed contact, including a genuine zero | `contact_observed` |
| $D(d)>0$ | Expected contact in compatible measurement units | `contact_expected` |
| $\lambda(E,t)\in[0,1]$ | Independently assessed reliability | `contact_reliability` |
| Source | Matched tissue, declared surrogate, distance prior or unknown | `contact_source` |

For $P=1$ and $H/D=4$, adjusted contact is 1 at $\lambda=0$, 1.75 at $\lambda=0.25$, and 4 at $\lambda=1$. Weakly supported observations therefore change the prior less. A missing observation uses the prior. A measured zero with $\lambda=1$ yields zero contact.

This is useful when a livestock study has a limited or borrowed contact map. Missing expectations, missing reliability or undeclared provenance trigger prior fallback and an explanatory state. PACE does not infer reliability from contact strength or correct tissue mismatch automatically. Expected contact must match the sample, normalization, resolution and distance of the observation.

The [original ABC predictor](https://github.com/EngreitzLab/ABC-Enhancer-Gene-Prediction-20250314-archive/blob/NG2019/src/predictor.py) already scales contact and adds a distance-derived pseudocount. Its $C_{\mathrm{ABC}}$ should therefore be understood as a processed contact estimate. The PACE distinction is the explicit observation/expectation ratio, local reliability and provenance rule.

### 2. Activity: which measurements were actually observed?

The original two-assay ABC profile uses an unshifted geometric mean. PACE uses scaled signals $x(E,i)=S(E,i)/a_i$ and, in the primary configuration, computes

$$
A(E)=\exp\left[\frac{\sum_i m(E,i)q(E,i)w_i\log(1+x(E,i))}{\sum_i m(E,i)q(E,i)w_i}\right]-1.
$$

The mask $m$ identifies available measurements; $w_i$ gives assay importance; $q$ supplies local quality; $a_i$ fixes the signal scale. The [parameter reference](PARAMETERS.md#1-activity-integration) explains their separate roles.

With equal weights, known unit quality and unit scales:

| Accessibility | H3K27ac | PACE activity | Interpretation |
| --- | --- | --- | --- |
| 8 | `NA` | 8 | Use the available measurement; the planned missing assay reduces evidence coverage |
| 8 | 0 | 2 | Both assays were measured; zero participates in the shifted geometric mean |
| 0 | 0 | 0 | All effective measurements are observed zero |
| `NA` | `NA` | `NA` | Activity cannot be calculated |

This supports studies with incomplete multi-omics coverage without silently imputing absent assays. Unknown assay quality uses unit weight for provisional activity and stays unknown in evidence reporting. The shift also changes behavior for observed low signals: the original geometric mean of 8 and 0 is zero. Its effect on sensitivity and false positives needs evaluation.

Single-assay inputs and distance contact are not exclusive PACE capabilities; the comparison concerns the full missing-value and quality treatment. Removing an assay from the planned JSON changes the input regime; leaving it declared with missing values records incomplete coverage.

### 3. Target allocation: how broadly is an enhancer's contact distributed?

$$
B(E,G)=\frac{\overline C(E,G)}{\sum_{g\in\mathcal G(E)}\overline C(E,g)},
\qquad R(E,G)=A(E)\overline C(E,G)B(E,G)^\eta.
$$

For one enhancer with contacts 3 and 1 to two genes, the allocation factors are 0.75 and 0.25. At $A=1$ and the default $\eta=1$, raw support changes from 3 and 1 to 2.25 and 0.25. Gene-centred normalization is then applied to these adjusted supports.

Allocation reduces raw support for broadly shared enhancers and can help distinguish targets in dense genomic regions. Final scores depend on all candidates for each gene: allocation does not guarantee a lower normalized score for every link. If each gene has only one scored candidate and $U=0$, both gene-normalized scores can still be 1.

The idea is attributed to [generalized ABC / STARE](https://doi.org/10.1093/bioinformatics/btad062). PACE exposes $\eta$ through `--competition-power`: 1 uses allocation and 0 removes it for an ablation. Incomplete gene annotation changes the allocation denominator and must be considered when comparing catalogues.

### 4. Promoters: how should multiple transcript records contribute?

$$
\overline C(E,G)=\sum_{t\in\mathcal T(G)}\pi(G,t)C_{\mathrm{adj}}(E,t),
\qquad\sum_{t\in\mathcal T(G)}\pi(G,t)=1.
$$

Suppose a gene has two distinct TSSs with adjusted contacts 2 and 6. Uniform weights give gene contact 4. Repeating a transcript record at either TSS leaves that result unchanged. Independent use weights 0.75 and 0.25 instead give contact 3.

This retains annotated promoter alternatives while avoiding support inflation from duplicate transcript records. Set weights over the full available TSS catalogue before applying candidate windows; do not renormalize them for each enhancer. Without supplied weights, the table interface can only infer uniform weights from the TSSs visible in its input.

The multiple-TSS principle also draws on generalized ABC. PACE implements a weighted average. Missing promoters remain missing, and gene-boundary fallback is labelled by `prepare_tss.py`.

### 5. Evidence: does a high score have adequate input support?

$$
Q(E,G)=\min\{Q_A(E),Q_C(E,G),Q_T(G),Q_{\mathrm{cat}}(G)\}.
$$

The four components describe activity, contact, TSS and catalogue evidence. Unknown components remain unknown. The resulting $Q$ is reported separately, with no additional $Q$ multiplier in the structural score. Assay quality and contact reliability still affect their respective activity and contact calculations.

At the default threshold 0.5, `sufficient_input_evidence` requires all four components to meet the threshold and no supplied gene candidates to remain unscored. Undefined scores or zero activity/TSS quality are `insufficient`; other cases are `provisional`. Prior-only predictions have contact quality zero and ordinarily remain provisional.

The distinction matters when normalized scores look strong despite sparse inputs. For example, one observed candidate can score 1 even when another candidate's activity is missing. Its evidence state exposes that limitation. $Q$ is an input-quality index, not a calibrated probability of functional regulation.

### 6. Unscored candidates and residual support: what is known about the denominator?

$\mathcal E^{\mathrm{obs}}(G)$ contains supplied candidates with finite raw support. Candidates outside that set remain in the output with missing scores and contribute to `unscored_candidates`.

For raw supports 2, 1 and `NA`, and unknown $U$, the finite scores are 2/3 and 1/3. The third score remains `NA`; the output reports one unscored candidate. If independent evidence justifies $U=1$ in the same raw-support units, the finite scores become 1/2 and 1/4. The unscored row remains missing.

Unknown $U$ is computationally zero, with `score_scope=observed_candidates_only` and reason `unassigned_mass_unknown`. Explicit zero declares a zero-residual assumption. PACE provides no automatic estimator of omitted support and does not turn a missing-peak fraction into $U$.

Two forms of incompleteness therefore remain distinct: unscorable rows in the supplied table are counted; undiscovered candidates outside that table require independent coverage information. Neither is treated as a measured absence of regulation.

The [NG2019 output code](https://github.com/EngreitzLab/ABC-Enhancer-Gene-Prediction-20250314-archive/blob/NG2019/src/predict.py) also writes missing values and failed-gene records. PACE adds the explicit finite-support denominator, per-gene unscored-candidate count, residual-support input and linked evidence state described above.

## Which parameters changed, and which are inherited starting settings?

| Category | Settings | Reason |
| --- | --- | --- |
| Activity extension | $m,q,w,a$ and fixed `log1p` shift | Represent observed assays, their quality and their scales explicitly |
| Contact extension | $H,D,\lambda$ and source state | Modulate the distance prior with qualified observations |
| Promoter/target extension | $\pi$ and $\eta$ | Integrate TSS use and enhancer-side allocation |
| Evidence/coverage extension | Four QC components, evidence threshold, $\mathcal E^{\mathrm{obs}}$ and $U$ | Keep ranking, evidence sufficiency and missing support distinguishable |
| Retained candidate choices | 5-Mb cis window, approximately 500-bp candidates, peak cap | Define the search background; these are not livestock-specific innovations |
| Configurable starting values | Distance decay/scale/offset and score cutoff 0.02 | Require a declared input profile and independent calibration when used for decisions |
| Disabled extension | Inhibitory strength $\kappa=0$ | Avoid requiring methylation or repressive assays for primary predictions |

[Parameters](PARAMETERS.md) lists exact values, units, input columns, CLI options and active YAML entries. RNA remains an expression-context annotation and optional downstream filter; it does not multiply the primary PACE score.

### Numerical defaults that require care in a comparison

The original command-line values come from [NG2019 `predict.py`](https://github.com/EngreitzLab/ABC-Enhancer-Gene-Prediction-20250314-archive/blob/NG2019/src/predict.py); its contact operations are defined in `predictor.py` above.

| Setting | Original ABC implementation | Current PACE file adapters | Reason for the PACE choice |
| --- | --- | --- | --- |
| Cis window | `--window 5000000` | `--max_distance 5000000` | Retain a broad cis search background |
| Distance exponent | `--hic_gamma 1`, reference exponent 1 | `--hic_gamma 1.024238616787792` | Retain the current adapter's starting decay profile; this is not a fitted livestock optimum |
| Near-distance treatment | `hic_pseudocount_distance=1000000` controls an additive contact pseudocount | A 5,000-bp offset enters $P(d)=\mathrm{Scale}(d+5000)^{-\gamma}$ | Keep the prior finite at zero distance; the two parameters have different meanings |
| Contact scale | Row-maximum rescaling to 100, pseudocount addition and clipping at 100 | Prior scale 5.9594510043736655, modulated by qualified $H/D$ | Define explicit prior/support units; a common scale cancels from normalized scores when $U=0$ |
| Selection threshold | Required CLI argument; the original README discusses 0.02 and shows a 0.022 command example | 0.02 in file predictors and filters | Provide an illustrative output selection; calibrate each model on independent validation |

PACE's retained YAML key `hic_pseudocount_distance` is not forwarded by the current file-adapter prediction path. It does not configure the 5,000-bp prior offset. Use a custom `contact_prior` in the table interface when another prior is required.

## Applying the model to livestock and other species

Species suitability comes from the data conditions, not a hard-coded species preset. Use the appropriate assembly, chromosome sizes, stable gene IDs and promoter catalogue. Apply consistent candidate and annotation policies when comparing tissues or breeds. For chicken, review microchromosome/scaffold inclusion and the effective genome size used for peak calling rather than reusing a mammalian placeholder.

Pig, cattle and chicken motivate practical use cases for these extensions. They do not establish universal superiority over ABC. The [validation record](../VALIDATION.md) describes the biological evidence and its scope.

## Designing a direct ABC comparison

Identify the ABC branch and commit. Match assembly, available data, candidate background and eligibility rules wherever the declared comparison permits. Report each method's preprocessing and input profile. Select thresholds using independent validation, then retain omitted predictions in coverage and recall calculations.

Separate a full native-pipeline comparison from a component ablation. With one TSS, $\eta=0$, $U=0$ and exactly the same $A$ and $C$, PACE recovers the ABC normalization rule. It does not reproduce original ABC preprocessing, contact handling or candidate construction. The [tutorial](TUTORIAL.md#compare-predictions-and-annotate-genetic-support) provides the evaluation command and label schema.
