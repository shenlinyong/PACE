# Validation under limited experimental data

This is a study protocol, not a report of completed biological experiments. No performance increase is asserted by the software release.

## Primary question

Does contact regularization and explicit handling of missing support preserve useful regulatory-link ranking and coverage as experimental data become sparse? This is a measurable claim. Human downsampling experiments provide a controlled stress test; they cannot by themselves establish transfer accuracy in chicken, pig or cattle.

## Functional reference and split

Use genuinely tested K562 enhancer–gene pairs from the [EngreitzLab CRISPR benchmark](https://github.com/EngreitzLab/CRISPR_comparison), recording release, genome assembly, assay and effect/power definitions. Fulco/Gasperini datasets overlap across compilations: deduplicate perturbation regions, genes and assays before splitting. Never count a reused benchmark as an independent replication. Negatives should be tested, sufficiently powered non-effects; untested pairs and low-power experiments are not verified negatives.

Freeze chromosome or connected region/gene groups into training, calibration and final test sets. All tuning, allocation selection, promoter weighting choices and contact preprocessing decisions must precede final test evaluation. Keep repeated guide/region/gene entities in one split. Record both the original benchmark's screening selection and exclusions caused by candidate mapping.

## Controlled degradation

Use a full-quality K562 reference and predeclare contact fractions such as 100%, 50%, 25%, 10% and 5%, with several random seeds. These are experimental design choices, not PACE defaults.

1. Thin raw read pairs or unique raw Hi-C pixel counts using binomial sampling, preserving symmetry and counting each pair once. Do not thin a balanced matrix or separately resample duplicated enhancer–TSS queries from the same bin pair.
2. Rebuild/balance the sparse map at each depth and refit its prior on training chromosomes. A full-depth fitted prior reused at low depth is an extra-information condition; label it separately.
3. Remove a whole activity assay and explicitly choose the remaining fixed single-assay panel. Separately simulate missing regions while retaining the fixed panel, so assay ablation is not confused with selective missingness.
4. Reduce independent donors/biological replicates. Technical replicates do not increase independent animal count. Keep population summaries distinct from individual estimates.
5. Freeze the candidate universe for the primary comparison. Separately examine candidate-construction sensitivity, including peak-centered windows and fixed-grid offsets. Report mapping/coverage differences rather than selecting a favorable catalog post hoc.

## Baselines and ablations

Run the pinned original ABC implementation with its own documented contact processing; export its scores as an external benchmark method with version/configuration provenance. PACE's `ABC_style_single_TSS` is a mathematical control using PACE contacts, not a substitute for running original ABC.

Compare: PACE eta=0; no contact pseudocount; no near-diagonal correction; prior-only; explicit missing-contact fallback; single versus multiple TSSs; equal versus measured promoter-signal weights; B allocation versus a pure contact-power control. `pace benchmark` includes `contact_power_2` (A×Cbar²) to compare with eta=1 allocation. If tuning contact power or eta, tune each with the same independent validation budget. Report annotation and biotype sensitivity for any B-related claim.

## Metrics and uncertainty

Report average precision (AP, with its non-trapezoidal definition), precision–recall curves, recall at an independently selected precision target, and candidate/positive-label coverage. Because default PACE withholds incomplete primary scores, a higher AP among fewer available predictions is not sufficient. Report AP on a prespecified common evaluable set **and** abstention/coverage on all tested pairs. Show how many positives remain unscored and how coverage changes with depth.

Use paired bootstrap over independent gene/region groups, with the same sampled groups for every method. Chromosome holdout and the number of independent genes constrain precision; random edge bootstrap is not an independent biological replication. Report per-distance, promoter-count, gene-density and chromosome strata to expose confounding. Do not select a seed, cutoff or degradation level after inspecting final test performance.

For score intervals, use complete high-quality supports to check containment after deliberately masking candidates. Evaluate interval width and any assumptions used for finite support bounds. The default [0,infinity) missing-support bounds can be very wide: this is honest lack of information. These mathematical ranges are not 95% confidence intervals and are not automatically a novel method simply because they are implemented here.

## Livestock evidence

Use species- and tissue-matched FarmGTEx/eQTL resources as orthogonal association evidence. Match or stratify control pairs by enhancer–gene distance, allele frequency, LD, gene expression, variant and candidate density, and assay accessibility. Prefer fine-mapped/colocalized signals where available; do not treat each correlated SNP as an independent regulatory event. eQTL association is not a direct enhancer perturbation label.

Use available independent chicken, pig and cattle perturbations or targeted experiments as final functional evidence. Keep the target of each conclusion precise: human sparse-data robustness, livestock association enrichment, or livestock functional-link accuracy. Public gene annotations and a more elaborate QC system alone do not establish improved accuracy in livestock.
