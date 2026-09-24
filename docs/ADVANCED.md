# Advanced use

The [README](../README.md) covers installation, input preparation and ordinary
scoring. This guide covers optional evidence, calibration, independent validation
and every setting with its default. Every command takes command-line options
(`pace COMMAND --help`); no configuration file is needed. Examples with biological filenames are templates:
replace their paths and context before running them. Bundled training data are
synthetic.

- [Evidence preparation](#evidence-preparation)
- [Sparse activity and alternative TSSs](#sparse-activity-and-alternative-tsss)
- [Contact policies and priors](#contact-policies-and-priors)
- [Boundary priors and weak calibration](#boundary-priors-and-weak-calibration)
- [Incomplete scores and support bounds](#incomplete-scores-and-support-bounds)
- [Additional omics](#additional-omics)
- [Experimental target allocation](#experimental-target-allocation)
- [Classifier training and inference](#classifier-training-and-inference)
- [Comparisons and benchmarks](#comparisons-and-benchmarks)
- [Worked calculations](#worked-calculations)
- [Validation](#validation)
- [Limitations](#limitations)
- [Complete configuration defaults](#complete-configuration-defaults)

## Evidence preparation

Keep one sample identifier per experimental replicate and register its donor,
biological replicate, technical replicate and source in `samples.tsv`. Assays from
different animals do not become individual-level measurements merely because their
tissue names match. Use `--target-level population_mean` when pooling donors.

Create empty input-table templates with:

```bash
pace init --species chicken --assembly GRCg7w --tissue liver \
  --panel H3K27ac --out chicken_liver
```

Fill the tables before validation. `--catalog-dir prepared/catalog` imports the
catalog tables and settings exported by `pace catalog`.

Catalog preparation maps supplied peaks and required promoters onto unique grid
cells; it does not tile the whole genome. GTF input needs transcript records, with
positive-strand TSS=start−1 and negative-strand TSS=end−1. Physical TSSs are
deduplicated, and versioned gene/transcript IDs are preserved. Convert annotations
without transcript records to `promoters.tsv` first. Chromosome aliases must be
explicit (`alias`, `canonical` columns); no prefix guessing is performed.
Incomplete terminal grid cells are reported and excluded. Keep the exported
`run_catalog_config.yaml` settings so reference-boundary exclusions remain
verifiable. Requantify the resulting cells: copying a whole peak's read count to
every overlapping cell inflates support.

bigWig input must contain nonnegative quantitative signals. Unstored bases remain
missing unless the source explicitly represents them as measured zero and
`--missing-as-zero` is given. The callable fraction measures coverage, not
activity. Normalization identifiers describe the actual processing; changing a
name does not make different protocols comparable.

Contact input must retain the distance background. Convert O/E to contact using its
matched distance expectation, and undo any log transform first. P values and
correlation coefficients are not contact values. Convert `.hic` upstream to
cool/mcool and record the converter and settings. Invalid balancing weights remain
unavailable; unstored pixels count as zero only when requested and both bins are
valid.

### Raw activity counts

If your window table contains **raw filtered fragment counts**, prepare these TSV headers:

```text
# counts.tsv
element_id sample_id assay count
# library_sizes.tsv
sample_id library_size
```

Use actual tabs, not the spaces shown above. One count per element/sample/assay; one full-library total per sample. Missing counts are NA; measured zero counts are 0.

```bash
pace normalize-activity --counts counts.tsv --library-sizes library_sizes.tsv \
  --units prepared/catalog/units.tsv --out prepared/normalized
```

Output signal = count × 10^6 / library_size / window_length, in CPM_per_bp. The library total counts all retained fragments, not only peak-overlapping fragments. Paired-end data must use the same fragment-counting convention in numerator and denominator. Variable-width regions require `--window-id`. The default identifier is `CPM_density_v1`; a project can supply a more specific `--normalization-id`.

This command does not infer total library size from a bigWig, renormalize TPM, or normalize fold-enrichment tracks. Use the bigWig adapter for already normalized quantitative tracks and name their actual processing protocol. Changing a protocol string does not make incompatible assays comparable. No default cross-animal quantile normalization is imposed: it could remove real global biological differences.

### Measured promoter weights

Prepare `promoter_signal.tsv` with `promoter_id` and a nonnegative `signal` measured over the same specified promoter window. Supply one already aggregated value per physical promoter.

```bash
pace prepare-promoter-weights --promoters prepared/catalog/promoters.tsv \
  --signals promoter_signal.tsv --assay H3K4me3 \
  --normalization-id H3K4me3_CPM_TSS500 --out prepared/promoter_weights
```

Use the resulting promoters.tsv in the run. Weights are fixed before enhancer scoring. An unavailable TSS signal errors. All-zero genes error unless `--zero-policy equal` is explicitly chosen. ATAC, DNase and CAGE are also supported; these measurements are promoter-activity proxies and require their own biological validation.

### Imported evidence and extra output fields

The README lists ordinary input tables. The following tables support imports and
advanced analyses; headers below use tabs.

```text
# resolved_activity.tsv
element_id	assay	observed_value	resolved_value	evidence_id	evidence_type	observation_sample_id	parent_evidence_ids	resolution_status	reason	unit	normalization_id	window_id
# resolved_contacts.tsv
element_id	promoter_id	resolved_value	evidence_id	evidence_type	observation_sample_id	prior_id	reliability	resolved_mode	bin_pair_id	resolution_status	reason	scale
# evidence.tsv
evidence_id	evidence_type	source_id	parent_evidence_ids	model_id	unit	processing_method	checksum
```

Activity imports accept only observed or experimentally aggregated evidence and
must identify the correct assay/sample, unit, normalization and measurement
window. Contact evidence can be `observed`, `aggregate`, `contact_prior`,
`regularized` or `fused`; the last refers only to explicit observation/prior
shrinkage. Imported contacts require a resolution from the table or compatible
run inputs. Scale alone is insufficient; normalization, balancing and window
identifiers must agree when specified. `resolution_status` is `resolved`,
`unresolved` or `invalid`.

`resolved_activity.tsv` preserves measured values and records any
`activity_pseudocount` separately; reimporting it does not apply an offset twice.
`promoter_weights.tsv` records `pi_original`, effective `pi`,
`tss_selection_reason` and `n_missing_candidate_contacts`. Score and gene-summary
rows include `tss_retained_weight`, `tss_dropped_ids`, `tss_policy_status` and
`tss_contact_scope`. `n_tss_used` counts positive effective weights; `n_tss` counts
all original TSSs.

`region_scores.tsv` sums each source/region/gene's unique member cells, without a
new denominator. Do not add overlapping region scores as independent units.
`evidence.tsv` and `sources.tsv` preserve input and derived provenance.
`ml_feature_contract.json` records measurement definitions for the optional
classifier; `pace_ml_score` and `pace_ml_probability` remain separate from the
primary score and stay NA when inference is unavailable.

## Sparse activity and alternative TSSs

These options are off by default. Add only the options needed for the analysis:

```bash
pace run ... --panel ATAC H3K27ac \
  --activity-pseudocount H3K27ac 0.1 \
  --tss-weights provided --minimum-tss-weight 0.01 \
  --missing-tss-policy drop_missing --minimum-retained-tss-weight 0.9
```

The 0.1 offset is an example only; use the units of your normalized signal.

The pseudocount is added after replicate aggregation; NA remains NA. Compare
results with zero offsets because positive offsets can increase background support.
There is no universally appropriate 0.1 offset across assays or normalizations.

TSS filtering is performed once per gene. Every candidate for that gene uses the
same retained set; a TSS missing any planned contact is removed globally under
`drop_missing`. At least 90% of the original weight must remain in this example.
The original candidate catalog is preserved. Inspect `promoter_weights.tsv`,
`tss_contact_scope` and `tss_retained_weight` before interpreting scores.
Different selected promoter sets do not define the same biological comparison.

```math
A_\star(E)=\left[\prod_{m\in\mathcal M}x_m(E)\right]^{1/|\mathcal M|}.
```

Supported panels are ATAC, DNase, H3K27ac, ATAC+H3K27ac and DNase+H3K27ac. The panel is fixed across elements. A missing required assay gives NA; a measured zero gives zero activity by default. Technical replicates are averaged within biological replicates, biological replicates within donors, and donors equally. Signals must be normalized comparably before averaging. `individual` requires one donor; `population_mean` pools donors explicitly.

Multiplying every element's signal for assay m by the same positive constant multiplies all activities by a common factor. It cancels in a gene's normalized score. This does not remove differences in signal-to-noise, measurement windows, missingness or replicate-specific scaling.

For shallow measurements, users may set an assay-specific nonnegative pseudocount:

```math
A_{\star,\epsilon}(E)=
\left[\prod_{m\in\mathcal M}(x_m(E)+\epsilon_m)\right]^{1/|\mathcal M|}.
```

`activity.pseudocounts` defaults to `{}` (all epsilon values zero). Pseudocounts are in the units of the normalized assay signal and are applied after replicate aggregation. Raw resolved measurements are preserved; `activity_pseudocount` records the addition. NA remains NA. Rescaling an assay preserves scores only if its pseudocount is rescaled too. A fixed positive offset can increase background support; choose it independently of test outcomes and report sensitivity to zero.

`pace normalize-activity` converts raw counts to count × 10^6 / filtered-library-size / window-length. Already normalized bigWigs must not be normalized twice. Library-size normalization is not batch correction.

### TSS selection

```math
\overline C(E,G)=\sum_{t\in\mathcal T(G)}\pi(t\mid G)\widetilde C(E,t).
```

Identical physical TSSs are deduplicated. Distinct TSSs in the same measured bin may reuse the same element/sample contact query without adding replicates. The distance prior may differ at their exact coordinates. By default all positive-weight TSSs are required; missing one gives NA. Zero-weight TSSs do not require contact.

`pace prepare-promoter-weights` calculates proportional weights from measured ATAC, DNase, H3K4me3 or CAGE promoter signals. These weights are proxies for promoter usage. Gene TPM alone cannot determine TSS usage.

Optional filtering uses **one retained set per gene for every candidate**:

- `promoters.minimum_weight` removes TSSs below the original weight threshold (default 0).
- `promoters.missing_policy: drop_missing` removes a TSS if any planned candidate lacks a resolved contact to it. The default `strict` keeps it.
- The retained original weight must reach `minimum_retained_weight` (default 0.9); otherwise the gene is left unscoreable.

```math
\pi_{\mathrm{used}}(t\mid G)=
\frac{\pi(t\mid G)}{\sum_{u\in\mathcal T_{\mathrm{keep}}(G)}\pi(u\mid G)},
\qquad t\in\mathcal T_{\mathrm{keep}}(G).
```

The original candidates and promoter units remain in the denominator catalog. `promoter_weights.tsv` records original/effective weights, dropped TSSs, missing-contact counts and retained weight. `tss_contact_scope=selected_tss_set` identifies scores for a filtered promoter definition. They do not recover regulation through excluded promoters. Bounds refer to this selected definition too. Different retained sets cannot be treated as full biological score changes in `pace compare`. With sparse maps, dropping every TSS touched by a missing contact may leave too little weight; use a justified contact prior or retain the default NA result.

## Contact policies and priors

```math
P(d)=a\left[\frac{\max(d,d_{\min})}{d_{\mathrm{ref}}}\right]^{-\gamma},\qquad
p(d)=\kappa\min\{P(d),P(d_0)\}.
```

| Situation | Contact used |
|---|---|
| Same bin or configured near-distance range, compatible prior available | P(d) |
| Same bin, no prior, default prior_or_neighbor policy | Maximum of valid neighboring contacts recorded by the cooler adapter |
| Near diagonal with no eligible replacement | NA |
| Finite off-diagonal observation, including zero | Observed contact plus the configured pseudocount |
| Missing observation with allow_prior_fallback=true and compatible prior | P(d), labeled as a prior |
| Missing observation without allowed fallback | NA |
| Explicit prior_only mode | P(d) |
| Explicit shrinkage mode | r C_obs + (1-r) P(d) |

Near-diagonal handling runs first. Same-bin membership depends on bin boundaries and resolution; it is not equivalent to every pair within a fixed ±5 kb distance. `near_diagonal_policy: unresolved` leaves these contacts unavailable. `prior_or_unresolved` permits only a matching prior. The default `prior_or_neighbor` also accepts the recorded neighbor correction. A local maximum is a contact surrogate, not a measured regulatory loop.

`contact.pseudocount: auto` adds p(d) only to finite, off-diagonal observations in observed mode when a compatible prior is available. `powerlaw` requires that prior; `none` disables the addition. The defaults are kappa=1 and d0=5000 bp; these have not been optimized for livestock. Pseudocounts are not added again to near-diagonal replacements, prior-only contacts or shrinkage. Invalid balanced bins remain unavailable unless an explicitly allowed prior replaces them.

Outputs retain `observed_value`, `prior_value`, `pseudocount_value`, evidence type, prior identity and processing reason. A regularized contact is labeled `regularized`; its retained observation coefficient is not a confidence probability.

### Fitting and transferring a prior

`pace fit-prior --cooler ...` estimates decay from distance-bin means, including all callable zero pixels. The default lower fitting distance is the matrix resolution. Zero-mean bins cannot enter the log fit and are reported. Optional held-out chromosomes diagnose contact-decay fit; they do not validate enhancer–gene predictions. A prior mixed with observed contacts must match their scale, normalization, resolution, balancing and measurement window.

A different tissue normally fails the context check. `contact.allow_cross_context_prior: true` permits a prior from the same species and assembly, with the same target level. Its original context is retained, the target context is recorded, and the transfer is marked unvalidated. It cannot be used in the `validated` profile. Renaming a scale or tissue does not calibrate a prior; assess transferred decay on target data when possible.

The optional `prior_preset: abc_human` requires `mode: prior_only` and the research profile. It uses the [ABC reference configuration](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction/blob/main/config/config.yaml) gamma=1.024238616787792, with a=1 and d_ref=d_min=5000 on a relative scale. This amplitude cancels in prior-only scores. It is a transferred human baseline, not a fitted livestock parameter, and cannot be mixed with measured contact tables. Measured activity is still required.

### Fit from a contact map

```bash
pace fit-prior --cooler data/animal1.mcool::/resolutions/5000 \
  --species chicken --assembly GRCg7w --tissue liver \
  --scale balanced_contact_protocol_1 --normalization-id hic_norm_protocol_1 \
  --test-chromosomes chr2 --out models/animal1_contact
```

Use chromosome names actually present in the file. `--test-chromosomes` is optional; omitting it reports that no holdout was performed. It is not correct to rename a liver fit as a brain-validated asset. The default lower fitting distance is the matrix resolution, the upper distance is 5 Mb, and the fit uses 30 logarithmic bins. An explicit `--min-distance` cannot be smaller than the resolution. `--no-balanced` uses raw counts and requires matching raw-count preparation. The command only accepts fixed-bin, symmetric-upper coolers.

The output includes `manifest.json`, `distance_bins.tsv`, and `fit_report.json`. All callable bin opportunities enter distance means, including zeros absent from sparse storage; invalid balancing weights are excluded. Fitting is streamed over sparse pixels. The optional chromosome holdout assesses decay fit, not enhancer function or generalization across tissues.

Use the same scale/normalization/balancing/resolution when preparing contact observations. For observations extracted from this map:

```bash
pace contacts -i data/animal1.mcool -r 5000 -d prepared/catalog \
  --scale balanced_contact_protocol_1 --normalization-id hic_norm_protocol_1 \
  -o prepared/animal1_hic
pace run ... --contacts prepared/animal1_hic/observed_contacts.tsv \
  --contact-prior models/animal1_contact --allow-prior-fallback
```

`pace predict --hic` performs both steps automatically.

The matching prior supports near-diagonal replacement, sparse-count regularization and explicitly requested missing-contact fallback. Rows retain their individual evidence sources. Multiple samples require comparable scales before pooling; fitting one sample's raw-count scale does not calibrate another sample.

`pace fit-contact-prior` remains available for explicitly prepared distance/contact tables.

### Fit from a prepared contact table

```bash
pace fit-contact-prior --contacts contact_fit.tsv --model-id contact_prior_research \
  --species chicken --assembly assembly_identifier --tissue liver \
  --scale balanced_contact_protocol_1 --resolution 5000 \
  --normalization-id hic_norm_protocol_1 --balancing balanced \
  --bin-edges 5000 10000 20000 50000 100000 500000 5000000 \
  --d-ref 10000 --d-min 5000 -o models/prior
```

The TSV has `bin_pair_id,distance_bp,contact_value,split,region_id`. Each measured pair appears
once, including real zeros. Train/test regions cannot overlap. Distance bins use all valid
values; positive mean bins enter a log-linear fit weighted by pair count. Empty and zero-mean
bins are reported separately. A nondecreasing fitted curve is rejected. The d_min distance floor is
an explicit near-distance rule. Parameters, bins, training regions and held-out residuals are
saved. This is a simple empirical prior, not an optimal count-noise model.

### Transfer between tissues

Fit a prior using the source tissue name. For 25 kb matrices the fitting lower
bound defaults to 25 kb:

```bash
pace fit-prior --cooler reference_liver.mcool::/resolutions/25000 \
  --species chicken --assembly GRCg7w --tissue liver \
  --out priors/liver
```

Score the target tissue (muscle) with it:

```bash
pace run -d prepared/catalog --activity prepared/muscle_activity.tsv \
  --species chicken --assembly GRCg7w --tissue muscle \
  --contact-mode prior_only --contact-prior priors/liver \
  --allow-cross-context-prior -o results/muscle
```

or `pace predict ... --tissue muscle --prior priors/liver` after adding the same
cross-context permission through `pace run`. Keep the
prior's source tissue name intact. The output records source and target contexts
and marks the transfer unvalidated; it cannot use the `validated` profile. Measured
muscle activity remains necessary. If combining a prior with measured target
contacts, first establish matching normalization, scale, resolution and windows;
a source-tissue amplitude is not automatically calibrated to a target library.

## Incomplete scores and support bounds

Writing S=A_star×Cbar (or its explicitly selected experimental extension),

```math
\mathrm{PACE}_{\mathrm{conditional}}(E,G)=
\frac{S(E,G)}{\sum_{e\in\mathcal E^{\mathrm{score}}(G)}S(e,G)}.
```

`scoring.partial_policy: conditional` also writes this quantity into the primary column for compatibility; `score_scope` still marks it conditional. The default is `withhold`. An incomplete gene with observed supports 2 and 1 has a conditional first score of 2/3. An unseen support of 97 would make its full score 0.02.

For nonnegative support ranges L_i <= S_i <= U_i:

```math
\mathrm{PACE}_{i,\mathrm{lo}}=\frac{L_i}{L_i+\sum_{j\ne i}U_j},\qquad
\mathrm{PACE}_{i,\mathrm{hi}}=\frac{U_i}{U_i+\sum_{j\ne i}L_j}.
```

Resolved supports are fixed. Unresolved supports default to [0,infinity); `inputs.support_bounds` accepts justified limits and a `bound_source`. Bounds do not impute a primary score. For supports (2,1,unknown), the first full-score range is [0,2/3]; bounding the unknown support to [2,4] narrows it to [2/7,2/5]. These are sensitivity ranges conditional on positive total support, **not confidence intervals**. A zero total yields NA. Correlated unknowns can make these ranges conservative.

Supply optional bounds as `inputs.support_bounds`:

```text
element_id	gene_id	support_lower	support_upper	bound_source
```

These bounds apply to final unnormalized support, including any selected
allocation exponent, and only to unresolved edges. An NA upper bound is unbounded.
Record the basis for each bound. Reassess it after changing activity/contact scales
or the promoter definition, and do not estimate it using held-out test outcomes.

## Additional omics

RNA, auxiliary histone marks, CTCF and DNA methylation are annotations by default.
They do not change the primary score. Register their samples and sources, and keep
measurement units, normalization, assay and windows explicit. There is no automatic
penalty for methylation or multiplier for CTCF binding.

### Histone and CTCF tracks

The bigWig adapter can quantify additional marks using the same cells:

```bash
pace features -b data/H3K4me1_peaks.bed -d prepared/catalog \
  --feature-prefix H3K4me1 -o prepared/H3K4me1
```

Signal (rather than peak overlap) from other marks can be written in the activity
table format and merged with `pace merge -t observed_activity`; register their
sample metadata with `--samples`. Keep `--panel` unchanged: additional assays remain named
annotations. For promoter H3K4me3, quantify the actual promoter window with a
distinct identifier and import it as a promoter feature. A ±2 kb promoter window
must not be labeled `grid:500:mean`. Gene-level summaries of promoter features use
the run's promoter weights; optional TSS filtering changes that promoter definition.

For interval overlap and motif orientation:

```bash
pace features -b data/CTCF_motifs.bed -d prepared/catalog \
  --source-id CTCF_motif_source --evidence-id CTCF_motif_overlap \
  --feature-prefix CTCF --motif-strands -o prepared/CTCF
```

Add the resulting table to `pace run` with `--features`. `--motif-strands`
requires strand-bearing motif intervals. ChIP-seq peaks alone provide occupancy,
not motif orientation; leave the option out for those peaks. Neither overlap nor orientation
is a measured loop probability. Keep element, promoter, gene and edge features
separate.

### RNA

Supply `pace run --expression` with `gene_id`, `sample_id`, `tpm`, `status`. Genes must
exist in the catalog and samples must be registered as RNA. A non-observed status
remains unavailable even if a numeric value appears in the row.

Transcript TPM can be prepared using the catalog's mapping:

```bash
pace expression --tpm data/transcript_tpm.tsv -d prepared/catalog -o prepared/rna
```

TPM is summed only through the supplied mapping. Gene TPM does not determine
transcript or TSS usage. Extra RNA measurements from other animals require an
explicit population-level interpretation.

### CpG methylation

`pace run --methylation` takes raw dyad counts, not percentages or summary tables:

```text
chrom	dyad_start0	methylated_count	total_count	sample_id	assay
```

`dyad_start0` is the zero-based C position in the forward-reference CpG. Keep one
row per dyad/sample/assay, with methylated count no greater than total count. Merge
complementary strands before importing; they are not independent CpGs. The Python
`merge_stranded_cpg` adapter checks explicit-strand cytosine calls against a
reference FASTA and rejects duplicate strands or non-CG bases.

Add these options to `pace run`:

```bash
pace run ... --methylation data/cpg_dyad_counts.tsv \
  --methylation-min-coverage 5 --promoter-upstream 2000 --promoter-downstream 500 \
  --reference-cpg data/reference_cpg_counts.tsv
```

Coverage 5 is an example; the default is 1. Choose the threshold before evaluating
predictions. Upstream/downstream follow transcription direction: positive-strand
windows are `[tss0-upstream, tss0+downstream+1)` and negative-strand windows exchange
the two sides. Starts are clipped to zero. Reference CpG counts must match the
actual window and reference boundaries.

The reference table has `entity_type`, `entity_id`, `n_cpg`; types are `element` or
`promoter`, and IDs match the corresponding catalog. `n_cpg` counts all reference
CpGs in that window, not only observed CpGs. The older `element_id`, `n_cpg` form is
also accepted for elements. Without reference counts, coverage fraction is NA.

The run produces both element and promoter features. WGBS and RRBS remain separate,
for example `DNA_methylation:WGBS:M_site` and `DNA_methylation:RRBS:M_site`. Missing
RRBS coverage is unknown, not zero methylation.

| Field | Calculation |
|---|---|
| `M_site` | Mean of individual CpG methylation fractions after the coverage filter |
| `M_pooled` | Sum of methylated counts divided by sum of total counts; deeper CpGs weigh more |
| `covered_cpg` | Number of CpGs meeting the coverage threshold |
| `cpg_coverage_fraction` | Qualified CpGs divided by reference CpGs in the same window |
| Status | Distinguishes measured zero, no CpG, no coverage and low coverage |

Ordinary bisulfite assays do not independently distinguish 5mC from 5hmC.

For a standalone element summary:

```bash
pace methylation --counts data/cpg_dyad_counts.tsv -d prepared/catalog \
  --min-coverage 5 --reference-cpg data/element_reference_cpg_counts.tsv \
  -o prepared/methylation_qc
```

Here `reference_cpg` has `element_id`, `n_cpg` columns. The output
`methylation_summary.tsv` is for inspection; it cannot replace raw counts in
`--methylation`. To import a summary as custom features, convert it explicitly
and register the evidence.

### Other feature tables

`pace run --features` accepts:

```text
entity_type	entity_id	feature_name	value	evidence_id	status
```

Types are `element`, `promoter`, `gene` or `edge`; an edge ID is
`element_id|gene_id`. Keep one row per entity type, entity ID and feature name.
Aggregate replicate measurements by a stated rule before importing them; repeated
rows cannot serve as weights. Unknown, invalid and unmeasured values remain
missing. Also specify `assay`, `unit`, `normalization_id` and `window_id` so training
and inference can check the same measurement definitions.

Feature names and entity levels must match the trained model. Accepting a feature
format does not show that it improves prediction accuracy. With independent
functional labels, a separate model can be trained as described below and used via:

```bash
pace run ... --ml-model models/liver_classifier
```

Compare the base-only classifier with the added features on independent data.
Coverage, label power and the training population affect that comparison.

## Experimental target allocation

The optional extension adds a cross-gene contact share:

```math
B(E,G)=\frac{\overline C(E,G)}{\sum_{H\in\mathcal G(E)}\overline C(E,H)},\qquad
\mathrm{PACE}_{\eta}(E,G)=
\frac{A_\star(E)\overline C(E,G)B(E,G)^\eta}
{\sum_{e\in\mathcal E(G)}A_\star(e)\overline C(e,G)B(e,G)^\eta}.
```

At eta=0 the allocation term and its data requirements are omitted, giving the main formula. `--eta auto` (the default) uses zero without functional calibration or an explicitly supplied eQTL weak model. A fixed nonzero value is an explicit experimental choice. At nonzero eta, missing contact to any candidate target prevents resolving B for that element; candidate genes are not silently removed.

Cbar × B^eta equals Cbar^(1+eta) / (sum_H Cbar)^eta. This strengthens contact contrasts and depends on gene annotation density. It is not an established biological competition law. Compare it against eta=0 and the separate `contact_power_2` control using held-out functional data. The calibration procedure below uses independent functional labels.

```bash
pace fit-eta -d prepared/catalog --activity ... --contacts ... \
  --species pig --assembly Sscrofa11.1 --tissue liver \
  --eta-labels functional_labels.tsv -o results/eta
```

The same fitting is available during a run with `--eta-labels`. Reuse
`results/eta/eta_calibration.json` with `--eta-model` for compatible later runs.
Any finite exponent in [0,1] can be specified manually; an explicit number never
triggers fitting. True-zero support stays zero and zero normalization totals give
NA. Fitting eta does not train the separate classifier.

### Calibration objective

The estimator uses functional enhancer inhibition/deletion experiments with
direction-qualified positives and powered negatives **within each gene**. Define
`a = log(A*Cbar)` and `b = log(B)`. For a positive p and negative n for the same gene,
their log PACE ratio equals their log support ratio; the gene denominator cancels.
The specified ranking objective is

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

### Required functional-label table

`--eta-labels` takes a tab-delimited file with these columns:

| Columns | Requirements |
|---|---|
| `label_id` | Unique nonempty assay-label identity |
| `element_id`, `gene_id` | One fixed canonical element–gene edge; do not copy a multi-element perturbation onto individual tiles |
| `species`, `assembly`, `context_id`, `target_level` | Must match the scored scientific context |
| `perturbation_type` | `CRISPRi`, `deletion`, or `enhancer_inhibition`; `synthetic_inhibition` is demonstration-only |
| `label_status`, `effect_direction` | `enhancing_positive` with `down`, or `powered_negative` with `none` |
| `split` | train, calibration or test; train fits when present, calibration then confirms without refitting; calibration-only data require grouped validation |
| `group_id` | Prespecified independent locus/experimental group; no group may cross splits |
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
subset must also meet this requirement. The minimum may be prespecified as an integer
at least 2. These checks do not supply a biological sample-size calculation.

### Deployment validation after candidate fitting

The convex fit above proposes a candidate; it does not by itself authorize a
nonzero deployed exponent. Shared group_id, gene_id or element_id joins records
into connected independent components. At least `minimum_groups` (default 3)
components are required. Up to `validation_folds` (default 5) grouped folds refit
the continuous candidate using only each training fold.

The selection procedure tests 21 prespecified shrinkage multipliers from 0 to 1 in
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

### Options and reuse

```bash
pace run ... --eta-labels functional_labels.tsv \
  --eta-min-genes 3 --eta-min-groups 3 --eta-validation-folds 5
```

Alternatively reuse a saved calibration with `--eta-model results/eta/eta_calibration.json`.
The two sources are mutually exclusive. For a prespecified ablation give
`--eta 0`, `--eta 1` or another finite value in the interval without either source.
`--eta-labels`/`--eta-model` select automatic mode; an explicit `--eta NUMBER`
selects manual mode.

Every run exports `eta_calibration.json`; the run manifest repeats its metadata and
stores the numeric exponent in `comparison_contract.eta`. Fitted artifacts record
the objective, fitted/boundary/fallback status, label and fitting-evidence hashes,
used genes/elements/groups, excluded-label reasons, baseline/fitted loss and scope.
Fixed reuse checks scope hashes, context, target level, regime, panel/scales,
catalog/candidate/promoter identities, contact policy and quantitative asset hashes.
Synthetic calibration cannot enter a research/validated run. An unchanged fixed
calibrator may be reused across compatible individuals, not across incompatible
tissues, candidate universes or measurement/model compatibility checks.

For between-animal comparisons, fit once and reuse the same fixed artifact.
Refitting separately can confound genetic/evidence changes with parameter changes.
The benchmark command reports calibrated PACE in addition to its fixed baselines
when a calibration source is configured, and rejects evaluation labels that
overlap fitting genes, elements or groups. Calibrated in-sample scores are training
outputs; reporting them as independent performance is invalid.

## Classifier training and inference

Saved-model inference never refits weights, scales, medians, quantiles or thresholds.

```bash
pace train --data labelled_edge_features.tsv --model-id regulatory_link_classifier \
  --species chicken --assembly assembly_identifier --tissue liver \
  --extra-features H3K4me1 promoter_H3K4me3 methylation_M_site \
  --penalty 0.01 0.01 --penalty 0.1 0.01 --folds 3 --seed 17 --calibrate \
  --feature-contract results/measured/ml_feature_contract.json \
  -o models/classifier
```

For a non-synthetic model the full feature definition is required: pass the
`ml_feature_contract.json` of a compatible scoring run with `--feature-contract`.

The definition includes activity panels/scales, contact resolution and normalization,
candidate definitions, target/estimand and auxiliary feature preprocessing.
[The bundled training fixture](../examples/training/learning.yaml) demonstrates a
complete synthetic definition.

Prepare one-to-one label mappings before assembling this table. Required columns:
`element_id,gene_id,assayed_region_id,mapping_count,group_id,split,label_status,effect_direction,
A_used,Cbar,distance_bp,pace_score,regime,activity_sources,contact_sources`, plus named extra
features and optional `sample_weight`. Splits are train/calibration/test. Enhancing positives
have `label_status=enhancing_positive,effect_direction=down`; negatives must be explicitly
`powered_negative`. Upregulation, low power, untested regions and ambiguous multi-tile mappings
are excluded and exported. These status labels assert the user's experiment-specific effect,
significance and power rules; PACE does not invent those rules from a p value.

Group, edge and perturbation-region identities cannot cross splits. Within training,
connected repeated entities also stay together across tuning folds. Grouped tuning refits
preprocessing inside each training fold; infeasible groups/classes cause an error. A single
prespecified penalty pair permits research fitting without pretending cross-validation occurred.
Core-missing edges are excluded; only extra features can be imputed. Training medians/IQR,
missing indicators and removed all-missing columns are fixed in JSON. The objective is mean
sample-weighted logistic loss + lambda1·L1 + lambda2·L2²/2, with unpenalized intercept, implemented
by proximal gradient. Convergence is reported. A base-only classifier is saved and evaluated.

Optional sigmoid calibration uses the independent calibration split, with a specified 1e-6
slope penalty for numerical stability. Without it, probability is NA and `pace_ml_score` is an
uncalibrated classifier score. Neither is averaged with PACE. Real calibration probabilities
remain specific to the recorded perturbation and candidate sampling design.

```bash
pace predict-ml --run results/measured --model models/classifier \
  --features results/measured/multiomics_features.tsv.gz -o results/ml
```

Inference checks the complete feature definition, regime, evidence-source and context scope.
Mismatches return out_of_scope with unavailable scores/probabilities. Old unbound models
are research-ineligible; a demonstration-only legacy score is explicitly unverified.
Synthetic probabilities, where present, are labelled synthetic_demonstration_only. Extra features must have unambiguous
entity-qualified names, and repeated annotation measurements need specified aggregation before
pivoting into a learning feature. Gene TPM is off by default; use it only as an explicit ablation.

### Held-out reports

The test report uses the same inference path as deployment. `test` evaluates `pace_ml_score`; `test_probability` separately evaluates calibrated probabilities. `test_status_counts` reports out-of-scope or unavailable predictions, and coverage includes their effect. Inspect `optimization_converged` and `calibration_optimization_converged` before interpreting a fit. Training and inference are restricted to measured activity evidence.

## Comparisons and benchmarks

PACE preserves each original run and recomputes conditional support shares on the intersection
of measurable units. A simple join-and-subtract of existing normalized scores is not valid when
their denominators differ.

### Compare two runs

```bash
pace compare --left results/animal_a --right results/animal_b \
  --min-common-units 2 -o results/comparison
```

The comparison checks catalog, candidate and promoter hashes, target level, context, panel,
units, normalization and formula parameters. At least two common units and positive common
denominators are required. Complete Delta also requires complete planned backgrounds and
compatible structural status. Technical absence leaves complete Delta as NA, even when a
conditional Delta can be reported. Eta comparisons can be explicitly enabled, but then no
complete-background score difference is claimed. B retains its original candidate-gene set in each run.
Eta differences are enabled with `--allow-eta-difference`. Contact-policy and measurement-quality comparisons can similarly enable `--allow-evidence-difference`; these
produce conditional comparisons, with full Delta withheld because the inference policy differs.

The table contains original/common scores, full/conditional Delta, delta A, delta support,
gene totals, common denominators, comparison IDs and reasons. All deltas are right minus left.
Do not infer enhancer activity or expression direction from the sign of a support-share change.

### Replicate stability

```bash
pace stability --replicates replicate_runs.tsv -o results/stability
```

The TSV contains `run_path,donor_id,replicate_type`, where paths resolve relative to the TSV and
replicate type is biological or technical. Every pair is compared on recomputed common
denominators. Pearson correlation is NA for constant or insufficient vectors. Distinct-donor
counts are reported; technical replicates are never counted as independent animals. No
confidence interval is manufactured from a small number of animals.

### Functional-label benchmark

When the run config supplies eta calibration labels or a fixed artifact, the benchmark
adds `PACE_calibrated` alongside the endpoint baselines. It rejects functional evaluation
labels sharing fitting genes, elements or groups. A fixed continuous exponent adds
`PACE_fixed_eta`. Freeze eta before between-animal comparisons; separate refits can
confound evidence changes with parameter changes. See [eta calibration](#experimental-target-allocation).

```bash
pace benchmark --run results/measured --labels functional_labels.tsv \
  --membership prepared/catalog/region_membership.tsv --stratify gene_id \
  --external gABC gabc_scores.tsv recorded_external_version recorded_external_configuration \
  -o results/benchmark
```

`--external` is optional and repeatable. Decision thresholds are omitted by default. If a
threshold has been independently fixed from appropriate training/calibration data, give
`--threshold METHOD VALUE SOURCE_SPLIT SOURCE_ID`; omit it when no such evidence exists. The benchmark never selects
a deployment threshold using the test set. Labels use the [standard dictionary](../README.md#input-table-reference)
and preserve every tested positive in the evaluation universe, including absent predictions.
The current implementation trains/evaluates only unambiguous one-to-one region mappings; ambiguous
regions are exported rather than copied to several tiles.

Built-ins are negative distance ranking, `ABC_style_single_TSS`, PACE eta 0/1, and
`contact_power_2` (A×Cbar²). The single TSS
is chosen by smallest coordinate then promoter ID before evaluation, and saved in the report.
All configured formula methods are additionally normalized on their common scoreable sets. External
scores require version/configuration provenance and columns `element_id,gene_id,score`; absent
files become `not_available`. Without external raw support, their common-denominator status is
explicitly not assessed. PACE does not claim to reimplement the external gABC package.

Average precision is the tie-aware sum of recall increments times precision, not trapezoidal PR
area. AUROC is reported only with both classes. Coverage, missing tested positives, candidate
recall and end-to-end recall are separate. Missing scores are not filled with zero. Requested
strata are explicit label columns and include their sample counts. QTL association and input
Hi-C are not substituted for independent functional perturbation labels.

## Worked calculations

These small analytical examples test the current formula; they are not biological
performance results. The numerical tests cite these same hand-derived quantities.

### Fixed-panel activity and TSS contact

ATAC=4 and H3K27ac=9 give activity sqrt(4*9)=6. If H3K27ac is measured zero,
activity is zero. If it is unavailable, activity is NA under the two-assay panel.
Contacts (2,6) with promoter weights (3/4,1/4) give Cbar=3. A missing required TSS
contact remains unavailable by default. An explicit gene-wide TSS filter may
renormalize a retained set, subject to a minimum retained-weight check.

### Endpoint and continuous allocation

Let activities be (4,2,1), G1 contacts (3,2,2), and G2 contacts (1,2,6).
Then B for G1 is (3/4,1/2,1/4), and B for G2 is (1/4,1/2,3/4).

| Exponent | G1 support | G1 PACE | G2 PACE |
|---|---|---|---|
| 0 | (12,4,2) | (2/3,2/9,1/9) | (2/7,2/7,3/7) |
| 1 | (9,2,1/2) | (18/23,4/23,1/23) | (2/15,4/15,9/15) |

For eta=1/2, G1 support is (6sqrt(3),2sqrt(2),1) and G2 support is
(2,2sqrt(2),3sqrt(3)). Divide each vector by its own sum. Every listed denominator
is exactly the sum of those supports. A completely zero vector produces NA scores.

### Bulk aggregation and composition

Two equally weighted measured replicates of E1 have assay pairs (9,1) and (1,9). Two equally weighted measured replicates of E2 each have (4,4).
The bulk calculation first averages assays, producing activities 5 and 4. With
identical contact and eta zero, E1 receives 5/9. Computing support separately per
replicate and summing gives 3/7 instead; that is a different quantity.

Supports (1,1,2) and (1,1,NA) give conditional E1 shares 1/4 and 1/2. The second
primary score is withheld because its full denominator is unknown; its default
sensitivity interval is [0,1/2]. On the common first two elements, both recomputed
shares are 1/2; a complete-background difference remains unavailable.
Changing supports (1,1) to (1,2) reduces the first share from
1/2 to 1/3 even though its own support is unchanged.

### Continuous calibration check

Two equally weighted pairwise margins eta-1/2 and 1/2-eta yield the convex loss
`[softplus(1/2-eta)+softplus(eta-1/2)]/2`. Its unique minimum is eta=1/2 and its value
is log(2). This independently checks an interior solution, not just the two endpoints.

[Scoring tests](../tests/test_canonical_core.py) ·
[Continuous calibration tests](../tests/test_canonical_allocation.py) ·
[Main formula](../README.md#how-it-works).

## Validation

Software verification checks calculations and input handling. The tests use
synthetic data and small genomic file fixtures; they do not establish predictive
accuracy for a species or tissue. See [CONTRIBUTING](../CONTRIBUTING.md) for test,
lint and packaging commands, and inspect the
[CI run for the exact commit](https://github.com/shenlinyong/PACE/actions/workflows/ci.yml).
Optional IO tests can be skipped when their dependencies are absent.

| Area | Invariant |
|---|---|
| Formula | Independent hand calculations; fixed assay and candidate sets; exact zeros and NA |
| Measured scope | Retired modes/options/configurations fail; unavailable activity is never imputed |
| Contact | Compatible resolutions and sources; same-bin imports use the same policy as raw inputs |
| Imports | Real sample identities and correct assay; measured activity only |
| Allocation | Zero fallback, grouped validation, bounded fit and test-set isolation |
| Comparison | Common denominators; underflow-safe reloads; partial results remain conditional |
| Multiomics | Experimental annotations, measurement definitions, CpG coverage and RNA status |
| Classifier | Grouped splits, inference scope, separate score and probability evaluation |
| Packaging | Supported Python versions, wheel outside checkout, measured demo, Conda and Docker |
| Documentation | README equation, current defaults, working links and executable interfaces |

This is a study protocol, not a report of completed biological experiments. No performance increase is asserted by the software release.

### Primary question

Does contact regularization and explicit handling of missing support preserve useful regulatory-link ranking and coverage as experimental data become sparse? This is a measurable claim. Human downsampling experiments provide a controlled stress test; they cannot by themselves establish transfer accuracy in chicken, pig or cattle.

### Functional reference and split

Use experimentally tested K562 enhancer–gene pairs from the [EngreitzLab CRISPR benchmark](https://github.com/EngreitzLab/CRISPR_comparison), recording release, genome assembly, assay and effect/power definitions. Fulco/Gasperini datasets overlap across compilations: deduplicate perturbation regions, genes and assays before splitting. Never count a reused benchmark as an independent replication. Negatives should be tested, sufficiently powered non-effects; untested pairs and low-power experiments are not verified negatives.

Assign chromosome or connected region/gene groups to training, calibration and final test sets. All tuning, allocation selection, promoter weighting choices and contact preprocessing decisions must precede final test evaluation. Keep repeated guide/region/gene entities in one split. Record both the original benchmark's screening selection and exclusions caused by candidate mapping.

### Controlled degradation

Use a full-quality K562 reference and predeclare contact fractions such as 100%, 50%, 25%, 10% and 5%, with several random seeds. These are experimental design choices, not PACE defaults.

1. Thin raw read pairs or unique raw Hi-C pixel counts using binomial sampling, preserving symmetry and counting each pair once. Do not thin a balanced matrix or separately resample duplicated enhancer–TSS queries from the same bin pair.
2. Rebuild/balance the sparse map at each depth and refit its prior on training chromosomes. A full-depth fitted prior reused at low depth is an extra-information condition; label it separately.
3. Remove a whole activity assay and explicitly choose the remaining fixed single-assay panel. Separately simulate missing regions while retaining the fixed panel, so assay ablation is not confused with selective missingness.
4. Reduce independent donors/biological replicates. Technical replicates do not increase independent animal count. Keep population summaries distinct from individual estimates.
5. Freeze the candidate universe for the primary comparison. Separately examine candidate-construction sensitivity, including peak-centered windows and fixed-grid offsets. Report mapping/coverage differences rather than selecting a favorable catalog post hoc.

### Baselines and ablations

Run the pinned original ABC implementation with its own documented contact processing; export its scores as an external benchmark method with version/configuration provenance. PACE's `ABC_style_single_TSS` is a mathematical control using PACE contacts, not a substitute for running original ABC.

Compare: PACE eta=0; no contact pseudocount; no near-diagonal correction; prior-only; explicit missing-contact fallback; single versus multiple TSSs; equal versus measured promoter-signal weights; B allocation versus a pure contact-power control. `pace benchmark` includes `contact_power_2` (A×Cbar²) to compare with eta=1 allocation. If tuning contact power or eta, tune each with the same independent validation budget. Report annotation and biotype sensitivity for any B-related claim.

### Metrics and uncertainty

Report average precision (AP, with its non-trapezoidal definition), precision–recall curves, recall at an independently selected precision target, and candidate/positive-label coverage. Because default PACE withholds incomplete primary scores, a higher AP among fewer available predictions is not sufficient. Report AP on a prespecified common evaluable set **and** abstention/coverage on all tested pairs. Show how many positives remain unscored and how coverage changes with depth.

Use paired bootstrap over independent gene/region groups, with the same sampled groups for every method. Chromosome holdout and the number of independent genes constrain precision; random edge bootstrap is not an independent biological replication. Report per-distance, promoter-count, gene-density and chromosome strata to expose confounding. Do not select a seed, cutoff or degradation level after inspecting final test performance.

For score intervals, use complete high-quality supports to check containment after deliberately masking candidates. Evaluate interval width and any assumptions used for finite support bounds. The default [0,infinity) missing-support bounds can be very wide: this is honest lack of information. These mathematical ranges are not 95% confidence intervals and are not automatically a novel method simply because they are implemented here.

### Livestock evidence

Use species- and tissue-matched FarmGTEx/eQTL resources as orthogonal association evidence. Match or stratify control pairs by enhancer–gene distance, allele frequency, LD, gene expression, variant and candidate density, and assay accessibility. Prefer fine-mapped/colocalized signals where available; do not treat each correlated SNP as an independent regulatory event. eQTL association is not a direct enhancer perturbation label.

Use available independent chicken, pig and cattle perturbations or targeted experiments as final functional evidence. Keep the target of each conclusion precise: human sparse-data robustness, livestock association enrichment, or livestock functional-link accuracy. Public gene annotations and a more elaborate QC system alone do not establish improved accuracy in livestock.

## Limitations

PACE is research software for relative regulatory support from measured activity and specified contact evidence. Software tests check calculations and input handling; they do not establish predictive accuracy for a livestock species, tissue, breed or experimental condition.

| Capability | Boundary |
|---|---|
| Activity | Requires measured ATAC, DNase or H3K27ac; no replacement of unavailable activity |
| Candidate catalog | Supplied/experimental regions and trusted promoter annotations; completeness refers to this planned catalog |
| Contact | Measured contacts or an explicitly supplied applicable distance prior; prior evidence is not an observed regulatory loop |
| Contact shrinkage | Fixed weights or per-pair Gamma–Poisson posteriors; the count model and matching measurement scale must be appropriate |
| Multiple TSSs | Fixed distinct physical promoters and weights; incorrect annotations remain a limitation |
| Replicates | Specified technical/biological/donor aggregation; no automatic batch correction or increase in independent sample size |
| Additional omics | Named experimental annotations; optional separate classifier with its own validation |
| Allocation | Separate functional and eQTL weak calibration; neither forces a positive eta, and independent evaluation remains necessary |
| Calibrated probability | Specific to supplied functional labels and sampling design; not universal causality |
| Between-animal comparison | Comparable measured states on common denominators; not an isolated genetic or causal effect |
| Cross-species/tissue use | Explicit background checks; any accuracy claim requires target-scope evaluation |
| Scalability | Sparse candidate contact processing; training tables reside in memory |

No universal livestock weights or independent biological benchmark results are bundled. A high conditional fraction can occur when only a small subset is measurable; the default primary score is withheld in that case. Always report candidate construction, assay availability, exclusions and coverage. Do not interpret partial normalization as complete regulatory discovery.

## Boundary priors and weak calibration

### Contact prior

The boundary-aware prior is an optional contact model:

```math
C_0(E,t)=a\left(\frac{\max(d(E,t),d_{\min})}{d_{\mathrm{ref}}}\right)^{-\gamma}
\exp\!\left[-\beta\sum_{b\;\mathrm{between}\;E,t}s_b\right].
```

Boundary positions strictly between the two anchors contribute; boundary points
at either endpoint do not. Beta is nonnegative. At beta=0 the implementation
returns the existing distance-prior function exactly. This nests the distance
shape, without asserting that all ABC preprocessing is identical. The standard
reference distance is 5,000 bp; fitted assets record their actual reference and
minimum distances.

`pace boundaries` accepts either motif hits (`--chip` is optional; omit it without occupancy data):

```bash
pace boundaries --motifs ctcf_motifs.tsv --chip ctcf_peaks.tsv --max-gap 1000000 -o boundary_asset
```

or a FASTA scan (12.0 is an example log2-odds threshold; choose it for your PWM):

```bash
pace boundaries --fasta reference.fa --pwm ctcf_pwm.tsv --threshold 12.0 --max-gap 1000000 -o boundary_asset
```

Motif tables require `chrom`, `start`, `end`, `strand` (+ or -), and `strength`
(0–1). Coordinates are zero-based, half-open. PWM tables have one row per motif
position and columns `A`, `C`, `G`, `T` containing counts or probabilities. A 0.01
pseudocount is added to each entry before row normalization; scores use uniform
background frequencies and both strands. Ambiguous bases cannot match. This
scanner does not calculate motif p-values or learn a species-specific motif.
ChIP tables require `chrom`, `start`, `end`; only overlapping motifs are retained.

Overlapping motif hits collapse to their strongest representative. Each adjacent
minus-strand then plus-strand pair within `max_gap_bp` contributes a point at the
midpoint of the intervening gap. Its strength is the smaller normalized motif
strength. This is a sequence/occupancy heuristic, not a demonstrated boundary.
Outputs retain `candidate_unvalidated` status. Supply a calibrated threshold and
check chromosome naming before interpreting a whole-genome scan.

`pace prior` packages supplied parameters (`a` is illustrative and must match the
intended contact scale; `gamma` is an explicit assumption until fitted; `--kappa`
and `--pairs` are optional):

```bash
pace prior --model-id pig_liver_boundary_prior \
  --species "Sus scrofa" --assembly Sscrofa11.1 --tissue liver --target-level population_mean \
  --scale cooler_native --resolution 10000 --normalization-id liver_cooler_v1 --balancing balanced \
  --boundaries boundary_asset/boundaries.tsv \
  --a 1.0 --gamma 1.0 --beta 0.0 --d-ref 5000 --d-min 10000 -o prior_asset
```

For prior-only ranking an arbitrary positive amplitude cancels within a gene.
For fusion it does not: fit amplitude on the same measurement scale as the Hi-C.
The boundary table is copied into the asset and checksum-bound. A positive beta
requires an explicit table. Empty tables are allowed if they have headers and
represent an intentional absence of candidate boundaries.

### Fitting raw Hi-C and estimating reliability

`pace fit-hic` uses the same metadata, `--boundaries`, `--d-ref` and `--d-min`,
replacing `--a`, `--gamma`, `--beta` and `--kappa` with grids. Measurement labels
(scale, resolution, normalization, balancing, window) are read from the contact
table when not given. `--test-chromosomes` is optional here (mandatory for fit-labels):

```bash
pace fit-hic --contacts prepared/pig1_hic/observed_contacts.tsv \
  --gamma-grid 0.6 0.8 1.0 1.2 --beta-grid 0 0.5 1 2 --test-chromosomes 18 \
  --boundaries boundary_asset/boundaries.tsv --d-ref 5000 --d-min 10000 \
  --species "Sus scrofa" --assembly Sscrofa11.1 --tissue liver -o fitted_prior
```

Input columns include `chrom`, `anchor0`, `tss0`, `sample_id`, `bin_pair_id`,
`resolution`, `measurement_status`, `raw_count`, `count_to_contact`,
`contact_value`, `scale`, `normalization_id`, `balancing`, and `window_id`.
The cooler adapter exports these with the two balance weights. Raw counts must
be integers; `raw_count * count_to_contact` must reproduce the observed value.
For native balanced coolers the conversion is the product of the two weights;
for native unbalanced coolers it is one. If contacts undergo further depth
normalization, the conversion factor must incorporate that exact transformation.
Do not feed O/E, log contacts, correlations, or normalized noninteger values as
Poisson counts. This model treats balancing factors as fixed exposures; it does
not model uncertainty introduced by balancing itself.

Fit input must retain measured zero pixels. With sparse coolers, use the
default of `pace contacts` (absent pixels are measured zeros) only when that is
true; otherwise use `--missing-as-missing`. Missing or masked rows are excluded
from fitting. Shared sample/bin pairs count once. Diagonal pairs are excluded;
other pairs use bin centers, matching the resolution of the measurement. At
least three distinct callable pairs with some positive counts are required.

For each gamma/beta grid point, the Poisson profile likelihood estimates the
amplitude analytically. The best grid point defines the expected raw counts.
The Gamma shape is then estimated by moments, bounded to [1e-6, 1e6]:

```math
\widehat{\kappa}^{-1}=
\frac{\sum_i[(y_i-\mu_i)^2-y_i]}{\sum_i\mu_i^2}.
```

Nonpositive excess variation uses the upper bound, recorded in `fit_report.json`.
Held-out chromosomes do not enter amplitude, grid selection or dispersion.
`held_out_residuals.tsv` reports their observed and expected raw counts. The
report also includes the correlation between log distance and boundary strength;
large correlation makes individual gamma/beta estimates difficult to interpret.
Select/query fitting pairs independently of their positive contact counts.
A candidate-enriched fitting set estimates its own background, not necessarily
the genome-wide background.

Enable the posterior with `pace fuse`, which takes the same options as `pace run`
and sets `--contact-mode shrinkage --contact-reliability per_pair`:

```bash
pace fuse -d prepared/catalog --activity prepared/activity/observed_activity.tsv \
  --contacts prepared/pig1_hic/observed_contacts.tsv --contact-prior fitted_prior \
  --species "Sus scrofa" --assembly Sscrofa11.1 --tissue liver -o fused
``` With conversion factor f and prior contact C0:

```math
C\sim\mathrm{Gamma}(\kappa,\mathrm{rate}=\kappa/C_0),\qquad
Y\mid C\sim\mathrm{Poisson}(C/f),\qquad
r=\frac{C_0/f}{C_0/f+\kappa},\qquad
\mathrm{E}(C\mid Y)=r\,fY+(1-r)C_0.
```

Reliability depends on expected information under the prior, not the observed
count alone. True zero counts remain observations and receive a positive
posterior mean. Masked/unavailable pairs use the prior with reliability zero and
no invented observation sample. Diagonal pairs retain the configured
near-diagonal policy. Each valid sample is shrunk separately, then aggregated
using the existing technical/biological/donor hierarchy. Missing samples do not
become extra prior-only replicates.

The default (`kappa: auto`) first uses the fitted asset's kappa; without it, the run estimates
kappa from unique callable bin pairs. A positive numeric kappa overrides this.
The resolved table records the actual kappa, its source, the bin-center prior
policy, and `posterior_samples` containing per-sample means, variances, raw
counts, factors and reliabilities. The variance is conditional on fitted
parameters; it is not a posterior for the final normalized gene score.
Imported per-pair resolved tables require their raw observations and are checked
against recomputation. This prevents reusing posteriors after changing the prior
or count conversion.

### Learning gamma, beta and eta from eQTLs

`pace fit-labels` accepts fine-mapping results
converted to a simple table. Downloading or harmonizing FarmGTEx releases is not
part of this command. The table requires `variant_id`, `chrom`, `pos0`, `gene_id`,
and `pip`; `independent_signals` also requires `signal_id`. Use one row per
variant/gene with PIP in [0,1], already harmonized to the reference assembly and
tissue. Restrict genes to the cis candidate catalog. Duplicate gene/variant rows
are rejected rather than silently counted twice.

```bash
pace fit-labels --run results/pig_liver_baseline --eqtl eqtl_finemapping.tsv \
  --aggregation independent_signals \
  --gamma-grid 0.8 1.0 1.2 --beta-grid 0 0.5 1 --eta-grid 0 0.5 1 \
  --test-chromosomes 18 -o weak_fit
```

`--run` is a finished `pace run` result folder; species, assembly and tissue are
taken from it.

Two aggregation choices are explicit:

- `independent_variants`: element/gene mass is `1 - product(1 - PIP)`. This is an
  independence approximation, not the exact probability for mutually exclusive
  variants in one credible set.
- `independent_signals`: sum PIPs of variants inside the element within each
  signal, then take `1 - product(1 - signal_mass)` across signals. Each gene/signal
  must have total PIP at most one. This still assumes the supplied signals can be
  treated as independent.

All candidate elements for genes present in the fine-mapping table enter the
assessed universe. Elements without PIP mass are marked `unlabelled_background`,
not validated negatives; genes without fine-mapping rows are excluded. Omitted
credible-set mass, eQTL power, LD, expression and variant frequency all affect
this weak target. The model does not infer LD or perform fine-mapping itself.

The run must use `prior_only` or `shrinkage` contact with a frozen source prior,
`--eta auto` (the default), and no existing functional or weak calibrator. For
per-pair fusion, freeze kappa in the prior or configuration before fitting.
Boundary checksums, candidate coordinates, promoter weights, activity scales and
processing policies bind the model to its source setup.

At least three training chromosomes and one final test chromosome are required.
The command performs nested leave-one-chromosome-out selection on the training
set, selects final parameters on all training chromosomes, and evaluates the
untouched test chromosomes once. Candidate gene normalization always uses the
full candidate set, including unlabelled elements. Exact selection ties prefer
eta=0 and beta=0; both values must be in the grid. Eta is allowed to remain zero.
At most 50 values per parameter and 2,000 combinations are accepted; start small
because contact resolution is repeated for each gamma/beta combination.

The metric is tie-aware fractional AP of PIP mass, averaged equally across
chromosomes. It is a descriptive soft-label ranking score, not CRISPR AP, a
calibrated probability, or an unbiased estimate of causal accuracy. Baselines
include nearest TSS (ties retained) and a pure power-law ABC-style score with
matched activity, catalog and TSS weights; this is not a full external ABC
software reproduction. The source contact prior remains fixed across folds;
the held-out claim concerns eQTL labels, not every earlier Hi-C training input.

Apply the saved parameters in the same compatible run:

```bash
pace run <same options as the baseline run> --weak-model weak_fit/weak_model.json -o calibrated
```

`pace run` loads gamma/beta before resolving contacts and eta before final
scoring. Output marks `allocation_evidence: eqtl_weak` and retains chromosome
membership, source hashes and weak validation status in `eta_calibration.json`.
Functional-label calibrators continue to use `--eta-model`; the two cannot
be combined. Weak calibration is unavailable in the `validated` execution
profile. `benchmark` keeps `PACE_eqtl_weak` separate from fixed-eta ablations.

The four-chromosome [offline example](../examples/contact/run.sh) exercises
all five commands. It verifies the implementation only. Evaluate real livestock
data with distance/LD-aware controls and independent functional evidence before
claiming an improvement. Cross-tissue prior transfer remains explicit;
automatic tissue-hierarchical contact learning is not implemented.

## Complete configuration defaults

These defaults match `pace_livestock.config.DEFAULTS`; the documentation check
compares the complete mapping during tests. The same structure is written to
`resolved_config.yaml` in every result folder; `null` means not provided. Ordinary runs need only the settings shown in the
[README command reference](../README.md#command-reference).

<!-- configuration-defaults:start -->
```yaml
schema_version: pace-1
run_id: pace
regime: measured
execution_profile: research
estimand: bulk_proxy
target_level: individual
context:
  species: null
  assembly: null
  context_id: null
inputs:
  units: null
  region_membership: null
  promoters: null
  candidates: null
  samples: null
  observed_activity: null
  observed_contacts: null
  resolved_activity: null
  resolved_contacts: null
  features: null
  methylation: null
  expression: null
  labels: null
  evidence: null
  sources: null
  support_bounds: null
catalog:
  profile: canonical_grid
  width_bp: 500
  offset_bp: 0
  include_promoter_units: true
  chrom_sizes_path: null
  candidate_radius_bp: 5000000
activity:
  panel:
  - ATAC
  - H3K27ac
  combine: geometric_equal
  missing_policy: unresolved
  minimum_callable_fraction: 0.0
  replicate_aggregation: equal_donor_mean
  pseudocounts: {}
contact:
  mode: observed
  scale: depth_normalized_contact
  prior_path: null
  prior_preset: null
  near_diagonal_policy: prior_or_neighbor
  near_diagonal_bp: 0
  allow_prior_fallback: false
  allow_cross_context_prior: false
  reliability: null
  kappa: auto
  reliability_source: null
  resolution: null
  normalization_id: null
  balancing: null
  window_id: null
  pseudocount: auto
  pseudocount_distance_bp: 5000
  pseudocount_strength: 1.0
scoring:
  partial_policy: withhold
promoters:
  weights: provided
  minimum_weight: 0.0
  missing_policy: strict
  minimum_retained_weight: 0.9
allocation:
  eta: auto
  missing_policy: fixed_gene_set
  labels_path: null
  calibrator_path: null
  weak_model_path: null
  minimum_genes: 3
  minimum_groups: 3
  validation_folds: 5
  minimum_positive_fraction: 0.8
multiomics:
  mode: annotate
  model_path: null
methylation:
  minimum_coverage: 1
  promoter_upstream_bp: 2000
  promoter_downstream_bp: 500
  reference_cpg_path: null
comparison:
  full_delta_requires_complete: true
  allow_conditional_intersection: true
  minimum_common_units: 2
output:
  format: tsv_gz
  retain_all_candidates: true
seed: 17
```
<!-- configuration-defaults:end -->
