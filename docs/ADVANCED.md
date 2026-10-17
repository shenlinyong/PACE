# Advanced use

The [README](../README.md) covers installation, inputs and ordinary scoring. This guide documents optional evidence, the contact models, calibration, comparisons and every setting with its default. Every command takes command-line options (`pace COMMAND --help`).

- [Evidence preparation](#evidence-preparation)
- [Activity and alternative TSSs](#activity-and-alternative-tsss)
- [Contact policies and priors](#contact-policies-and-priors)
- [Boundary priors and weak calibration](#boundary-priors-and-weak-calibration)
- [Incomplete scores and support bounds](#incomplete-scores-and-support-bounds)
- [Additional omics](#additional-omics)
- [Experimental target allocation](#experimental-target-allocation)
- [Classifier training and inference](#classifier-training-and-inference)
- [Comparisons and benchmarks](#comparisons-and-benchmarks)
- [Worked calculations](#worked-calculations)
- [Evaluation design](#evaluation-design)
- [Complete configuration defaults](#complete-configuration-defaults)

## Evidence preparation

Register one sample identifier per experimental replicate, with its donor, biological replicate, technical replicate and source, in `samples.tsv`. Use `--target-level population_mean` when pooling donors.

Empty input-table templates:

```bash
pace init --species chicken --assembly GRCg7w --tissue liver --panel H3K27ac --out chicken_liver
```

`--catalog-dir prepared/catalog` imports the catalog tables and settings exported by `pace catalog`.

**Catalog.** Supplied peaks and all TSSs are mapped onto unique cells of a fixed grid (`--width`, `--offset`); the genome is not tiled. TSSs are read from GTF transcript records, or from exon records when a file has none: TSS = start−1 on the plus strand and end−1 on the minus strand. Identical physical TSSs are merged, and versioned gene and transcript IDs are kept. Chromosome aliases are given explicitly (`alias`, `canonical` columns). Incomplete terminal grid cells are excluded and reported; keep the exported `run_catalog_config.yaml` with the catalog. Quantify signal per cell: copying a whole peak's read count to every overlapping cell inflates support.

**bigWig.** Signals must be nonnegative. Unstored bases are missing unless `--missing-as-zero` is given. The callable fraction measures coverage; `--min-callable-fraction` sets the minimum.

**Contact.** Contact values must keep the distance background: convert observed/expected values back with their distance expectation and undo log transforms; p-values and correlations are not contact values. Convert `.hic` to cool/mcool. Bins with invalid balancing weights are unavailable; unstored pixels count as zero when both bins are valid (`pace contacts` default; `--missing-as-missing` turns this off).

### Raw activity counts

For a window table of raw filtered fragment counts:

```text
# counts.tsv
element_id	sample_id	assay	count
# library_sizes.tsv
sample_id	library_size
```

One count per element, sample and assay; one full-library total per sample; NA for missing counts, 0 for measured zeros.

```bash
pace normalize-activity --counts counts.tsv --library-sizes library_sizes.tsv \
  --units prepared/catalog/units.tsv --out prepared/normalized
```

Signal = count × 10⁶ / library size / window length (CPM per bp). The library total counts all retained fragments, and paired-end data use the same fragment convention in numerator and denominator. Variable-width regions need `--window-id`; the default normalization identifier is `CPM_density_v1`. Already normalized bigWigs are read with `pace activity` instead. No cross-animal quantile normalization is applied, so global biological differences are kept.

### Measured promoter weights

`promoter_signal.tsv` has `promoter_id` and a nonnegative `signal` over the same promoter window, one value per physical promoter (CAGE, ATAC, DNase or H3K4me3):

```bash
pace promoter-weights --promoters prepared/catalog/promoters.tsv \
  --signals promoter_signal.tsv --assay H3K4me3 \
  --normalization-id H3K4me3_CPM_TSS500 --out prepared/promoter_weights
```

Use the resulting `promoters.tsv` with `--tss-weights provided`. Every TSS needs a signal; a gene whose signals are all zero is an error unless `--zero-policy equal` is chosen.

### Imported evidence and extra output fields

```text
# resolved_activity.tsv
element_id	assay	observed_value	resolved_value	evidence_id	evidence_type	observation_sample_id	parent_evidence_ids	resolution_status	reason	unit	normalization_id	window_id
# resolved_contacts.tsv
element_id	promoter_id	resolved_value	evidence_id	evidence_type	observation_sample_id	prior_id	reliability	resolved_mode	bin_pair_id	resolution_status	reason	scale
# evidence.tsv
evidence_id	evidence_type	source_id	parent_evidence_ids	model_id	unit	processing_method	checksum
```

Activity imports accept observed or aggregated evidence with sample, assay, unit, normalization and window. Contact evidence is `observed`, `aggregate`, `contact_prior`, `regularized` or `fused` (observation/prior shrinkage). Imported contacts need a resolution, and their normalization, balancing and window identifiers must agree when given. `resolution_status` is `resolved`, `unresolved` or `invalid`.

`resolved_activity.tsv` keeps measured values and records `activity_pseudocount` separately, so reimporting it does not add the offset twice. `promoter_weights.tsv` records `pi_original`, the effective `pi`, `tss_selection_reason` and `n_missing_candidate_contacts`. Score and gene-summary rows include `tss_retained_weight`, `tss_dropped_ids`, `tss_policy_status` and `tss_contact_scope`; `n_tss` counts all TSSs and `n_tss_used` those with positive weight.

`region_scores.tsv` sums the unique member cells of each source region and gene, without a new denominator. `evidence.tsv` and `sources.tsv` record input and derived provenance. `ml_feature_contract.json` records the measurement definitions used by the optional classifier.

## Activity and alternative TSSs

```math
A_\star(E)=\left[\prod_{m\in\mathcal M}x_m(E)\right]^{1/|\mathcal M|}
```

Panels: ATAC, DNase, H3K27ac, ATAC + H3K27ac or DNase + H3K27ac, fixed across elements. A missing assay gives NA; a measured zero gives zero activity. Technical replicates are averaged within biological replicates, biological replicates within donors, and donors with equal weight; `individual` requires one donor, `population_mean` pools donors. Multiplying all signals of an assay by one constant multiplies every activity by a common factor, which cancels in the score.

An assay-specific pseudocount, in units of the normalized signal, can be added after replicate aggregation:

```math
A_{\star,\epsilon}(E)=\left[\prod_{m\in\mathcal M}(x_m(E)+\epsilon_m)\right]^{1/|\mathcal M|}
```

```bash
pace run ... --panel ATAC H3K27ac --activity-pseudocount H3K27ac 0.1 \
  --minimum-tss-weight 0.01 --missing-tss-policy drop_missing --minimum-retained-tss-weight 0.9
```

`activity.pseudocounts` defaults to `{}` (all zero). NA stays NA. When an assay is rescaled, its pseudocount must be rescaled with it.

### TSS selection

```math
\overline C(E,G)=\sum_{t\in\mathcal T(G)}\pi(t\mid G)\,\widetilde C(E,t)
```

Distinct TSSs in the same Hi-C bin share one contact query per element and sample, without counting as replicates. All positive-weight TSSs are required by default; a missing one gives NA. Zero-weight TSSs need no contact.

Optional filtering keeps one TSS set per gene, shared by all its candidates:

- `--minimum-tss-weight` removes TSSs below an original weight (default 0).
- `--missing-tss-policy drop_missing` removes a TSS if any candidate lacks a resolved contact to it (default `strict` keeps it).
- The retained original weight must reach `--minimum-retained-tss-weight` (default 0.9); otherwise the gene is unscoreable.

```math
\pi_{\mathrm{used}}(t\mid G)=\frac{\pi(t\mid G)}{\sum_{u\in\mathcal T_{\mathrm{keep}}(G)}\pi(u\mid G)},\qquad t\in\mathcal T_{\mathrm{keep}}(G)
```

The candidate catalog is unchanged. Scores of filtered genes carry `tss_contact_scope=selected_tss_set`.

## Contact policies and priors

```math
P(d)=a\left[\frac{\max(d,d_{\min})}{d_{\mathrm{ref}}}\right]^{-\gamma},\qquad
p(d)=\kappa\min\{P(d),P(d_0)\}
```

A prior fitted by `pace fit-prior` with CTCF peaks and the genome FASTA also carries the boundary term, $P(d)\,e^{-\beta\min(S,1)}$ ([Fit from a contact map](#fit-from-a-contact-map)); every use of $P(d)$ below then includes it.

| Situation | Contact used |
|---|---|
| Same bin or configured near-distance range, prior available | $P(d)$ |
| Same bin, no prior, default `prior_or_neighbor` | Maximum of the valid neighbouring contacts |
| Near-diagonal pair with no replacement | NA |
| Finite off-diagonal observation, including zero | Observed contact + pseudocount $p(d)$ |
| Missing observation, `allow_prior_fallback` and a prior | $P(d)$, labelled as prior |
| Missing observation without fallback | NA |
| `prior_only` mode | $P(d)$ at the element–TSS distance, $d_{\min}$ at most 5 kb |
| `shrinkage` mode | $r\,C_{\mathrm{obs}} + (1-r)\,P(d)$ |

Near-diagonal handling runs first; same-bin membership follows the bin boundaries. `near_diagonal_policy: unresolved` leaves these pairs NA, `prior_or_unresolved` accepts only a prior, and the default `prior_or_neighbor` also accepts the neighbour maximum.

`contact.pseudocount: auto` adds $p(d)$ to finite off-diagonal observations in observed mode when a prior is available; `powerlaw` requires the prior; `none` disables it. Defaults: $\kappa=1$, $d_0=5$ kb. The pseudocount is not added to near-diagonal replacements, prior-only contacts or shrinkage.

Outputs keep `observed_value`, `prior_value`, `pseudocount_value`, the evidence type, the prior identity and the processing reason.

### Fitting and transferring a prior

`pace fit-prior --cooler ...` estimates the decay from distance-class means over all callable bin pairs, including zero pixels. Zero-mean classes are reported and left out of the log fit. Optional `--test-chromosomes` are held out and their decay fit is reported. A prior combined with observed contacts must match their scale, normalization, resolution, balancing and window.

`pace predict --prior` accepts a prior fitted on another tissue of the same species and assembly. With `pace run`, add `--allow-cross-context-prior`. Source and target contexts are recorded; the `validated` profile does not accept transferred priors.

A prior fitted on bins is evaluated at bin centres when it is combined with measured contacts (observed and shrinkage modes), on the scale of those bins. In `prior_only` mode it replaces Hi-C and is evaluated at the element and TSS, with distances floored at $\min(d_{\min}, 5\,\mathrm{kb})$; `prior_coordinate_policy` records which applies.

The distance-only estimate (`pace predict --distance-only`, or `pace run --prior-preset abc_human --contact-mode prior_only`) sets $\gamma=1$, $a=1$ and $d_{\mathrm{ref}}=d_{\min}=5$ kb, so contact is $1/\max(d,5\,\mathrm{kb})$. It runs in the research profile and cannot be combined with measured contact tables.

### Fit from a contact map

```bash
pace fit-prior --cooler data/animal1.mcool::/resolutions/5000 \
  --ctcf data/animal1_CTCF.narrowPeak --genome data/GRCg7w.fa \
  --species chicken --assembly GRCg7w --tissue liver \
  --scale balanced_contact_protocol_1 --normalization-id hic_norm_protocol_1 \
  --test-chromosomes chr2 --out models/animal1_contact
```

The fit uses 30 logarithmic distance classes from the bin size (or `--min-distance`, not below it) to 5 Mb. `--no-balanced` uses raw counts. Fixed-bin, symmetric-upper coolers are required.

Boundary positions come from CTCF and their effect from the Hi-C. With `--ctcf peaks.bed` (CTCF ChIP-seq or CUT&Tag peaks, BED or narrowPeak) and `--genome genome.fa`, the strongest CTCF motif (JASPAR MA0139.1) in each peak orients the site, and its strength is the peak signal scaled by its 95th percentile (the relative motif score if the peaks carry no signal). A boundary lies midway between a reverse site and the next forward site within 1 Mb, with the smaller of the two strengths. One of the two files alone is an error; with neither, the prior has no boundary term. The prior is

```math
P(d,S)=a\left(\frac{\max(d,d_{\min})}{d_{\mathrm{ref}}}\right)^{-\gamma}e^{-\beta\min(S,1)},
```

with $d_{\mathrm{ref}}$ = `--reference-distance` (5 kb), $d_{\min}$ the larger of the bin size and `--minimum-distance` (5 kb), and $S$ the summed boundary strength strictly between the two bin centres. $\beta$ is estimated from pairs up to 1 Mb with one free level per distance class, comparing pairs at the same distance; $\gamma$ is then fitted on the boundary-corrected distance levels. A nonpositive $\beta$ is set to 0 (`beta_status: no_boundary_effect`). `--no-boundaries` fits the distance decay only.

Outputs: `manifest.json`, `distance_bins.tsv`, `fit_report.json`, and with boundaries `boundaries.tsv` and `boundary_cells.tsv` (pairs and contact by distance and boundary class).

Contacts extracted from the same map use the same labels:

```bash
pace contacts -i data/animal1.mcool -r 5000 -d prepared/catalog \
  --scale balanced_contact_protocol_1 --normalization-id hic_norm_protocol_1 -o prepared/animal1_hic
pace run ... --contacts prepared/animal1_hic/observed_contacts.tsv \
  --contact-prior models/animal1_contact --allow-prior-fallback
```

`pace predict --hic` performs both steps.

### Fit from a prepared contact table

```bash
pace fit-contact-prior --contacts contact_fit.tsv --model-id contact_prior_research \
  --species chicken --assembly GRCg7w --tissue liver \
  --scale balanced_contact_protocol_1 --resolution 5000 \
  --normalization-id hic_norm_protocol_1 --balancing balanced \
  --bin-edges 5000 10000 20000 50000 100000 500000 5000000 \
  --d-ref 10000 --d-min 5000 -o models/prior
```

The table has `bin_pair_id`, `distance_bp`, `contact_value`, `split`, `region_id`, with each measured pair once, zeros included. Training and test regions do not overlap. Positive-mean distance bins enter a log-linear fit weighted by pair count; empty and zero-mean bins are reported. A nondecreasing fitted curve is rejected.

### Transfer between tissues

```bash
pace fit-prior --cooler reference_liver.mcool::/resolutions/25000 \
  --ctcf muscle_CTCF.narrowPeak --genome GRCg7w.fa \
  --species chicken --assembly GRCg7w --tissue liver --out priors/liver

pace predict ... --species chicken --assembly GRCg7w --tissue muscle --prior priors/liver
# or
pace run -d prepared/catalog --activity prepared/muscle_activity.tsv \
  --species chicken --assembly GRCg7w --tissue muscle \
  --contact-mode prior_only --contact-prior priors/liver --allow-cross-context-prior -o results/muscle
```

For a 25 kb matrix the lower fitting distance is 25 kb. Keep the prior's source tissue name; the output records source and target contexts. To combine a prior with measured target contacts, their normalization, scale, resolution and windows must match.

## Boundary priors and weak calibration

### Contact prior from a boundary table

```math
C_0(E,t)=a\left(\frac{\max(d(E,t),d_{\min})}{d_{\mathrm{ref}}}\right)^{-\gamma}
\exp\!\left[-\beta\sum_{b\;\mathrm{between}\;E,t}s_b\right]
```

Boundaries strictly between the two anchors contribute. Priors fitted by `pace fit-prior` cap the summed strength at 1 (`boundary_cap`); priors packaged by `pace prior` use the uncapped sum. At $\beta=0$ the prior is the distance decay.

`pace boundaries` builds a boundary table from motif hits (optionally restricted to CTCF ChIP peaks):

```bash
pace boundaries --motifs ctcf_motifs.tsv --chip ctcf_peaks.tsv --max-gap 1000000 -o boundary_asset
```

or from a FASTA scan with a supplied matrix (12.0 is an example log2-odds threshold):

```bash
pace boundaries --fasta reference.fa --pwm ctcf_pwm.tsv --threshold 12.0 --max-gap 1000000 -o boundary_asset
```

Motif tables have `chrom`, `start`, `end`, `strand` (+ or −) and `strength` (0–1), zero-based half-open. PWM tables have one row per position with columns `A`, `C`, `G`, `T` (counts or probabilities); 0.01 is added to each entry before row normalization, and scores are log2 odds against a uniform background on both strands. ChIP tables have `chrom`, `start`, `end`; only overlapping motifs are kept. Overlapping hits collapse to the strongest; each adjacent reverse-then-forward pair within `max_gap_bp` gives a boundary at the midpoint of the gap, with the smaller of the two strengths.

`pace prior` packages supplied parameters (`--kappa` and `--pairs` are optional):

```bash
pace prior --model-id pig_liver_boundary_prior \
  --species "Sus scrofa" --assembly Sscrofa11.1 --tissue liver --target-level population_mean \
  --scale cooler_native --resolution 10000 --normalization-id liver_cooler_v1 --balancing balanced \
  --boundaries boundary_asset/boundaries.tsv \
  --a 1.0 --gamma 1.0 --beta 0.0 --d-ref 5000 --d-min 10000 -o prior_asset
```

The amplitude cancels in prior-only scores; for shrinkage it must be on the scale of the Hi-C. The boundary table is copied into the asset and bound by checksum; a positive $\beta$ requires a table.

### Fitting raw Hi-C and estimating reliability

`pace fit-hic` takes the same metadata, `--boundaries`, `--d-ref` and `--d-min`, with grids in place of `--a`, `--gamma`, `--beta` and `--kappa`. Measurement labels are read from the contact table when not given.

```bash
pace fit-hic --contacts prepared/pig1_hic/observed_contacts.tsv \
  --gamma-grid 0.6 0.8 1.0 1.2 --beta-grid 0 0.5 1 2 --test-chromosomes 18 \
  --boundaries boundary_asset/boundaries.tsv --d-ref 5000 --d-min 10000 \
  --species "Sus scrofa" --assembly Sscrofa11.1 --tissue liver -o fitted_prior
```

Input columns: `chrom`, `anchor0`, `tss0`, `sample_id`, `bin_pair_id`, `resolution`, `measurement_status`, `raw_count`, `count_to_contact`, `contact_value`, `scale`, `normalization_id`, `balancing`, `window_id` (exported by `pace contacts`). Raw counts are integers and `raw_count × count_to_contact` reproduces the observed value: the product of the two balancing weights for balanced coolers, 1 for unbalanced ones. Measured zero pixels are kept; missing or masked rows, diagonal pairs and repeated sample/bin pairs are excluded; pairs are placed at bin centres. At least three callable pairs with some positive counts are needed.

For each $\gamma$/$\beta$ grid point the amplitude is profiled analytically under a Poisson likelihood; the best point defines the expected counts $\mu_i$. The Gamma shape is the moment estimate, bounded to [10⁻⁶, 10⁶]:

```math
\widehat{\kappa}^{-1}=\frac{\sum_i[(y_i-\mu_i)^2-y_i]}{\sum_i\mu_i^2}
```

Nonpositive excess variation uses the upper bound. Held-out chromosomes enter neither the fit nor $\kappa$; `held_out_residuals.tsv` reports their observed and expected counts, and the report gives the correlation between log distance and boundary strength.

`pace predict --hic` uses this shrinkage by default (`--no-shrinkage` turns it off); with `--chunk-pairs`, $\kappa$ is first fitted on every chromosome and then fixed for all chunks. `pace fuse` takes the options of `pace run` and sets `--contact-mode shrinkage --contact-reliability per_pair`:

```bash
pace fuse -d prepared/catalog --activity prepared/activity/observed_activity.tsv \
  --contacts prepared/pig1_hic/observed_contacts.tsv --contact-prior fitted_prior \
  --species "Sus scrofa" --assembly Sscrofa11.1 --tissue liver -o fused
```

With conversion factor $f$ and prior contact $C_0$:

```math
C\sim\mathrm{Gamma}(\kappa,\ \mathrm{rate}=\kappa/C_0),\qquad
Y\mid C\sim\mathrm{Poisson}(C/f),\qquad
r=\frac{C_0/f}{C_0/f+\kappa},\qquad
\mathrm{E}(C\mid Y)=r\,fY+(1-r)\,C_0
```

Zero counts are observations and receive a positive posterior mean. Masked pairs use the prior with $r=0$. Diagonal pairs follow the near-diagonal policy. Each sample is shrunk separately and then aggregated over technical replicates, biological replicates and donors. `kappa: auto` takes the prior's $\kappa$, or estimates it from the run's unique callable bin pairs; a positive number overrides it. `resolved_contacts.tsv` records $\kappa$, its source and `posterior_samples` (per-sample means, variances, raw counts, factors and reliabilities). Imported per-pair tables are recomputed from their raw observations and must agree.

### Learning gamma, beta and eta from eQTLs

`pace fit-labels` takes fine-mapping results as a table with `variant_id`, `chrom`, `pos0`, `gene_id` and `pip` (plus `signal_id` for `independent_signals`), one row per variant and gene, PIP in [0, 1], on the run's assembly and tissue. Genes outside the candidate catalog are counted in `n_label_genes_outside_catalog`; duplicate gene/variant rows are rejected.

```bash
pace fit-labels --run results/pig_liver_baseline --eqtl eqtl_finemapping.tsv \
  --aggregation independent_signals \
  --gamma-grid 0.8 1.0 1.2 --beta-grid 0 0.5 1 --eta-grid 0 0.5 1 \
  --test-chromosomes 18 -o weak_fit
```

`--run` is a finished `pace run` folder, from which species, assembly and tissue are taken. Aggregation:

- `independent_variants`: element/gene mass is $1-\prod(1-\mathrm{PIP})$.
- `independent_signals`: PIPs of variants in the element are summed within each signal (total at most 1 per gene and signal), then combined as $1-\prod(1-\text{signal mass})$.

All candidates of genes in the table enter the evaluation; elements without PIP mass are `unlabelled_background`. The run must use `prior_only` or `shrinkage` with a fixed prior, `--eta auto` and no other calibrator; for per-pair shrinkage, $\kappa$ is fixed beforehand. At least three training chromosomes and one test chromosome are required. Parameters are selected by nested leave-one-chromosome-out on the training chromosomes, refitted on all of them, and evaluated once on the test chromosomes. Ties prefer $\eta=0$ and $\beta=0$, which must be in the grids; at most 50 values per parameter and 2,000 combinations are accepted. The metric is tie-aware fractional average precision of PIP mass, averaged over chromosomes; baselines are the nearest TSS and a power-law score with the run's activity, catalog and TSS weights.

```bash
pace run <options of the baseline run> --weak-model weak_fit/weak_model.json -o calibrated
```

Output carries `allocation_evidence: eqtl_weak`, and `eta_calibration.json` records chromosomes, source hashes and the held-out result. `--weak-model` and `--eta-model` are mutually exclusive.

## Incomplete scores and support bounds

With $S = A_\star\overline C$ (or the allocation-extended support),

```math
\mathrm{PACE}_{\mathrm{conditional}}(E,G)=\frac{S(E,G)}{\sum_{e\in\mathcal E^{\mathrm{score}}(G)}S(e,G)}
```

`--partial-policy conditional` writes this into `pace_score`; `score_scope` marks it conditional. Supports 2 and 1 with a third unknown give a conditional first score of 2/3; an unseen support of 97 would make its full score 0.02.

For support ranges $L_i \le S_i \le U_i$:

```math
\mathrm{PACE}_{i,\mathrm{lo}}=\frac{L_i}{L_i+\sum_{j\ne i}U_j},\qquad
\mathrm{PACE}_{i,\mathrm{hi}}=\frac{U_i}{U_i+\sum_{j\ne i}L_j}
```

Resolved supports are fixed; unresolved supports default to [0, ∞). For supports (2, 1, unknown) the first score lies in [0, 2/3]; bounding the unknown to [2, 4] narrows it to [2/7, 2/5]. A zero total gives NA. Bounds for unresolved edges are supplied as `inputs.support_bounds`:

```text
element_id	gene_id	support_lower	support_upper	bound_source
```

They apply to the final unnormalized support, including any allocation exponent; an NA upper bound is unbounded.

## Additional omics

Data added with `--expression`, `--features` or `--methylation` are annotations and do not change `pace_score`. (CTCF peaks given to `pace fit-prior`/`pace predict --hic` with `--ctcf --genome` place the domain boundaries of the contact prior.)

### Histone and CTCF tracks

`pace features` records the overlap of each element with peaks of other marks:

```bash
pace features -b data/H3K4me1_peaks.bed -d prepared/catalog --feature-prefix H3K4me1 -o prepared/H3K4me1
pace features -b data/CTCF_motifs.bed -d prepared/catalog --source-id CTCF_motif_source \
  --evidence-id CTCF_motif_overlap --feature-prefix CTCF --motif-strands -o prepared/CTCF
```

`--motif-strands` needs strand-bearing motif intervals. Add the tables to `pace run` with `--features`. Signal from other marks can also be written in the activity-table format, merged with `pace merge -t observed_activity` and registered in `samples.tsv`; with `--panel` unchanged they remain annotations. Promoter windows (e.g. H3K4me3 over ±2 kb) need their own window identifier; gene-level promoter summaries use the run's TSS weights.

### RNA

`pace run --expression` takes `gene_id`, `sample_id`, `tpm`, `status`, for catalog genes and samples registered as RNA. Transcript TPM is summed to genes through the catalog's mapping:

```bash
pace expression --tpm data/transcript_tpm.tsv -d prepared/catalog -o prepared/rna
```

### CpG methylation

`pace run --methylation` takes raw dyad counts:

```text
chrom	dyad_start0	methylated_count	total_count	sample_id	assay
```

`dyad_start0` is the zero-based C of the forward-strand CpG; complementary strands are merged first (the `merge_stranded_cpg` adapter checks them against the reference FASTA).

```bash
pace run ... --methylation data/cpg_dyad_counts.tsv \
  --methylation-min-coverage 5 --promoter-upstream 2000 --promoter-downstream 500 \
  --reference-cpg data/reference_cpg_counts.tsv
```

Minimum coverage defaults to 1. Promoter windows follow the strand: `[tss0−upstream, tss0+downstream+1)` on the plus strand, sides exchanged on the minus strand, clipped at zero. The reference table has `entity_type` (`element` or `promoter`), `entity_id` and `n_cpg`, the number of reference CpGs in the window; without it the coverage fraction is NA. WGBS and RRBS are kept separate.

| Field | Calculation |
|---|---|
| `M_site` | Mean CpG methylation fraction after the coverage filter |
| `M_pooled` | Methylated counts / total counts |
| `covered_cpg` | CpGs meeting the coverage threshold |
| `cpg_coverage_fraction` | Covered CpGs / reference CpGs in the window |
| Status | Measured zero, no CpG, no coverage or low coverage |

`pace methylation --counts ... -d prepared/catalog --min-coverage 5 --reference-cpg ... -o prepared/methylation_qc` writes a per-element `methylation_summary.tsv` for inspection.

### Other feature tables

`pace run --features` accepts

```text
entity_type	entity_id	feature_name	value	evidence_id	status
```

with types `element`, `promoter`, `gene` or `edge` (ID `element_id|gene_id`), one row per entity and feature, and optional `assay`, `unit`, `normalization_id`, `window_id`. Unknown, invalid and unmeasured values are missing. A trained classifier is applied with `pace run ... --ml-model models/liver_classifier`.

## Experimental target allocation

```math
B(E,G)=\frac{\overline C(E,G)}{\sum_{H\in\mathcal G(E)}\overline C(E,H)},\qquad
\mathrm{PACE}_{\eta}(E,G)=\frac{A_\star(E)\,\overline C(E,G)\,B(E,G)^\eta}{\sum_{e\in\mathcal E(G)}A_\star(e)\,\overline C(e,G)\,B(e,G)^\eta}
```

$\eta=0$ gives the main formula; `--eta auto` (default) is 0 unless a calibration source is given. At $\eta>0$, a missing contact to any candidate gene leaves $B$ unresolved for that element. Since $\overline C\,B^\eta=\overline C^{\,1+\eta}/(\sum_H\overline C)^\eta$, the term sharpens contact contrasts and depends on gene annotation density; `pace benchmark` compares it with $\eta=0$ and with `contact_power_2` (A×C̄²).

```bash
pace fit-eta -d prepared/catalog --activity ... --contacts ... \
  --species pig --assembly Sscrofa11.1 --tissue liver --eta-labels functional_labels.tsv -o results/eta
```

The same fit runs inside `pace run --eta-labels`; `--eta-model results/eta/eta_calibration.json` reuses it. A number in [0, 1] given with `--eta` is used as is.

### Calibration objective

Functional inhibition or deletion experiments provide positives and powered negatives within each gene. With $a=\log(A_\star\overline C)$ and $b=\log B$, the gene denominator cancels in a positive–negative pair, and

```math
\widehat\eta=\arg\min_{0\le\eta\le1}
\frac1{|\mathcal G_{fit}|}\sum_{G\in\mathcal G_{fit}}
\frac1{|P_G||N_G|}\sum_{p\in P_G,n\in N_G}
\log\left[1+\exp\{-[(a_p-a_n)+\eta(b_p-b_n)]\}\right]
```

Genes and, within a gene, positive–negative pairs have equal weight. The objective is convex; endpoint derivative checks and 40 bisection iterations locate the minimum on [0, 1].

### Required functional-label table

| Columns | Requirements |
|---|---|
| `label_id` | Unique label |
| `element_id`, `gene_id` | One canonical element–gene edge |
| `species`, `assembly`, `context_id`, `target_level` | Match the scored context |
| `perturbation_type` | `CRISPRi`, `deletion` or `enhancer_inhibition` (`synthetic_inhibition` in the demonstration profile) |
| `label_status`, `effect_direction` | `enhancing_positive` with `down`, or `powered_negative` with `none` |
| `split` | `train`, `calibration` or `test` |
| `group_id` | Independent locus or experimental group; no group crosses splits |
| `assay_id`, `source_id` | Experimental assay and source |

Genes, elements and groups cannot cross splits; duplicates, missing fields and leaking splits are errors. Context mismatches, unsupported perturbations, uncertain or repressive labels and zero or unresolved support are excluded with recorded reasons, and test rows are set aside before fitting. Each fitting gene needs an eligible positive, an eligible negative and a contrast in $\log B$ above 10⁻¹²; with fewer than `minimum_genes` (default 3, at least 2) informative genes, $\eta=0$ is returned with the reason.

### Deployment validation after candidate fitting

Records sharing a group, gene or element form independent components; at least `minimum_groups` (default 3) are required, split into up to `validation_folds` (default 5) grouped folds. Each fold refits the candidate on its training part, and 21 shrinkage multipliers from 0 to 1 (step 0.05) are scored by held-out per-gene average precision. The smallest multiplier within one standard error of the best is chosen; it is accepted when the paired gain over $\eta=0$ exceeds its standard error and is positive in at least 80% of held-out groups. The candidate is then refitted on all fitting labels and multiplied by the chosen shrinkage. Calibration rows confirm a candidate fitted on train rows; test rows never enter fitting, selection or confirmation. Otherwise $\eta=0$ is returned with reasons.

Artifacts (schema `pace-eta-2`) record `candidate_eta`, the applied `eta`, `cross_validation`, `estimation_label_ids`, `confirmation_label_ids` and `deployment_validated`.

### Options and reuse

```bash
pace run ... --eta-labels functional_labels.tsv --eta-min-genes 3 --eta-min-groups 3 --eta-validation-folds 5
```

`--eta-labels` and `--eta-model` are mutually exclusive; an explicit `--eta NUMBER` uses neither. Every run writes `eta_calibration.json`, and the run manifest stores the exponent in `comparison_contract.eta`. A saved calibration is reused only when context, target level, panel and scales, catalog, candidates, promoters, contact policy and asset hashes match. For between-animal comparisons, fit once and reuse the same artifact. `pace benchmark` reports calibrated PACE next to the fixed baselines and rejects evaluation labels that share genes, elements or groups with the fit.

## Classifier training and inference

```bash
pace train --data labelled_edge_features.tsv --model-id regulatory_link_classifier \
  --species chicken --assembly GRCg7w --tissue liver \
  --extra-features H3K4me1 promoter_H3K4me3 methylation_M_site \
  --penalty 0.01 0.01 --penalty 0.1 0.01 --folds 3 --seed 17 --calibrate \
  --feature-contract results/measured/ml_feature_contract.json -o models/classifier
```

`--feature-contract` takes the `ml_feature_contract.json` of a compatible scoring run (activity panel and scales, contact resolution and normalization, candidates, target level and auxiliary preprocessing). Required columns: `element_id`, `gene_id`, `assayed_region_id`, `mapping_count`, `group_id`, `split`, `label_status`, `effect_direction`, `A_used`, `Cbar`, `distance_bp`, `pace_score`, `regime`, `activity_sources`, `contact_sources`, the extra features and optional `sample_weight`. Positives are `enhancing_positive`/`down`, negatives `powered_negative`; upregulated, low-power, untested and multi-tile labels are excluded and exported.

Groups, edges and perturbed regions stay within one split and one tuning fold; preprocessing is refitted inside each training fold. Only extra features can be imputed; medians, IQRs, missing indicators and dropped columns are stored in JSON. The objective is sample-weighted logistic loss + λ₁·L1 + λ₂·L2²/2 with an unpenalized intercept, solved by proximal gradient. Optional sigmoid calibration uses the calibration split (slope penalty 10⁻⁶); without it `pace_ml_probability` is NA. A base-only classifier is saved and evaluated alongside.

```bash
pace predict-ml --run results/measured --model models/classifier \
  --features results/measured/multiomics_features.tsv.gz -o results/ml
```

Inference checks the feature definition, regime, evidence sources and context; mismatches return `out_of_scope`. Saved models are applied without refitting. The test report evaluates `pace_ml_score` (`test`) and the calibrated probability (`test_probability`) separately and counts unavailable predictions in coverage.

## Comparisons and benchmarks

### Compare two runs

```bash
pace compare --left results/animal_a --right results/animal_b --min-common-units 2 -o results/comparison
```

Both runs are rescored on their common measurable units. Catalog, candidate and promoter hashes, target level, context, panel, units, normalization and formula parameters must match; at least two common units and positive common denominators are needed. The full difference also requires complete backgrounds; otherwise only the conditional difference is reported. `--allow-eta-difference` and `--allow-evidence-difference` permit different exponents or contact policies and report conditional differences. The table holds original and common scores, full and conditional differences, activity and support differences, gene totals and reasons, all right minus left.

### Replicate stability

```bash
pace stability --replicates replicate_runs.tsv -o results/stability
```

The table has `run_path`, `donor_id`, `replicate_type` (`biological` or `technical`); paths are relative to it. Every pair is compared on recomputed common denominators; Pearson correlation is NA for constant or short vectors, and distinct donors are counted.

### Functional-label benchmark

```bash
pace benchmark --run results/measured --labels functional_labels.tsv \
  --membership prepared/catalog/region_membership.tsv --stratify gene_id \
  --external gABC gabc_scores.tsv recorded_external_version recorded_external_configuration \
  -o results/benchmark
```

Built-in methods are negative distance, `ABC_style_single_TSS` (the TSS with the smallest coordinate, then promoter ID), PACE with $\eta$ = 0 and 1, and `contact_power_2`; a calibration source adds `PACE_calibrated`, a fixed exponent `PACE_fixed_eta`. Formula methods are also renormalized on their common scoreable sets. External scores (`--external`, repeatable) need version and configuration and columns `element_id`, `gene_id`, `score`. A threshold is given only with `--threshold METHOD VALUE SOURCE_SPLIT SOURCE_ID`. Labels use the [label table](../README.md#reference); only one-to-one region mappings are evaluated, and every tested positive stays in the evaluation set. Average precision is the tie-aware sum of recall increments times precision; AUROC needs both classes. Coverage, missing positives, candidate recall and end-to-end recall are reported separately, and missing scores are not set to zero.

## Worked calculations

These cases are computed by hand and checked in [scoring tests](../tests/test_core.py) and [calibration tests](../tests/test_allocation.py); the main formula is in the [README](../README.md#the-model).

**Activity and TSS contact.** ATAC = 4 and H3K27ac = 9 give activity √(4·9) = 6; measured H3K27ac = 0 gives 0; unavailable H3K27ac gives NA. Contacts (2, 6) with weights (3/4, 1/4) give C̄ = 3.

**Allocation.** Activities (4, 2, 1), G1 contacts (3, 2, 2) and G2 contacts (1, 2, 6) give B for G1 = (3/4, 1/2, 1/4) and for G2 = (1/4, 1/2, 3/4).

| Exponent | G1 support | G1 PACE | G2 PACE |
|---|---|---|---|
| 0 | (12, 4, 2) | (2/3, 2/9, 1/9) | (2/7, 2/7, 3/7) |
| 1 | (9, 2, 1/2) | (18/23, 4/23, 1/23) | (2/15, 4/15, 9/15) |

For $\eta=1/2$, G1 support is (6√3, 2√2, 1) and G2 support (2, 2√2, 3√3), each divided by its sum.

**Bulk aggregation.** Two replicates of E1 with assay pairs (9, 1) and (1, 9) and two of E2 with (4, 4) average to activities 5 and 4, so with equal contact E1 receives 5/9. Summing per-replicate supports instead gives 3/7, a different quantity.

**Partial genes.** Supports (1, 1, 2) and (1, 1, NA) give conditional E1 shares 1/4 and 1/2; the second primary score is withheld, with range [0, 1/2]. Changing supports (1, 1) to (1, 2) lowers the first share from 1/2 to 1/3 with its own support unchanged.

**Calibration.** Two equally weighted margins η − 1/2 and 1/2 − η give the loss [softplus(1/2 − η) + softplus(η − 1/2)]/2, minimized at η = 1/2 with value log 2.

## Evaluation design

A protocol for measuring how PACE performs as data become sparse:

1. **Reference.** Experimentally tested K562 enhancer–gene pairs from the [EngreitzLab CRISPR benchmark](https://github.com/EngreitzLab/CRISPR_comparison), with release, assembly, assay and effect and power definitions recorded. Deduplicate regions, genes and assays shared by Fulco and Gasperini compilations. Negatives are powered non-effects.
2. **Splits.** Chromosomes or connected region/gene groups are assigned to training, calibration and test before any tuning.
3. **Degradation.** Predeclared Hi-C depths (e.g. 100, 50, 25, 10, 5%) by binomial thinning of raw pairs, with the map rebuilt, rebalanced and the prior refitted at each depth; removal of one activity assay; fewer donors; the candidate universe fixed.
4. **Methods.** The original ABC implementation (exported with `--external`), `ABC_style_single_TSS`, PACE with $\eta=0$, without pseudocount, without near-diagonal correction, prior-only, with prior fallback, single versus multiple TSSs, equal versus activity-weighted TSSs, and allocation versus `contact_power_2`.
5. **Metrics.** Average precision on a common evaluable set together with coverage of all tested pairs, recall at a fixed precision, and paired bootstrap over gene or region groups; results stratified by distance, promoter count, gene density and chromosome.
6. **Livestock.** Species- and tissue-matched fine-mapped eQTL (e.g. FarmGTEx), with controls matched for distance, allele frequency, LD, expression and candidate density, and perturbation experiments in chicken, pig or cattle where available.

## Complete configuration defaults

These defaults match `pace_livestock.config.DEFAULTS`; the same structure is written to `resolved_config.yaml` in every result folder (`null` means not provided).

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
  weights: usage
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
