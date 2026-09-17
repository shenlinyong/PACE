# PACE: current Methods

Mathematical symbols follow [the unified notation](NOTATION.md).

These methods describe the current model and the associated analysis profile. See [the input contract](INPUTS.md) for software fields and [validation scope](../VALIDATION.md) for the distinction between software tests and biological evaluation.


### M1. Study resources and input preparation

We audited the original Supplementary Tables S1 and S2 against public metadata and validated 334 livestock source files by checksum and gzip integrity. The revised tables retain the original entries alongside accession-level corrections. Of 385 S2 accession occurrences, 223 disagreed with GEO species or tissue labels, 15 had assay-only discrepancies identifying DNase-seq rather than the reported ATAC-seq assay, and 51 fell outside GSE158430. The 15 assay discrepancies had concordant species and tissue labels; 96 occurrences matched all three fields. Analysis sample identities followed GEO records [1].

Livestock analyses used Ensembl release 99 annotations on GRCg6a, ARS-UCD1.2 and Sscrofa11.1. Accessibility and H3K27ac summits were each expanded to 500-bp windows; overlapping windows from available replicate pairs were retained and merged within assays. Accessibility windows overlapping the supported H3K27ac windows were combined with 500-bp windows around every distinct eligible TSS. Overlapping intervals were merged. Gene eligibility required RNA counts per million of at least one in both available source replicates. For each assay and candidate, the feature was the maximum overlapping called-peak enrichment, averaged across available replicates. Absence of a called peak was encoded as zero in this detection-limited profile; it does not establish zero biological activity. These features were not continuous read-coverage measurements. Individual replicate availability is documented in the run manifests.

Human analyses used GRCh38 with GENCODE v29 and cell-matched ENCODE quantitative signal tracks. Accessibility defined candidate summits; blacklist regions were removed and expressed-gene promoter windows were added before scoring. K562 used DNase-seq and H3K27ac; GM12878 used ATAC-seq and H3K27ac. Gene eligibility followed the public K562 expressed-gene resource or GM12878 TPM of at least one in both RNA replicates. Signal was integrated over each region and divided by its length. Unreported bases on available chromosomes were zero; unavailable chromosomes remained missing. The full accession, reference and input-profile inventory accompanies Supplementary Tables S1–S2.

### M2. Scoring definition and version change

The intended previous expression-weighted score was

$$
\mathit{PACE}(E,G)=W_{\mathrm{expr}}(G)
\frac{A(E)C(E,G)}{\sum_{e\in\mathcal E(G)}A(e)C(e,G)}. \qquad (1)
$$

When positive, this expression weight rescales every candidate link to the same gene equally. It can affect pooled cross-gene ranking and threshold selection, but cannot change the ordering of enhancers within a gene. The previous documentation and prediction interfaces did not consistently represent this ordering of operations. PACE uses one shared scoring kernel and removes expression weighting from the primary score. Historical implementation differences and their implications for the previous ablation are documented in Supplementary Methods SM1.

The revised score is

$$
\mathit{PACE}(E,G)=
\frac{A(E)C(E,G)B(E,G)^{\eta}}
{\sum_{e\in\mathcal E^{\mathrm{obs}}(G)}A(e)C(e,G)B(e,G)^{\eta}+U(G)}. \qquad (2)
$$

The score $\mathit{PACE}(E,G)$ represents PACE support for enhancer $E$ and gene $G$ (Fig. 1a). Here, $A(E)$ is enhancer activity, $C(E,G)$ is TSS-aggregated contact, $B(E,G)$ is enhancer-centric target allocation, and $\mathcal E^{\mathrm{obs}}(G)$ contains available candidates with finite raw support for gene $G$. Candidates with undefined support remain in the output but do not contribute a numerical term to this sum. The optional nonnegative quantity $U(G)$ represents independently estimated residual support not included in that sum, in the same units as its terms. No estimator of $U(G)$ is implemented in this revision. In the primary configuration it is zero for computation and explicitly marked unknown. Thus, the score remains conditional on the observed candidate set and is not a calibrated probability, causal effect size or false discovery rate.

### M3. Candidate elements and promoter annotation

Enhancers use zero-based, half-open genomic intervals. Candidate pairs are generated before score filtering when the absolute distance from the interval midpoint, (start + end)/2, to a TSS is strictly less than the configurable cis window, initially 5 Mb. Distinct transcript TSSs are retained and shared TSSs are deduplicated. The negative-strand TSS is the GTF end coordinate minus one. Gene-boundary fallback is labelled when transcript records are unavailable. Stable gene identifiers, rather than gene symbols, define normalization groups.

TSS-use weights $\pi(G,t)$ are normalized over the full available set of distinct TSSs for each gene before window filtering. Independent promoter-use measurements can supply these weights; otherwise, equal weights are used. Out-of-window TSSs do not trigger per-enhancer renormalization. Conflicting duplicate weights or quality annotations are rejected. The method does not discover genes or promoters absent from the input annotation.

### M4. Missing-aware activity integration

Nonnegative assay signals from the declared preprocessing profile undergo fixed assay-specific scaling, $x(E,i)=S(E,i)/a_i$, where $a_i>0$. Each scale was the median positive signal across the complete unique candidate catalogue for that assay and sample, calculated without functional labels. With observation indicator $m(E,i)$, quality weight $q(E,i)$ and assay prior $w_i$, activity is

$$
A(E)=\left\{\exp\left[
\frac{\sum_i m(E,i)q(E,i)w_i\log(1+x(E,i))}
{\sum_i m(E,i)q(E,i)w_i}
\right]-1\right\}\exp(-\kappa I(E)). \qquad (3)
$$

Missing measurements are omitted from the mean; measured zeros remain observations. No effective observation produces an undefined activity, retained as missing. Unknown quality uses a unit weight only for provisional ranking and remains unknown in the evidence output. Equal activity weights are declared design priors, not values learned from a GM12878 benchmark. The shifted geometric mean allows one zero-valued assay to coexist with positive aggregate activity, which may also increase false positives and therefore requires empirical evaluation.

The attenuation strength is $\kappa\geq0$, with $\kappa=0$ in all primary analyses, so $\exp(-\kappa I(E))=1$. If the optional inhibitory extension is explicitly used, $I(E)$ is a weighted mean of observed inhibitory values on a defined [0,1] scale. Methylation fractions are not min–max rescaled across regions. The core API supports this extension, whereas the primary raw-file adapters reject unspecified inhibitory inputs. No primary performance claim is attributed to WGBS or a repressive-mark module.

### M5. Contact estimation and provenance

A positive distance prior is defined as $P(d)=\mathrm{Scale}(d+\mathrm{Pseudocount})^{-\gamma}$. The current starting values, inherited for traceability, are $\mathrm{Scale}=5.9594510043736655$, $\mathrm{Pseudocount}=5000$ bp and $\gamma=1.024238616787792$. They are not described as independently optimized livestock parameters. A compatible observed contact is combined with this prior as

$$
C_{\mathrm{adj}}(E,t)=P(d(E,t))
\left[(1-\lambda(E,t))+\lambda(E,t)\frac{H(E,t)}{D(d(E,t))}\right]. \qquad (4)
$$

$H(E,t)$ and its positive distance expectation $D(d(E,t))$ must have the same sample, resolution and normalization. The expectation is estimated from eligible genomic bin pairs, without reference to functional labels. The reliability $\lambda(E,t)\in[0,1]$ describes measurement quality rather than the strength of the observed edge. Local marginal coverage, valid normalization and replicate agreement are potential inputs to a fixed reliability procedure; these values are not estimated automatically by the scoring kernel.

Missing observation, expected contact or reliability causes fallback to the prior with a reason code. A valid measured zero can reduce contact support and is distinct from an unobserved bin. Tissue-matched contacts, surrogate contacts and prior-only predictions have separate output states. The source and target tissues, species and assemblies must be recorded in the analysis manifest. A liver map used for a non-liver tissue is surrogate evidence. Across-species comparison of independently generated maps is not cross-species transfer of a contact matrix.

Gene-level contact is

$$
C(E,G)=\sum_{t\in\mathcal T(G)}\pi(G,t)C_{\mathrm{adj}}(E,t),
\qquad \sum_{t\in\mathcal T(G)}\pi(G,t)=1. \qquad (5)
$$

This weighted-average adaptation of the multiple-TSS principle avoids increasing contact support solely because more duplicate transcript records are present [3]. It does not establish tissue equivalence or correct missing promoters.

### M6. Allocation, evidence annotation and implementation

Enhancer-centric target allocation is defined as [3]

$$
B(E,G)=\frac{C(E,G)}{\sum_{g\in\mathcal G(E)}C(E,g)},
\qquad R(E,G)=A(E)C(E,G)B(E,G)^{\eta}. \qquad (6)
$$

The default is $\eta=1$. When all candidate contacts for an enhancer are zero, allocation is undefined; finite activity with zero contact is explicitly assigned zero raw support. Missing activity remains unknown. A zero gene-level denominator remains unscorable. With one TSS, $\eta=0$, $U(G)=0$ and identical supplied activity/contact values, Eq. (2) reduces to ABC normalization. This is a score-level relationship, not equivalence of preprocessing pipelines.

Evidence is reported separately as

$$
Q(E,G)=\min\{Q_A(E),Q_C(E,G),Q_T(G),Q_{\mathrm{cat}}(G)\}. \qquad (7)
$$

These components describe activity measurement quality, contact-observation support, TSS annotation quality and candidate-catalogue quality. Unknown components yield an unknown index. The initial sufficient-input-evidence rule requires all four components to be known and at least 0.5, with no unscored candidates for that gene. Unscorable links or explicitly zero activity/TSS quality are insufficient; other links are provisional. This operational rule is uncalibrated and does not identify functional true positives by itself.

All primary current predictions are formula-based. The legacy supervised module is archived and does not contribute to the primary model. Both prediction adapters delegate to `pace_core.py`. Outputs retain scores, unknown values, quality reasons and contact provenance before any optional filtering. Detailed field definitions, quality calculations and resampling semantics are provided in Supplementary Methods SM3–SM4.

### M7. Human CRISPRi benchmarking and comparator scope

We used the filtered GRCh38 experimental labels distributed with the public enhancer–gene benchmark [11]. Its K562 development table contained 10,356 tested pairs, including 471 regulated pairs. Context-matched external records comprised 1,918 K562 pairs and 68 GM12878 pairs, with 118 and 16 regulated pairs, respectively. Labels were used unchanged; no method-specific negatives were sampled. Other cells in the five-cell resource were excluded because this analysis lacked matching processed inputs. These external records had been inspected historically and do not constitute an untouched prospective test. The current formula, activity scales and candidate construction were not fitted to these labels.

Predictions were generated from the complete candidate and gene background before label matching. We selected the candidate with greatest genomic overlap with each tested interval, breaking ties by midpoint distance and lower start coordinate. Matching also required the same stable target-gene identifier. This mapping did not use prediction scores or experimental labels. Candidate omission and missing scores remained visible in the evaluation table.

ABC was the primary comparator [2]. We invoked the official scoring implementation at commit `92ac50360231a6bcfd654f3147839d078be445d9`, using the same processed activity features and its supported distance-contact profile. ABC used one TSS at the GTF gene boundary, unshifted geometric activity and the declared self-promoter rule. PACE retained all distinct TSSs and its declared shifted activity and distance prior. This was a controlled comparison of scoring on processed inputs; neither model used matched Hi-C in the primary human benchmark. It was not a comparison of complete native raw-read pipelines. Distance ranking, removal of target allocation and each single-assay PACE profile provided additional controls. EpiTensor, TargetFinder, JEME, EPIPDLF and GATv2EPI were methodological context [4–8]; they were not run in this revision. The current results therefore establish a restricted ABC comparison and do not rank PACE against these five methods.

The primary metric was end-to-end average precision (AP), computed as a stepwise precision–recall integral with tied scores grouped. All experimental positives remained in its recall denominator, including unscored pairs. Missing scores did not add a terminal prediction group. We reported scoring and positive coverage alongside AP; AUROC was conditional on scored pairs. Uncertainty used 2,000 paired target-gene bootstrap resamples with seed 20260915. The same gene formed one cluster across cells. Exploratory external subsets excluded development-region overlaps, then additionally excluded development genes. These restrictions cannot remove prior inspection history or all dependence between neighboring genes.

### M8. Livestock atlas and genetic support

All 24 primary atlases used distance-prior contact and retained evidence limitations. A score cutoff of 0.02 summarized distal pairs at least 2 kb from the gene-boundary TSS; it was descriptive and not a validated false-discovery threshold. All candidates, including promoter-proximal elements, remained in score denominators. Within-species tissue overlap was the Jaccard index of gene–interval-component sets. Overlapping or abutting retained intervals were joined into connected components across all eight tissues of each species; links shared a component and target gene. It did not test cross-species enhancer conservation.

Chicken-liver genetic support used significant ChickenGTEx associations [12], mapped with the tissue-matched genotype BIM coordinates on GRCg6a. A supported pair required both variant overlap with its candidate interval and agreement of the associated and predicted gene. Absence of eQTL support does not define a false positive. We ranked all distal candidates within each target gene, using average ranks for ties, and averaged supported-pair percentile ranks within genes. The paired endpoint was ABC rank minus PACE rank; positive values favored PACE. Recovery at the top 1%, 5% and 10% of candidates used the same supported-pair denominator.

We also rescored fixed, previously constructed matches retaining the same target gene and covariate-caliper requirements. Both interval coordinates and same-gene support membership had to remain consistent in the new atlas. One match whose comparison region also had current same-gene eQTL support was excluded without using prediction scores. Matching controlled distance, length, GC, central-50-mer mappability, accessibility, H3K27ac and local tested-SNP density. Same-gene matching held target expression fixed. We audited local LD separately; LD was not a matching variable. Supported-pair wins counted ties as one-half and were averaged within genes. Both endpoints used 2,000 gene bootstrap resamples. These comparisons concern association support, not functional precision, and do not fully remove cross-gene LD dependence. Supplementary Methods describe calipers, retained coverage and balance diagnostics.

### M9. Measured-contact sensitivity

We rescored chicken-liver candidates using the available FR-AgENCODE 40-kb raw-contact map. Source galGal5 bin midpoints were lifted to GRCg6a; nearest mapped anchors within 40 kb were used. Expected contact was the mean raw count across all cis bin pairs at the corresponding source-bin separation, including zeros. Residual coverage bias and coarse resolution remained. Assumed contact reliabilities of 0, 0.1, 0.25, 0.5 and 1 defined a fixed sensitivity grid. They were not measured quality estimates or optimized parameters. All TSS-level candidates were rescored through the shared kernel, and reliability zero reproduced the primary liver scores. The endpoint was paired change in gene-averaged eQTL-supported rank, with evaluable-gene counts recorded for each setting.

### M10. Published functional-locus reanalysis

We recovered reporter measurements and cloning primers from the published chicken abdominal-fat study [9]. Its reporter assays used the chicken preadipocyte line ICP2. Construct activity does not establish endogenous enhancer dependence in the source adipose tissue. Source spreadsheet cells, construct identities and available replicate counts were preserved. Reporter activity was calculated from raw Fluc/Rluc values relative to the within-assay vector mean. The source described nine biologically independent samples per group; one forward IGFBP2 G/In raw value was blank, leaving eight available values. We did not impute it.

The four forward promoter contrasts compared G/In with G/Del, holding the second variant background fixed. Effects were log2 ratios of arithmetic mean raw activities, with 2,000 independent-group bootstrap intervals. Two-sided Welch tests on log-transformed raw activities were adjusted across four genes by Holm's method. Cloning arms were removed from primers before exact matching to GRCg6a. The resulting 1,330-bp interval was compared with the frozen adipose catalogue without adding candidates after outcome inspection. It locates the primer-bounded region; the source reported a 1,300-bp construct on GRCg7b, so complete construct equivalence is not established. This was reanalysis of published experiments, with no new animal or cell experiments.
