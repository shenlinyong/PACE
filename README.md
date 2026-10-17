# PACE

**Enhancer–gene prediction for domestic animals from chromatin activity and 3D contact.**

PACE ranks the candidate regulatory elements of each gene in a tissue. It uses the Activity-by-Contact (ABC) score — an element's activity times its contact with the gene, as a share of the total over all candidates of that gene — and changes how activity and contact are estimated for the data that pig, cattle, sheep, goat and chicken projects usually have. Any species with a reference genome, a gene annotation and measured chromatin data can be analysed.

```bash
pace predict -b peaks.bed -g genes.gtf --atac atac.bw --h3k27ac h3k27ac.bw \
  --hic liver.mcool --hic-resolution 10000 --ctcf ctcf.narrowPeak --genome genome.fa \
  --species pig --assembly Sscrofa11.1 --tissue liver -t 4 -o pig1_liver
```

## What PACE improves

PACE keeps the ABC score and changes how its terms are estimated where domestic-animal data differ from human data.

**1. Sparse Hi-C.**

*Problem:* human cell lines have deeply sequenced, kilobase-resolution Hi-C for many cell types; in domestic animals Hi-C exists for few tissues, mostly at 10–40 kb and low depth, and most tissues have none.

*PACE:*

- *Low depth:* each element–TSS contact is the Gamma–Poisson posterior mean $\widetilde C = w\,C_{\mathrm{obs}} + (1-w)\,P(d)$ with $w = \frac{P(d)/f}{P(d)/f+\kappa}$, so pairs supported by few reads are drawn toward the distance expectation (default with `--hic`).
- *Coarse bins:* one bin pair holds element–TSS pairs at very different distances; each is placed within it by $P(d)/P(d_{\mathrm{bin}})$, the prior at its exact distance over the prior at the bin centres. An element in the TSS's own bin, whose diagonal value is dominated by self-ligation, takes $P(d)$ fitted on the same map; masked bins also take $P(d)$.
- *Hi-C only from another tissue:* $P(d)$ is fitted on it (Case B). *No Hi-C:* $\widetilde C = 1/\max(d,\,5\,\mathrm{kb})$ (Case C).

**2. Predicted TSSs and alternative promoters.**

*Problem:* gene annotations of domestic-animal genomes (Ensembl, NCBI RefSeq) are largely computational predictions. Many genes carry several TSSs or alternative promoters, few of them supported by CAGE or other TSS data, so the annotation does not say which promoter a tissue uses. ABC assigns each gene a single TSS; if it is the wrong one, the distance and Hi-C contact of every candidate element are measured from the wrong position.

*PACE:* keeps every distinct TSS of a gene (transcripts starting at the same base are merged) and weights each by its promoter activity in the same tissue, $\pi_t = A_t / \sum_{t'} A_{t'}$, where $A_t$ is the ATAC/H3K27ac activity of the 500 bp cell containing the TSS. The element–gene contact is the weighted mean, $\overline C(E,G) = \sum_t \pi_t\,\widetilde C(E,t)$: the promoter that is open and acetylated in this tissue dominates, and predicted but inactive TSSs contribute little.

**3. Domain boundaries.**

*Problem:* enhancer–promoter contact drops across CTCF-anchored domain boundaries, but domestic animals lack the boundary maps available for human cell types, and ABC has no boundary term: with distance-based contact, an element beyond a boundary counts as much as one at the same distance inside the gene's domain.

*PACE:* multiplies $P(d)$ by $e^{-\beta\min(S,1)}$, where $S$ is the strength of the boundaries between element and TSS. Boundaries lie between a reverse and the next forward CTCF site; sites are CTCF ChIP-seq or CUT&Tag peaks of the tissue, oriented by the built-in CTCF motif (`--ctcf peaks.bed --genome genome.fa`). $\beta$ is fitted on the Hi-C (Cases A and B) by comparing pairs at the same distance.

**4. One or a few animals.**

*Problem:* a tissue is usually profiled in one or a few animals, often with several libraries or technical replicates per animal.

*PACE:* averages technical replicates, then biological replicates, then animals, with equal weight per animal, so repeated sequencing of one animal is not counted as several animals.

---

## Contents

1. [Installation](#installation)
2. [Input data](#input-data)
3. [Running PACE: choose your case](#running-pace-choose-your-case)
4. [Reading the results](#reading-the-results)
5. [The model](#the-model)
6. [Step by step, several animals and comparisons](#step-by-step-several-animals-and-comparisons)
7. [RNA-seq, methylation and other data](#rna-seq-methylation-and-other-data)
8. [Reference](#reference)
9. [PACE and ABC](#pace-and-abc)
10. [Troubleshooting](#troubleshooting)
11. [Citation](#citation)

---

## Installation

Linux or macOS, Python 3.11 or newer; no GPU.

```bash
git clone https://github.com/shenlinyong/PACE.git && cd PACE
conda env create -f environment.yml && conda activate pace     # recommended
# or: python3 -m venv .venv && source .venv/bin/activate && pip install '.[io,ml]'
pace --version
```

The `io` extra (pyBigWig, cooler, pandas) is needed to read bigWig and Hi-C files; `ml` (scikit-learn, scipy) only for the optional classifier. With Docker, `docker build -t pace .` and mount your data folder (`-v "$PWD":/work -w /work`); every path you pass must be inside it.

`pace` lists all commands; `pace COMMAND --help` shows each option with its default. Record `pace --version` and `git rev-parse HEAD` for publications.

## Input data

PACE starts from processed files. All files must use the **same genome assembly and the same chromosome names** (`chr1` and `1` do not match; see `--aliases`).

| Data | Needed | Format | Option | Used for |
|---|---|---|---|---|
| Candidate regions | Required | BED, e.g. MACS2/3 peaks of ATAC/DNase or H3K27ac; several files are united | `-b` | Candidate elements |
| Gene annotation | Required | GTF with `transcript` or `exon` lines (Ensembl, NCBI) | `-g` | TSSs |
| Chromatin accessibility | Recommended | ATAC-seq or DNase-seq bigWig | `--atac` / `--dnase` | Activity |
| Active-enhancer mark | Recommended | H3K27ac ChIP-seq or CUT&Tag bigWig | `--h3k27ac` | Activity |
| Hi-C | Optional | `.cool`/`.mcool` contact frequencies, balanced | `--hic` (A) or `pace fit-prior` (B) | Contact |
| CTCF + genome | Optional, with Hi-C | CTCF ChIP-seq or CUT&Tag peaks (BED or narrowPeak) **and** the uncompressed genome FASTA | `--ctcf` + `--genome` | Domain boundaries |
| Chromosome sizes | Optional | `chrom.sizes`, `genome.fa.fai` or a bigWig | `-c` | Grid bounds (default: first bigWig) |

- **Activity.** At least one of ATAC, DNase or H3K27ac is required; accessibility (ATAC or DNase) plus H3K27ac is recommended, because their geometric mean requires an element to be both open and acetylated. Make bigWigs with, for example, deepTools `bamCoverage --normalizeUsing CPM`. Several bigWigs of one assay are averaged as replicates of the same animal.
- **Hi-C.** Give contact frequencies, not observed/expected matrices, p-values or loop calls: these remove the distance dependence PACE models. Convert `.hic` files to `.mcool` (for example with `hic2cool`) and balance them with `cooler balance`.
- **CTCF.** Both files are needed: the peaks say where CTCF is bound, and the built-in CTCF motif (JASPAR MA0139.1) found in each peak on the genome sequence gives the site's orientation. The peak signal (narrowPeak column 7) is the site strength. Without them, the model has no boundary term.
- **Preprocessing.** Read alignment, peak calling and Hi-C processing are done beforehand.

## Running PACE: choose your case

Choose by the Hi-C you have. All three cases use the same score; only the element–TSS contact $\widetilde C(E,t)$ differs.

| | A. Hi-C of the same tissue | B. Hi-C of another tissue | C. No Hi-C |
|---|---|---|---|
| Hi-C | Of the sample itself | Same species and assembly: another tissue, animal or public dataset | None |
| Contact $\widetilde C(E,t)$ | Measured contact + distance term | $P(d)$ fitted on that Hi-C | $1/\max(d,\,5\,\mathrm{kb})$ |
| Domain boundaries | With `--ctcf --genome`; $\beta$ fitted on the sample's Hi-C | With `--ctcf --genome`; $\beta$ fitted on that Hi-C | Not used |
| Option | `--hic` | `--prior` from `pace fit-prior` | `--distance-only` |

Cases A and B use the distance–boundary prior

```math
P(d)=a\left(\frac{\max(d,\,d_{\min})}{5\,\mathrm{kb}}\right)^{-\gamma}\,e^{-\beta\,\min(S,\,1)},
```

where $d$ is the element–TSS distance, $d_{\min}$ the larger of the Hi-C bin size and 5 kb where $P(d)$ is combined with measured bins (below one bin, distances are not resolved) and 5 kb where it replaces Hi-C (Case B), $\gamma$ the rate at which contact falls with distance, $S$ the summed strength (0–1 each) of the boundaries between element and TSS, and $\beta$ the drop in log contact across one full-strength boundary. $S$ is capped at 1 because contact falls at the first boundary and not again at each further one. The amplitude $a$ cancels in the score; $\gamma$ and $\beta$ are fitted on the Hi-C, so nothing is set by hand. Without `--ctcf` and `--genome` the prior is $P(d)$ without the boundary term.

### Case A: Hi-C of the same tissue

```bash
pace predict -b peaks.bed -g genes.gtf --atac pig1_ATAC.bw --h3k27ac pig1_H3K27ac.bw \
  --hic pig1.mcool --hic-resolution 10000 --ctcf pig1_CTCF.narrowPeak --genome Sscrofa11.1.fa \
  --species pig --assembly Sscrofa11.1 --tissue liver -o results/pig1
```

$\widetilde C(E,t)$ is the Gamma–Poisson posterior mean of the contact between the bins of element and TSS, with prior $P(d)$ fitted on the same map ([Contact from sparse Hi-C](#contact-from-sparse-hi-c)): bin pairs with many reads keep their measured value, sparsely sampled ones move toward $P(d)$, so low-depth maps add information instead of noise. A bin pair holds many element–TSS pairs at different distances; each is placed within it by the prior's ratio $P(d)/P(d_{\mathrm{bin}})$ between the exact distance (floored at 5 kb) and the distance of the bin centres, so the bin's measured level is kept and spread by distance (`within_bin_factor` in `resolved_contacts.tsv`; `--bin-level-contacts` turns this off). This matters most for coarse maps, where one 50–100 kb bin pair holds near and distant pairs alike. An element in the TSS's own bin, and a pair in a masked bin, takes $P(d)$; such pairs are marked in `resolved_contacts.tsv`. `--no-shrinkage` uses the balanced contact plus $P(\max(d,5\,\mathrm{kb}))$ instead (needed for coolers with non-integer counts), and `--strict-contacts` uses measured bin-pair contacts only, without pseudocount or prior, and leaves same-bin and masked pairs NA. `--hic-resolution` selects the bin size in an `.mcool`. This is the most direct case: loops between distal enhancers and promoters are measured.

### Case B: Hi-C of another tissue of the same species

```bash
# 1. fit the prior once on the available Hi-C, with CTCF peaks of your target tissue
pace fit-prior --cooler public_pig_liver.mcool::/resolutions/10000 \
  --ctcf muscle_CTCF.narrowPeak --genome Sscrofa11.1.fa \
  --species pig --assembly Sscrofa11.1 --tissue liver -o priors/pig_liver

# 2. score any sample of the species
pace predict -b peaks.bed -g genes.gtf --atac pig2_ATAC.bw --h3k27ac pig2_H3K27ac.bw \
  --prior priors/pig_liver --species pig --assembly Sscrofa11.1 --tissue muscle -o results/pig2
```

`pace fit-prior`:

1. **Boundaries** are placed between each reverse CTCF site and the next forward site, as CTCF loop anchors face into their own domain; strengths are scaled to 0–1.
2. **Valid bin pairs** are grouped by distance (30 logarithmic classes from $d_{\min}$ to 1 Mb, where enhancer–promoter links lie; closer pairs lie on the flat part of the prior) and by the crossed boundary strength $S$. Pixels absent from the sparse matrix count as zeros.
3. **$\beta$** is the drop in log contact per unit of $S$ between pairs at the same distance.
4. **$a$ and $\gamma$** describe the boundary-corrected decay with distance.

Both steps treat the summed contact of each class as Poisson with the number of bin pairs as exposure, so classes without contacts stay in the fit and the overall scale of the map changes only $a$. Beyond 1 Mb the prior extrapolates the fitted power law.

The fitted values are in `manifest.json`, the boundaries in `boundaries.tsv`. Without Hi-C of the sample, $P(d)$ is evaluated at the exact element–TSS distance with $d_{\min}=5$ kb, as in Case C. The prior's tissue may differ from the sample's; the transfer is recorded in the outputs. One prior serves every sample of the species and keeps them comparable.

### Case C: no Hi-C

```bash
pace predict -b peaks.bed -g genes.gtf --atac pig3_ATAC.bw --h3k27ac pig3_H3K27ac.bw \
  --distance-only --species pig --assembly Sscrofa11.1 --tissue liver -o results/pig3
```

Contact is estimated from distance alone, $\widetilde C(E,t)=1/\max(d,\,5\,\mathrm{kb})$; there is nothing to fit. Without Hi-C, $\beta$ cannot be estimated, so the boundary term is not used.

### Useful options

| Option | Default | Meaning |
|---|---|---|
| `-w/--width` | 500 | Element width on a fixed genome grid (bp) |
| `-r/--radius` | 5000000 | Pair each gene with elements up to this distance from its TSSs; 5 Mb matches ABC |
| `--gene-types` | all | Keep only these biotypes, e.g. `protein_coding` (Ensembl) or `mRNA` (NCBI) |
| `--skip-unlisted-chroms` | off | Ignore peaks and genes on contigs missing from the chromosome sizes |
| `--min-callable-fraction` | 0 | Elements with less bigWig coverage become NA; 0.8 is a common choice |
| `--missing-as-zero` | off | Read bases absent from the bigWig as measured zero signal |
| `--tss-weights` | `usage` | TSS weights: promoter activity share, `equal`, or `provided` (`promoters.tsv`) |
| `--no-boundaries` | off | Fit the distance decay only, even with CTCF data |
| `--no-shrinkage` | off | With `--hic`: use measured contacts as they are, without per-pair shrinkage |
| `--partial-policy` | `withhold` | Genes with unscored candidates: `withhold` or `conditional` |
| `-t/--threads` | 1 | Chromosomes scored in parallel; memory grows with each thread |

The output folder is self-contained: scores, every intermediate table (`prepared/`) and `resolved_config.yaml`, which records every setting with paths relative to the folder. `pace validate` checks the inputs without writing results. `-o` must be a new folder; `--force` replaces a PACE result after keeping the old one as a backup.

## Reading the results

| File | Content |
|---|---|
| `scores.tsv.gz` | **Main result**, one row per element–gene pair |
| `gene_summary.tsv` | Per gene: candidates, how many could be scored, denominator |
| `qc_report.json` | Overall quality summary; read it first |
| `resolved_contacts.tsv` | The contact used for every element–TSS pair and its source (measured, prior, neighbour, shrinkage) |
| `resolved_activity.tsv` | Activity values used, per element and assay |
| `promoter_weights.tsv` | TSS weights per gene |
| `run_manifest.json` | Input hashes, software version and environment, for the methods section |
| `resolved_config.yaml` | Every setting, defaults included |

Key columns of `scores.tsv.gz`:

| Column | Meaning |
|---|---|
| `element_id`, `gene_id` | The pair |
| `pace_score` | Share of the gene's support, 0–1; withheld (NA) for genes with unscored candidates |
| `A_used`, `Cbar` | Activity of the element; its contact with the gene |
| `support` | Activity × contact before dividing by the gene total |
| `normalization_status` | `complete`, `partial`, `zero_support` or `empty` |
| `reason` | Why a value is NA |

```bash
zcat results/pig1/scores.tsv.gz | head -5 | column -t
```

How to read them:

- **Check coverage first.** Only `complete` genes have a primary score; a share computed from an incomplete total would inflate the remaining candidates.
- **NA is not zero.** NA means the value could not be computed; zero means a measured zero signal. Never replace NA with 0.
- **A score of 1** means that one candidate carries all the gene's support.
- **Thresholds.** Set a cutoff on functional data of your species where available, and filter after scoring so that gene totals are unchanged.
- **The score ranks** the candidate elements of each gene.

## The model

```math
\boxed{
\mathrm{PACE}(E,G)
=\frac{A_\star(E)\,\overline C(E,G)}{\displaystyle\sum_{e\in\mathcal E(G)} A_\star(e)\,\overline C(e,G)}
=\frac{\left[\prod_{m\in\mathcal M}x_m(E)\right]^{1/|\mathcal M|}\sum_{t\in\mathcal T(G)}\pi(t\mid G)\,\widetilde C(E,t)}{\displaystyle\sum_{e\in\mathcal E(G)}\left[\prod_{m\in\mathcal M}x_m(e)\right]^{1/|\mathcal M|}\sum_{t\in\mathcal T(G)}\pi(t\mid G)\,\widetilde C(e,t)}
}
```

An enhancer regulates a gene when it is active in the tissue and chromatin folding brings it to the gene's promoter. The product of activity and contact is the element's *support* for the gene, and the score is its share of the gene's total support, so the scores of one gene sum to 1.

- **Activity $A_\star(E)$:** the geometric mean of the assays $\mathcal M$ used, e.g. $\sqrt{\text{ATAC}\times\text{H3K27ac}}$, where $x_m(E)$ is the replicate-averaged signal of assay $m$. The same assays are used for every element of a run. A constant scale factor of a track (library size, normalization) cancels in the share.
- **Contact $\overline C(E,G)$:** the weighted mean of the contacts with the gene's distinct TSSs $\mathcal T(G)$; transcripts that start at the same position count once. $\pi(t\mid G)$ is the TSS's share of the gene's promoter activity (the activity of the 500 bp cell containing it). A TSS whose cell is unmeasured takes the mean activity of the gene's measured TSSs (`pi_source: gene_mean_promoter_activity`); genes without measured TSSs, or with all-zero promoter activity, use equal weights. `--tss-weights equal` or `provided` (the `pi` column of `promoters.tsv`, e.g. from CAGE with `pace promoter-weights`) change this.
- **TSS contact $\widetilde C(E,t)$:** from Hi-C and the distance–boundary prior, depending on the case (above and below).
- **Denominator:** all candidates $\mathcal E(G)$ within the radius of the gene's TSSs, including promoters of the gene and its neighbours. A gene with one strong enhancer gives it a large share; a gene with many shares the credit.

*Worked example.* ATAC = 4 and H3K27ac = 9 give $A_\star=\sqrt{36}=6$. Two TSSs with weights 3/4 and 1/4 and contacts 2 and 6 give $\overline C=\tfrac34\cdot2+\tfrac14\cdot6=3$, so the support is 18. If the gene's total support is 60, the score is $18/60=0.3$.

### Contact from sparse Hi-C

| Situation | Contact used |
|---|---|
| Measured pair, `pace predict --hic` | Gamma–Poisson posterior mean, below |
| Measured pair, `--no-shrinkage` or `pace run` | Balanced contact + $P(\max(d,5\,\mathrm{kb}))$ (`--pseudocount auto`, when a fitted prior is available) |
| Element and TSS in the same bin | $P(d)$, or without a prior the maximum of the neighbouring bins, instead of the self-ligation diagonal |
| Masked or unbalanceable bin | $P(d)$ (`--allow-prior-fallback` with `pace run`; automatic in `pace predict`) |

**Shrinkage.** `pace predict --hic` (and `pace fuse` in step-by-step runs) treats the true contact as Gamma-distributed with mean $P(d)$ and shape $\kappa$, and the raw count $y$ as Poisson with mean $C/f$, where $f$ converts counts to balanced contact. The posterior mean is

```math
\widetilde C = w\,f\,y + (1-w)\,P(d),\qquad w=\frac{P(d)/f}{P(d)/f+\kappa}.
```

A pair whose expected count is small gets a small $w$ and stays near the prior; a pair with a large expected count follows the data. A zero count is kept as an observation. $\kappa$ comes from the prior or is estimated once from all unique bin pairs of the run, also when chromosomes are scored in chunks, and each replicate is shrunk before averaging. A fixed weight can be given instead (`pace run --contact-mode shrinkage --contact-reliability 0.5 --reliability-source "<how it was chosen>"`).

### Zero and NA

A measured zero and an unmeasured value are kept apart.

| Situation | Result |
|---|---|
| An assay is 0 at an element | Activity 0, support 0, score 0; the pair stays in the gene's total |
| An assay is missing (NA) at an element | Activity NA, even if another assay is 0 (`reason: unresolved_activity`) |
| Contact is NA for a TSS with positive weight | Contact NA (`reason: unresolved_contact`) |
| Hi-C pixel absent from the sparse matrix | Measured zero |
| bigWig without data for part of an element | Covered bases are averaged; elements below `--min-callable-fraction` or without coverage are NA |
| Some candidates of a gene are NA | `pace_score` withheld; `pace_score_conditional` gives the share among the scored candidates |
| All candidates of a gene have support 0 | No share can be computed; scores NA (`zero_support`) |

### Optional allocation term (experimental)

An enhancer near several genes may share its effect among them. The optional term multiplies support by $B(E,G)^\eta$, where $B(E,G)=\overline C(E,G)/\sum_{H}\overline C(E,H)$ is the share of the element's contact that goes to gene $G$ among all genes it is a candidate for. $\eta=0$ (default) leaves the score unchanged; $\eta=1$ favours elements whose contact is concentrated on one gene. $\eta$ can be fitted on functional perturbation data; see [equations and calibration](docs/ADVANCED.md#experimental-target-allocation).

## Step by step, several animals and comparisons

`pace predict` runs four commands, which can be run separately to reuse a catalog across animals or inspect each step:

```bash
# 1. candidate elements (500 bp grid cells overlapping peaks and TSSs) and distinct TSSs
pace catalog -b data/peaks.bed -g data/genes.gtf -c data/genome.fa.fai -o prepared/catalog

# 2. activity per bigWig, then one table
pace activity -i data/pig1_ATAC.bw    -a ATAC    -d prepared/catalog -o prepared/pig1_atac
pace activity -i data/pig1_H3K27ac.bw -a H3K27ac -d prepared/catalog -o prepared/pig1_k27ac
pace merge -t observed_activity \
  -i prepared/pig1_atac/observed_activity.tsv prepared/pig1_k27ac/observed_activity.tsv \
  -o prepared/activity

# 3. Hi-C contacts of every element-TSS pair, and the prior of the same map
pace contacts -i data/pig1.mcool -r 10000 -d prepared/catalog -o prepared/pig1_hic
pace fit-prior --cooler data/pig1.mcool::/resolutions/10000 --ctcf data/pig1_CTCF.narrowPeak \
  --genome data/Sscrofa11.1.fa --species pig --assembly Sscrofa11.1 --tissue liver -o prepared/pig1_prior

# 4. score
pace run -d prepared/catalog --activity prepared/activity/observed_activity.tsv \
  --contacts prepared/pig1_hic/observed_contacts.tsv \
  --contact-prior prepared/pig1_prior --allow-prior-fallback \
  --species pig --assembly Sscrofa11.1 --tissue liver --by-chromosome -t 4 -o results/pig1
```

`pace run` without `--samples` treats every sample as a replicate of one animal. bigWigs often omit regions without reads, which may mean "zero" or "not measured"; use `--missing-as-zero` only if you know which.

**Several animals.** Score each animal first and check that the top predictions agree. For a tissue-level map, run once with all animals and a sample table naming each library's animal:

```
sample_id	donor_id	assay	biological_replicate	technical_replicate	species	assembly	context_id	source_id
pig1_ATAC	pig1	ATAC	1	1	pig	Sscrofa11.1	liver	pig1_ATAC
pig2_ATAC	pig2	ATAC	1	1	pig	Sscrofa11.1	liver	pig2_ATAC
```

```bash
pace run -d prepared/catalog --activity prepared/all_activity.tsv --contacts prepared/all_contacts.tsv \
  --samples prepared/samples.tsv --target-level population_mean \
  --species pig --assembly Sscrofa11.1 --tissue liver -o results/liver_population
```

Technical replicates are averaged within each biological replicate, biological replicates within each animal, and every animal has equal weight. `species`, `assembly` and `context_id` must match the command line.

**Comparing animals, tissues or treatments.** Scores are shares, so a change in one element shifts the others of the same gene; do not subtract two score tables. Score both runs on the same candidates (the same peak files, e.g. the union of both animals' peaks, GTF, width, radius and assays) and the same contact prior, then:

```bash
pace compare --left results/pig1 --right results/pig2 -o results/pig2_vs_pig1
```

`pace compare` rescores both runs over their common candidates and reports right minus left. Runs with their own Hi-C each have their own contact model and need `--allow-evidence-difference`, which gives conditional differences only. Interpret a difference through activity, contact and the gene total, not `pace_score` alone.

**Whole genomes.** Every term is computed within a chromosome, so `pace predict` scores chromosomes separately (`--by-chromosome` in `pace run`) and in parallel (`-t`) with identical results. Memory is set by the largest chromosome; a smaller `--radius` reduces it.

## RNA-seq, methylation and other data

The score uses only activity and contact. Other data are added as annotations next to each prediction:

| Data | Prepare with | Add to `pace run` with |
|---|---|---|
| RNA-seq | `pace expression --tpm transcripts.tsv -d prepared/catalog` | `--expression` |
| DNA methylation (WGBS, RRBS) | `pace methylation --counts cpg.tsv -d prepared/catalog` | `--methylation` |
| Other marks or peaks (H3K4me1, H3K4me3, H3K27me3) | `pace features -b marks.bed -d prepared/catalog` | `--features` |

Expression is not multiplied into the score because it would scale numerator and denominator alike and cancel; use it to set aside genes not expressed in the tissue. Uncovered CpGs are reported as missing, not unmethylated.

**Optional classifier.** With functional labels (e.g. CRISPRi), `pace train` fits a separate classifier, applied with `pace run --ml-model` or `pace predict-ml`. Its `pace_ml_score` and `pace_ml_probability` are kept separate from `pace_score`.

## Reference

| Command | Purpose | Key options |
|---|---|---|
| `predict` | Peaks + GTF + bigWig (+ Hi-C, CTCF) to scores in one step | `-b -g --atac/--dnase/--h3k27ac --hic/--prior/--distance-only --ctcf --genome -o` |
| `catalog` | Candidate elements, TSSs and pairs | `-b -g -c -w -r --gene-types -o` |
| `activity` | Signal per element from a bigWig | `-i -a -d -o` |
| `contacts` | Element–TSS contacts from cool/mcool | `-i -r -d --no-balance -o` |
| `fit-prior` | Distance–boundary prior from cool/mcool (+ CTCF) | `--cooler --ctcf --genome --species --assembly --tissue -o` |
| `merge` | Combine tables with duplicate checks | `-t -i -o` |
| `run` / `validate` | Score prepared tables / check them | `-d --activity --contacts --contact-prior -o` |
| `fuse` | Per-pair shrinkage of sparse contacts | as `run` |
| `compare`, `stability`, `benchmark` | Downstream evaluation | see `--help` |
| `boundaries`, `prior`, `fit-hic`, `fit-labels`, `fit-eta`, `train`, `predict-ml` | Advanced models | see [Advanced use](docs/ADVANCED.md) |
| `features`, `methylation`, `expression`, `pairs`, `normalize-activity`, `promoter-weights`, `init` | Annotation tables and helpers | see `--help` |

Main `pace run` options:

| Option | Default | Meaning |
|---|---|---|
| `--target-level` | `individual` | `individual`: one animal; `population_mean`: animals averaged with equal weight |
| `--panel` | assays present | Activity assays: ATAC, DNase or H3K27ac alone, or accessibility + H3K27ac |
| `--contact-mode` | `observed` | `observed`, `shrinkage` or `prior_only` |
| `--contact-prior` | none | Prior folder from `fit-prior`; needed for `shrinkage` and `prior_only` |
| `--allow-prior-fallback` | off | Masked or missing contacts take the prior |
| `--allow-cross-context-prior` | off | Accept a prior fitted on another tissue |
| `--partial-policy` | `withhold` | Genes with unscored candidates get no primary score |
| `--eta` | `auto` | Allocation exponent; `auto` is 0 without a calibrator |
| `--by-chromosome`, `-t` | off, 1 | Chunked and parallel scoring |

All settings and defaults are listed in [Advanced use](docs/ADVANCED.md#complete-configuration-defaults).

**Input tables** (written by the commands above; usually only `samples.tsv` is written by hand). Tab-separated with a header, optionally gzipped; coordinates 0-based, half-open as in BED; `NA` for missing values and `0` only for measured zeros.

| Table | Columns |
|---|---|
| `units` | `element_id chrom start end anchor0 element_roles canonical_catalog_id` |
| `promoters` | `gene_id promoter_id chrom tss0 strand pi pi_source` |
| `candidates` | `element_id gene_id candidate_universe_id` |
| `samples` | `sample_id donor_id assay biological_replicate technical_replicate species assembly context_id source_id` |
| `observed_activity` | `element_id sample_id assay signal measurement_status callable_fraction unit normalization_id window_id` |
| `observed_contacts` | `element_id promoter_id sample_id contact_value measurement_status bin_pair_id scale resolution source_id` |
| `region_membership` | `region_id element_id source_id membership_rule` |
| `sources` (optional) | `source_id path_or_accession source_type assembly processing_method normalization_id checksum` |
| `expression` (optional) | `gene_id sample_id tpm status` |
| `methylation` (optional) | `chrom dyad_start0 methylated_count total_count sample_id assay` |
| `features` (optional) | `entity_type entity_id feature_name value evidence_id status` |
| `labels` (optional) | `label_id assayed_region_id gene_id context_id perturbation_type effect_direction effect_size label_status assay_id group_id source_id` |

`measurement_status` is one of `observed`, `unmeasured`, `low_coverage`, `unmappable`, `invalid`, `not_applicable`.

## PACE and ABC

The score is the ABC score: with identical activity, contact, one TSS and a complete candidate set, PACE and ABC give the same value. They differ in how the inputs are obtained:

| | ABC | PACE |
|---|---|---|
| TSSs of a gene | One TSS per gene | All distinct TSSs, weighted by promoter activity |
| Low-depth Hi-C | Distance-based pseudocount | Per-pair Gamma–Poisson shrinkage toward the distance prior |
| Domain boundaries | Not modelled | CTCF boundaries, effect $\beta$ fitted on Hi-C (`--ctcf --genome`) |
| Contact without the sample's Hi-C | Fixed power law fitted on human Hi-C | $\gamma$ (and $\beta$) fitted on Hi-C of the same species (Case B), or $1/d$ (Case C) |
| Missing contact | Set to zero | Kept as NA; the gene's score is withheld and a conditional share is reported |

Scores from different runs are comparable only when candidates, assays, normalization and TSS definitions are the same; see [comparison methods](docs/ADVANCED.md#comparisons-and-benchmarks).

## Troubleshooting

**`pace: command not found`** — activate the environment (`conda activate pace` or `source .venv/bin/activate`).

**Chromosome name mismatch** — BED, GTF, bigWig, Hi-C and FASTA disagree on names (`chr1` vs `1`, `chrMT` vs `MT`). Rename in one place or give `--aliases`; drop unneeded contigs with `--skip-unlisted-chroms`.

**"needs both --ctcf peaks and --genome FASTA"** — the boundary term needs both files. The FASTA must be uncompressed (`gunzip genome.fa.gz`).

**Most genes are `partial` or NA** — usually Hi-C is missing for some pairs. See the `reason` column. `pace predict --hic` fills masked bins from the prior; with `pace run`, add `--contact-prior` and `--allow-prior-fallback`. Do not replace NA with 0.

**"requires tables for every declared activity assay"** — `--panel` names an assay absent from the activity table; check the `assay` column or omit `--panel`.

**`scale` or `normalization_id` mismatch** — all contact tables and the prior must describe the same measurement. Reprocess the files the same way.

**`Output already exists`** — choose a new `-o` folder, or `--force` to replace it after a backup.

**Out of memory** — fewer `--threads`, a smaller `--radius`, or chunked scoring (default in `predict`, `--by-chromosome` in `run`).

**Only one animal?** — supported (`--target-level individual`).

Still stuck? Open an [issue](https://github.com/shenlinyong/PACE/issues) with `pace --version`, the command and `qc_report.json`.

## Citation

PACE is not yet published. Please cite this repository with the commit used (`git rev-parse HEAD`) and the ABC papers it builds on:

> Fulco CP, Nasser J, Jones TR, et al. Activity-by-contact model of enhancer–promoter regulation from thousands of CRISPR perturbations. *Nature Genetics* 51, 1664–1669 (2019).
>
> Nasser J, Bergman DT, Fulco CP, et al. Genome-wide enhancer maps link risk variants to disease genes. *Nature* 593, 238–243 (2021).

The built-in CTCF motif is JASPAR MA0139.1 (Castro-Mondragon JA, et al. JASPAR 2022. *Nucleic Acids Research* 50, D165–D173 (2022)).

## License and contact

MIT License. Written and maintained by Linyong Shen, Northwest A&F University. Questions and bug reports: [GitHub Issues](https://github.com/shenlinyong/PACE/issues). Methods, multiomics, classifier training, calibration and complete defaults: [Advanced use](docs/ADVANCED.md).
