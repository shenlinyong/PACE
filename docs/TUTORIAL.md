> **Interface scope:** This retained guide documents the legacy region-based scripts/workflow.
> For the installable canonical-grid package and its September 2026 defaults, see [software.md](software.md).

# Running an analysis

PACE (Prediction of Activity-based regulatory Connections for Enhancers) offers a numerical-table interface and genomic-file adapters. This tutorial starts with files shipped in the repository, then shows how to use a dataset from your species and tissue. Run commands from the repository root with `conda activate pace`.

## Step-by-step genomic-file example

### 1. Prepare the promoter catalogue

```bash
mkdir -p results/tutorial
python scripts/prepare_tss.py \
  --gtf example/reference/example_annotation.gtf.gz \
  --output results/tutorial/tss.tsv
```

This converts GTF coordinates to BED0 TSS positions, deduplicates transcript TSSs and assigns uniform weights across distinct TSSs per gene. When no transcript record exists, the gene boundary is used and labelled `gene_boundary_fallback`. It does not discover unannotated promoters.

### 2. Construct candidate intervals

```bash
python workflow/scripts/pace_candidate_regions.py \
  --narrowPeak example/data/example_peaks.narrowPeak \
  --chrom_sizes example/reference/example.chrom.sizes \
  --output results/tutorial/candidates.bed \
  --nStrongestPeaks 150000 \
  --peakExtendFromSummit 250
```

The bundled peaks were prepared in advance. For a biological run, call peaks on suitably processed accessibility data first, or use the optional Snakemake route below. The candidate tool ranks the narrowPeak `signalValue`, extends summits and retains the resulting intervals. It also accepts `--blacklist` and `--tss_regions` BED files; the latter marks overlaps internally and does not add promoter candidates. If needed, construct an enhancer/promoter union upstream and pass that candidate BED to the next step.

### 3. Quantify activity

```bash
python workflow/scripts/pace_neighborhoods.py \
  --candidate_regions results/tutorial/candidates.bed \
  --genes results/tutorial/tss.tsv \
  --chrom_sizes example/reference/example.chrom.sizes \
  --accessibility_file example/data/example_ATAC.tagAlign.gz \
  --accessibility_type ATAC \
  --H3K27ac example/data/example_H3K27ac.tagAlign.gz \
  --activity_method missing_geometric \
  --output_dir results/tutorial/neighborhoods
```

Outputs: `EnhancerList.txt` and `GeneList.txt`. BED/tagAlign/BAM signals are counted per kb of candidate width; bigWig tracks use exact regional summaries. Read density is not library-depth normalization. The example demonstrates execution with synthetic signals. For biological comparisons, specify assay-appropriate preprocessing; use normalized tracks or the table route when fixed scales and local quality need explicit control.

For an accessibility-only run, omit `--H3K27ac`. If a second assay was planned but is missing at some loci, represent this with `NA` in a table and declare both modalities in the activity JSON. Omitting an assay by design and failing to observe a planned assay are different input regimes.

### 4. Score every candidate pair

```bash
python workflow/scripts/pace_predict.py \
  --enhancers results/tutorial/neighborhoods/EnhancerList.txt \
  --genes results/tutorial/neighborhoods/GeneList.txt \
  --output results/tutorial/all_predictions.tsv.gz \
  --max_distance 5000000
```

Without qualified contact measurements, this is a labelled distance-prior run. The main file retains all generated gene-level links, including missing scores. The script also writes a convenience `.filtered.tsv` companion. Explicit filtering below makes the chosen threshold visible.

### 5. Filter and summarize

```bash
python workflow/scripts/pace_filter.py \
  --predictions results/tutorial/all_predictions.tsv.gz \
  --output results/tutorial/selected.tsv \
  --threshold 0.02

python workflow/scripts/pace_metrics.py \
  --predictions results/tutorial/all_predictions.tsv.gz \
  --output_dir results/tutorial/metrics \
  --sample_name tutorial
```

Retain the unfiltered file. The QC directory contains `QCSummary_tutorial.tsv` and `QCPlots_tutorial.pdf`. QC plots describe predictions; they do not validate enhancer function. See [output fields](IO_FORMATS.md).

The convenience command `bash example/run_example_direct.sh` performs the corresponding five-step read example using the bundled gene BED. It writes to `example/results/Example_Sample/` and produces 12,000 unfiltered links.

## A single-command signal-file calculator

The [example YAML](../config/config_direct.example.yaml) and [sample sheet](../config/biosamples_direct.example.tsv) select the same synthetic signals:

```bash
python scripts/calculate_pace_score.py \
  --config config/config_direct.example.yaml \
  --sample Example_Sample \
  --candidates results/tutorial/candidates.bed \
  --genes results/tutorial/tss.tsv \
  --output results/tutorial/standalone.tsv
```

The main file contains all predictions; `standalone.tsv.filtered.tsv` contains selected rows. Each biosample name must occur exactly once in the sample sheet. The standalone calculator and the two-stage neighborhood/predictor route delegate to the same scoring kernel. [Parameters](PARAMETERS.md) documents which controls each wrapper exposes.

## Replace the example with your species and tissue

1. **Choose the assembly and annotation together.** Use a matched FASTA, GTF and chromosome-size file. Keep the same identifiers (`1` versus `chr1`, scaffold names and stable gene IDs) in reads, signals, peaks and annotation. PACE does not lift coordinates between assemblies.
2. **Prepare processed assay inputs.** Use aligned, quality-controlled accessibility reads or an appropriately normalized bigWig. Add matched H3K27ac if available. Keep replicate handling and normalization in the run record; they are not inferred from filenames.
3. **Prepare candidate intervals.** Use same-assembly peaks, artifact filters and a declared promoter-candidate policy. An effective genome size used for peak calling must match the species/assembly and mapping protocol; do not use the generic mammalian placeholder for chicken.
4. **Run `prepare_tss.py` on the matched GTF.** It keeps all supplied gene biotypes. Apply any justified eligibility rule before analysis and record it; missing genes cannot be recovered by a score formula.
5. **Replace the input paths in the commands or copy and edit the direct sample/YAML templates.** Leave absent optional files blank. Process one sample/assembly per invocation.
6. **Inspect coverage before selecting links.** Review missing-score counts, contact states, evidence reasons and candidate coverage. Calibrate decision thresholds using independent validation relevant to the species and tissue.

To make chromosome sizes from an uncompressed reference FASTA:

```bash
samtools faidx reference/animal.fa
cut -f1,2 reference/animal.fa.fai > reference/animal.chrom.sizes
```

`animal.fa` is a path placeholder for your real reference, not an included dataset. Record its assembly accession and checksum. A size file generated this way gives chromosome lengths; it is not an estimate of MACS2 effective genome size.

## Add measured contacts or RNA context

For qualified measured contact, run the complete command below after preparing the enhancer and gene lists. Replace the two `data/` paths with your own contact file and metadata; these files are not bundled inputs:

```bash
python workflow/scripts/pace_predict.py \
  --enhancers results/tutorial/neighborhoods/EnhancerList.txt \
  --genes results/tutorial/neighborhoods/GeneList.txt \
  --output results/tutorial/measured_contact.tsv.gz \
  --max_distance 5000000 \
  --hic_file data/tissue_contacts.bedpe \
  --hic_type bedpe --hic_resolution 5000 \
  --contact_metadata data/contact_metadata.tsv
```

The metadata is keyed by enhancer coordinates, gene ID and TSS and supplies compatible expected contacts, reliability and `matched` or `surrogate` provenance. If required metadata are absent, the observation remains visible but the structural score falls back to the distance prior. Setting every reliability to 1 simply because a file exists is not a valid QC procedure. [Contact schemas](INPUTS.md#contact-metadata)

The smoke workflow provides a complete, independently checked BEDPE example in `results/smoke/inputs/` after running `python scripts/smoke_test.py --output-dir results/smoke`. Its numerical QC values are synthetic fixtures.

For RNA annotation, add `--expression data/rna.tsv` to `pace_predict.py`. The table must contain stable `gene_id` and `TPM` columns. Low expression, zero and absent identifiers remain distinguishable. Expression does not multiply the primary score. Use `pace_filter.py --only_expressed` only if expression-based output selection is part of the declared analysis.

## Prepared tables and model ablations

Use `scripts/pace.py` when activity scales, quality or contact expectations have already been quantified. It expects one enhancer–gene–TSS row and scores the supplied candidate set without generating or distance-filtering pairs. The [input specification](INPUTS.md) lists required fields.

The [worked examples](WORKED_EXAMPLES.md) provide runnable cases for a missing assay, missing contact and measured zero contact. To remove enhancer-centred allocation while holding the input table fixed:

```bash
python scripts/pace.py \
  --pairs example_quantified/candidates.tsv \
  --activity-config example_quantified/activity.json \
  --competition-power 0 \
  --output results/quickstart/no_allocation.tsv
```

This is a PACE component ablation, not a command for running the complete original ABC pipeline. Evidence threshold changes are available through `--evidence-threshold`; they change evidence designation rather than structural score.

## Optional Snakemake workflow

Use `conda activate pace-workflow` after creating the [optional environment](INSTALLATION.md#optional-peak-calling-and-snakemake-environment). This route starts with aligned reads and includes MACS2. It does not align FASTQ or prepare Hi-C matrices.

First inspect the bundled example DAG, then execute it:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile example/config.yaml --cores 2 --dry-run all

snakemake --snakefile workflow/Snakefile \
  --configfile example/config.yaml --cores 2 all
```

The explicit target `all` requests predictions and QC. Run inside the prepared environment; `--use-conda` is unnecessary for these commands. The example YAML controls the output path. Do not run the direct example concurrently into the same output directory.

For your own project, copy `config/config.yaml` and `config/config_biosamples.tsv` to new files. Set the biosample-table path, output directory and every reference path. Keep the full sample-table header: the wrapper uses `biosample`, `DHS`, `ATAC`, `default_accessibility_feature`, `HiC_file`, `HiC_type`, `HiC_resolution`, `alt_TSS` and `alt_genes` even when optional cells are empty. Provide exactly one of DHS/ATAC per sample. Add activating marks only when their corresponding configuration entries are enabled.

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile config/my_project.yaml --cores 8 --dry-run all

snakemake --snakefile workflow/Snakefile \
  --configfile config/my_project.yaml --cores 8 all
```

`my_project.yaml` must be created and populated with real paths first. The wrapper currently restricts binary/contact sample metadata to 5-kb resolution. Its exposed controls are narrower than the direct/table interfaces; consult [the parameter mapping](PARAMETERS.md#5-which-configuration-entries-are-active) before changing YAML keys. The small synthetic workflow is not a validation of all species, sequencing protocols or binary contact files.

## Compare predictions and annotate genetic support

For experimentally tested enhancer–gene pairs, the comparison utility accepts original-ABC and PACE predictions:

```bash
python scripts/benchmark_compare.py \
  --labels tested_pairs.tsv --abc abc.tsv --pace pace.tsv \
  --abc-threshold 0.02 --pace-threshold 0.02 \
  --output results/benchmark
```

These filenames stand for real experimental tables. Labels must be 0/1 or missing; absent eQTL support is not a negative functional label. Thresholds shown illustrate syntax; freeze each on independent validation. Supply unfiltered predictions so omitted positives remain visible in recall and coverage.

To annotate independent eQTL support:

```bash
python scripts/eqtl_pair_support.py \
  --predictions pace.tsv --eqtls independently_filtered_eqtls.tsv \
  --output results/eqtl_supported.tsv
```

The eQTL table needs `chr`, zero-based `pos0`, and `gene_id`. Match tissue/assembly and filter associations independently. Pair support requires both variant overlap and target-gene agreement. This is association annotation, not an FDR estimate.

Finally, save the [environment and commit](INSTALLATION.md#save-the-exact-software-environment), configuration, input manifest, normalization choices and complete prediction table.
