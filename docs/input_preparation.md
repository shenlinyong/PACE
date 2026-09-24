# Preparing PACE inputs

PACE starts from processed genomic measurements and annotations. It does not align FASTQ,
call peaks, call variants, or perform general liftover. Keep species, assembly, chromosome
names, normalization and quantitative windows consistent across all files.

## A project can have many assays and many replicates

Keep one canonical table per evidence role, even when the source data arrive as many files. For example, merge the prepared rows from `animal_001_rep1_liver_ATAC.bw`, `animal_001_rep2_liver_ATAC.bw` and `animal_001_rep3_liver_ATAC.bw` into an activity table with one `sample_id` per replicate; do the same for contact, RNA-seq and methylation tables. Register the biological and technical replicate relationships in `samples.tsv` rather than encoding them only in filenames.

Use `observed_activity.tsv` for the fixed primary panel (ATAC, DNase and/or H3K27ac), `observed_contacts.tsv` for Hi-C/Prom-Hi-C or a compatible contact assay, `expression.tsv` for RNA-seq, `methylation.tsv` for WGBS/RRBS counts, and `features.tsv` for additional histone marks or CTCF. Additional assays remain separately named evidence; they do not become extra factors in the PACE formula automatically. The [multi-omics guide](MULTIOMICS.md) explains the role and definition of each table.

## Canonical catalog from BED and GTF

```yaml
kind: catalog
bed: atlas.bed
gtf: annotation.gtf
chrom_sizes: chrom_sizes.tsv
source_id: atlas_source
width: 500
offset: 0
radius: 5000000
include_promoters: true
# aliases: aliases.tsv  # columns alias, canonical; no automatic chr-prefix guessing
```

`chrom_sizes.tsv` has header `chrom<TAB>length`. BED is zero-based half-open; GTF is one-based
closed. The adapter uses transcript records, preserves versioned gene/transcript IDs, computes
positive-strand TSS as start−1 and negative-strand TSS as end−1, and deduplicates physical TSSs.
Transcript annotations without transcript records must first be converted to promoters.tsv.
The atlas is mapped by at least one base overlap to fixed grid cells. Duplicate sources never
add another scoring unit. Incomplete edge cells are reported and excluded. Default promoter
cells are included once. Candidate generation binary-searches local cis neighborhoods.

```bash
pace prepare --config prepare_catalog.yaml --out prepared/catalog
```

This writes units, promoters, candidates, region membership, transcript mappings,
chrom_sizes.tsv and run_catalog_config.yaml. Merge that exported catalog block
into the run configuration (or use --catalog-dir, which imports the metadata),
so excluded terminal cells retain their reference-boundary justification. It does
not invent quantitative values for the new windows: extract tracks again or provide independently
precomputed matching quantities. Add samples.tsv and sources.tsv with the scientific context
before scoring. `provided_regions` is a distinct measured-only profile for user-defined regions;
its candidate background must be matched explicitly when comparing results.

## bigWig activity and mark annotations

```yaml
kind: bigwig
track: H3K27ac.bw
units: catalog/units.tsv
sample_id: animal1_H3K27ac
assay: H3K27ac
unit: normalized_signal
normalization_id: library_norm_protocol_1
window_id: grid:500:mean
missing_is_measured_zero: false
minimum_callable_fraction: 0.8
```

The adapter computes a nonnegative mean and callable fraction. Unstored signal is missing
unless the source format explicitly represents measured zero and the corresponding switch
is true. Negative tracks are rejected as activity. The callable threshold above is an example,
not a universal recommendation. The same adapter accepts ATAC/DNase and auxiliary histone/CTCF
tracks. Supply promoter-window intervals with a distinct target identifier for promoter marks;
map their values to `entity_type=promoter` in features.tsv. Do not reuse the central element
window label for a different promoter quantification protocol.

The `kind: bed_features` preparation command (and Python `interval_features` adapter)
extracts peak overlap and CTCF interval occupancy from BED sources; motif orientation requires an explicit strand-bearing motif table and `motif_strands: true`. Its output is
annotation evidence, never an extra contact multiplier. Gene-level promoter features use fixed
pi weights by default. Optional gene-wide filtering and its retained-weight report are described in [the model](FORMULA.md#multiple-tsss).

## cool/mcool contact

```yaml
kind: cooler
contact: contacts.mcool::/resolutions/5000
pairs: element_promoter_pairs.tsv
resolution: 5000
balanced: true
missing_pixels_are_zero: true
sample_id: animal1_HiC
source_id: contact_source
scale: balanced_contact_protocol_1
normalization_id: hic_norm_protocol_1
```

The pairs table contains `element_id,promoter_id,chrom,anchor0,tss0`. The unit anchor is
floor((start+end−1)/2); TSS uses its exact base. The reader queries sparse bin rows and preserves
shared bin-pair identifiers. Invalid balanced bins stay NA. Unstored pixels are zero only when
explicitly specified and both bins are valid. Choose a matching resolution and state the scale;
raw counts and balanced values are not interchangeable. Near-diagonal handling is configured
separately in the run. The adapter records same-bin neighbor maxima. A human reference contact shape is available only as an explicitly selected, unvalidated prior-only baseline. See [practical workflow](PRACTICAL_WORKFLOW.md) for direct cooler prior fitting.

PACE expects contact retaining distance background. An O/E track must first be multiplied by
its matched distance expectation; log-O/E requires the corresponding inverse transform before
multiplication. P values and correlation coefficients are not contact. Convert `.hic` upstream
to cool/mcool with an independently verified converter and record its version and settings.

## Methylation and RNA

```yaml
kind: methylation
counts: cpg_dyad_counts.tsv
units: catalog/units.tsv
minimum_coverage: 1
# reference_cpg: reference_cpg_counts.tsv  # element_id,n_cpg
```

CpG counts must already refer to the forward-strand dyad, with methylated ≤ total reads and
one row per dyad/sample/assay. The Python `merge_stranded_cpg` adapter converts explicit-strand
cytosine calls using a reference FASTA. It rejects duplicate strands or non-CG reference bases.
`M_site` averages site fractions; `M_pooled` pools reads. Zero methylation, no CpG, no coverage
and low coverage remain different. Without reference CpG counts the coverage fraction is NA.
Ordinary bisulfite assays do not distinguish 5mC and 5hmC.

```yaml
kind: rna
expression: transcript_tpm.tsv
transcript_mapping: catalog/transcript_mapping.tsv
```

Transcript TPM is summed only through the supplied transcript-to-gene/TSS map. Gene TPM can
be supplied directly as expression.tsv; it is not distributed back to transcripts or used to
invent TSS usage. Prepared GTF TSS weights are equal across distinct physical TSSs. `pace prepare-promoter-weights` can freeze explicit weights from measured ATAC/DNase/H3K4me3/CAGE promoter signals; these are proxies requiring validation.

For automatic element **and promoter** methylation annotations during scoring,
pass the raw CpG table as `inputs.methylation` and use the run's `methylation` block
for minimum coverage, strand-aware promoter windows and reference CpG denominators.
The standalone summary above is not itself a valid raw counts input. See the
[complete multiomics guide](MULTIOMICS.md) for exact formats.
