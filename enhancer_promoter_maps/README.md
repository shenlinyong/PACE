# Enhancer-gene maps of pig, cattle and chicken

This directory contains the enhancer-gene maps of eight tissues of pig, cattle and chicken: adipose, cerebellum, cortex, hypothalamus, liver, lung, muscle and spleen. Each species and tissue has one tab-separated file, named `<species>_<tissue>_enhancer_promoter_links.txt`.

The maps were built with PACE 1.0.0 (commit 5c36ea7) in Case B. The distance-boundary prior of each species was fitted on fibroblast Hi-C, and domain boundaries were placed from the CTCF ChIP-seq peaks of each tissue. Activity was computed from ATAC-seq (pig and cattle) or DNase-seq (chicken) and H3K27ac ChIP-seq of the FAANG pilot project. Each file lists the element-gene links with a PACE score of at least 0.02.

| Column | Description |
|---|---|
| `element_id` | Identifier of the regulatory element |
| `target_gene_id` | Symbol of the target gene |
| `chromosome_index` | Chromosome |
| `element_position` | Position of the element (bp) |
| `target_tss_position` | Position of the TSS of the target gene (bp) |
| `element_tss` | Distance between the element and the TSS (bp) |

Coordinates refer to Sscrofa11.1 for pig and ARS-UCD1.2 for cattle (Ensembl release 110), and to GRCg6a for chicken (Ensembl release 105).
