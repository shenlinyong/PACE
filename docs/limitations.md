# Scope and limitations

PACE is research software for relative regulatory support from measured activity and specified contact evidence. Software tests check calculations and input handling; they do not establish predictive accuracy for a livestock species, tissue, breed or experimental condition.

| Capability | Boundary |
|---|---|
| Activity | Requires measured ATAC, DNase or H3K27ac; no replacement of unavailable activity |
| Candidate catalog | Supplied/experimental regions and trusted promoter annotations; completeness refers to this planned catalog |
| Contact | Measured contacts or an explicitly supplied applicable distance prior; prior evidence is not an observed regulatory loop |
| Contact shrinkage | Specified observed/prior point estimate and reliability source; no automatically fitted edge-wise posterior |
| Multiple TSSs | Fixed distinct physical promoters and weights; incorrect annotations remain a limitation |
| Replicates | Specified technical/biological/donor aggregation; no automatic batch correction or increase in independent sample size |
| Additional omics | Named experimental annotations; optional separate classifier with its own validation |
| Allocation | Grouped functional-label selection with zero fallback; final independent evaluation remains necessary |
| Calibrated probability | Specific to supplied functional labels and sampling design; not universal causality |
| Between-animal comparison | Comparable measured states on common denominators; not an isolated genetic or causal effect |
| Cross-species/tissue use | Explicit background checks; any accuracy claim requires target-scope evaluation |
| Scalability | Sparse candidate contact processing; training tables reside in memory |

No universal livestock weights or independent biological benchmark results are bundled. A high conditional fraction can occur when only a small subset is measurable; the default primary score is withheld in that case. Always report candidate construction, assay availability, exclusions and coverage. Do not interpret partial normalization as complete regulatory discovery.

[Formula](FORMULA.md) · [Verification](validation.md) · [Migration](migration.md)
