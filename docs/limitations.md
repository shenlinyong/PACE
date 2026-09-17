# Scope and limitations

The new package is research software implementing a specified mathematical contract. Its tests
verify software behavior, not biological predictive performance. No real livestock model weights,
independent individual-effect validation or enhancer–gene perturbation validation is bundled.

| Capability | Current boundary |
|---|---|
| Three evidence regimes | Implemented with shared bulk-proxy order and asset checks; synthetic examples only |
| Quantitative CNN | Real CPU training/inference and safe weights; a baseline, not an optimized livestock model |
| Genome-only candidates | Existing atlas and trusted promoter cells; no whole-genome novel-element discovery |
| SNV/short indel | Normalized SNVs and fixed-target flanking indels of at most 50 bp; verified additional context when needed |
| Target-changing indel | Unresolved if the central training target changes length or correspondence |
| E–TSS geometry changes | Conservatively unresolved; no haplotype-specific distance reconstruction |
| Phase and ploidy | Explicit haploid/diploid; unphased heterozygotes and disconnected local PS blocks are unresolved |
| Callability | Callable BED or explicit reference-assumption policy; no automatic gVCF confidence-block conversion |
| SV/CNV | Flags reported input SV effects; does not discover unreported SVs or reconstruct complex dosage |
| Nine omics | Standard tables and track/count/interval adapters; annotation by default, named features in independent ML |
| Functional labels | First implementation uses unambiguous one-to-one region/unit mappings; multi-tile labels are excluded and reported |
| Continuous eta | Bounded within-gene ranking fit with explicit information/split checks; no claim of an optimal biological exponent, confidence interval or held-out improvement |
| Probability calibration | Specific to independent supplied labels and their sampling design; no universal causal probability |
| Contact uncertainty | Explicit observed/prior/shrinkage point estimates; no Poisson–Gamma or posterior sampling |
| Statistical uncertainty | Replicate consistency and coverage; no invented small-sample confidence intervals |
| External methods | Provenance-tagged score ingestion, not a claimed complete reimplementation of gABC |
| Scalability | Sparse candidate/bin-pair algorithms; current training TSV and model input batches reside in memory |
| Cross-species transfer | Scope mismatch is rejected; train/provide a target-scope asset rather than silently reusing another species |

Additional-reference context is conservative: unknown sequence, callability, target correspondence
or structural relationships can reduce coverage. A complete denominator only means every planned
candidate was handled, not that every biological enhancer was discovered. Conditional comparisons
describe their common measurable subset and must not be represented as complete genetic effects.

The supplied documents describe future research extensions; they do not mean those extensions
already exist. Unsupported extensions have no fake-success command. Full posterior inference,
copy-resolved support, arbitrary assembly mapping and de novo contact-network training remain
explicit research extensions.
