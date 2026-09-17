# Canonical implementation decisions

Updated 2026-09-17. The development and review documents at the repository root are the primary
scientific specification; the existing scripts implement a different historical region profile.

| Decision | Source and rationale |
|---|---|
| Separate src package | Supplied contract A/B: preserve existing APIs and analyses while implementing revised rules |
| eta=0 default, no residual mass | Model §§2/6 and review: allocation is optional; unknown support is not a fitted constant |
| float64 log normalization | Development E: protect high-dynamic-range shares without pseudocounts; preserve support=0 |
| Shared bulk marginal average | Review §3.1: 5/9 is distinct from copy-sum 3/7 |
| Fixed panel and promoter pi | Model §§4/6: missing required layers/TSSs do not silently redefine the estimand |
| Same-bin contacts use diagonal policy | Adapter/near-diagonal contract: these are not measured promoter self-support by default |
| Sparse candidate search | Development I.13: index local anchors and emit actual pairs, no genome-wide dense E×G array |
| Conservative indel/SV handling | Development F: preserve target identity, verify newly exposed context, mark changed geometry unavailable |
| Safe model storage | Development G: JSON/NPZ for simple models and safetensors for CNN; no pickle loader |
| Transparent elastic-net solver | Direct mean-weighted objective avoids ambiguous library C scaling; intercept unpenalized |
| Training-only preprocessing | Development G/I.9: medians/IQR/tuning cannot observe external calibration/test data |
| Explicit capability/validation boundaries | Development H: synthetic/research/validated assets and three biological tasks stay distinct |

Core reference tolerances are 1e-10 absolute for hand-derived fractions and 1e-12 where comparing
the same deterministic kernel under reordered input. Optimizer tests use a 2e-6 intercept
tolerance against the analytical Bernoulli intercept; classifier optimization records convergence.
CNN tests check gradients, updates, masks and persistence, not an invented accuracy target.

Documentation audit: README has a new canonical entry section; all original substantial legacy
content remains below its explicit interface boundary. Existing guides received scope notices,
with their scientific formulas left unchanged for historical reproducibility. New model/data/
parameter/training/comparison guides cover the new package. CITATION and AUTHORS use only the
confirmed author and institution, with no invented DOI or coauthor.

Skills applied: research-software-engineering and its resource/documentation companions,
installed from a-attia/scicomp-research-skills commit
8435b16d91972c4f31b006de7bcacf1f5eb47e8e. The workflow emphasizes independent numerical references,
small tested changes, safe artifacts and honest limits. The skill installation is separate
from the software distribution and is not vendored into this repository.
