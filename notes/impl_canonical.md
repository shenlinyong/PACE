# Canonical implementation decisions

Updated 2026-09-17. The authoritative scientific definition is [FORMULA.md](../docs/FORMULA.md).
All current entry points use the src package; independent hand calculations remain
in the root review document.

| Decision | Source and rationale |
|---|---|
| Single src package | All current entry points must implement the same scientific contract; history retains retired APIs |
| Automatic eta with zero fallback | Formula and calibration contract: fit eligible functional training/calibration data within [0,1]; no added support term |
| float64 log normalization | Development E: protect high-dynamic-range shares without pseudocounts; preserve support=0 |
| Shared bulk marginal average | Review §3.1: 5/9 is distinct from copy-sum 3/7 |
| Fixed panel and promoter pi | Model §§4/6: missing required layers/TSSs do not silently redefine the estimand |
| Same-bin contacts use diagonal policy | Adapter/near-diagonal contract: these are not measured promoter self-support by default |
| Sparse candidate search | Development I.13: index local anchors and emit actual pairs, no genome-wide dense E×G array |
| Conservative indel/SV handling | Development F: preserve target identity, verify newly exposed context, mark changed geometry unavailable |
| Safe model storage | Development G: JSON contact priors and classifiers; no pickle loader |
| Transparent elastic-net solver | Direct mean-weighted objective avoids ambiguous library C scaling; intercept unpenalized |
| Training-only preprocessing | Development G/I.9: medians/IQR/tuning cannot observe external calibration/test data |
| Explicit capability/validation boundaries | Development H: synthetic/research/validated assets and three biological tasks stay distinct |

Core reference tolerances are 1e-10 absolute for hand-derived fractions and 1e-12 where comparing
the same deterministic kernel under reordered input. Optimizer tests use a 2e-6 intercept
tolerance against the analytical Bernoulli intercept; classifier optimization records convergence.
Classifier tests check optimization, held-out inference and persistence; no accuracy target is claimed for animal data.

Public documentation is consolidated around the current model. Incompatible retained
formulas, workflows and outputs were removed from the current branch after the user
identified the conflicting descriptions. The [public audit](public_docs_audit.md)
records the content disposition and verification. CITATION and AUTHORS use confirmed
authorship, without an invented DOI or coauthor.

Skills applied: research-software-engineering and its resource/documentation companions,
installed from a-attia/scicomp-research-skills commit
8435b16d91972c4f31b006de7bcacf1f5eb47e8e. The workflow emphasizes independent numerical references,
small tested changes, safe artifacts and honest limits. The skill installation is separate
from the software distribution and is not vendored into this repository.
