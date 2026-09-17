# Verification and release acceptance

Software tests check implementation behavior with synthetic inputs and small real
file formats. They do not establish biological accuracy for an animal species,
tissue, breed or individual. The package includes no independently validated
livestock weights or claimed functional benchmark result.

## Reproduce the software checks

From an editable installation with `dev,io,sequence,ml` dependencies:

```bash
python -m pytest -q
python scripts/check_public_docs.py
python -m ruff check .
```

The exact passing count depends on the checked-out revision. Use the test log for
that revision and the [CI run](https://github.com/shenlinyong/PACE/actions/workflows/ci.yml),
not an earlier release's count. Required extras should be installed when evaluating
the complete suite; inspect skips as well as failures.

## Regression coverage

| Area | Important invariant |
|---|---|
| Core formula | Hand-derived activity, TSS, allocation and fractional exponent results; fixed candidate sets |
| Numerical comparison | Finite log support survives raw product underflow/overflow and round trips |
| Contact integrity | Mixed resolution or declared measurement contracts cannot be silently pooled or compared as full biological deltas |
| Automatic allocation | Zero fallback, grouped held-out selection, stability and calibration confirmation; test-label exclusion |
| Classifier | Connected repeated entities cannot leak across tuning folds; feature scope must match inference |
| Genomic reconstruction | Only carried alleles matter; missing GT, unsupported geometry and both BND endpoints remain conservative |
| External predictions | Genomic binding is checked and imports cannot restore invalid windows |
| Omics | Invalid RNA cannot enter ML; extra marks and methylation have resolvable evidence identities |
| CpG features | Element/promoter windows, coverage thresholds, reference CpG denominators and WGBS/RRBS separation |
| Catalog preparation | End-of-chromosome partial windows are explicitly handled by the generated configuration |
| User interface | Three modes, CLI/YAML paths, transactional outputs, consistent defaults and resolvable documentation links |

## Installation checks

CI configurations cover supported Python environments, optional IO/sequence
adapters, source/wheel distribution, Conda installation and the three Docker
example modes. A configured job is not evidence that its latest execution passed:
inspect the workflow attached to the exact commit. Local environments without
Conda or Docker cannot substitute a source-file inspection for an actual build.

## Biological validation

Quantitative assay prediction, within-locus individual effects and functional
regulatory-link prediction require separate appropriate validation. Use independent
animals, loci or experimental groups according to the intended claim. Include
coverage and unavailable tested positives; avoid selecting parameters with the
final test set. Association/QTL and contact data used as inputs are not independent
functional ground truth. More details: [training](training.md),
[comparison and benchmark](comparison.md), [limitations](limitations.md).
