# Worked examples: incomplete inputs and score components

These examples use the bundled synthetic table to show how PACE handles missing assays, missing contacts and measured zero. They run entirely from repository files. Activate the `pace` Conda environment and run all commands below from the repository root.

For the underlying equations and hand calculations, see [the ABC comparison](ABC_COMPARISON.md). The example quality values demonstrate the interface; they are not a recipe for estimating biological QC.

## 1. Score the original example

```bash
python scripts/pace.py \
  --pairs example_quantified/candidates.tsv \
  --activity-config example_quantified/activity.json \
  --output results/worked/baseline.tsv
```

There are three enhancers and two genes. Gene `g1` has two TSSs; `g2` has one. Nine enhancer–TSS rows collapse to six enhancer–gene rows. Enhancer `e3` has no measured activity, so its two scores remain missing.

## 2. Prepare three input scenarios

```bash
python - <<'PYCODE'
from pathlib import Path
import numpy as np
import pandas as pd

out = Path('results/worked/inputs')
out.mkdir(parents=True, exist_ok=True)
pairs = pd.read_csv('example_quantified/candidates.tsv', sep='\t')

# A planned H3K27ac assay is unavailable; ATAC remains as supplied.
missing_mark = pairs.copy()
missing_mark['H3K27ac'] = np.nan
missing_mark['H3K27ac_quality'] = np.nan
missing_mark.to_csv(out / 'missing_mark.tsv', sep='\t', index=False, na_rep='NA')

# No contact measurements are available; positive priors are retained.
missing_contact = pairs.copy()
for column in ['contact_observed', 'contact_expected', 'contact_reliability']:
    missing_contact[column] = np.nan
missing_contact['contact_source'] = 'distance_prior'
missing_contact.to_csv(out / 'missing_contact.tsv', sep='\t', index=False, na_rep='NA')

# e1 has observed zero contact at every candidate TSS, with known reliability.
zero_contact = pairs.copy()
e1 = zero_contact['name'].eq('e1')
zero_contact.loc[e1, 'contact_observed'] = 0.0
zero_contact.loc[e1, 'contact_expected'] = 1.0
zero_contact.loc[e1, 'contact_reliability'] = 1.0
zero_contact.loc[e1, 'contact_source'] = 'matched'
zero_contact.to_csv(out / 'zero_contact.tsv', sep='\t', index=False, na_rep='NA')
PYCODE
```

Keep both assays in the activity JSON for the `missing_mark` scenario. This records a missing planned assay. An intentional ATAC-only design would instead declare only ATAC in the JSON, changing the planned evidence denominator.

## 3. Score and inspect the scenarios

```bash
for scenario in missing_mark missing_contact zero_contact; do
  python scripts/pace.py \
    --pairs "results/worked/inputs/${scenario}.tsv" \
    --activity-config example_quantified/activity.json \
    --output "results/worked/${scenario}.tsv"
done
```

```bash
python - <<'PYCODE'
import pandas as pd

columns = ['name', 'TargetGeneEnsemblID', 'activity', 'activity_observed_fraction',
           'contact_gene', 'contact_quality',
           'PACE.Score', 'contact_state', 'evidence_status', 'unscored_candidates']
for scenario in ['baseline', 'missing_mark', 'missing_contact', 'zero_contact']:
    result = pd.read_csv(f'results/worked/{scenario}.tsv', sep='\t')
    print(f'\n{scenario}: {len(result)} gene-level rows')
    print(result.sort_values(['name', 'TargetGeneEnsemblID'])[columns].to_string(index=False))
PYCODE
```

Expected changes:

| Scenario | What to inspect | Expected result |
| --- | --- | --- |
| Baseline | `e1` activity from ATAC 3 and H3K27ac 0 | Activity 1; both measurements participate |
| Missing H3K27ac | `e1` activity and `activity_observed_fraction` | Activity 3; observed fraction 0.5 because only one of two planned assays is available |
| Missing contact | `contact_state`, `contact_quality` | Every row is `prior_only`, with contact quality 0; finite scores are provisional |
| Measured zero contact | `e1` contact and score | Gene contact and score are 0; this differs from missing-contact fallback |
| Every scenario | `e3` scores, output row count | Two missing scores remain in six output rows; each gene reports one unscored candidate |

`e1` has undefined `target_share` in the zero-contact scenario because all its target contacts are zero. Its finite activity and known zero contact still establish zero raw support. Missing `e3` activity remains undefined in every scenario.

The baseline itself is incomplete: an unscored candidate prevents either gene from having sufficient input evidence under the default rule. A high score for `e2` does not fill that evidence gap.

## 4. Remove target allocation while keeping the data fixed

```bash
python scripts/pace.py \
  --pairs example_quantified/candidates.tsv \
  --activity-config example_quantified/activity.json \
  --competition-power 0 \
  --output results/worked/no_allocation.tsv
```

```bash
python - <<'PYCODE'
import pandas as pd

key = ['chr', 'start', 'end', 'TargetGeneEnsemblID']
full = pd.read_csv('results/worked/baseline.tsv', sep='\t')
ablation = pd.read_csv('results/worked/no_allocation.tsv', sep='\t')
comparison = full[key + ['PACE.Score']].merge(
    ablation[key + ['PACE.Score']], on=key,
    suffixes=('_allocation', '_no_allocation'), validate='one_to_one')
print(comparison.to_string(index=False))
PYCODE
```

This changes $\eta$ from 1 to 0. Activity, contact processing, TSS integration and missing-value reporting are unchanged. It is a component ablation; a native ABC comparison must run the specified ABC implementation separately.

## Apply the examples to a biological analysis

Choose signal scales before comparing samples. Supply contact expectations and reliability from independent, compatible measurements. Preserve the complete candidate background and TSS weights when testing a component. Keep `NA` in inputs and unfiltered outputs, and calibrate selection thresholds on independent validation. [Input schemas](INPUTS.md) · [Parameter controls](PARAMETERS.md) · [Real-data tutorial](TUTORIAL.md#replace-the-example-with-your-species-and-tissue)
