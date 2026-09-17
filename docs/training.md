# Training and calibration

Frozen-model inference never refits weights, scales, medians, quantiles or thresholds.
Eta is the explicit exception: providing `--eta-labels` requests its fitting during a run;
`--eta-model` instead reuses a frozen estimate. All distributed weights are synthetic.

## Continuous allocation eta

`PACE fit-eta --config run.yaml --eta-labels functional_labels.tsv --out results/eta`
fits a bounded continuous exponent from applicable training/calibration perturbations.
The same behavior is available in a normal run with `--eta-labels`. Without suitable
data, automatic eta remains zero with an explicit fallback reason. Freeze
`results/eta/eta_calibration.json` with `--eta-model` for subsequent compatible runs.
The [calibration specification](eta_calibration.md) defines the label schema, convex
ranking objective, conservative information requirements and split-isolation rules.

## Quantitative sequence CNN

```yaml
data: sequence_training.tsv
model_id: chicken_liver_quantitative_research
species: chicken
assembly: assembly_identifier
context_id: liver
target_level: individual
assays: [ATAC, H3K27ac]
signal_unit: normalized_signal
normalization_id: frozen_normalization_protocol
is_synthetic: false
input_length: 8192
output_window: 500
seed: 17
epochs: 50
batch_size: 64
learning_rate: 0.001
patience: 5
channels: [64, 128, 128]
threads: 1
max_n_fraction: 0.05
```

The training TSV contains `sequence_id,sequence,split,group_id,chrom,start,end` and one signal
column per assay. Each row contains an actual sequence of input_length bases; start/end describe
the central quantitative target. Labels are nonnegative signals or NA. `split` is train,
validation or test. Group IDs and overlapping input contexts cannot cross splits. If a design
also needs donor/family isolation, encode it in the groups before calling the trainer.

The architecture is three same-padded convolutions (kernels 15/7/5), GELU and factor-four average
pooling, followed by a local target-overlap pool plus a separate global-context pool. Softplus
heads predict nonnegative signal. Input length must be divisible by 64. N is encoded as zeros;
excess N is rejected. Training uses reverse complements and masked Huber loss in log1p space.
Each scale is the median positive **training** label; missing/all-zero heads are rejected.
AdamW and early stopping select a validation-loss checkpoint, not a biological accuracy claim.
Test rows do not influence fitting. Weights use safetensors, never pickle.

```bash
pace-livestock train-sequence --config sequence_training.yaml --out models/sequence
pace-livestock predict-sequence --config individual_run.yaml --out predictions/individual
```

External quantitative adapters can import `inputs.predictions` with an `external_table`
manifest; source-model hash, output type, species/assembly/context, target level, assay, signal
unit, normalization and target window are required. Imported predictions with a supplied
reference or individual genome must retain the matching genomic input configuration
and genome_binding_id; the software independently rechecks windows and binding. Peak probabilities and logits are rejected as quantitative intensity.
The manifest schema and fully working synthetic assets are visible in `examples/*/models/`.

## Contact prior

```yaml
data: contact_fit.tsv
model_id: contact_prior_research
species: chicken
assembly: assembly_identifier
context_id: liver
target_level: individual
is_synthetic: false
scale: balanced_contact_protocol_1
resolution: 5000
normalization_id: hic_norm_protocol_1
balancing: balanced
window_id: bin_pair
bin_edges: [5000, 10000, 20000, 50000, 100000, 500000, 5000000]
d_ref: 10000
d_min: 5000
```

The TSV has `bin_pair_id,distance_bp,contact_value,split,region_id`. Each measured pair appears
once, including real zeros. Train/test regions cannot overlap. Distance bins use all valid
values; positive mean bins enter a log-linear fit weighted by pair count. Empty and zero-mean
bins are reported separately. A nondecreasing fitted curve is rejected. The d_min platform is
an explicit near-distance rule. Parameters, bins, training regions and held-out residuals are
saved. This is a simple empirical prior, not an optimal count-noise model.

## Fusion calibration

```yaml
data: fusion_measurements.tsv
model_id: individual_fusion_research
species: chicken
assembly: assembly_identifier
context_id: liver
target_level: individual
is_synthetic: false
calibration_target: individual_state
signal_unit: normalized_signal
normalization_id: frozen_normalization_protocol
output_window: 500
scales: {ATAC: 1.0, H3K27ac: 1.0}
minimum_samples: 20
measurement_design: independent_measurements_from_the_same_animal
```

Scales and sample-count rules must come from a recorded training protocol; example values are
not generally calibrated choices. The TSV columns are `assay,quality_stratum,observed,predicted,
target,split,group_id,input_donor,target_donor,input_measurement_id,target_measurement_id`.
Splits are calibration/test. Source and target measurement IDs must differ. Individual-state
calibration requires matching donors. A population-mean calibrator cannot be used to restore
an individual state. These checks cannot prove laboratory independence: the experiment design
must establish it. The constrained least-squares fit is closed form. Identical observed and
predicted vectors have unidentifiable weight and an explicit single-source fallback.

## Independent elastic-net classifier

```yaml
data: labelled_edge_features.tsv
model_id: regulatory_link_classifier
is_synthetic: false
context: {species: chicken, assembly: assembly_identifier, context_id: liver}
extra_features: [H3K4me1, promoter_H3K4me3, methylation_M_site]
penalties: [[0.01, 0.01], [0.1, 0.01]]
folds: 3
seed: 17
calibrate: true
# feature_contract: insert the full mapping exported by the representative run
```

The YAML block is a template, not runnable until its paths and `feature_contract`
are supplied. For a non-synthetic model the full feature contract is required.
Copy the mapping from a compatible scoring run's `ml_feature_contract.json` into
training YAML; do not write a path string where a mapping is expected. For example:

```bash
python - <<'PYTHON'
import json, yaml
from pathlib import Path
config = yaml.safe_load(Path('classifier.yaml').read_text())
config['feature_contract'] = json.loads(Path('results/measured/ml_feature_contract.json').read_text())
Path('classifier_with_contract.yaml').write_text(yaml.safe_dump(config, sort_keys=False))
PYTHON
PACE train --config classifier_with_contract.yaml --out models/classifier
```

The contract includes activity panels/scales, contact resolution and normalization,
candidate definitions, target/estimand and auxiliary feature preprocessing.
[The bundled training config](../examples/training/learning.yaml) demonstrates a
complete synthetic contract.

Prepare one-to-one label mappings before assembling this table. Required columns:
`element_id,gene_id,assayed_region_id,mapping_count,group_id,split,label_status,effect_direction,
A_used,Cbar,distance_bp,pace_score,regime,activity_sources,contact_sources`, plus named extra
features and optional `sample_weight`. Splits are train/calibration/test. Enhancing positives
have `label_status=enhancing_positive,effect_direction=down`; negatives must be explicitly
`powered_negative`. Upregulation, low power, untested regions and ambiguous multi-tile mappings
are excluded and exported. These status labels assert the user's experiment-specific effect,
significance and power rules; PACE does not invent those rules from a p value.

Group, edge and perturbation-region identities cannot cross splits. Within training,
connected repeated entities also stay together across tuning folds. Grouped tuning refits
preprocessing inside each training fold; infeasible groups/classes cause an error. A single
predeclared penalty pair permits research fitting without pretending cross-validation occurred.
Core-missing edges are excluded; only extra features can be imputed. Training medians/IQR,
missing indicators and removed all-missing columns are frozen in JSON. The objective is mean
sample-weighted logistic loss + lambda1·L1 + lambda2·L2²/2, with unpenalized intercept, implemented
by proximal gradient. Convergence is reported. A base-only classifier is saved and evaluated.

Optional sigmoid calibration uses the independent calibration split, with a declared 1e-6
slope penalty for numerical stability. Without it, probability is NA and `pace_ml_score` is an
uncalibrated classifier score. Neither is averaged with PACE. Real calibration probabilities
remain specific to the recorded perturbation and candidate sampling design.

```yaml
# predict_ml.yaml
run: results/measured
model: models/classifier
features: results/measured/multiomics_features.tsv.gz
```

Inference checks the complete feature contract, regime, evidence-source and context scope.
Mismatches return out_of_scope with unavailable scores/probabilities. Old unbound models
are research-ineligible; a demonstration-only legacy score is explicitly unverified.
Synthetic probabilities, where present, are labelled synthetic_demonstration_only. Extra features must have unambiguous
entity-qualified names, and repeated annotation measurements need declared aggregation before
pivoting into a learning feature. Gene TPM is off by default; use it only as an explicit ablation.
