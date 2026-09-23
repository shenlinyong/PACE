# Migration to measured activity

Software 0.4.0 restricts PACE to experimental activity. The removed `hybrid`, `genome` and `genome_only` modes, sequence activity prediction, activity-fusion calibration and variant-effect scenarios are no longer available.

Removed commands: `train-sequence`, `predict-sequence`, `prepare-genome`, `fit-fusion`, `variant-effects`.
Removed run sections: `sequence`, `fusion`, `genome`; removed input: `inputs.predictions`.
Removed extras: `sequence`. PyTorch, safetensors and VCF/BCF readers are no longer required.

Old configurations using these entries fail validation. Do not convert a sequence-predicted activity column into an observed table or merely change a mode name. A new measured run requires actual qualified activity measurements and sample provenance.

Existing measured configurations remain supported. Use `pace run` or `pace measured`; mode defaults to measured. Imported activity now requires normalization_id, unit and window_id and is limited to observed/aggregate evidence. Contact observations, optional distance priors, multiple promoters, automatic eta and measured multiomics interfaces are retained. Reference FASTA access for CpG annotation remains available as `pace_livestock.io.reference.Reference`.

Rerun analyses with the current version before preparing new comparisons or benchmarks. Refit auxiliary classifier contracts if their measurement definitions changed. Record the exact software commit. Older outputs retain their original meaning and should not be relabelled as experimental evidence.

Historical source: [0.3.0-era snapshot](https://github.com/shenlinyong/PACE/tree/8c6f51e2ced192d26c17384893e4dd67305bdb25). Git history preserves it without exposing parallel current implementations.

## Experimental-contact and score reporting update

Only the lowercase `pace` executable is installed. Reinstall the package and update scripts that invoked `PACE` or `pace-livestock`. The Python module remains `pace_livestock`; package name remains `pace-livestock`.

Default incomplete primary scores are now NA. Use `pace_score_conditional` for the old measurable-subset quantity, or explicitly request `scoring.partial_policy: conditional`. Neither should be compared across different missing backgrounds without recomputing a common denominator. New sensitivity intervals are not confidence intervals.

Near-diagonal default is `prior_or_neighbor`; the previous `prior_or_unresolved` remains supported. Contact pseudocount defaults to auto and changes observed-mode scores when a compatible supplied prior exists. Set pseudocount:none for a predeclared unsmoothed control. Rerun analyses and retrain/revalidate frozen calibrators/classifiers when their evidence contracts differ.

`--force` only replaces identifiable PACE results, keeping the previous directory as a sibling backup. Root-level temporary planning/development documents have been replaced by maintained documentation.
