# PACE workflow scripts

The shared kernel is `pace_core.py`. `predictor.py` and the standalone calculator use the same model definitions. Current usage and required fields are documented in [the tutorial](../../docs/TUTORIAL.md) and [the input contract](../../docs/INPUTS.md).

The ordered command-line stages are `pace_candidate_regions.py`, `pace_neighborhoods.py`, `pace_predict.py`, `pace_filter.py` and `pace_metrics.py`. Run `bash example/run_example_direct.sh` from the repository root for a bundled read example. Activity uses `missing_geometric`; RNA may annotate eligibility but does not multiply scores.
