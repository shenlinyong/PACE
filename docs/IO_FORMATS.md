# Input and output formats

See the maintained [data dictionary](data_dictionary.md) for every supported table
and the [command guide](cli.md) for file arguments. Main outputs are `scores.tsv.gz`,
`gene_summary.tsv`, resolved evidence tables, `eta_calibration.json`, QC, manifest
and resolved configuration.

The primary field is `pace_score`, computed by [the current formula](FORMULA.md).
`pace_ml_score` and `pace_ml_probability` belong to a separate classifier.
Missing values are NA, measured zeros are zero, and partial normalization is explicit.
