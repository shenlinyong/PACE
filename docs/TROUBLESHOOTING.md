# Troubleshooting

| Symptom | Check and correction |
|---|---|
| `PACE: command not found` | Activate the environment where the package was installed or add the installer prefix's bin directory to PATH |
| Python version error | Select Python >=3.11; the prefix installer accepts `--python` |
| Missing optional reader/model library | Install the matching `io` extra using the same environment's Python |
| Input path not found | CLI paths use the working directory; YAML paths use the YAML file's directory |
| Existing output directory | Select a new output path; outputs are not overwritten |
| Schema/unknown configuration error | Use the [current schemas](data_dictionary.md) and [parameters](parameters.md), not a retired configuration |
| Missing activity | Supply every declared assay for the fixed panel, or explicitly choose an appropriate separate single-layer run |
| Missing contact or B | Check required TSSs and the fixed candidate-gene set; do not remove missing genes to change allocation |
| All gene scores are NA | Inspect unresolved inputs and `gene_summary.tsv`; zero total support cannot be normalized |
| Primary scores are NA but conditional scores sum to one | Planned support is incomplete; inspect coverage, intervals and contact resolution reasons |
| Automatic eta remains zero | Read `eta_calibration.json`: no labels, inapplicable context, insufficient informative genes, or a valid boundary optimum can explain zero |
| Eta calibration scope mismatch | Use the same scientific context, catalog and evidence/model policy, or perform a new justified calibration |
| Synthetic asset rejected | Use demonstration mode for bundled fixtures; provide real appropriate assets for research |

Report reproducible bugs with the command, commit, minimal non-sensitive example and
QC reason. Software checks and their limitations are described in [validation](validation.md).
