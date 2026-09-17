# Troubleshooting

PACE (Prediction of Activity-based regulatory Connections for Enhancers) reports invalid numerical inputs and preserves missing evidence. Start by checking the command, environment, assembly and unfiltered output.

| Symptom | Check / resolution |
| --- | --- |
| `conda: command not found` | Install Miniforge/Conda and initialize the shell; reopen it before creating the environment |
| Conda solve fails | Use the supplied recipe with strict priority; check network access and platform package availability. Avoid mixing default channels into this environment |
| `conda activate` is unavailable in a batch job | Initialize Conda for the job shell, or use `conda run -n pace python ...` |
| Existing environment name | Activate it and verify versions, or create a distinct prefix with `conda env create --prefix /path/to/env -f environment.yml` |
| Missing NumPy/pandas/pyBigWig/bedtools | Confirm `which python` and `which bedtools`; use the active `pace` environment. Re-run the import and test commands in [Installation](INSTALLATION.md) |
| Snakemake / PuLP API error | Use the supplied `pace-workflow` recipe with its compatible versions |
| Snakemake reports wildcard target errors | Include the explicit target `all` in the documented command; verify the repository commit |
| No candidate links | Match contig names, stable IDs and assembly; check TSS coordinates and the strict cis-distance window |
| Missing-column error from `pace.py` | Use a headered TSV with every required [input column](INPUTS.md#quantified-candidate-table), not a BED file |
| Activity/TSS conflict error | Deduplicate consistent records; do not attach different enhancer signal or TSS-use weight values to the same biological entity |
| Every score is `NA` for a gene | Inspect missing activity and total raw support; zero denominator is unscorable and must not be replaced by a fabricated score |
| Many rows are `provisional` | Expected with distance-only contact or unknown input QC. Read `evidence_reasons`; do not set unknown quality to 1 to obtain a preferred label |
| Hi-C file supplied but scores still use the prior | Inspect expected contact, reliability and source metadata. An available file alone is insufficient for contact mixing |
| Binary contact reader unavailable | Install the relevant optional reader and check the logs; the BEDPE example does not validate all binary file variants |
| H3K27ac/methylation option rejected | Use `missing_geometric`; primary file adapters reject inhibitory inputs. Explicit fraction-scale inhibition is a core-API extension, disabled by default |
| Filtered table is empty | Inspect the unfiltered scores and missingness before changing the threshold; verify that your biological decision cutoff is independently justified |
| `.filtered.tsv` or `_Full.tsv` seems smaller than expected | These contain selected rows. The main predictor output is the complete record |
| Unexpected effect from editing YAML | Consult [the interface mapping](PARAMETERS.md#5-which-configuration-entries-are-active); some compatibility keys are not forwarded |
| Long genome-wide runtime or high memory | Run one chromosome as a pilot. Partition on boundaries that preserve each gene's and enhancer's candidate background; arbitrary row chunks alter normalization |

To report a reproducible problem, include the Git commit, exact command, environment export, relevant error log and the smallest non-sensitive input that reproduces it in [GitHub Issues](https://github.com/shenlinyong/PACE/issues). Do not substitute a filtered result for the original input when diagnosing missing candidates.
