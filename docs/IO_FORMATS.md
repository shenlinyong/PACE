# Output specification

PACE (Prediction of Activity-based regulatory Connections for Enhancers) writes one row per candidate enhancer–gene link after combining its TSS alternatives. Files are tab-separated; a `.gz` output suffix requests gzip compression through pandas.

## Files produced by each interface

| Command | Main output | Companion outputs |
| --- | --- | --- |
| `scripts/pace.py` | Path supplied to `--output`, all scored/unscored gene-level candidates | None |
| `scripts/calculate_pace_score.py` | `--output`, all candidates | The same path plus `.filtered.tsv` |
| `workflow/scripts/pace_predict.py` | `--output`, all candidates | A filtered companion when threshold filtering is requested |
| `workflow/scripts/pace_filter.py` | Selected rows at `--output`, compact columns | Selected rows with all fields at `<stem>_Full.tsv`, or the explicit `--full_output_file` path |
| `workflow/scripts/pace_metrics.py` | `QCSummary_<sample_name>.tsv` in `--output_dir` | `QCPlots_<sample_name>.pdf` |
| `scripts/smoke_test.py` | Independent-check report `smoke_report.json` | Inputs, per-step logs, predictions, QC and `SHA256SUMS.txt` |

**“Full” in the filter output means all columns for the selected rows. It does not mean all candidate rows.** Keep the predictor's unfiltered output as the primary record.

## Identity and score components

| Column | Meaning |
| --- | --- |
| `sample_id`, `chr`, `start`, `end`, `TargetGeneEnsemblID` | Link identity and normalization scope |
| `TargetGene` | Optional display name; not the normalization key |
| `PACE.Score` | `raw_support / (observed_mass + supplied residual support)`; `NA` when unscorable |
| `activity` | Aggregated enhancer activity |
| `contact_gene` | Contact averaged over the gene's contributing TSS weights |
| `target_share` | Enhancer-centred contact share over candidate target genes |
| `raw_support` | Activity × gene contact × target share raised to the allocation exponent |
| `observed_mass` | Sum of finite raw support for the gene |
| `unassigned_mass` | Supplied residual support, or missing if unknown |
| `score_scope` | `observed_candidates_only` or `residual_adjusted` |
| `unscored_candidates` | Count of supplied gene candidates with undefined raw support |
| `n_tss`, `TargetGeneTSSs` | Number and coordinates of contributing TSS alternatives |
| `model_version` | Kernel version string; also record the Git commit |

`TargetGeneTSS`, `distance`, `contact` and other single-TSS fields can come from one representative input row after aggregation. They are diagnostics, **not** a claim that this one promoter mediates a gene-level prediction. Use `contact_gene` for the score component and retain the full input TSS table for promoter-specific inspection.

When residual support is computationally zero, finite scores sum to one within the observed support of a gene. Omitted or unscored candidates can therefore inflate the apparent relative support of remaining candidates. A score of 1 is not proof of regulation.

## Contact states

| `contact_state` | Interpretation |
| --- | --- |
| `prior_only` | No usable measured-contact combination; distance prior used |
| `unscaled_observation` | Observed contact supplied without compatible expected contact |
| `quality_unknown` | Observed and expected contact supplied but reliability unknown |
| `matched_shrunk` | Qualified matched contact blended with the prior |
| `surrogate_shrunk` | Qualified surrogate contact blended with the prior |

Multiple TSS states can be combined using `|`. A `matched_shrunk` state with reliability zero still gives prior-only numerical contact; the label describes the supplied source and usable metadata, not the strength of its contribution. Inspect `contact_quality` as well.

## Evidence states

| Field or state | Interpretation |
| --- | --- |
| `activity_quality` | Quality over planned assay layers; unknown retained |
| `contact_quality` | TSS-weighted contribution of qualified contact reliability |
| `tss_quality` | Conservative promoter-annotation quality across contributing TSSs |
| `catalogue_quality` | Independently supplied candidate-catalogue quality |
| `evidence_quality` | Minimum of the four components; unknown if any is unknown |
| `sufficient_input_evidence` | All components meet the configured threshold and the gene has no unscored supplied candidates |
| `provisional` | Finite score with incomplete or subthreshold input evidence |
| `insufficient` | Unscorable link or explicitly zero activity/TSS quality |
| `evidence_reasons` | Semicolon-separated reasons such as `activity_quality_unknown`, `contact_quality_low`, `unscored_candidates`, `unassigned_mass_unknown` or `surrogate_contact` |

These states concern input evidence. They are not functional ground-truth labels. Review the score and the reasons together before experimental prioritization. The [quick-start fixture](QUICKSTART.md) deliberately contains missing activity and incomplete quality so these cases can be inspected.
