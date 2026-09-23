# Notation and output fields

All symbols refer to the [current formula](FORMULA.md).

| Symbol | Meaning | Software field or contract |
|---|---|---|
| E, G | Canonical scoring unit and target gene | `element_id`, `gene_id` |
| t | Deduplicated physical promoter/TSS | `promoter_id`, `tss0` |
| M | Fixed assay panel | `activity.panel` |
| x_star,m | Qualified measured assay signal | `resolved_activity.resolved_value` |
| A_star | Equal geometric activity of the fixed panel | `A_used` |
| r | Declared observed-contact mixing weight | Contact resolution policy |
| pi(t given G) | Frozen within-gene promoter weight | `promoters.pi` |
| Cbar | Weighted mean contact across required TSSs | `Cbar` |
| G(E) | Frozen candidate-gene set for E | `inputs.candidates` |
| B | Contact share across G(E) | `B`; skipped at eta zero |
| eta_used | Actual fixed or calibrated allocation exponent | `comparison_contract.eta`, `eta_calibration.json` |
| S | Unnormalized activity/contact/allocation support | `support`, `log_support` |
| E_score(G) | Actually scoreable candidates for G | Normalization identifiers and coverage |
| PACE(E,G) | Within-gene relative support | `pace_score` |

Auxiliary classifier outputs `pace_ml_score` and `pace_ml_probability` are separate.
QC and coverage are reported explicitly; no universal combined quality index is
multiplied into the main score. See the [data dictionary](data_dictionary.md).
