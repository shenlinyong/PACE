# Real results available for the manuscript placeholders (K562 / human only)

Computed by the scripts in this folder from public data; CIs: 1,000 gene-level bootstrap resamples. ABC = ENCODE K562 ABC predictions; GBM = in-house TargetFinder-style re-implementation.

## Benchmark (Fig. 2)

| Method | Combined AUPRC (95% CI) | Held-out AUPRC (95% CI) | Coverage (comb./held-out) |
|---|---|---|---|
| PACE | 0.592 (0.539–0.648) | 0.561 (0.455–0.664) | 0.956 / 0.990 |
| PACE, accessibility only | 0.559 (0.508–0.612) | 0.513 (0.402–0.628) | 0.960 / 0.991 |
| PACE, H3K27ac only | 0.538 (0.484–0.597) | 0.503 (0.403–0.610) | 0.956 / 0.990 |
| PACE, single canonical TSS | 0.600 (0.546–0.657) | 0.553 (0.448–0.655) | 0.956 / 0.990 |
| ABC (ENCODE K562 predictions) | 0.574 (0.523–0.625) | 0.519 (0.400–0.637) | 0.927 / 0.889 |
| ABC power-law contact (ENCODE) | 0.565 (0.512–0.618) | 0.479 (0.361–0.600) | 0.927 / 0.889 |
| Inverse distance | 0.438 (0.386–0.498) | 0.369 (0.270–0.487) | 1.000 / 1.000 |
| Distance prior P_s(d) | 0.367 (0.318–0.421) | 0.320 (0.228–0.422) | 1.000 / 1.000 |
| GBM (TargetFinder-style) | 0.541 (0.484–0.593) | 0.600 (0.494–0.705) | 1.000 / 1.000 |

PACE AUROC (combined) 0.923 (0.898–0.943); prevalence 0.045.

## Paired differences, PACE − method (pairs scored by both)

| set     | method                         |   n_common |   pos_common |   auprc_pace |   auprc_other |   delta |     lo |    hi |   p_boot |
|:--------|:-------------------------------|-----------:|-------------:|-------------:|--------------:|--------:|-------:|------:|---------:|
| heldout | ABC (ENCODE K562 predictions)  |       1688 |          104 |        0.596 |         0.564 |   0.033 | -0.023 | 0.092 |    0.276 |
| heldout | ABC power-law contact (ENCODE) |       1688 |          104 |        0.596 |         0.52  |   0.076 |  0.022 | 0.131 |    0.002 |
| heldout | Inverse distance               |       1899 |          114 |        0.579 |         0.377 |   0.202 |  0.111 | 0.29  |    0     |
| heldout | Distance prior P_s(d)          |       1899 |          114 |        0.579 |         0.326 |   0.253 |  0.181 | 0.331 |    0     |
| heldout | GBM (TargetFinder-style)       |       1899 |          114 |        0.579 |         0.597 |  -0.018 | -0.079 | 0.043 |    0.586 |
| main    | ABC (ENCODE K562 predictions)  |       9404 |          429 |        0.608 |         0.609 |  -0.001 | -0.037 | 0.036 |    0.932 |
| main    | ABC power-law contact (ENCODE) |       9404 |          429 |        0.608 |         0.599 |   0.009 | -0.023 | 0.04  |    0.566 |
| main    | Inverse distance               |       9904 |          455 |        0.611 |         0.446 |   0.166 |  0.114 | 0.214 |    0     |
| main    | Distance prior P_s(d)          |       9904 |          455 |        0.611 |         0.373 |   0.238 |  0.197 | 0.281 |    0     |
| main    | GBM (TargetFinder-style)       |       9904 |          455 |        0.611 |         0.549 |   0.062 |  0.021 | 0.107 |    0.006 |

## Degradation (Fig. 3a–c; pairs scored under all conditions, combined set)

| condition                 |   auprc |     lo |     hi |   n_common |   pos_common |
|:--------------------------|--------:|-------:|-------:|-----------:|-------------:|
| A                         |   0.611 |  0.559 |  0.665 |       9904 |          455 |
| B_same                    |   0.601 |  0.548 |  0.656 |       9904 |          455 |
| B_other                   |   0.601 |  0.548 |  0.656 |       9904 |          455 |
| C                         |   0.619 |  0.566 |  0.672 |       9904 |          455 |
| depth100_ice              |   0.608 |  0.556 |  0.662 |       9904 |          455 |
| res5000                   |   0.629 |  0.578 |  0.681 |       9904 |          455 |
| res10000                  |   0.611 |  0.559 |  0.665 |       9904 |          455 |
| res25000                  |   0.565 |  0.511 |  0.618 |       9904 |          455 |
| res40000                  |   0.516 |  0.461 |  0.571 |       9904 |          455 |
| rep1                      |   0.599 |  0.546 |  0.653 |       9904 |          455 |
| B_same_rep1               |   0.59  |  0.536 |  0.646 |       9904 |          455 |
| depth50_mean              |   0.608 |  0.557 |  0.662 |       9904 |          455 |
| depth25_mean              |   0.606 |  0.553 |  0.66  |       9904 |          455 |
| depth10_mean              |   0.606 |  0.555 |  0.66  |       9904 |          455 |
| depth5_mean               |   0.606 |  0.554 |  0.659 |       9904 |          455 |
| diff_A_minus_Bsame        |   0.01  | -0.002 |  0.022 |       9904 |          455 |
| diff_Bsame_minus_Bother   |  -0     | -0     |  0     |       9904 |          455 |
| diff_Bother_minus_C       |  -0.018 | -0.03  | -0.004 |       9904 |          455 |
| diff_A_minus_Bother       |   0.01  | -0.002 |  0.022 |       9904 |          455 |
| diff_2rep_minus_1rep      |   0.012 |  0.006 |  0.019 |       9904 |          455 |
| diff_depth50_minus_Bother |   0.007 | -0.005 |  0.019 |       9904 |          455 |
| diff_depth25_minus_Bother |   0.004 | -0.008 |  0.016 |       9904 |          455 |
| diff_depth10_minus_Bother |   0.005 | -0.007 |  0.017 |       9904 |          455 |
| diff_depth5_minus_Bother  |   0.005 | -0.008 |  0.017 |       9904 |          455 |
| mask5                     |   0.593 |  0.533 |  0.662 |       3586 |          142 |
| mask10                    |   0.605 |  0.551 |  0.662 |       4593 |          260 |
| mask20                    |   0.586 |  0.533 |  0.641 |       7763 |          346 |

## Missing activity (combined set)

| condition   |   auprc |    lo |    hi |   ref_auprc |   cov_default |   cov_cond |
|:------------|--------:|------:|------:|------------:|--------------:|-----------:|
| mask5       |   0.593 | 0.533 | 0.662 |       0.61  |         0.529 |      0.903 |
| mask10      |   0.605 | 0.551 | 0.662 |       0.618 |         0.353 |      0.865 |
| mask20      |   0.586 | 0.533 | 0.641 |       0.613 |         0.078 |      0.757 |

## Held-out repeat (Table S14 content)

| condition               |   auprc |     lo |    hi |
|:------------------------|--------:|-------:|------:|
| A                       |   0.579 |  0.473 | 0.684 |
| B_same                  |   0.558 |  0.444 | 0.666 |
| B_other                 |   0.558 |  0.444 | 0.665 |
| C                       |   0.563 |  0.451 | 0.673 |
| diff_A_minus_Bsame      |   0.021 |  0.002 | 0.04  |
| diff_Bsame_minus_Bother |  -0     | -0.002 | 0.001 |
| rep1                    |   0.565 |  0.457 | 0.674 |
| res40000                |   0.475 |  0.369 | 0.59  |
| mask5                   |   0.582 |  0.438 | 0.717 |
| mask10                  |   0.57  |  0.448 | 0.696 |
| mask20                  |   0.526 |  0.425 | 0.638 |

## Sensitivity (Table S11 content)

| parameter                    | value   | set     |   auprc |    lo |    hi |   coverage |   jaccard_top10 |
|:-----------------------------|:--------|:--------|--------:|------:|------:|-----------:|----------------:|
| default                      | default | main    |   0.592 | 0.539 | 0.648 |      0.956 |           1     |
| default                      | default | heldout |   0.561 | 0.455 | 0.664 |      0.99  |           1     |
| candidate window (Mb)        | 1.0     | main    |   0.591 | 0.537 | 0.648 |      0.91  |           0.933 |
| candidate window (Mb)        | 1.0     | heldout |   0.565 | 0.459 | 0.668 |      0.991 |           0.979 |
| candidate window (Mb)        | 2.0     | main    |   0.594 | 0.541 | 0.649 |      0.958 |           0.976 |
| candidate window (Mb)        | 2.0     | heldout |   0.566 | 0.462 | 0.669 |      0.991 |           0.979 |
| candidate window (Mb)        | 10.0    | main    |   0.578 | 0.525 | 0.637 |      0.954 |           0.968 |
| candidate window (Mb)        | 10.0    | heldout |   0.562 | 0.456 | 0.665 |      0.99  |           1     |
| candidate window (Mb)        | 5.0     | main    |   0.592 | 0.539 | 0.648 |      0.956 |           1     |
| candidate window (Mb)        | 5.0     | heldout |   0.561 | 0.455 | 0.664 |      0.99  |           1     |
| pseudocount distance d0 (kb) | 2.0     | main    |   0.592 | 0.539 | 0.648 |      0.956 |           1     |
| pseudocount distance d0 (kb) | 2.0     | heldout |   0.561 | 0.455 | 0.664 |      0.99  |           1     |
| pseudocount distance d0 (kb) | 10.0    | main    |   0.592 | 0.539 | 0.648 |      0.956 |           1     |
| pseudocount distance d0 (kb) | 10.0    | heldout |   0.561 | 0.455 | 0.664 |      0.99  |           1     |
| pseudocount distance d0 (kb) | 5.0     | main    |   0.592 | 0.539 | 0.648 |      0.956 |           1     |
| pseudocount distance d0 (kb) | 5.0     | heldout |   0.561 | 0.455 | 0.664 |      0.99  |           1     |
| element width (bp)           | 250     | main    |   0.593 | 0.536 | 0.65  |      0.952 |           0.891 |
| element width (bp)           | 250     | heldout |   0.548 | 0.442 | 0.651 |      0.99  |           0.854 |
| element width (bp)           | 1000    | main    |   0.582 | 0.529 | 0.642 |      0.961 |           0.857 |
| element width (bp)           | 1000    | heldout |   0.522 | 0.416 | 0.628 |      0.99  |           0.818 |
| element width (bp)           | 500     | main    |   0.592 | 0.539 | 0.648 |      0.956 |           1     |
| element width (bp)           | 500     | heldout |   0.561 | 0.455 | 0.664 |      0.99  |           1     |
| expression threshold (TPM)   | 0.5     | main    |   0.6   | 0.548 | 0.654 |      0.957 |           0.99  |
| expression threshold (TPM)   | 0.5     | heldout |   0.561 | 0.455 | 0.664 |      0.991 |           0.99  |
| expression threshold (TPM)   | 2.0     | main    |   0.57  | 0.515 | 0.626 |      0.955 |           0.979 |
| expression threshold (TPM)   | 2.0     | heldout |   0.558 | 0.448 | 0.66  |      0.989 |           0.99  |
| expression threshold (TPM)   | 5.0     | main    |   0.531 | 0.479 | 0.586 |      0.951 |           0.94  |
| expression threshold (TPM)   | 5.0     | heldout |   0.536 | 0.425 | 0.643 |      0.986 |           0.984 |
| expression threshold (TPM)   | 1.0     | main    |   0.592 | 0.539 | 0.648 |      0.956 |           1     |
| expression threshold (TPM)   | 1.0     | heldout |   0.561 | 0.455 | 0.664 |      0.99  |           1     |

## Not computed
All livestock-dependent placeholders (γ_pig/cattle/chicken, N_links, eQTL OR/concordance, IMVGI recovery, Fig. 3d–e, 4, 5, 6, S1, S3) and the EpiTensor, JEME, EPIPDLF and GATv2EPI comparators.
