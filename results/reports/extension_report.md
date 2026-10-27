# Extension results (v1.4.0; 2026-10-13 17:43)

> Cell-line sensitivity is in-vitro sensitivity, not a clinical dose or response. Associations are not mechanism.

## 1. Ferroptosis subtypes: associations
| cohort | n | n_high | endpoint | n_surv | events | logrank_p | logHR | HR_high_vs_low | HR_lo95 | HR_hi95 | cox_p | HR_adj_age_stage | HR_adj_lo95 | HR_adj_hi95 | p_adj | n_adj | pct_stage34_high | pct_stage34_low | stage34_fisher_p | age_median_high | age_median_low | age_mwu_p | pct_dMMR_high | pct_dMMR_low | mmr_fisher_p | panel_version |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| TCGA-COAD | 458 | 72 | OS | 434 | 95 | 0.6425 | -0.1384 | 0.8707 | 0.4849 | 1.564 | 0.643 | 0.83 | 0.3516 | 1.96 | 0.6708 | 350 | 42.62 | 40.32 | 0.7764 | 68.38 | 68.9 | 0.6859 |  |  |  | 03d5812616df |
| GSE39582 | 585 | 72 | OS | 573 | 190 | 0.6374 | 0.1015 | 1.107 | 0.7267 | 1.686 | 0.6365 | 1.003 | 0.6578 | 1.53 | 0.9882 | 568 | 52.78 | 45.76 | 0.3126 | 69 | 68 | 0.4391 | 3.077 | 15.92 | 0.003893 | 03d5812616df |
| GSE17536 | 177 | 26 | OS | 177 | 73 | 0.09884 | -0.6486 | 0.5228 | 0.2389 | 1.144 | 0.1046 | 0.4679 | 0.2127 | 1.029 | 0.05901 | 177 | 57.69 | 53.64 | 0.8319 | 66.5 | 66 | 0.8212 |  |  |  | 03d5812616df |
| GSE14333 | 195 | 16 | DFS | 153 | 29 | 0.4245 | -0.7912 | 0.4533 | 0.06178 | 3.326 | 0.4365 | 0.3916 | 0.05314 | 2.886 | 0.3576 | 153 | 56.25 | 48.04 | 0.6068 | 71 | 67 | 0.2884 |  |  |  | 03d5812616df |

Stability:
| k_prespecified | bootstrap_ARI_mean | bootstrap_ARI_p05 | n_genes | panel_version | silhouette_k2 | silhouette_k3 | silhouette_k4 | silhouette_k5 |
|---|---|---|---|---|---|---|---|---|
| 2 | 0.5531 | 0.151 | 36 | 03d5812616df | 0.1885 | 0.08906 | 0.07595 | 0.07124 |

## 2. Subtypes vs random gene sets (identical procedure)
| cohort | n_controls | control_set_size | seed | obs_logHR | obs_HR | emp_p_abs_logHR | random_median_abs_logHR | obs_dMMR_diff | emp_p_abs_dMMR_diff | random_median_abs_dMMR_diff |
|---|---|---|---|---|---|---|---|---|---|---|
| TCGA-COAD | 200 | 36 | 42 | -0.1384 | 0.8707 | 0.4229 | 0.1085 |  |  |  |
| GSE39582 | 200 | 36 | 42 | 0.1015 | 1.107 | 0.6368 | 0.1479 | -0.1285 | 0.2587 | 0.09625 |
| GSE17536 | 200 | 36 | 42 | -0.6486 | 0.5228 | 0.3682 | 0.5449 |  |  |  |
| GSE14333 | 200 | 36 | 42 | -0.7912 | 0.4533 | 0.204 | 0.3932 |  |  |  |

## 3. Cell-line GPX4-inhibitor sensitivity classifier
| label_source | expression_source | drugs | metric | label_rule | n_lines | n_sensitive | n_resistant | n_genes | cv | auc_oof | auc_lo95 | auc_hi95 | permutation_p | n_perm | random_set_median_auc | empirical_p_vs_random_sets | n_random | ferroscore_auc_high_means_sensitive | panel_version | panel_status | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| crispr | DepMap (OmicsExpressionProteinCodingGenesTPMLogp1.csv) | gpx4crispr | AUC | sensitive = composite within-drug AUC/IC50 percentile below the median across CRC lines (label only) | 59 | 29 | 30 | 36 | 5-fold x3 repeats, cell-line level | 0.7598 | 0.6452 | 0.8867 | 0.00995 | 200 | 0.7724 | 0.6119 | 200 | 0.3195 | 03d5812616df | PROVISIONAL (not confirmed by the user) | in-vitro sensitivity, not a clinical dose or response |

Drug-label consistency (Spearman):
_not available_

Coefficients:
| gene | coefficient | gene_role |
|---|---|---|
| DHODH | -0.4812 | suppressor |
| FANCD2 | -0.3832 | suppressor |
| ALOX12 | -0.2473 | driver |
| GLS2 | -0.227 | driver |
| SLC7A11 | -0.2247 | suppressor |
| NQO1 | -0.1978 | suppressor |
| CISD1 | -0.1948 | suppressor |
| FTL | -0.1648 | suppressor |
| GCH1 | -0.162 | suppressor |
| DPP4 | -0.1494 | driver |
| NOX1 | -0.138 | driver |
| SCD | -0.1307 | suppressor |
| TP53 | -0.1075 | driver |
| GPX4 | -0.09673 | suppressor |
| ACSL4 | -0.09531 | driver |
| GCLC | -0.07891 | suppressor |
| AIFM2 | -0.04896 | suppressor |
| GCLM | -0.02108 | suppressor |
| LPCAT3 | -0.01853 | driver |
| TXNRD1 | -0.0162 | suppressor |
| NFE2L2 | 0.02466 | suppressor |
| SLC3A2 | 0.04313 | suppressor |
| PEBP1 | 0.06697 | driver |
| POR | 0.07529 | driver |
| NCOA4 | 0.09128 | driver |
| FTH1 | 0.09747 | suppressor |
| BECN1 | 0.1006 | driver |
| HSPB1 | 0.1056 | suppressor |
| TFRC | 0.1065 | driver |
| ALOX15 | 0.1226 | driver |
| SLC40A1 | 0.1305 | suppressor |
| ACSL3 | 0.1476 | suppressor |
| SAT1 | 0.1663 | driver |
| AKR1C1 | 0.1765 | suppressor |
| ALOX5 | 0.2253 | driver |
| HMOX1 | 0.3798 | driver |

## 4. Projection of the frozen cell-line classifier onto tumours
| cohort | n | endpoint | n_surv | events | HR_per_SD | HR_lo95 | HR_hi95 | p | HR_per_SD_adj_age_stage | p_adj | n_adj | emp_p_abs_logHR_vs_random | random_median_abs_logHR | mean_score_high_subtype | mean_score_low_subtype | subtype_mwu_p | auc_for_dMMR | mmr_mwu_p |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| TCGA-COAD | 458 | OS | 434 | 95 | 0.9941 | 0.814 | 1.214 | 0.9536 | 0.9077 | 0.4539 | 350 | 0.9453 | 0.05316 | 0.6517 | 0.465 | 6.887e-09 |  |  |
| GSE39582 | 585 | OS | 573 | 190 | 1.117 | 0.9716 | 1.284 | 0.12 | 1.12 | 0.1145 | 568 | 0.592 | 0.1251 | 0.5575 | 0.475 | 0.01196 | 0.4058 | 0.0081 |
| GSE17536 | 177 | OS | 177 | 73 | 1.057 | 0.8421 | 1.326 | 0.6331 | 1.005 | 0.9647 | 177 | 0.6866 | 0.09892 | 0.5399 | 0.4832 | 0.3479 |  |  |
| GSE14333 | 195 | DFS | 153 | 29 | 1.236 | 0.8539 | 1.789 | 0.2616 | 1.143 | 0.4848 | 153 | 0.5423 | 0.2207 | 0.6208 | 0.4853 | 0.03986 |  |  |

## 5. How to read this
- An empirical p near or above 0.05 against random gene sets means the panel is not distinguishable from arbitrary genes for that task.
- Projection assumes cohort-wise z-scoring makes cell-line and tumour scales comparable; tumour stroma/purity is a confounder.
- Subtype and classifier definitions are frozen before any outcome of an external cohort is examined.

## 6. Run summary
- ext_ferroptosis_subtypes: completed - TCGA-COAD: HR(high vs low)=0.87, log-rank p=0.643; GSE39582: HR(high vs low)=1.11, log-rank p=0.637; GSE17536: HR(high vs low)=0.52, log-rank p=0.0988; GSE14333: HR(high vs low)=0.45, log-rank p=0.424
- ext_cellline_sensitivity_classifier: completed - crispr/gpx4crispr: OOF AUC=0.760 [0.645,0.887], perm p=0.010, emp p vs 200 random sets=0.612
- ext_positive_control: completed - panel (ferroptosis): AUC=0.760 (emp p vs random=0.637) | EMT (mesenchymal+epithelial): AUC=0.590 (emp p vs random=0.980) | proliferation: AUC=0.671 (emp p vs random=0.771) | ceiling AUC=0.818 | stage 
- ext_dose_prediction: completed - 5 drug models; 1 beat random sets (emp p<0.05); best gdsc/oxaliplatin: rho=0.48, emp p=0.045