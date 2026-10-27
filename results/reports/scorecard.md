# Scorecard (pipeline v1.5.0; 2026-10-13 17:43)

**Gene panel status: PROVISIONAL (not confirmed by the user)** (version 03d5812616df). Numbers from a PROVISIONAL panel must not be reported as the paper's results.

Reading rule: every ferroptosis result is judged against random gene sets of the same size run through the identical procedure.

## Output 1 - Ferroptosis classification
- `classification_dMMR`: completed - GSE39582: AUC=0.958 (OOF CV), empirical p vs 200 random sets=0.234
- `ext_ferroptosis_subtypes`: completed - TCGA-COAD: HR(high vs low)=0.87, log-rank p=0.643; GSE39582: HR(high vs low)=1.11, log-rank p=0.637; GSE17536: HR(high vs low)=0.52, log-rank p=0.0988; GSE14333: HR(high vs low)=0.45, log-rank p=0.424
- subtypes vs survival, panel vs random sets: TCGA-COAD emp p=0.42; GSE39582 emp p=0.64; GSE17536 emp p=0.37; GSE14333 emp p=0.20
- verdict: NOT distinguishable from random gene sets (smallest p across cohorts; not multiplicity-corrected)

## Output 2 - Prognosis prediction
- `prognosis_external_validation`: completed - GSE39582: C=0.539 [0.499,0.580]; GSE17536: C=0.486 [0.415,0.562]; GSE14333: C=0.581 [0.471,0.701]
- `novelty1_specificity_prognosis`: completed - GSE39582: C=0.539, empirical p=0.338; GSE17536: C=0.486, empirical p=0.831; GSE14333: C=0.581, empirical p=0.343

## Output 3 - Drug response (sensitive vs resistant, cell lines -> tumours)
- `ext_cellline_sensitivity_classifier`: completed - crispr/gpx4crispr: OOF AUC=0.760 [0.645,0.887], perm p=0.010, emp p vs 200 random sets=0.612
- `drug_sensitivity_classification`: not run
- verdict: NOT distinguishable from random gene sets (panel AUC 0.760 vs random-set median 0.772; label = gpx4crispr)

## Output 4 - IC50 / dose-response prediction
- `ext_dose_prediction`: completed - 5 drug models; 1 beat random sets (emp p<0.05); best gdsc/oxaliplatin: rho=0.48, emp p=0.045

| source | drug | response_unit | n_lines | spearman_oof | rho_lo95 | rho_hi95 | random_median_rho | empirical_p_vs_random_sets | q_bh |
|---|---|---|---|---|---|---|---|---|---|
| gdsc | fluorouracil | ln IC50 (uM) | 40 | 0.1265 | -0.198 | 0.4408 | 0.1833 | 0.6468 | 0.8085 |
| gdsc | irinotecan | ln IC50 (uM) | 40 | -0.03096 | -0.3379 | 0.3019 | -0.1083 | 0.3881 | 0.8085 |
| gdsc | oxaliplatin | ln IC50 (uM) | 40 | 0.4795 | 0.1834 | 0.7024 | 0.2009 | 0.04478 | 0.2239 |
| gdsc | sn38 | ln IC50 (uM) | 39 | -0.1711 | -0.4786 | 0.1644 | -0.1464 | 0.5572 | 0.8085 |
| crispr | gpx4crispr | gene effect (lower = more GPX4-dependent) | 59 | -0.05377 | -0.2779 | 0.2163 | 0.3337 | 0.9851 | 0.9851 |

- 1 of 5 drug models beat random gene sets at p<0.05 (0 after BH).

## Can the pipeline detect a real signal at all? (positive controls)
- `ext_positive_control`: completed - panel (ferroptosis): AUC=0.760 (emp p vs random=0.637) | EMT (mesenchymal+epithelial): AUC=0.590 (emp p vs random=0.980) | proliferation: AUC=0.671 (emp p vs random=0.771) | ceiling AUC=0.818 | stage HR/SD: TCGA-COAD=2.35, GSE39582=1.72, GSE17536=2.77, GSE14333=2.12

| target | signature | n_genes | n_lines | auc_oof | auc_lo95 | auc_hi95 | random_median_auc | empirical_p_vs_random_sets | n_random | single_score_auc_high_means_sensitive |
|---|---|---|---|---|---|---|---|---|---|---|
| crispr/gpx4crispr | panel (ferroptosis) | 36 | 59 | 0.7598 | 0.6452 | 0.8867 | 0.7759 | 0.6368 | 200 |  |
| crispr/gpx4crispr | EMT (mesenchymal+epithelial) | 24 | 59 | 0.5897 | 0.438 | 0.7451 | 0.7552 | 0.9801 | 200 | 0.3644 |
| crispr/gpx4crispr | proliferation | 12 | 59 | 0.6713 | 0.5308 | 0.7975 | 0.7218 | 0.7711 | 200 |  |
| crispr/gpx4crispr | data-driven ceiling (top 36 genes chosen inside CV folds) | 36 | 59 | 0.8184 |  |  |  |  |  |  |

| cohort | endpoint | n_surv | events | signature | n_genes | HR_per_SD | HR_lo95 | HR_hi95 | p | emp_p_abs_logHR_vs_random | random_median_abs_logHR |
|---|---|---|---|---|---|---|---|---|---|---|---|
| TCGA-COAD | OS | 434 | 95 | clinical stage (per SD) | 0 | 2.353 | 1.782 | 3.107 | 1.644e-09 |  |  |
| TCGA-COAD | OS | 434 | 95 | proliferation | 12 | 0.9001 | 0.7357 | 1.101 | 0.3067 | 0.2338 | 0.06145 |
| TCGA-COAD | OS | 434 | 95 | EMT (mesenchymal-epithelial) | 24 | 1.127 | 0.9372 | 1.355 | 0.2042 | 0.209 | 0.06565 |
| GSE39582 | OS | 573 | 190 | clinical stage (per SD) | 0 | 1.718 | 1.472 | 2.004 | 6.29e-12 |  |  |
| GSE39582 | OS | 573 | 190 | proliferation | 12 | 0.8691 | 0.7642 | 0.9883 | 0.03244 | 0.1841 | 0.07452 |
| GSE39582 | OS | 573 | 190 | EMT (mesenchymal-epithelial) | 24 | 1.215 | 1.065 | 1.387 | 0.003775 | 0.07463 | 0.09704 |
| GSE17536 | OS | 177 | 73 | clinical stage (per SD) | 0 | 2.772 | 2.068 | 3.717 | 9.251e-12 |  |  |
| GSE17536 | OS | 177 | 73 | proliferation | 12 | 1.047 | 0.8385 | 1.307 | 0.6863 | 0.7711 | 0.1268 |
| GSE17536 | OS | 177 | 73 | EMT (mesenchymal-epithelial) | 24 | 1.117 | 0.9226 | 1.353 | 0.2563 | 0.6269 | 0.1452 |
| GSE14333 | DFS | 153 | 29 | clinical stage (per SD) | 0 | 2.125 | 1.351 | 3.342 | 0.0011 |  |  |
| GSE14333 | DFS | 153 | 29 | proliferation | 12 | 0.8074 | 0.5747 | 1.134 | 0.2175 | 0.2985 | 0.1437 |
| GSE14333 | DFS | 153 | 29 | EMT (mesenchymal-epithelial) | 24 | 1.099 | 0.7651 | 1.579 | 0.6087 | 0.6667 | 0.1434 |

## Still missing / checklist
- [ ] Original gene panel supplied (--panel my_panel.csv)
- [ ] Drug-response file with a GPX4 inhibitor (ML162/ML210/RSL3)
- [x] IC50 (not only AUC/CRISPR) models exist
- [ ] Positive control recovered a real signal (EMT or ceiling AUC clearly above random)
- [x] All four tumour cohorts evaluated

In-vitro sensitivity is not a clinical dose or response. Associations are not mechanism.