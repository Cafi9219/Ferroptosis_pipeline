# Results summary (2026-10-13 17:43)

**Gene panel status: PROVISIONAL (not confirmed by the user)** (hash `03d5812616df`). All ferroptosis-specific statements below are provisional.

> This summary reports only quantities computed in this run. A finished script does not imply a successful analysis; see per-analysis status.

## 1. Data retrieval
| dataset | files | ok | failed |
|---|---|---|---|
| CTRPv2 | 1 | 0 | 1 |
| DepMap | 3 | 3 | 0 |
| GDSC | 1 | 0 | 1 |
| GDSC2 | 1 | 1 | 0 |
| GEO | 3 | 3 | 0 |
| GEO-platform | 1 | 1 | 0 |
| TCGA-COAD | 459 | 459 | 0 |

## 2. Dataset inventory
| dataset | data_type | platform | n_files | n_samples | n_patients | n_genes | pct_missing_expression | endpoint | events | note | n_probes | value_range | scale_note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| TCGA-COAD | RNA-seq STAR counts -> log2(TPM+1) | Illumina RNA-seq (GDC harmonized) | 524 | 458 | 458 | 19938 | 0 | OS | 95 | MSI not available; DFS not fetched |  |  |  |
| GSE39582 | microarray | GPL570 |  | 585 | 585 | 20846 | 0 | OS | 190 |  | 5.468e+04 | -0.40..16.11 | already log-scale (max<=100); no transform |
| GSE17536 | microarray | GPL570 |  | 177 | 177 | 20846 | 0 | OS | 73 |  | 5.468e+04 | 3.47..14.38 | already log-scale (max<=100); no transform |
| GSE14333 | microarray | GPL570 |  | 290 | 290 | 20846 | 0 | DFS | 50 | OS unavailable -> DFS/RFS used (endpoint differs from TCGA OS training; exploratory) | 5.468e+04 | -4.65..13.27 | already log-scale (max<=100); no transform |

## 3. Cohort inclusion / exclusion
| cohort | step | n_before | n_after | n_excluded | reason |
|---|---|---|---|---|---|
| TCGA-COAD | STAR-Counts open-access files returned | 524 | 524 | 0 | 458 cases, 514 samples, 524 aliquots |
| TCGA-COAD | primary tumour + one sample per patient | 524 | 458 | 66 | non-primary samples and extra primary samples removed |
| TCGA-COAD | expression file parsed | 458 | 458 | 0 | 0 files missing locally |
| TCGA-COAD | clinical record with survival fields | 458 | 455 | 3 | vital status/time missing |
| TCGA-COAD | valid OS (time>0, event known) | 455 | 434 | 21 | non-positive or missing follow-up |
| GSE39582 | samples in series matrix | 585 | 585 | 0 |  |
| GSE39582 | unique samples with expression | 585 | 585 | 0 | duplicate/all-NaN columns removed |
| GSE39582 | valid survival (OS) | 585 | 573 | 12 | time<=0 or event unknown |
| GSE17536 | samples in series matrix | 177 | 177 | 0 |  |
| GSE17536 | unique samples with expression | 177 | 177 | 0 | duplicate/all-NaN columns removed |
| GSE17536 | valid survival (OS) | 177 | 177 | 0 | time<=0 or event unknown |
| GSE14333 | samples in series matrix | 290 | 290 | 0 |  |
| GSE14333 | unique samples with expression | 290 | 290 | 0 | duplicate/all-NaN columns removed |
| GSE14333 | valid survival (none) | 290 | 0 | 290 | time<=0 or event unknown |
| GSE14333 | valid survival (DFS) | 290 | 226 | 64 | time<=0 or event unknown |

## 4. Analysis status
| analysis | status | n | key_result | note |
|---|---|---|---|---|
| prognosis_external_validation | completed | 903 | GSE39582: C=0.539 [0.499,0.580]; GSE17536: C=0.486 [0.415,0.562]; GSE14333: C=0.581 [0.471,0.701] | panel=PROVISIONAL (not confirmed by the user) |
| novelty1_specificity_prognosis | completed | 200 | GSE39582: C=0.539, empirical p=0.338; GSE17536: C=0.486, empirical p=0.831; GSE14333: C=0.581, empirical p=0.343 | 200 random sets of size 36, seed 42; panel=PROVISIONAL (not confirmed by the user) |
| classification_dMMR | completed | 536 | GSE39582: AUC=0.958 (OOF CV), empirical p vs 200 random sets=0.234 | labels from GEO characteristics; FerroScore coverage {'drivers_present': 16, 'drivers_total': 16, 'suppressors_present': 20, 'suppressors_total': 20, 'missing': ''}; panel=PROVISIONAL (not confirmed by the user) |
| ext_ferroptosis_subtypes | completed | 458 | TCGA-COAD: HR(high vs low)=0.87, log-rank p=0.643; GSE39582: HR(high vs low)=1.11, log-rank p=0.637; GSE17536: HR(high vs low)=0.52, log-rank p=0.0988; GSE14333: HR(high vs low)=0.45, log-rank p=0.424 | k=2 prespecified; stability ARI=0.55; 200 random sets; panel=PROVISIONAL (not confirmed by the user) |
| ext_cellline_sensitivity_classifier | completed | 59 | crispr/gpx4crispr: OOF AUC=0.760 [0.645,0.887], perm p=0.010, emp p vs 200 random sets=0.612 | projection onto 4 tumour cohorts; panel=PROVISIONAL (not confirmed by the user); in-vitro sensitivity only |
| ext_positive_control | completed | 200 | panel (ferroptosis): AUC=0.760 (emp p vs random=0.637) | EMT (mesenchymal+epithelial): AUC=0.590 (emp p vs random=0.980) | proliferation: AUC=0.671 (emp p vs random=0.771) | ceiling AUC=0.818 | stage HR/SD: TCGA-COAD=2.35, GSE39582=1.72, GSE17536=2.77, GSE14333=2.12 |  |
| ext_dose_prediction | completed | 218 | 5 drug models; 1 beat random sets (emp p<0.05); best gdsc/oxaliplatin: rho=0.48, emp p=0.045 | ridge on 36 panel genes; panel=PROVISIONAL (not confirmed by the user); in-vitro only |

## 5. Novelty 2: frozen external validation (prognosis)
| cohort | endpoint | model | scheme | n | events | cindex | cindex_lo95 | cindex_hi95 | hr_per_sd | hr_p | calibration_slope | logrank_p_high_vs_low | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| TCGA-COAD (CV, training-only selection) | OS | molecular | cross_validation | 350 | 61 | 0.5481 | 0.4816 | 0.6967 |  |  |  |  | mean over 10 CV test folds; interval = 2.5-97.5 pct of fold values (not a CI) |
| TCGA-COAD (CV, training-only selection) | OS | clinical | cross_validation | 350 | 61 | 0.7697 | 0.6611 | 0.8636 |  |  |  |  | mean over 10 CV test folds; interval = 2.5-97.5 pct of fold values (not a CI) |
| TCGA-COAD (CV, training-only selection) | OS | combined | cross_validation | 350 | 61 | 0.7736 | 0.639 | 0.8741 |  |  |  |  | mean over 10 CV test folds; interval = 2.5-97.5 pct of fold values (not a CI) |
| TCGA-COAD (training, in-sample; optimistic) | OS | molecular | in_sample | 434 | 95 | 0.6173 | 0.5494 | 0.6926 | 1.502 | 1.169e-05 | 18.62 | 0.0001668 |  |
| GSE39582 | OS | molecular | external_frozen | 573 | 190 | 0.539 | 0.4988 | 0.5804 | 1.129 | 0.07985 | 5.961 | 0.2174 |  |
| GSE39582 | OS | clinical | external_frozen | 568 | 189 | 0.675 | 0.6315 | 0.7175 | 1.923 | 4.401e-16 | 0.7139 |  | frozen TCGA coefficients; complete-case age+stage n=568 |
| GSE39582 | OS | combined | external_frozen | 568 | 189 | 0.6831 | 0.6415 | 0.7262 | 2.012 | 3.345e-17 | 0.7583 |  | frozen TCGA coefficients; complete-case age+stage n=568 |
| GSE17536 | OS | molecular | external_frozen | 177 | 73 | 0.4864 | 0.4147 | 0.5619 | 0.9666 | 0.7831 | -1.706 | 0.2928 |  |
| GSE17536 | OS | clinical | external_frozen | 177 | 73 | 0.7467 | 0.6881 | 0.8009 | 2.676 | 1.257e-11 | 0.9335 |  | frozen TCGA coefficients; complete-case age+stage n=177 |
| GSE17536 | OS | combined | external_frozen | 177 | 73 | 0.7248 | 0.6667 | 0.7815 | 2.372 | 5.085e-10 | 0.7815 |  | frozen TCGA coefficients; complete-case age+stage n=177 |
| GSE14333 | DFS | molecular | external_frozen | 153 | 29 | 0.5809 | 0.4706 | 0.7008 | 1.339 | 0.1198 | 14.01 | 0.1944 | ; endpoint mismatch (DFS vs OS-trained) - exploratory |
| GSE14333 | DFS | clinical | external_frozen | 153 | 29 | 0.624 | 0.5069 | 0.7193 | 1.499 | 0.04577 | 0.4524 |  | frozen TCGA coefficients; complete-case age+stage n=153 |
| GSE14333 | DFS | combined | external_frozen | 153 | 29 | 0.6416 | 0.5247 | 0.7402 | 1.582 | 0.0224 | 0.4769 |  | frozen TCGA coefficients; complete-case age+stage n=153 |

External cohorts evaluated: 3; cohorts whose bootstrap 95% CI for C-index excludes 0.5: 0. Association with survival is not evidence of clinical utility.

## 6. Novelty 1: ferroptosis specificity
| cohort | endpoint | signature | cindex | empirical_p_one_sided | percentile_vs_random | n_controls | note |
|---|---|---|---|---|---|---|---|
| GSE39582 | OS | ferroptosis_panel_frozen | 0.539 | 0.3383 | 66.5 | 200 | same procedure applied to controls |
| GSE39582 | OS | random_sets_median | 0.5278 |  |  | 200 | screening+ridge identical to panel |
| GSE39582 | OS | random_sets_95th_percentile | 0.5751 |  |  | 200 |  |
| GSE39582 | OS | proliferation_markers (author-compiled, exploratory) | 0.5698 |  |  | 200 | not a published validated signature |
| GSE39582 | OS | clinical_model_frozen | 0.675 |  |  | 200 | TCGA-frozen age+stage |
| GSE39582 | OS | combined_model_frozen | 0.6831 |  |  | 200 | TCGA-frozen age+stage+risk |
| GSE17536 | OS | ferroptosis_panel_frozen | 0.4864 | 0.8308 | 17 | 200 | same procedure applied to controls |
| GSE17536 | OS | random_sets_median | 0.5324 |  |  | 200 | screening+ridge identical to panel |
| GSE17536 | OS | random_sets_95th_percentile | 0.594 |  |  | 200 |  |
| GSE17536 | OS | proliferation_markers (author-compiled, exploratory) | 0.4978 |  |  | 200 | not a published validated signature |
| GSE17536 | OS | clinical_model_frozen | 0.7467 |  |  | 200 | TCGA-frozen age+stage |
| GSE17536 | OS | combined_model_frozen | 0.7248 |  |  | 200 | TCGA-frozen age+stage+risk |
| GSE14333 | DFS | ferroptosis_panel_frozen | 0.5809 | 0.3433 | 66.5 | 200 | same procedure applied to controls |
| GSE14333 | DFS | random_sets_median | 0.555 |  |  | 200 | screening+ridge identical to panel |
| GSE14333 | DFS | random_sets_95th_percentile | 0.6465 |  |  | 200 |  |
| GSE14333 | DFS | proliferation_markers (author-compiled, exploratory) | 0.5088 |  |  | 200 | not a published validated signature |
| GSE14333 | DFS | clinical_model_frozen | 0.624 |  |  | 200 | TCGA-frozen age+stage |
| GSE14333 | DFS | combined_model_frozen | 0.6416 |  |  | 200 | TCGA-frozen age+stage+risk |

- GSE39582: panel C=0.539; empirical p=0.338 vs 200 random sets (not distinguishable from random at 0.05)
- GSE17536: panel C=0.486; empirical p=0.831 vs 200 random sets (not distinguishable from random at 0.05)
- GSE14333: panel C=0.581; empirical p=0.343 vs 200 random sets (not distinguishable from random at 0.05)

## 7. Classification (dMMR vs pMMR)
| cohort | scheme | model | n | n_pos | auc | auc_lo95 | auc_hi95 | balanced_accuracy | f1 | empirical_p_vs_random | note |
|---|---|---|---|---|---|---|---|---|---|---|---|
| GSE39582 | oof_cv_patient_level | ferroptosis_panel_LR | 536 | 77 | 0.9575 | 0.9324 | 0.9757 | 0.8839 | 0.7174 | 0.2338 |  |
| GSE39582 | oof_cv_patient_level | proliferation_markers_LR | 536 | 77 | 0.9336 | 0.8994 | 0.962 | 0.887 | 0.6699 |  | exploratory |
| GSE39582 | unsupervised_score | FerroScore_raw (driver-suppressor), orientation not prespecified | 536 | 77 | 0.4526 |  |  |  |  |  | AUC<0.5 means higher FerroScore in pMMR; descriptive only |

## 8. Novelty 3: FerroScore vs drug response (GDSC2, CRC cell lines)
_not run_

### Drug models
_not available_

### Drug-sensitivity classification (prespecified drugs)
_not available_

### CTRP replication
_not available_

## 9. Deviations / limitations flagged during the run
- TCGA MSI/MMR labels and DFS not retrieved (no standard open-access GDC field); MSI classification uses GEO cohorts with MMR labels.
- GSE14333: no survival endpoint with enough events
- Overlap warning: possible duplicate samples GSE17536 vs GSE14333 (r>0.9 on cohort-z top-variance genes) = 95 -> external cohorts may not be fully independent of each other
- Gene panel is PROVISIONAL (not the original project panel); all ferroptosis findings are provisional.
- GSE14333: OS unavailable -> DFS/RFS used (endpoint differs from TCGA OS training; exploratory)
- 95 GEO samples duplicate an earlier cohort and are removed from analysis (data/processed/overlap_drop.csv); re-run prognosis/benchmark/classification after preprocess.

## 10. Stage execution status
| stage | status | seconds | note |
|---|---|---|---|
| fetch | ok | 91.9 |  |
| preprocess | ok | 160.9 |  |
| prognosis | ok | 106.1 |  |
| benchmark | ok | 760.7 |  |
| classification | ok | 57.8 |  |
| drug | skipped/unavailable | 0 | GDSC2 processed data missing; run fetch+preprocess (or place manual files) |
| figures | ok | 14.5 |  |
| report | ok | 0.1 |  |
| subtypes | ok | 45.5 |  |
| cellline | ok | 90.5 |  |
| ext_figures | ok | 8.1 |  |
| ext_report | ok | 0 |  |
| poscontrol | ok | 130.2 |  |
| dose | ok | 38.3 |  |
| scorecard | ok | 0.5 |  |

No dose prediction or treatment recommendation is made; cell-line IC50 is not a clinical dose.