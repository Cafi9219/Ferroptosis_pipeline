# Data limitations

- TCGA MSI/MMR and DFS were not retrieved; MSI classification relies on MMR labels parsed from GEO characteristics (GSE39582 typically).
- RNA-seq (TCGA) and microarray (GEO) differ in technology and scale; external evaluation uses cohort-specific z-scoring of selected genes, which assumes similar cohort case-mix.
- GEO clinical fields are parsed heuristically (regex over characteristics keys); check `data_dictionary.csv`. Time units are assumed months when not stated.
- GSE14333 typically lacks OS; if so, DFS is used and flagged (endpoint mismatch).
- GEO cohorts may share patients (see train_validation_overlap_checks.csv); independence is not guaranteed.
- The gene panel may be PROVISIONAL; replace it with the original project panel before drawing conclusions.
- GDSC2/CTRP: few CRC lines (tens); ML162/ML210/RSL3 may be absent from GDSC2 and/or CTRP; CTRP lines/compounds matched by name.
- Cell-line IC50/AUC is not a clinical dose; association is not mechanism.
- TCGA-COAD contains colon tumours only, whereas GDSC lineage COREAD also includes rectal lines.
- GEO samples duplicating an earlier cohort (r>0.9) are removed from analysis (data/processed/overlap_drop.csv); near-duplicates below the threshold are not detected.
- The drug-sensitivity classifier uses a median split of cell-line LN_IC50 (label only); n is small, so AUCs have wide intervals.

## Deviations recorded in this run
- TCGA MSI/MMR labels and DFS not retrieved (no standard open-access GDC field); MSI classification uses GEO cohorts with MMR labels.
- GSE14333: no survival endpoint with enough events
- Overlap warning: possible duplicate samples GSE17536 vs GSE14333 (r>0.9 on cohort-z top-variance genes) = 95 -> external cohorts may not be fully independent of each other
- Gene panel is PROVISIONAL (not the original project panel); all ferroptosis findings are provisional.
- GSE14333: OS unavailable -> DFS/RFS used (endpoint differs from TCGA OS training; exploratory)
- 95 GEO samples duplicate an earlier cohort and are removed from analysis (data/processed/overlap_drop.csv); re-run prognosis/benchmark/classification after preprocess.