# Ferroptosis colorectal-cancer pipeline (single-file edition)

## Setup (Windows 10/11, VS Code)
```
cd C:\Users\Mayan\OneDrive\Documents\paper_ferroptosis
py -3.11 -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
```
(Python 3.10-3.12 recommended. Put ferro_pipeline.py in the project root.)

## Run
```
python ferro_pipeline.py --dry-run                    # check that sources are reachable (small requests)
python ferro_pipeline.py --accept-provisional-panel   # full run with the provisional gene panel
python ferro_pipeline.py --stage fetch                # single stage
python ferro_pipeline.py --stage prognosis benchmark
python ferro_pipeline.py --refresh                    # re-download
python ferro_pipeline.py --test                       # unit tests
```
## Gene panel
`config/ferroptosis_genes.csv` is generated as a PROVISIONAL literature-based candidate list. It is NOT your
original 35-gene panel. Replace it with your panel (columns: gene_symbol, gene_role[driver|suppressor|other],
source, evidence_or_notes, confirmed=True) or set `panel_confirmed: true` in config/config.yaml.
Without confirmation, the analysis stages need `--accept-provisional-panel` and all reports are labelled PROVISIONAL.

## Manual download fallbacks
If an automatic download fails, the log prints the exact URL and destination. Place the file in
`data/raw/gdsc/` or `data/raw/ctrp/` using the same file name and re-run; valid files are never re-downloaded.
 - GDSC2 response: https://www.cancerrxgene.org/downloads/bulk_download  -> data/raw/gdsc/GDSC2_fitted_dose_response.xlsx
 - GDSC expression (Cell_line_RMA_proc_basalExp.txt.zip) -> data/raw/gdsc/Cell_line_RMA_proc_basalExp.txt.zip
 - CTRP v2: https://portals.broadinstitute.org/ctrp.v2.1/ -> data/raw/ctrp/CTRPv2_ExpandedDataset.zip
Outputs: results/csv, results/figures, results/reports, results/provenance, results/logs.
