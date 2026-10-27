#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
FERROPTOSIS END-TO-END PIPELINE (single-file edition, v1.3.2)
Integrated pipeline for ferroptosis classification, prognosis prediction and
drug-response analysis in colorectal cancer.

Usage (from the project root, in the VS Code terminal):
    python ferro_pipeline.py                       # everything
    python ferro_pipeline.py --dry-run             # validate sources only, no big downloads
    python ferro_pipeline.py --stage fetch         # one stage (fetch preprocess prognosis
                                            #   benchmark classification drug figures report)
    python ferro_pipeline.py --test                # unit tests
    python ferro_pipeline.py --refresh             # re-download everything
    python ferro_pipeline.py --accept-provisional-panel   # allow the provisional gene panel (see README)
    python ferro_pipeline.py --panel my_panel.csv  # import YOUR original panel (marks it CONFIRMED)
    python ferro_pipeline.py --inspect-geo GSE14333  # print GEO characteristics keys and exit

FILE NAME: this file is ferro_pipeline.py (renamed from run_all.py so it can never be confused with an older copy).
ONE COMMAND:   python ferro_pipeline.py --final
    (add  --panel my_panel.csv  to use your original panel; otherwise the provisional panel is used and labelled so)
    --final = remove stale penalizer_grid from config.yaml (backup kept) + clear cached random controls
              + run every stage + print the key result tables. No other file is needed.

Changes in v1.5.0:
  * new stages: poscontrol (positive controls: EMT/proliferation/stage + nested data-driven ceiling), dose (IC50/AUC regression for
    ferroptosis inducers + CRC chemotherapy, random-set controls, projection to tumours), scorecard (one-page mapping of the four outputs)
  * new flag --dose-drugs

Changes in v1.4.0:
  * new stages: subtypes (k=2 ferroptosis subtypes of tumours, frozen centroids, projected to GEO cohorts, tested vs
    survival/stage/age/MMR), cellline (GPX4-inhibitor sensitive vs resistant CRC cell lines -> panel classifier ->
    projection onto tumours), ext_figures (pipeline overview, subtype KM, classifier/projection benchmark), ext_report
  * every new test is benchmarked against random gene sets of the same size (identical procedure)
  * cell-line expression may come from GDSC RMA OR DepMap (data/raw/depmap/); new flags --check-data --drugs
    --response-source --metric --min-lines --n-perm

Changes in v1.3.2:
  * renamed to ferro_pipeline.py; --final now implies --reprocess (GEO clinical files are always rebuilt)
  * --final prints a DATA SANITY block (flags cohorts where every evaluated patient is an event, or where most
    samples were dropped) and saves all final tables to results/reports/final_tables.txt
  * startup banner prints the exact file that is running

Changes in v1.3.1:
  * BUG FIX: GEO event coding "no death" / "no recurrence" was not recognised (-> NaN), so GSE17536 lost all censored
    patients (73 of 177 kept, all events). Now any value starting with "no " is a non-event. After updating run:
        python ferro_pipeline.py --final --reprocess

Changes in v1.3.0:
  * --final flag (see above); run_final.py / fix_config.py are no longer needed

Changes in v1.2.2:
  * download fallbacks that are ALWAYS tried after the URLs stored in config/config.yaml (an old config.yaml keeps old,
    dead URLs and silently overrides code defaults): GDSC2 fitted dose-response from the Sanger FTP current release
    (24Jul22 file confirmed in the FTP listing). The GDSC RMA expression and CTRP fallbacks are unverified guesses;
    if they fail, place the files manually (see README).
  * preprocess skips cohorts whose processed files already exist (use --reprocess to force) - saves ~5 minutes per run
  * --reprocess flag

Changes in v1.2.1:
  * stale results: when a stage is skipped/failed in a run, its old entries in the "Analyses" summary are marked
    "not run in this run" instead of still showing "completed" from an earlier run
  * random-gene-set cache is fingerprinted (external cohort samples after overlap removal + panel + penalizer grid +
    seed + set size); a changed fingerprint discards cached prognosis controls automatically
  * warning when a GEO DFS/OS event rate looks implausible or a censoring flag was inverted (verify with --inspect-geo)

Changes in v1.2.0:
  * clinical CRC drugs (5-fluorouracil, oxaliplatin, irinotecan) added as a PRESPECIFIED 'clinical' tier
    (own BH family together with the primary ferroptosis inducers); user drug_panel.csv is auto-merged with defaults
  * new drug-sensitivity CLASSIFICATION (sensitive = LN_IC50 below the drug median) for prespecified drugs,
    out-of-fold AUC for panel-LR and FerroScore-LR  -> results/csv/drug_sensitivity_classification.csv
  * GEO cohorts: samples that duplicate a sample of an earlier cohort (r>0.9) are written to
    data/processed/overlap_drop.csv and removed from external validation (config: drop_overlapping_samples)
  * warning when config.yaml overrides the penalizer grid; wording: cell-line sensitivity, never a clinical dose
  * new unit test (tier merge)

Changes in v1.1.0:
  * penalizer grid extended to [0.05, 0.5, 5, 50, 500]
  * TCGA follow-up time also taken from follow_ups.days_to_follow_up (recovers alive patients)
  * --panel importer, --inspect-geo inspector
  * manual GDSC/CTRP file adoption (files with different names are copied to canonical names)

First run creates: config/, data/, results/ and the config/requirements/README files.
Nothing is simulated: if a dataset cannot be retrieved, the affected analysis is skipped
and the reason is reported.
"""
from __future__ import annotations

import argparse, copy, datetime as dt, gzip, hashlib, io, json, logging, math, os, re, shutil, sys, time
import traceback, unittest, warnings, zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

warnings.filterwarnings("ignore")
try:
    import numpy as np
    import pandas as pd
    import requests
    import yaml
    from scipy import stats
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.linear_model import LogisticRegressionCV, RidgeCV
    from sklearn.model_selection import StratifiedKFold, KFold, RepeatedStratifiedKFold
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.metrics import (roc_auc_score, roc_curve, confusion_matrix, mean_absolute_error,
                                 mean_squared_error, r2_score, precision_recall_fscore_support)
    from lifelines import CoxPHFitter, KaplanMeierFitter
    from lifelines.statistics import logrank_test
    from lifelines.utils import concordance_index
    from sklearn.cluster import KMeans
    from sklearn.metrics import adjusted_rand_score, silhouette_score
    from matplotlib.patches import FancyBboxPatch
except ImportError as _e:  # pragma: no cover
    sys.exit(f"Missing Python package '{_e.name}'. Run:  pip install -r requirements.txt")

# ----------------------------------------------------------------------------------------
# Paths, logging, state
# ----------------------------------------------------------------------------------------
ROOT = Path(os.environ.get("FERRO_ROOT", Path(__file__).resolve().parent))
D = {k: ROOT / v for k, v in dict(
    cfg="config", raw="data/raw", proc="data/processed", csv="results/csv", fig="results/figures",
    rep="results/reports", log="results/logs", prov="results/provenance", tests="tests").items()}
LOG = logging.getLogger("ferro")
VERSION = "1.5.0"
COHORT_TCGA = "TCGA-COAD"
GDC = "https://api.gdc.cancer.gov"
SESSION = requests.Session()
SESSION.headers.update({"User-Agent": f"ferroptosis-pipeline/{VERSION} (research; python-requests)"})


class DataUnavailable(Exception):
    """A dataset could not be retrieved/parsed safely. The dependent analysis is skipped."""


class ModelFitError(Exception):
    pass


def setup_logging():
    D["log"].mkdir(parents=True, exist_ok=True)
    LOG.setLevel(logging.INFO)
    if LOG.handlers:
        return
    fmt = logging.Formatter("%(asctime)s | %(levelname)-7s | %(message)s", "%Y-%m-%d %H:%M:%S")
    fh = logging.FileHandler(D["log"] / "run.log", mode="a", encoding="utf-8")
    sh = logging.StreamHandler(sys.stdout)
    for h in (fh, sh):
        h.setFormatter(fmt)
        LOG.addHandler(h)


def ensure_dirs():
    for k in ("cfg", "raw", "proc", "csv", "fig", "rep", "log", "prov", "tests"):
        D[k].mkdir(parents=True, exist_ok=True)
    for sub in ("tcga", "geo", "gdsc", "ctrp"):
        (D["raw"] / sub).mkdir(parents=True, exist_ok=True)


def read_state():
    p = D["prov"] / "state.json"
    base = {"summary": {}, "deviations": [], "stages": {}, "datasets": {}, "incl": [], "inventory": {}}
    if p.exists():
        try:
            base.update(json.loads(p.read_text(encoding="utf-8")))
        except Exception:
            pass
    return base


def write_state(s):
    D["prov"].mkdir(parents=True, exist_ok=True)
    (D["prov"] / "state.json").write_text(json.dumps(s, indent=1, default=str), encoding="utf-8")


def summ(name, status, n="", metric="", note=""):
    s = read_state()
    s["summary"][name] = dict(analysis=name, status=status, n=n, key_result=metric, note=note)
    write_state(s)


def deviation(msg):
    s = read_state()
    if msg not in s["deviations"]:
        s["deviations"].append(msg)
    write_state(s)


def incl(cohort, step, before, after, reason=""):
    s = read_state()
    row = dict(cohort=cohort, step=step, n_before=before, n_after=after, n_excluded=before - after, reason=reason)
    for i, r in enumerate(s["incl"]):
        if r["cohort"] == cohort and r["step"] == step:
            s["incl"][i] = row
            break
    else:
        s["incl"].append(row)
    write_state(s)


def inventory(name, **kw):
    s = read_state()
    s["inventory"][name] = dict(dataset=name, **kw)
    write_state(s)


def save_csv(df, name, index=False):
    D["csv"].mkdir(parents=True, exist_ok=True)
    df.to_csv(D["csv"] / name, index=index)
    LOG.info("saved results/csv/%s (%d rows)", name, len(df))


def read_csv_opt(name):
    p = D["csv"] / name
    return pd.read_csv(p) if p.exists() else None


def norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def deep_update(a, b):
    for k, v in b.items():
        if isinstance(v, dict) and isinstance(a.get(k), dict):
            deep_update(a[k], v)
        else:
            a[k] = v
    return a


# ----------------------------------------------------------------------------------------
# Default configuration & project files (written on first run, never overwritten)
# NOTE: config/config.yaml, once written, OVERRIDES these defaults. Edit that file too.
# ----------------------------------------------------------------------------------------
DEFAULT_CFG = dict(
    seed=42, tcga_project="TCGA-COAD", workflow_type="STAR - Counts",
    panel_confirmed=False,
    p_screen=0.10, min_features=3, penalizer_grid=[0.05, 0.5, 5.0, 50.0, 500.0], inner_folds=3,
    cv_folds=5, cv_repeats=2, min_events=10, n_boot=500,
    n_random_prognosis=200, n_random_classification=200,
    min_panel_coverage=0.8, drop_overlapping_samples=True,
    gdsc_lineage_field="TCGA_DESC", gdsc_lineage_value="COREAD",
    gdsc_min_cell_lines=15, model_min_cell_lines=30, max_model_drugs=10,
    ctrp_site="large_intestine", expressed_quantile_floor=0.5,
    geo=dict(cohorts=["GSE39582", "GSE17536", "GSE14333"], overrides={}),
    proliferation_markers=["MKI67", "PCNA", "TOP2A", "MCM2", "CCNB1", "CDK1", "AURKA", "BIRC5", "UBE2C",
                           "CDC20", "PLK1", "MYBL2", "CCNE1", "BUB1", "TYMS", "RRM2"],
    urls=dict(
        gdsc2_response=["https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/GDSC2_fitted_dose_response_27Oct23.xlsx"],
        gdsc_expression=["https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Data/preprocessed/Cell_line_RMA_proc_basalExp.txt.zip",
                         "https://cog.sanger.ac.uk/cancerrxgene/GDSC_release6.0/Cell_line_RMA_proc_basalExp.txt.zip"],
        ctrp=["https://ctd2-data.nci.nih.gov/Public/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset/CTRPv2.0_2015_ctd2_ExpandedDataset.zip"],
    ),
)

FALLBACK_URLS = dict(
    gdsc2_response=["https://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/GDSC2_fitted_dose_response_24Jul22.xlsx"],
    gdsc_expression=["https://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/Cell_line_RMA_proc_basalExp.txt.zip"],  # unverified
    ctrp=[],
)


def urls_for(cfg, key):
    """URLs from config.yaml first, then code fallbacks (de-duplicated, order preserved)."""
    out = list(cfg["urls"].get(key, []))
    out += [u for u in FALLBACK_URLS.get(key, []) if u not in out]
    return out


PANEL_ROWS = (
    [(g, "driver") for g in ["ACSL4", "LPCAT3", "ALOX15", "ALOX12", "ALOX5", "POR", "PEBP1", "NCOA4", "TFRC",
                              "HMOX1", "SAT1", "NOX1", "TP53", "GLS2", "DPP4", "BECN1"]] +
    [(g, "suppressor") for g in ["GPX4", "SLC7A11", "AIFM2", "GCH1", "NFE2L2", "FTH1", "FTL", "SLC3A2", "HSPB1",
                                  "ACSL3", "SCD", "DHODH", "NQO1", "GCLC", "GCLM", "TXNRD1", "SLC40A1", "CISD1",
                                  "FANCD2", "AKR1C1"]])

DRUG_PANEL = """drug_name,tier,note
ML162,primary,canonical GPX4-targeting ferroptosis inducer (prespecified)
ML210,primary,canonical GPX4-targeting ferroptosis inducer (prespecified)
RSL3,primary,canonical GPX4 inhibitor (prespecified)
Erastin,exploratory,system xc- inhibitor; exploratory only
Sulfasalazine,exploratory,system xc- inhibitor (weak); exploratory only
Sorafenib,exploratory,reported ferroptosis inducer in some contexts; exploratory only
5-Fluorouracil,clinical,standard CRC chemotherapy (prespecified clinical comparator)
Oxaliplatin,clinical,standard CRC chemotherapy (prespecified clinical comparator)
Irinotecan,clinical,standard CRC chemotherapy (prespecified clinical comparator)
"""

REQUIREMENTS = """numpy==1.26.4
pandas==2.2.2
scipy==1.13.1
scikit-learn==1.5.1
lifelines==0.29.0
matplotlib==3.9.2
requests==2.32.3
PyYAML==6.0.2
openpyxl==3.1.5
"""

README = """# Ferroptosis colorectal-cancer pipeline (single-file edition)

## Setup (Windows 10/11, VS Code)
```
cd C:\\Users\\Mayan\\OneDrive\\Documents\\paper_ferroptosis
py -3.11 -m venv .venv
.venv\\Scripts\\activate
pip install -r requirements.txt
```
(Python 3.10-3.12 recommended. Put ferro_pipeline.py in the project root.)

## Run
```
python ferro_pipeline.py --dry-run                    # check that sources are reachable (small requests)
python ferro_pipeline.py --accept-provisional-panel   # full run with the provisional gene panel
python ferro_pipeline.py --panel my_panel.csv         # import your original panel (CONFIRMED) and run
python ferro_pipeline.py --inspect-geo GSE14333       # list GEO characteristics keys (after fetch)
python ferro_pipeline.py --stage fetch                # single stage
python ferro_pipeline.py --stage prognosis benchmark
python ferro_pipeline.py --refresh                    # re-download
python ferro_pipeline.py --test                       # unit tests
```
## Gene panel
`config/ferroptosis_genes.csv` is generated as a PROVISIONAL literature-based candidate list. It is NOT your
original panel. Replace it with `--panel my_panel.csv` (columns: gene_symbol, gene_role[driver|suppressor|other],
optional source, evidence_or_notes) or set `panel_confirmed: true` in config/config.yaml.
Without confirmation, the analysis stages need `--accept-provisional-panel` and all reports are labelled PROVISIONAL.

## Manual download fallbacks
If an automatic download fails, the log prints the exact URL and destination. Place the file in
`data/raw/gdsc/` or `data/raw/ctrp/` (any reasonable file name is adopted) and re-run; valid files are never re-downloaded.
 - GDSC2 response: https://www.cancerrxgene.org/downloads/bulk_download  -> data/raw/gdsc/GDSC2_fitted_dose_response.xlsx
 - GDSC expression (Cell_line_RMA_proc_basalExp.txt.zip) -> data/raw/gdsc/Cell_line_RMA_proc_basalExp.txt.zip
 - CTRP v2: https://portals.broadinstitute.org/ctrp.v2.1/ -> data/raw/ctrp/CTRPv2_ExpandedDataset.zip
Outputs: results/csv, results/figures, results/reports, results/provenance, results/logs.
"""


def ensure_project():
    ensure_dirs()
    (D["cfg"] / "config.yaml").exists() or (D["cfg"] / "config.yaml").write_text(
        yaml.safe_dump(DEFAULT_CFG, sort_keys=False), encoding="utf-8")
    pg = D["cfg"] / "ferroptosis_genes.csv"
    if not pg.exists():
        src = ("PROVISIONAL candidate list compiled by the code author from general ferroptosis literature; "
               "NOT the original project panel; NOT verified against FerrDb; generated 2026-10-10")
        pd.DataFrame([dict(gene_symbol=g, gene_role=r, source=src,
                           evidence_or_notes="canonical driver/suppressor in ferroptosis reviews; verify", confirmed=False)
                      for g, r in PANEL_ROWS]).to_csv(pg, index=False)
    (D["cfg"] / "drug_panel.csv").exists() or (D["cfg"] / "drug_panel.csv").write_text(DRUG_PANEL, encoding="utf-8")
    (ROOT / "requirements.txt").exists() or (ROOT / "requirements.txt").write_text(REQUIREMENTS, encoding="utf-8")
    (ROOT / "README.md").exists() or (ROOT / "README.md").write_text(README, encoding="utf-8")


def load_cfg():
    ensure_project()
    user = yaml.safe_load((D["cfg"] / "config.yaml").read_text(encoding="utf-8")) or {}
    if "penalizer_grid" in user and list(user["penalizer_grid"]) != list(DEFAULT_CFG["penalizer_grid"]):
        LOG.warning("config.yaml overrides penalizer_grid with %s (code default %s); edit/delete that key if unintended",
                    user["penalizer_grid"], DEFAULT_CFG["penalizer_grid"])
    return deep_update(copy.deepcopy(DEFAULT_CFG), user)


ALIASES = {"FSP1": "AIFM2", "NRF2": "NFE2L2", "P53": "TP53", "FPN1": "SLC40A1", "CARS": "CARS1"}


def load_panel(cfg, args=None):
    df = pd.read_csv(D["cfg"] / "ferroptosis_genes.csv")
    df["gene_symbol"] = df["gene_symbol"].astype(str).str.strip().str.upper().replace(ALIASES)
    df["gene_role"] = df["gene_role"].astype(str).str.strip().str.lower()
    df = df.drop_duplicates("gene_symbol")
    conf = bool(cfg.get("panel_confirmed")) or ("confirmed" in df and df["confirmed"].astype(str).str.lower().isin(
        ["true", "1", "yes"]).all())
    status = "CONFIRMED" if conf else "PROVISIONAL (not confirmed by the user)"
    ver = hashlib.sha256(",".join(sorted(df.gene_symbol + ":" + df.gene_role)).encode()).hexdigest()[:12]
    return df, status, ver


def panel_gate(cfg, args):
    panel, status, ver = load_panel(cfg)
    if status.startswith("PROVISIONAL") and not (args and args.accept_provisional_panel):
        raise DataUnavailable("Gene panel is PROVISIONAL. Confirm/replace config/ferroptosis_genes.csv (set confirmed=True "
                              "or panel_confirmed: true in config.yaml) or re-run with --accept-provisional-panel.")
    if status.startswith("PROVISIONAL"):
        deviation("Gene panel is PROVISIONAL (not the original project panel); all ferroptosis findings are provisional.")
    return panel, status, ver


def import_panel(path):
    """Replace config/ferroptosis_genes.csv with the user's original panel and mark it CONFIRMED.
    Stale panel-dependent results (frozen model, random-control cache) are removed."""
    src = pd.read_csv(path)
    src.columns = [c.strip().lower() for c in src.columns]
    for need in ("gene_symbol", "gene_role"):
        if need not in src.columns:
            sys.exit(f"--panel file must contain column '{need}' (found {list(src.columns)})")
    src["gene_symbol"] = src["gene_symbol"].astype(str).str.strip().str.upper().replace(ALIASES)
    src["gene_role"] = src["gene_role"].astype(str).str.strip().str.lower()
    bad = src[~src.gene_role.isin(["driver", "suppressor", "other"])]
    if len(bad):
        sys.exit(f"invalid gene_role values (use driver/suppressor/other): {bad.gene_role.unique().tolist()}")
    dups = src.gene_symbol[src.gene_symbol.duplicated()].tolist()
    if dups:
        LOG.warning("duplicate genes dropped: %s", dups)
    src = src.drop_duplicates("gene_symbol")
    target = D["cfg"] / "ferroptosis_genes.csv"
    if target.exists():
        shutil.copy(target, target.with_name(f"ferroptosis_genes.bak_{dt.datetime.now():%Y%m%d_%H%M%S}.csv"))
    out = pd.DataFrame(dict(
        gene_symbol=src.gene_symbol, gene_role=src.gene_role,
        source=src["source"] if "source" in src else "user-supplied original panel",
        evidence_or_notes=src["evidence_or_notes"] if "evidence_or_notes" in src else "",
        confirmed=True))
    out.to_csv(target, index=False)
    for p in (D["csv"] / "random_gene_set_benchmark.csv", D["proc"] / "frozen_prognosis_model.json"):
        if p.exists():
            p.unlink()
    LOG.info("Imported %d genes (%d driver, %d suppressor, %d other) -> %s; panel marked CONFIRMED",
             len(out), (out.gene_role == "driver").sum(), (out.gene_role == "suppressor").sum(),
             (out.gene_role == "other").sum(), target)
    if not (out.gene_role == "driver").any() or not (out.gene_role == "suppressor").any():
        LOG.warning("panel lacks driver and/or suppressor roles: FerroScore, classification FerroScore and drug stage will fail")


# ----------------------------------------------------------------------------------------
# Statistics utilities
# ----------------------------------------------------------------------------------------
def bh_adjust(p):
    p = np.asarray(p, float)
    q = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    n = len(pv)
    if n == 0:
        return q
    order = np.argsort(pv)
    r = pv[order] * n / np.arange(1, n + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(r, 0, 1)
    q[ok] = out
    return q


def spearman_ci(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = len(x)
    if n < 4 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return np.nan, np.nan, np.nan, np.nan
    r, p = stats.spearmanr(x, y)
    if abs(r) >= 1:
        return r, p, np.nan, np.nan
    z, se = np.arctanh(r), 1.06 / math.sqrt(n - 3)
    return r, p, float(np.tanh(z - 1.96 * se)), float(np.tanh(z + 1.96 * se))


def cindex(t, e, risk):
    return float(concordance_index(np.asarray(t, float), -np.asarray(risk, float), np.asarray(e, int)))


def boot_ci(fn, n, n_boot, seed, *arrs):
    rng = np.random.default_rng(seed)
    vals = []
    arrs = [np.asarray(a) for a in arrs]
    for _ in range(n_boot):
        i = rng.integers(0, n, n)
        try:
            v = fn(*[a[i] for a in arrs])
            if np.isfinite(v):
                vals.append(v)
        except Exception:
            pass
    return (float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))) if len(vals) > 20 else (np.nan, np.nan)


def zfit(X):
    mu, sd = X.mean(), X.std(ddof=0)
    return mu, sd.where(sd > 1e-8)


def empirical_p(observed, controls):
    c = np.asarray([v for v in controls if np.isfinite(v)])
    return (1 + np.sum(c >= observed)) / (1 + len(c)) if len(c) else np.nan


# ----------------------------------------------------------------------------------------
# Download infrastructure + provenance manifest
# ----------------------------------------------------------------------------------------
MANIFEST = []


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for ch in iter(lambda: f.read(1 << 20), b""):
            h.update(ch)
    return h.hexdigest()


def md5_file(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for ch in iter(lambda: f.read(1 << 20), b""):
            h.update(ch)
    return h.hexdigest()


def record(dataset, accession, url, path, status, kind="raw", version="", http_status="", md5_ok="", error=""):
    p = Path(path) if path else None
    ok = p is not None and p.exists()
    MANIFEST.append(dict(dataset=dataset, accession=accession, source_url=url,
                         download_timestamp=dt.datetime.now().isoformat(timespec="seconds"),
                         file_name=p.name if p else "", file_size_bytes=p.stat().st_size if ok else "",
                         sha256=sha256_file(p) if ok else "", retrieval_status=status, data_version=version,
                         data_kind=kind, http_status=http_status, vendor_md5_verified=md5_ok, error=error))


def flush_manifest():
    D["prov"].mkdir(parents=True, exist_ok=True)
    p = D["prov"] / "download_manifest.csv"
    new = pd.DataFrame(MANIFEST)
    if p.exists() and len(new):
        old = pd.read_csv(p, dtype=str)
        new = pd.concat([old, new.astype(str)]).drop_duplicates(["dataset", "accession", "file_name"], keep="last")
    if len(new):
        new.to_csv(p, index=False)
        new.to_csv(D["csv"] / "download_manifest.csv", index=False)
    MANIFEST.clear()


def valid_file(path, magic):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return False
    with open(path, "rb") as f:
        head = f.read(512)
    if magic == "gzip":
        return head[:2] == b"\x1f\x8b"
    if magic == "zip":
        return head[:2] == b"PK"
    low = head.lower()
    return not (b"<!doctype html" in low or b"<html" in low)


def http(method, url, retries=4, timeout=(20, 120), **kw):
    last = None
    for a in range(retries):
        try:
            r = SESSION.request(method, url, timeout=timeout, **kw)
            if r.status_code >= 500 or r.status_code == 429:
                raise requests.HTTPError(f"HTTP {r.status_code}")
            r.raise_for_status()
            return r
        except Exception as ex:
            last = ex
            time.sleep(min(2 ** a, 20))
    raise DataUnavailable(f"request failed after {retries} tries: {method} {url}: {last}")


def download_file(urls, dest, dataset, accession, kind="raw", magic="text", md5=None, refresh=False, version=""):
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and not refresh and valid_file(dest, magic) and (not md5 or md5_file(dest) == md5):
        record(dataset, accession, "cache", dest, "cached", kind, version, md5_ok=bool(md5))
        return dest
    errs = []
    for url in urls:
        for attempt in range(1, 4):
            tmp = dest.with_suffix(dest.suffix + ".part")
            try:
                with SESSION.get(url, stream=True, timeout=(20, 180)) as r:
                    r.raise_for_status()
                    with open(tmp, "wb") as f:
                        for ch in r.iter_content(1 << 20):
                            f.write(ch)
                    code = r.status_code
                if not valid_file(tmp, magic):
                    raise ValueError("downloaded content failed format validation (HTML error page or wrong format)")
                if md5 and md5_file(tmp) != md5:
                    raise ValueError("md5 mismatch")
                tmp.replace(dest)
                record(dataset, accession, url, dest, "downloaded", kind, version, code, bool(md5))
                return dest
            except Exception as ex:
                errs.append(f"{url} [try {attempt}]: {ex}")
                tmp.exists() and tmp.unlink()
                time.sleep(2 ** attempt)
    record(dataset, accession, "; ".join(urls), None, "failed", kind, version, error=" | ".join(errs[-3:]))
    raise DataUnavailable(f"{dataset}/{accession} download failed. Tried: {urls}. Last errors: {errs[-2:]}. "
                          f"Manual fallback: download from the official site and place the file at {dest}")


def adopt_manual(folder, patterns, canonical):
    """If the user dropped a file with a different name into the folder, copy it to the canonical name."""
    canonical = Path(canonical)
    if canonical.exists():
        return canonical
    for pat in patterns:
        for p in sorted(Path(folder).glob(pat)):
            if p.resolve() != canonical.resolve() and not p.name.endswith(".part"):
                shutil.copy(p, canonical)
                LOG.info("using manually placed file %s as %s", p.name, canonical.name)
                return canonical
    return None


# ----------------------------------------------------------------------------------------
# DATA ACQUISITION: TCGA-COAD via GDC
# ----------------------------------------------------------------------------------------
def gdc_query(endpoint, filters, fields, size=1000):
    hits, frm, total = [], 0, 0
    while True:
        body = {"filters": filters, "fields": ",".join(fields), "size": size, "from": frm, "format": "JSON"}
        d = http("POST", f"{GDC}/{endpoint}", json=body).json()["data"]
        hits += d["hits"]
        total = d["pagination"]["total"]
        frm += size
        if frm >= total or not d["hits"]:
            break
    return hits, total


def _gdc_in(field, vals):
    return {"op": "in", "content": {"field": field, "value": vals}}


def tcga_file_inventory(cfg):
    filt = {"op": "and", "content": [
        _gdc_in("cases.project.project_id", [cfg["tcga_project"]]),
        _gdc_in("data_category", ["Transcriptome Profiling"]),
        _gdc_in("data_type", ["Gene Expression Quantification"]),
        _gdc_in("experimental_strategy", ["RNA-Seq"]),
        _gdc_in("analysis.workflow_type", [cfg["workflow_type"]]),
        _gdc_in("access", ["open"])]}
    fields = ["file_id", "file_name", "md5sum", "file_size", "access", "cases.case_id", "cases.submitter_id",
              "cases.samples.submitter_id", "cases.samples.sample_type", "associated_entities.entity_submitter_id"]
    hits, total = gdc_query("files", filt, fields)
    rows = []
    for h in hits:
        cases = h.get("cases") or [{}]
        case = cases[0]
        ents = h.get("associated_entities") or [{}]
        aliq = ents[0].get("entity_submitter_id", "") or ""
        samp = aliq[:16]
        st = None
        for s in case.get("samples") or []:
            if s.get("submitter_id") == samp:
                st = s.get("sample_type")
        rows.append(dict(file_id=h["file_id"], file_name=h.get("file_name"), md5=h.get("md5sum"), size=h.get("file_size"),
                         case_id=case.get("case_id"), patient_id=case.get("submitter_id") or aliq[:12], sample_id=samp,
                         aliquot_id=aliq, sample_type_code=aliq[13:15], sample_type_api=st, n_cases_in_file=len(cases)))
    return pd.DataFrame(rows), total


def select_tcga_files(inv):
    """Deterministic rule: primary tumour (code 01 / 'Primary Tumor'); per patient keep the lexicographically first
    sample barcode, then first file name."""
    audit = []
    prim_mask = (inv.sample_type_code == "01") & (inv.sample_type_api.fillna("Primary Tumor") == "Primary Tumor")
    for _, r in inv[~prim_mask].iterrows():
        audit.append(dict(cohort=COHORT_TCGA, id=r.file_id, patient_id=r.patient_id, reason="not_primary_tumour_sample"))
    prim = inv[prim_mask].sort_values(["patient_id", "sample_id", "file_name"])
    sel = prim.drop_duplicates("patient_id", keep="first")
    for _, r in prim[~prim.file_id.isin(sel.file_id)].iterrows():
        audit.append(dict(cohort=COHORT_TCGA, id=r.file_id, patient_id=r.patient_id,
                          reason="additional_primary_sample_for_same_patient (deterministic rule)"))
    return sel, pd.DataFrame(audit)


TCGA_CASE_FIELDS = ["submitter_id", "case_id", "demographic.vital_status", "demographic.days_to_death",
                    "demographic.age_at_index", "demographic.gender", "diagnoses.age_at_diagnosis",
                    "diagnoses.days_to_last_follow_up", "diagnoses.ajcc_pathologic_stage", "diagnoses.primary_diagnosis",
                    "follow_ups.days_to_follow_up"]


def fetch_tcga(cfg, refresh, dry):
    st = http("GET", f"{GDC}/status").json()
    release = f"GDC {st.get('data_release', '?')} (api {st.get('version', '?')})"
    LOG.info("GDC status: %s", release)
    inv, total = tcga_file_inventory(cfg)
    if inv.empty:
        raise DataUnavailable("GDC returned no open-access STAR-Counts files for the requested filters")
    LOG.info("GDC: %d files, %d cases, %d samples, %d aliquots", len(inv), inv.case_id.nunique(),
             inv.sample_id.nunique(), inv.aliquot_id.nunique())
    sel, audit = select_tcga_files(inv)
    inv.to_csv(D["proc"] / "tcga_file_inventory.csv", index=False)
    sel.to_csv(D["proc"] / "tcga_selected_files.csv", index=False)
    audit.to_csv(D["proc"] / "tcga_audit.csv", index=False)
    LOG.info("TCGA primary-tumour files selected: %d (one per patient)", len(sel))
    if dry:
        LOG.info("[dry-run] skipping %d expression downloads and clinical download", len(sel))
        return
    filt = {"op": "in", "content": {"field": "project.project_id", "value": [cfg["tcga_project"]]}}
    hits, _ = gdc_query("cases", filt, TCGA_CASE_FIELDS)
    cp = D["raw"] / "tcga" / "clinical_cases.json"
    cp.write_text(json.dumps(hits), encoding="utf-8")
    record("TCGA-COAD", "clinical(cases endpoint)", f"{GDC}/cases", cp, "downloaded", "raw", release)

    def one(r):
        return download_file([f"{GDC}/data/{r.file_id}"], D["raw"] / "tcga" / f"{r.file_id}.tsv", "TCGA-COAD",
                             r.file_id, "raw", "text", md5=r.md5, refresh=refresh, version=release)
    fails = 0
    with ThreadPoolExecutor(4) as ex:
        futs = {ex.submit(one, r): r.file_id for r in sel.itertuples()}
        for i, f in enumerate(as_completed(futs), 1):
            try:
                f.result()
            except DataUnavailable as e:
                fails += 1
                LOG.warning("file %s failed: %s", futs[f], str(e)[:200])
            if i % 50 == 0:
                LOG.info("TCGA downloads: %d/%d", i, len(sel))
    LOG.info("TCGA expression download complete: %d failed of %d", fails, len(sel))


# ----------------------------------------------------------------------------------------
# DATA ACQUISITION: GEO (series matrix + platform annotation), GDSC, CTRP
# ----------------------------------------------------------------------------------------
def geo_stub(acc):
    return acc[:-3] + "nnn"


def fetch_geo(acc, cfg, refresh):
    raw = D["raw"] / "geo" / acc
    sm = download_file([f"https://ftp.ncbi.nlm.nih.gov/geo/series/{geo_stub(acc)}/{acc}/matrix/{acc}_series_matrix.txt.gz"],
                       raw / f"{acc}_series_matrix.txt.gz", "GEO", acc, "raw", "gzip", refresh=refresh)
    parsed = parse_series_matrix(sm, header_only=True)
    gpl = parsed["platform"]
    if not gpl:
        raise DataUnavailable(f"{acc}: platform id not found in series matrix")
    try:
        download_file([f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{geo_stub(gpl)}/{gpl}/annot/{gpl}.annot.gz"],
                      raw / f"{gpl}.annot.gz", "GEO-platform", gpl, "processed", "gzip", refresh=refresh)
    except DataUnavailable:
        LOG.warning("%s annot.gz unavailable; trying GEO text-export fallback", gpl)
        download_file([f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gpl}&targ=self&form=text&view=full"],
                      raw / f"{gpl}.soft.txt", "GEO-platform", gpl, "processed", "text", refresh=refresh)


def _open_text(path):
    path = Path(path)
    with open(path, "rb") as f:
        gz = f.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if gz else open(path, "rt", encoding="utf-8", errors="replace")


def parse_series_matrix(path, header_only=False):
    meta = {}
    chars = []
    with _open_text(path) as f:
        expr = None
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                if not header_only:
                    expr = pd.read_csv(f, sep="\t", index_col=0, low_memory=False)
                    expr = expr[~expr.index.astype(str).str.startswith("!")]
                break
            if not line.startswith("!"):
                continue
            parts = line.rstrip("\n").split("\t")
            key, vals = parts[0], [v.strip().strip('"') for v in parts[1:]]
            if key == "!Sample_characteristics_ch1":
                chars.append(vals)
            else:
                meta.setdefault(key, vals)
    gsm = meta.get("!Sample_geo_accession", [])
    ch = {}
    for vals in chars:
        for i, v in enumerate(vals):
            if ":" in v and i < len(gsm):
                k, val = v.split(":", 1)
                if re.search(r";\s*[A-Za-z_][\w .()-]*:", val):
                    # packed characteristics, e.g. "Right; DukesStage: A; DFS_Time: 3.64; DFS_Cens: 1"
                    parts = val.split(";")
                    ch.setdefault(k.strip().lower(), {}).setdefault(gsm[i], parts[0].strip())
                    for part in parts[1:]:
                        if ":" in part:
                            k2, v2 = part.split(":", 1)
                            ch.setdefault(k2.strip().lower(), {}).setdefault(gsm[i], v2.strip())
                else:
                    ch.setdefault(k.strip().lower(), {}).setdefault(gsm[i], val.strip())
    chdf = pd.DataFrame(ch, index=gsm) if gsm else pd.DataFrame()
    plat = (meta.get("!Series_platform_id") or [""])[0]
    return dict(samples=gsm, titles=dict(zip(gsm, meta.get("!Sample_title", []))), chars=chdf, platform=plat, expr=expr,
                n_platforms=len(set(meta.get("!Sample_platform_id", [plat]))))


def inspect_geo(acc):
    """Print every characteristics key (lower-cased, as used in config geo.overrides) with example values."""
    sm = D["raw"] / "geo" / acc / f"{acc}_series_matrix.txt.gz"
    if not sm.exists():
        sys.exit(f"{sm} not found; run:  python ferro_pipeline.py --stage fetch")
    P = parse_series_matrix(sm, header_only=True)
    print(f"\n{acc}: platform {P['platform']}, {len(P['samples'])} samples")
    if P["chars"].empty:
        print("  (no 'key: value' characteristics found)")
    for k in P["chars"].columns:
        vals = P["chars"][k].dropna().unique()
        print(f"  key={k!r:45} n_unique={len(vals):4}  e.g. {list(vals[:5])}")
    print("\nUse these keys in config/config.yaml, e.g.:\n"
          "geo:\n  overrides:\n    " + acc + ":\n      dfs_time: \"<key>\"\n      dfs_event: \"<key>\"\n"
          "(fields: os_time os_event dfs_time dfs_event age stage mmr). Event keys containing 'cens' are inverted; "
          "verify coding on the GEO page.")


def parse_platform_table(path):
    with _open_text(path) as f:
        for line in f:
            if line.lower().startswith("!platform_table_begin"):
                break
        else:
            raise DataUnavailable("platform table marker not found in annotation file")
        t = pd.read_csv(f, sep="\t", dtype=str, low_memory=False)
    t = t[~t.iloc[:, 0].astype(str).str.startswith("!")]
    cols = {c.lower().strip(): c for c in t.columns}
    sym = next((cols[k] for k in ("gene symbol", "gene_symbol", "symbol", "gene.symbol") if k in cols), None)
    if sym is None:
        raise DataUnavailable(f"no gene-symbol column in platform table; columns={list(t.columns)[:15]}")
    return t.rename(columns={t.columns[0]: "ID", sym: "SYMBOL"})[["ID", "SYMBOL"]]


def probes_to_genes(expr_probe, ann):
    """Rule: drop probes without symbol or mapping to >1 gene ('///'); per gene keep the probe with highest
    mean expression (outcome-free). Returns genes x samples."""
    a = ann.copy()
    a["SYMBOL"] = a["SYMBOL"].astype(str).str.strip().str.upper()
    a = a[(a.SYMBOL != "NAN") & (a.SYMBOL != "") & (a.SYMBOL != "---") & (~a.SYMBOL.str.contains("///|//|;"))]
    a["ID"] = a["ID"].astype(str)
    e = expr_probe.copy()
    e.index = e.index.astype(str)
    e = e[~e.index.duplicated(keep="first")]
    a = a[a.ID.isin(e.index)].drop_duplicates("ID")
    mean = e.loc[a.ID].mean(axis=1)
    a["m"] = mean.values
    best = a.sort_values("m", ascending=False).drop_duplicates("SYMBOL")
    g = e.loc[best.ID]
    g.index = best.SYMBOL.values
    return g.replace([np.inf, -np.inf], np.nan), dict(n_probes_total=len(e), n_probes_mapped=len(a), n_genes=len(g))


GEO_RULES = {
    "os_time": (r"(^|[^a-z])(os|overall)[^a-z]*.*(delay|time|month|year|day|follow)", r"event|status|dfs|rfs|dss|disease|relapse|recurr|free|cens"),
    "os_event": (r"(^|[^a-z])(os|overall).*(event|status|death|dead|vital)", r"dfs|rfs|dss|disease|relapse|recurr|free"),
    "dfs_time": (r"(dfs|rfs|disease[^a-z]*free|relapse[^a-z]*free|recurrence[^a-z]*free).*(delay|time|month|year|day|follow)", r"event|status|cens"),
    "dfs_event": (r"(dfs|rfs|disease[^a-z]*free|relapse[^a-z]*free).*(event|status|relapse|recurr|cens)", r"time|delay|month"),
    "age": (r"age", r"stage|dosage|passage|image|page"),
    "stage": (r"stage", r"grade|substage"),
    "mmr": (r"mmr|msi|microsatellite", r"^$"),
}


def pick_key(keys, field):
    inc, exc = GEO_RULES[field]
    cand = [k for k in keys if re.search(inc, k) and not re.search(exc, k)]
    if field == "stage":
        cand.sort(key=lambda k: (0 if ("tnm" in k or "ajcc" in k) else 1))
    return cand[0] if cand else None


def parse_event(s):
    pos = {"1", "1.0", "dead", "death", "deceased", "yes", "true", "event", "relapse", "recurrence", "died"}
    neg = {"0", "0.0", "alive", "censored", "no", "false", "none", "disease-free", "living"}

    def f(v):
        v = str(v).strip().lower()
        if v in pos or v.startswith("dead") or v.startswith("deceas"):
            return 1.0
        if v in neg or v.startswith("alive") or v.startswith("no ") or v.startswith("non-") or v.startswith("not "):
            return 0.0
        return np.nan
    return pd.Series(s).map(f).astype(float)


def parse_stage(v):
    t = str(v).strip().upper().replace("STAGE", "").replace("TNM", "").strip(" :=")
    if t in ("A", "B", "C", "D"):  # Dukes A-D -> 1-4
        return float("ABCD".index(t) + 1)
    m = re.match(r"^(IV|III|II|I|[1-4])[ABC]?\d?$", t)
    if not m:
        return np.nan
    return float({"I": 1, "II": 2, "III": 3, "IV": 4}.get(m.group(1), m.group(1)))


def parse_mmr(v):
    v = str(v).strip().lower()
    if any(k in v for k in ("dmmr", "deficient", "unstable", "msi-h", "msi-high")) or v in ("msi", "msi-h"):
        return 1.0
    if any(k in v for k in ("pmmr", "proficient", "mss", "stable", "msi-l", "msi-low")):
        return 0.0
    return np.nan


def _num(s):
    return pd.to_numeric(pd.Series(s).astype(str).str.extract(r"(-?\d+\.?\d*)")[0], errors="coerce")


def _to_months(vals, key):
    k = key.lower()
    if "year" in k:
        return vals * 12, "years->months"
    if "day" in k:
        return vals / 30.4375, "days->months"
    return vals, "assumed months (unit not in field name)"


def map_geo_clinical(acc, chars, titles, cfg):
    over = (cfg["geo"].get("overrides") or {}).get(acc, {})
    keys = list(chars.columns)
    mapping, dd = {}, []
    out = pd.DataFrame(index=chars.index)
    out["cohort"], out["sample_id"], out["patient_id"] = acc, chars.index, chars.index
    out["title"] = [titles.get(i, "") for i in chars.index]
    for fld in ("os_time", "os_event", "dfs_time", "dfs_event", "age", "stage", "mmr"):
        k = over.get(fld) or pick_key(keys, fld)
        if k is not None and k not in chars.columns:
            LOG.warning("%s: override key %r for %s not found in characteristics (available: %s)", acc, k, fld, keys)
            k = None
        mapping[fld] = k
        if k is None:
            out[fld] = np.nan
            dd.append(dict(cohort=acc, field=fld, source_key="(not found)", n_nonmissing=0, pct_missing=100.0, unit_note=""))
            continue
        raw = chars[k]
        note = ""
        if fld.endswith("_time"):
            v, note = _to_months(_num(raw), k)
        elif fld.endswith("_event"):
            v = parse_event(raw)
            if "cens" in k.lower():
                v, note = 1 - v, "censoring flag inverted"
        elif fld == "age":
            v = _num(raw)
            note = "years assumed" if "month" not in k.lower() else "months!"
        elif fld == "stage":
            v = raw.map(parse_stage)
        else:
            v = raw.map(parse_mmr)
        out[fld] = np.asarray(v, float)
        dd.append(dict(cohort=acc, field=fld, source_key=k, n_nonmissing=int(out[fld].notna().sum()),
                       pct_missing=round(100 * out[fld].isna().mean(), 2), unit_note=note))
    return out, pd.DataFrame(dd)


def choose_endpoint(c, cfg):
    def ok(t, e):
        v = c[t].gt(0) & c[e].isin([0, 1])
        return int(v.sum()), int(c.loc[v, e].sum())
    n_os, ev_os = ok("os_time", "os_event")
    n_dfs, ev_dfs = ok("dfs_time", "dfs_event")
    if ev_os >= cfg["min_events"]:
        return "OS", "os_time", "os_event", ""
    if ev_dfs >= cfg["min_events"]:
        return "DFS", "dfs_time", "dfs_event", "OS unavailable -> DFS/RFS used (endpoint differs from TCGA OS training; exploratory)"
    return "none", None, None, "no survival endpoint with enough events"


def process_geo(acc, cfg):
    raw = D["raw"] / "geo" / acc
    sm = raw / f"{acc}_series_matrix.txt.gz"
    if not sm.exists():
        raise DataUnavailable(f"{acc}: series matrix not downloaded")
    P = parse_series_matrix(sm)
    if P["expr"] is None or P["expr"].empty:
        raise DataUnavailable(f"{acc}: expression table missing/unparseable")
    gpl = P["platform"]
    if P["n_platforms"] > 1:
        deviation(f"{acc}: multiple platforms in one series; only first platform {gpl} handled -> verify")
    ex = P["expr"].apply(pd.to_numeric, errors="coerce")
    n_before = ex.shape[1]
    dup = ex.columns.duplicated()
    ex = ex.loc[:, ~dup]
    ex = ex.loc[:, ex.notna().any()]
    mn, mx, med = float(np.nanmin(ex.values)), float(np.nanmax(ex.values)), float(np.nanmedian(ex.values))
    if mx > 100:
        ex = np.log2(ex.clip(lower=0) + 1)
        scale = "raw intensities -> log2(x+1) applied"
    else:
        scale = "already log-scale (max<=100); no transform"
    LOG.info("%s: expression range [%.2f, %.2f] median %.2f -> %s", acc, mn, mx, med, scale)
    ann_path = next((p for p in [raw / f"{gpl}.annot.gz", raw / f"{gpl}.soft.txt"] if p.exists()), None)
    if ann_path is None:
        raise DataUnavailable(f"{acc}: platform annotation for {gpl} not downloaded")
    genes, pm = probes_to_genes(ex, parse_platform_table(ann_path))
    Xg = genes.T
    Xg = Xg.loc[:, Xg.notna().mean() > 0.9]
    clin, dd = map_geo_clinical(acc, P["chars"].reindex(Xg.index), P["titles"], cfg)
    ep, tcol, ecol, note = choose_endpoint(clin, cfg)
    LOG.info("%s clinical mapping: %s", acc, "; ".join(f"{r.field}->{r.source_key}" for r in dd.itertuples()))
    if ep == "none":
        LOG.warning("%s: no survival endpoint found. Characteristics keys available: %s", acc, list(P["chars"].columns))
    clin["endpoint"] = ep
    clin["endpoint_note"] = note
    clin["time"] = clin[tcol] if tcol else np.nan
    clin["event"] = clin[ecol] if ecol else np.nan
    valid = clin.time.gt(0) & clin.event.isin([0, 1])
    audit = [dict(cohort=acc, id=i, patient_id=i, reason="invalid/missing survival (time<=0 or event unknown)")
             for i in clin.index[~valid]]
    if dup.any():
        audit += [dict(cohort=acc, id=c, patient_id=c, reason="duplicate sample column") for c in P["expr"].columns[dup]]
    incl(acc, "samples in series matrix", n_before, n_before, "")
    incl(acc, "unique samples with expression", n_before, Xg.shape[0], "duplicate/all-NaN columns removed")
    incl(acc, f"valid survival ({ep})", Xg.shape[0], int(valid.sum()), "time<=0 or event unknown")
    Xg.astype("float32").round(4).to_csv(D["proc"] / f"{acc}_expr.csv.gz")
    clin.to_csv(D["proc"] / f"{acc}_clin.csv", index=False)
    pd.DataFrame(audit).to_csv(D["proc"] / f"{acc}_audit.csv", index=False)
    dd.to_csv(D["proc"] / f"{acc}_dd.csv", index=False)
    if valid.sum() >= 20:
        rate = float(clin.loc[valid, "event"].mean())
        inv_flag = any(r.unit_note == "censoring flag inverted" for r in dd.itertuples())
        if inv_flag or rate < 0.15 or rate > 0.80:
            LOG.warning("%s: %s event rate = %.1f%% (%d/%d)%s -> verify coding with:  python ferro_pipeline.py --inspect-geo %s",
                        acc, ep, 100 * rate, int(clin.loc[valid, "event"].sum()), int(valid.sum()),
                        "; a 'cens' flag was inverted" if inv_flag else "", acc)
    inventory(acc, data_type="microarray", platform=gpl, n_samples=Xg.shape[0], n_patients=clin.patient_id.nunique(),
              n_genes=Xg.shape[1], n_probes=pm["n_probes_total"], pct_missing_expression=round(100 * float(Xg.isna().mean().mean()), 3),
              value_range=f"{mn:.2f}..{mx:.2f}", scale_note=scale, endpoint=ep, events=int(clin.loc[valid, "event"].sum()),
              note=note)
    if note:
        deviation(f"{acc}: {note}")
    LOG.info("%s processed: %d samples, %d genes, endpoint=%s", acc, *Xg.shape, ep)


def fetch_gdsc(cfg, refresh):
    raw = D["raw"] / "gdsc"
    specs = [("GDSC2_fitted_dose_response.xlsx", ["*GDSC2*dose*.xlsx", "*fitted*dose*.xlsx"], "gdsc2_response", "GDSC2",
              "fitted_dose_response", "raw", "release 8.5 (config)"),
             ("Cell_line_RMA_proc_basalExp.txt.zip", ["*RMA*basalExp*.zip"], "gdsc_expression", "GDSC",
              "RMA_basal_expression", "processed", "")]
    errs = []
    for fname, pats, ukey, ds, acc, kind, ver in specs:
        dest = raw / fname
        adopt_manual(raw, pats, dest)
        try:
            download_file(urls_for(cfg, ukey), dest, ds, acc, kind, "zip", refresh=refresh, version=ver)
        except DataUnavailable as e:
            errs.append(str(e)[:250])
    if errs:
        raise DataUnavailable(" || ".join(errs))


def fetch_ctrp(cfg, refresh):
    raw = D["raw"] / "ctrp"
    dest = raw / "CTRPv2_ExpandedDataset.zip"
    adopt_manual(raw, ["*CTRP*Expanded*.zip", "*CTRPv2*.zip"], dest)
    return download_file(urls_for(cfg, "ctrp"), dest, "CTRPv2", "ExpandedDataset", "processed", "zip", refresh=refresh)


# ----------------------------------------------------------------------------------------
# PREPROCESSING
# ----------------------------------------------------------------------------------------
def process_tcga(cfg):
    sel_p = D["proc"] / "tcga_selected_files.csv"
    if not sel_p.exists():
        raise DataUnavailable("TCGA file inventory missing; run the fetch stage")
    inv = pd.read_csv(D["proc"] / "tcga_file_inventory.csv")
    sel = pd.read_csv(sel_p)
    audit = pd.read_csv(D["proc"] / "tcga_audit.csv")
    incl(COHORT_TCGA, "STAR-Counts open-access files returned", len(inv), len(inv), f"{inv.case_id.nunique()} cases, {inv.sample_id.nunique()} samples, {inv.aliquot_id.nunique()} aliquots")
    incl(COHORT_TCGA, "primary tumour + one sample per patient", len(inv), len(sel), "non-primary samples and extra primary samples removed")
    cols, miss = {}, []
    for r in sel.itertuples():
        p = D["raw"] / "tcga" / f"{r.file_id}.tsv"
        if not p.exists():
            miss.append(r.sample_id)
            continue
        t = pd.read_csv(p, sep="\t", comment="#", usecols=["gene_name", "gene_type", "tpm_unstranded"])
        t = t[t.gene_type == "protein_coding"].groupby("gene_name")["tpm_unstranded"].mean()
        cols[r.sample_id] = t
    if not cols:
        raise DataUnavailable("no TCGA expression files present locally")
    X = pd.DataFrame(cols).T
    X = np.log2(X + 1)
    X.index.name = "sample_id"
    incl(COHORT_TCGA, "expression file parsed", len(sel), X.shape[0], f"{len(miss)} files missing locally")
    cl_p = D["raw"] / "tcga" / "clinical_cases.json"
    if not cl_p.exists():
        raise DataUnavailable("TCGA clinical JSON missing")
    rows = []
    for h in json.loads(cl_p.read_text(encoding="utf-8")):
        dem = h.get("demographic") or {}
        dx = (h.get("diagnoses") or [{}])[0]
        vs, dtd, dlf = dem.get("vital_status"), dem.get("days_to_death"), dx.get("days_to_last_follow_up")
        fu = [f.get("days_to_follow_up") for f in (h.get("follow_ups") or []) if f.get("days_to_follow_up") is not None]
        days, ev = np.nan, np.nan
        if vs == "Dead":
            days, ev = dtd, 1
        elif vs == "Alive":
            days, ev = (dlf if dlf is not None else (max(fu) if fu else None)), 0
        age = dx.get("age_at_diagnosis")
        rows.append(dict(patient_id=h.get("submitter_id"), time=(float(days) / 30.4375) if days is not None else np.nan,
                         event=ev, age=(age / 365.25) if age else dem.get("age_at_index"),
                         stage=parse_stage(dx.get("ajcc_pathologic_stage") or ""), gender=dem.get("gender"),
                         stage_raw=dx.get("ajcc_pathologic_stage")))
    cl = pd.DataFrame(rows).drop_duplicates("patient_id")
    s2p = sel.set_index("sample_id").patient_id
    clin = pd.DataFrame({"sample_id": X.index, "patient_id": s2p.reindex(X.index).values}).merge(cl, on="patient_id", how="left")
    clin.insert(0, "cohort", COHORT_TCGA)
    clin["endpoint"], clin["mmr"] = "OS", np.nan
    clin["os_time"], clin["os_event"], clin["dfs_time"], clin["dfs_event"] = clin.time, clin.event, np.nan, np.nan
    n_have = int(clin.time.notna().sum())
    valid = clin.time.gt(0) & clin.event.isin([0, 1])
    incl(COHORT_TCGA, "clinical record with survival fields", X.shape[0], n_have, "vital status/time missing")
    incl(COHORT_TCGA, "valid OS (time>0, event known)", n_have, int(valid.sum()), "non-positive or missing follow-up")
    a2 = [dict(cohort=COHORT_TCGA, id=i, patient_id=i, reason="invalid/missing OS") for i in clin.loc[~valid, "sample_id"]]
    audit = audit[audit.reason != "invalid/missing OS"] if "reason" in audit.columns else audit
    pd.concat([audit, pd.DataFrame(a2)]).to_csv(D["proc"] / "tcga_audit.csv", index=False)
    X.astype("float32").round(4).to_csv(D["proc"] / f"{COHORT_TCGA}_expr.csv.gz")
    clin.to_csv(D["proc"] / f"{COHORT_TCGA}_clin.csv", index=False)
    pd.DataFrame([dict(cohort=COHORT_TCGA, field="mmr/MSI", source_key="(no open-access GDC field retrieved)", n_nonmissing=0,
                       pct_missing=100.0, unit_note="TCGA MSI not available through this pipeline"),
                  dict(cohort=COHORT_TCGA, field="dfs", source_key="(not defined)", n_nonmissing=0, pct_missing=100.0,
                       unit_note="TCGA DFS/PFI requires TCGA-CDR; not fetched"),
                  dict(cohort=COHORT_TCGA, field="os_time", source_key="demographic.days_to_death / diagnoses.days_to_last_follow_up / follow_ups.days_to_follow_up",
                       n_nonmissing=n_have, pct_missing=round(100 * (1 - n_have / max(len(clin), 1)), 2), unit_note="days -> months /30.4375"),
                  dict(cohort=COHORT_TCGA, field="stage", source_key="diagnoses.ajcc_pathologic_stage", n_nonmissing=int(clin.stage.notna().sum()),
                       pct_missing=round(100 * clin.stage.isna().mean(), 2), unit_note="I-IV -> 1-4 (sub-stages merged)"),
                  dict(cohort=COHORT_TCGA, field="age", source_key="diagnoses.age_at_diagnosis (days)", n_nonmissing=int(clin.age.notna().sum()),
                       pct_missing=round(100 * clin.age.isna().mean(), 2), unit_note="years")]).to_csv(D["proc"] / f"{COHORT_TCGA}_dd.csv", index=False)
    inventory(COHORT_TCGA, data_type="RNA-seq STAR counts -> log2(TPM+1)", platform="Illumina RNA-seq (GDC harmonized)",
              n_files=len(inv), n_samples=X.shape[0], n_patients=clin.patient_id.nunique(), n_genes=X.shape[1],
              pct_missing_expression=round(100 * float(X.isna().mean().mean()), 3), endpoint="OS", events=int(clin.loc[valid, "event"].sum()),
              note="MSI not available; DFS not fetched")
    deviation("TCGA MSI/MMR labels and DFS not retrieved (no standard open-access GDC field); MSI classification uses GEO cohorts with MMR labels.")
    LOG.info("TCGA processed: %d samples, %d genes, %d valid OS", X.shape[0], X.shape[1], int(valid.sum()))


def _read_zip_member(zpath, pattern):
    with zipfile.ZipFile(zpath) as z:
        names = [n for n in z.namelist() if re.search(pattern, n)]
        if not names:
            raise DataUnavailable(f"{Path(zpath).name}: no member matching {pattern}; contents={z.namelist()[:10]}")
        return z.read(names[0]), names[0]


def process_gdsc(cfg):
    rp = D["raw"] / "gdsc" / "GDSC2_fitted_dose_response.xlsx"
    ep = D["raw"] / "gdsc" / "Cell_line_RMA_proc_basalExp.txt.zip"
    if not (rp.exists() and ep.exists()):
        raise DataUnavailable("GDSC files not present (see README manual fallback)")
    resp = pd.read_excel(rp)
    need = {"COSMIC_ID", "CELL_LINE_NAME", "DRUG_ID", "DRUG_NAME", "LN_IC50", "AUC", cfg["gdsc_lineage_field"]}
    if not need.issubset(resp.columns):
        raise DataUnavailable(f"GDSC2 response schema unexpected; missing {need - set(resp.columns)}; columns={list(resp.columns)}")
    resp["COSMIC_ID"] = resp["COSMIC_ID"].astype("Int64").astype(str)
    n_all = resp.COSMIC_ID.nunique()
    crc = resp[resp[cfg["gdsc_lineage_field"]].astype(str) == cfg["gdsc_lineage_value"]].copy()
    dups = int(crc.duplicated(["COSMIC_ID", "DRUG_ID"]).sum())
    crc = crc.groupby(["COSMIC_ID", "CELL_LINE_NAME", "DRUG_ID", "DRUG_NAME"], as_index=False).agg(
        LN_IC50=("LN_IC50", "mean"), AUC=("AUC", "mean"))
    b, nm = _read_zip_member(ep, r"\.txt$")
    ex = pd.read_csv(io.BytesIO(b), sep="\t", low_memory=False)
    gcol = "GENE_SYMBOLS" if "GENE_SYMBOLS" in ex.columns else ex.columns[0]
    cell_cols = [c for c in ex.columns if str(c).startswith("DATA.")]
    if not cell_cols:
        raise DataUnavailable(f"GDSC expression: no DATA.<COSMIC> columns in {nm}")
    ex = ex[[gcol] + cell_cols].dropna(subset=[gcol])
    ex[gcol] = ex[gcol].astype(str).str.upper()
    ex = ex.groupby(gcol).mean()
    ex.columns = [c.replace("DATA.", "") for c in ex.columns]
    have = sorted(set(crc.COSMIC_ID) & set(ex.columns))
    incl("GDSC2-CRC", "CRC cell lines with drug response", n_all, crc.COSMIC_ID.nunique(), f"lineage rule {cfg['gdsc_lineage_field']}=={cfg['gdsc_lineage_value']}")
    incl("GDSC2-CRC", "CRC cell lines with expression (verified COSMIC ID join)", crc.COSMIC_ID.nunique(), len(have), "no expression for these lines")
    crc = crc[crc.COSMIC_ID.isin(have)]
    E = ex[have].T
    E.index.name = "COSMIC_ID"
    E.to_csv(D["proc"] / "gdsc2_expr_crc.csv.gz")
    crc.to_csv(D["proc"] / "gdsc2_response_crc.csv", index=False)
    inventory("GDSC2-CRC", data_type="RMA expression + LN_IC50/AUC", n_samples=len(have), n_patients=len(have),
              n_genes=E.shape[1], n_response_rows=len(crc), n_drugs=crc.DRUG_ID.nunique(),
              pct_missing_expression=round(100 * float(E.isna().mean().mean()), 3),
              note=f"{dups} duplicate (cell line, drug) pairs averaged; response units: LN_IC50 natural-log uM (source-defined), AUC")
    LOG.info("GDSC2 CRC: %d lines, %d drugs", len(have), crc.DRUG_ID.nunique())


def process_ctrp(cfg):
    zp = D["raw"] / "ctrp" / "CTRPv2_ExpandedDataset.zip"
    if not zp.exists():
        raise DataUnavailable("CTRP zip not present (see README manual fallback)")
    rd = lambda pat: pd.read_csv(io.BytesIO(_read_zip_member(zp, pat)[0]), sep="\t", low_memory=False)
    cur, exp, cl, cpd = rd(r"curves_post_qc"), rd(r"per_experiment"), rd(r"per_cell_line"), rd(r"per_compound")
    for df, need, nm in [(cur, {"experiment_id", "area_under_curve", "master_cpd_id"}, "curves"),
                         (exp, {"experiment_id", "master_ccl_id"}, "experiment"), (cl, {"master_ccl_id", "ccl_name"}, "cell_line"),
                         (cpd, {"master_cpd_id", "cpd_name"}, "compound")]:
        if not need.issubset(df.columns):
            raise DataUnavailable(f"CTRP {nm} schema unexpected: missing {need - set(df.columns)}")
    site = cl["ccle_primary_site"] if "ccle_primary_site" in cl.columns else pd.Series("", index=cl.index)
    cl = cl[site.astype(str) == cfg["ctrp_site"]]
    m = cur.merge(exp[["experiment_id", "master_ccl_id"]], on="experiment_id").merge(cl[["master_ccl_id", "ccl_name"]], on="master_ccl_id") \
        .merge(cpd[["master_cpd_id", "cpd_name"]], on="master_cpd_id")
    out = m.groupby(["ccl_name", "cpd_name"], as_index=False).agg(AUC=("area_under_curve", "mean"), n_experiments=("experiment_id", "nunique"))
    out["ccl_norm"], out["cpd_norm"] = out.ccl_name.map(norm), out.cpd_name.map(norm)
    out.to_csv(D["proc"] / "ctrp_response_crc.csv", index=False)
    inventory("CTRPv2-CRC", data_type="AUC (area under dose-response curve)", n_samples=out.ccl_name.nunique(),
              n_patients=out.ccl_name.nunique(), n_response_rows=len(out), n_drugs=out.cpd_name.nunique(),
              note="experiments averaged per (cell line, compound); site=" + cfg["ctrp_site"])
    LOG.info("CTRP CRC: %d lines, %d compounds", out.ccl_name.nunique(), out.cpd_name.nunique())


def load_cohort(name):
    ep, cp = D["proc"] / f"{name}_expr.csv.gz", D["proc"] / f"{name}_clin.csv"
    if not (ep.exists() and cp.exists()):
        raise DataUnavailable(f"processed data for {name} not found (run fetch+preprocess)")
    X = pd.read_csv(ep, index_col=0)
    X.index = X.index.astype(str)
    c = pd.read_csv(cp, dtype={"sample_id": str, "patient_id": str})
    c = c.set_index("sample_id", drop=False)
    c.index = c.index.astype(str)
    X = X.loc[X.index.intersection(c.index)]
    return X, c.loc[X.index]


def available_cohorts(cfg):
    out = {}
    for n in [COHORT_TCGA] + cfg["geo"]["cohorts"]:
        try:
            out[n] = load_cohort(n)
        except DataUnavailable as e:
            LOG.warning("%s", e)
    dp = D["proc"] / "overlap_drop.csv"
    if cfg.get("drop_overlapping_samples", True) and dp.exists() and dp.stat().st_size > 5:
        drop = pd.read_csv(dp, dtype=str)
        for n in list(out):
            ids = set(drop.loc[drop.cohort == n, "sample_id"])
            if n != COHORT_TCGA and ids:
                X, c = out[n]
                keep = ~X.index.isin(ids)
                LOG.info("%s: %d samples duplicate an earlier GEO cohort and are excluded from analysis (overlap_drop.csv)",
                         n, int((~keep).sum()))
                out[n] = (X.loc[keep], c.loc[keep])
    return out


def valid_surv(c):
    return c[c.time.gt(0) & c.event.isin([0, 1])]


def preprocess_finalize(cfg, args):
    coh = {}
    for n in [COHORT_TCGA] + cfg["geo"]["cohorts"]:
        try:
            coh[n] = load_cohort(n)
        except DataUnavailable:
            pass
    s = read_state()
    clin_all = pd.concat([c[["cohort", "sample_id", "patient_id", "time", "event", "endpoint", "os_time", "os_event", "dfs_time",
                              "dfs_event", "age", "stage", "mmr"]] for _, c in coh.values()]) if coh else pd.DataFrame()
    if len(clin_all):
        save_csv(clin_all, "clinical_data_clean.csv")
    save_csv(pd.DataFrame(s["incl"]), "cohort_inclusion_exclusion.csv")
    inv = pd.DataFrame(list(s["inventory"].values()))
    if len(inv):
        save_csv(inv, "dataset_inventory.csv")
    dd = [pd.read_csv(p) for p in D["proc"].glob("*_dd.csv")]
    if dd:
        save_csv(pd.concat(dd), "data_dictionary.csv")
    au = [pd.read_csv(p) for p in D["proc"].glob("*_audit.csv") if p.stat().st_size > 5]
    if au:
        save_csv(pd.concat(au), "sample_exclusion_audit.csv")
    panel, status, ver = load_panel(cfg)
    save_csv(panel.assign(panel_version=ver, panel_status=status), "ferroptosis_gene_panel.csv")
    rows = []
    allg = None
    for n, (X, _) in coh.items():
        pres = [g for g in panel.gene_symbol if g in X.columns]
        rows.append(dict(cohort=n, n_genes=X.shape[1], n_panel_genes=len(panel), n_panel_present=len(pres),
                         panel_missing=";".join(sorted(set(panel.gene_symbol) - set(pres)))))
        allg = set(X.columns) if allg is None else allg & set(X.columns)
    if allg is not None:
        rows.append(dict(cohort="INTERSECTION_ALL_COHORTS", n_genes=len(allg), n_panel_genes=len(panel),
                         n_panel_present=len(set(panel.gene_symbol) & allg), panel_missing=";".join(sorted(set(panel.gene_symbol) - allg))))
        save_csv(pd.DataFrame(rows), "gene_intersection_report.csv")
    if COHORT_TCGA in coh:
        X = coh[COHORT_TCGA][0]
        save_csv(X[[g for g in panel.gene_symbol if g in X.columns]], "tcga_expression_processed.csv", index=True)
    for n in cfg["geo"]["cohorts"]:
        if n in coh:
            X = coh[n][0]
            save_csv(X[[g for g in panel.gene_symbol if g in X.columns]], f"geo_{n}_expression_processed.csv", index=True)
    # overlap checks
    ov = []
    dup_rows = []
    for n, (X, c) in coh.items():
        if n != COHORT_TCGA:
            ov.append(dict(check=f"{n} sample titles containing 'TCGA'", result=int(c["title"].astype(str).str.contains("TCGA").sum()) if "title" in c else "n/a"))
    gl = [n for n in coh if n != COHORT_TCGA]
    for i in range(len(gl)):
        for j in range(i + 1, len(gl)):
            a, b = coh[gl[i]][0], coh[gl[j]][0]
            cg = list(set(a.columns) & set(b.columns))
            if len(cg) > 500:
                top = a[cg].var().nlargest(min(3000, len(cg))).index
                za = (a[top] - a[top].mean()) / a[top].std(ddof=0).replace(0, np.nan)
                zb = (b[top] - b[top].mean()) / b[top].std(ddof=0).replace(0, np.nan)
                za, zb = za.fillna(0).values, zb.fillna(0).values
                R = (za / np.linalg.norm(za, axis=1, keepdims=True)) @ (zb / np.linalg.norm(zb, axis=1, keepdims=True)).T
                ov.append(dict(check=f"possible duplicate samples {gl[i]} vs {gl[j]} (r>0.9 on cohort-z top-variance genes)", result=int((R > 0.9).sum())))
                for sid in b.index[(R > 0.9).any(axis=0)]:  # later cohort loses the duplicated sample
                    dup_rows.append(dict(cohort=gl[j], sample_id=str(sid), matched_with=gl[i]))
    pd.DataFrame(dup_rows, columns=["cohort", "sample_id", "matched_with"]).drop_duplicates(["cohort", "sample_id"]) \
        .to_csv(D["proc"] / "overlap_drop.csv", index=False)
    if dup_rows:
        deviation(f"{len(set((r['cohort'], r['sample_id']) for r in dup_rows))} GEO samples duplicate an earlier cohort and are removed from "
                  "analysis (data/processed/overlap_drop.csv); re-run prognosis/benchmark/classification after preprocess.")
    ov.append(dict(check="TCGA patients vs GEO: ID namespaces disjoint (TCGA barcodes vs GSM)", result="no ID overlap by construction"))
    save_csv(pd.DataFrame(ov), "train_validation_overlap_checks.csv")
    for r in ov:
        if isinstance(r["result"], (int, np.integer)) and r["result"] > 0 and "duplicate" in r["check"]:
            deviation(f"Overlap warning: {r['check']} = {r['result']} -> external cohorts may not be fully independent of each other")


# ----------------------------------------------------------------------------------------
# PROGNOSIS: models
# ----------------------------------------------------------------------------------------
def univariate_cox(Z, t, e):
    rows = []
    for g in Z.columns:
        d = pd.DataFrame({"x": Z[g].values, "T": t, "E": e})
        try:
            c = CoxPHFitter().fit(d, "T", "E")
            s = c.summary.loc["x"]
            rows.append(dict(gene=g, coef=s["coef"], HR=s["exp(coef)"], HR_lo95=s["exp(coef) lower 95%"],
                             HR_hi95=s["exp(coef) upper 95%"], p=s["p"], n=len(d), events=int(np.sum(e))))
        except Exception:
            rows.append(dict(gene=g, coef=np.nan, HR=np.nan, HR_lo95=np.nan, HR_hi95=np.nan, p=np.nan, n=len(d), events=int(np.sum(e))))
    u = pd.DataFrame(rows)
    u["q_BH"] = bh_adjust(u.p.values)
    return u


def choose_penalizer(Z, t, e, cfg, seed):
    res = {}
    skf = StratifiedKFold(cfg["inner_folds"], shuffle=True, random_state=seed)
    for pen in cfg["penalizer_grid"]:
        cs = []
        for tr, va in skf.split(Z, e):
            try:
                d = Z.iloc[tr].copy()
                d["T"], d["E"] = t[tr], e[tr]
                c = CoxPHFitter(penalizer=pen, l1_ratio=0.0).fit(d, "T", "E")
                cs.append(cindex(t[va], e[va], Z.iloc[va].values @ c.params_[Z.columns].values))
            except Exception:
                pass
        if cs:
            res[pen] = float(np.mean(cs))
    return max(res, key=lambda k: (res[k], k)) if res else max(cfg["penalizer_grid"])


def fit_prog_model(X, t, e, genes, cfg, seed):
    t, e = np.asarray(t, float), np.asarray(e, int)
    genes = [g for g in genes if g in X.columns]
    mu, sd = zfit(X[genes])
    keep = sd.dropna().index.tolist()
    Z = (X[keep] - mu[keep]) / sd[keep]
    Z = Z.fillna(0)
    uni = univariate_cox(Z, t, e)
    sel = uni.loc[uni.p < cfg["p_screen"], "gene"].tolist()
    fb = False
    if len(sel) < cfg["min_features"]:
        sel = uni.sort_values("p").head(cfg["min_features"]).gene.tolist()
        fb = True
    if not sel:
        raise ModelFitError("no features available")
    pen = choose_penalizer(Z[sel], t, e, cfg, seed)
    d = Z[sel].copy()
    d["T"], d["E"] = t, e
    try:
        c = CoxPHFitter(penalizer=pen, l1_ratio=0.0).fit(d, "T", "E")
    except Exception as ex:
        raise ModelFitError(str(ex))
    coef = c.params_[sel]
    risk = Z[sel].values @ coef.values
    return dict(selected=sel, coef=coef.to_dict(), mean=mu[sel].to_dict(), sd=sd[sel].to_dict(), penalizer=pen, fallback=fb,
                uni=uni, train_risk=risk, cut_median=float(np.median(risk)),
                cut_tertiles=[float(np.quantile(risk, 1 / 3)), float(np.quantile(risk, 2 / 3))], n_train=len(t), n_events=int(e.sum()))


def model_risk(m, X, mode):
    sel = m["selected"]
    A = X[sel]
    if mode == "train":
        mu, sd = pd.Series(m["mean"]), pd.Series(m["sd"])
    else:
        mu, sd = A.mean(), A.std(ddof=0).where(lambda s: s > 1e-8)
    Z = ((A - mu[sel]) / sd[sel]).fillna(0)
    return Z.values @ np.array([m["coef"][g] for g in sel])


def surv_metrics(t, e, lp, cfg, seed, cut=None):
    t, e, lp = np.asarray(t, float), np.asarray(e, int), np.asarray(lp, float)
    o = dict(n=len(t), events=int(e.sum()), cindex=np.nan, cindex_lo95=np.nan, cindex_hi95=np.nan, hr_per_sd=np.nan,
             hr_lo95=np.nan, hr_hi95=np.nan, hr_p=np.nan, calibration_slope=np.nan, calib_slope_lo95=np.nan,
             calib_slope_hi95=np.nan, logrank_p_high_vs_low=np.nan, n_high=np.nan)
    if o["events"] < 3:
        o["note"] = "too few events"
        return o
    try:
        o["cindex"] = cindex(t, e, lp)
        o["cindex_lo95"], o["cindex_hi95"] = boot_ci(lambda a, b, c: cindex(a, b, c), len(t), cfg["n_boot"], seed, t, e, lp)
        d = pd.DataFrame({"x": (lp - lp.mean()) / (lp.std() or 1), "T": t, "E": e})
        s = CoxPHFitter().fit(d, "T", "E").summary.loc["x"]
        o.update(hr_per_sd=s["exp(coef)"], hr_lo95=s["exp(coef) lower 95%"], hr_hi95=s["exp(coef) upper 95%"], hr_p=s["p"])
        d["x"] = lp
        s = CoxPHFitter().fit(d, "T", "E").summary.loc["x"]
        o.update(calibration_slope=s["coef"], calib_slope_lo95=s["coef lower 95%"], calib_slope_hi95=s["coef upper 95%"])
        if cut is not None:
            g = lp > cut
            o["n_high"] = int(g.sum())
            if g.sum() >= 5 and (~g).sum() >= 5:
                o["logrank_p_high_vs_low"] = logrank_test(t[g], t[~g], e[g], e[~g]).p_value
    except Exception as ex:
        o["note"] = f"estimation error: {str(ex)[:80]}"
    return o


def fit_clin(df, cols):
    d = df[cols + ["time", "event"]].dropna()
    c = CoxPHFitter(penalizer=1e-3).fit(d.rename(columns={"time": "T", "event": "E"}), "T", "E")
    return c.params_.to_dict()


def lp_from(params, df):
    return sum(df[k].values * v for k, v in params.items())


def eligible_features(panel, coh):
    g = set(panel.gene_symbol)
    for X, _ in coh.values():
        g &= set(X.columns)
    return sorted(g)


def prognosis_stage(cfg, args):
    panel, pstatus, pver = panel_gate(cfg, args)
    coh = available_cohorts(cfg)
    if COHORT_TCGA not in coh:
        raise DataUnavailable("TCGA-COAD processed data unavailable; cannot train the prognosis model")
    ext = {k: v for k, v in coh.items() if k != COHORT_TCGA}
    feats = eligible_features(panel, coh)
    miss = sorted(set(panel.gene_symbol) - set(feats))
    LOG.info("Eligible features (panel ∩ all loaded cohorts): %d of %d; excluded: %s", len(feats), len(panel), miss)
    if len(feats) < 5:
        raise DataUnavailable("fewer than 5 eligible panel genes across cohorts")
    X, c = coh[COHORT_TCGA]
    c = valid_surv(c)
    X = X.loc[c.index]
    seed = cfg["seed"]
    if int(c.event.sum()) < cfg["min_events"]:
        raise DataUnavailable(f"TCGA events {int(c.event.sum())} < min_events")
    # descriptive stats
    desc = pd.DataFrame([dict(cohort=COHORT_TCGA, n=len(c), events=int(c.event.sum()), median_followup_months=float(c.time.median()),
                              median_age=float(c.age.median()), pct_stage_missing=float(100 * c.stage.isna().mean()),
                              pct_age_missing=float(100 * c.age.isna().mean()))])
    save_csv(desc, "tcga_descriptive_stats.csv")
    # ---- training-only fit + freeze
    m = fit_prog_model(X, c.time.values, c.event.values, feats, cfg, seed)
    uni = m["uni"].copy()
    uni["selected_for_model"] = uni.gene.isin(m["selected"])
    uni.insert(0, "cohort", COHORT_TCGA)
    uni["panel_version"] = pver
    uni = uni.merge(panel[["gene_symbol", "gene_role"]], left_on="gene", right_on="gene_symbol", how="left").drop(columns="gene_symbol")
    save_csv(uni, "univariate_cox_results.csv")
    save_csv(pd.DataFrame(dict(gene=m["selected"], screening_threshold_p=cfg["p_screen"], fallback_used=m["fallback"],
                               training_cohort=COHORT_TCGA, panel_version=pver)), "selected_features.csv")
    cc = c.dropna(subset=["age", "stage"])
    cc = cc.assign(risk=model_risk(m, X.loc[cc.index], "train"))
    clin_p = fit_clin(cc, ["age", "stage"]) if len(cc) > 30 else {}
    comb_p = fit_clin(cc, ["age", "stage", "risk"]) if len(cc) > 30 else {}
    frozen = dict(version=VERSION, panel_version=pver, panel_status=pstatus, seed=seed, endpoint="OS (months)",
                  features=m["selected"], coef=m["coef"], train_mean=m["mean"], train_sd=m["sd"], penalizer=m["penalizer"],
                  risk_median=m["cut_median"], risk_tertiles=m["cut_tertiles"], clinical_params=clin_p, combined_params=comb_p,
                  eligible_features=feats, n_train=m["n_train"], events_train=m["n_events"], fallback_used=m["fallback"],
                  frozen_at=dt.datetime.now().isoformat(timespec="seconds"))
    (D["proc"] / "frozen_prognosis_model.json").write_text(json.dumps(frozen, indent=1), encoding="utf-8")
    coefdf = pd.DataFrame([dict(model="molecular_ridge_cox", term=g, coefficient=m["coef"][g], train_mean=m["mean"][g], train_sd=m["sd"][g],
                                penalizer=m["penalizer"], n_train=m["n_train"], panel_version=pver) for g in m["selected"]] +
                          [dict(model="clinical_cox", term=k, coefficient=v, n_train=len(cc)) for k, v in clin_p.items()] +
                          [dict(model="combined_cox", term=k, coefficient=v, n_train=len(cc)) for k, v in comb_p.items()])
    save_csv(coefdf, "prognosis_model_coefficients.csv")
    LOG.info("Frozen model: %d features, penalizer=%s from grid %s (frozen BEFORE external validation)",
             len(m["selected"]), m["penalizer"], cfg["penalizer_grid"])
    # ---- TCGA cross-validation (selection + tuning inside every training fold)
    metrics, scores = [], []
    if len(cc) >= 50:
        rskf = RepeatedStratifiedKFold(n_splits=cfg["cv_folds"], n_repeats=cfg["cv_repeats"], random_state=seed)
        rows = []
        Xc = X.loc[cc.index]
        for k, (tr, te) in enumerate(rskf.split(np.zeros(len(cc)), cc.event.values)):
            try:
                mf = fit_prog_model(Xc.iloc[tr], cc.time.values[tr], cc.event.values[tr], feats, cfg, seed + k)
                rte, rtr = model_risk(mf, Xc.iloc[te], "train"), model_risk(mf, Xc.iloc[tr], "train")
                tr_df, te_df = cc.iloc[tr].assign(risk=rtr), cc.iloc[te].assign(risk=rte)
                cp, bp = fit_clin(tr_df, ["age", "stage"]), fit_clin(tr_df, ["age", "stage", "risk"])
                tt, ee = cc.time.values[te], cc.event.values[te]
                rows.append(dict(split=k, molecular=cindex(tt, ee, rte), clinical=cindex(tt, ee, lp_from(cp, te_df)),
                                 combined=cindex(tt, ee, lp_from(bp, te_df))))
            except Exception as ex:
                LOG.warning("CV split %d failed: %s", k, str(ex)[:100])
        cvd = pd.DataFrame(rows)
        save_csv(cvd, "tcga_cv_splits.csv")
        for mod in ("molecular", "clinical", "combined"):
            if len(cvd):
                v = cvd[mod].values
                metrics.append(dict(cohort=f"{COHORT_TCGA} (CV, training-only selection)", endpoint="OS", model=mod, scheme="cross_validation",
                                    n=len(cc), events=int(cc.event.sum()), cindex=v.mean(), cindex_lo95=np.percentile(v, 2.5),
                                    cindex_hi95=np.percentile(v, 97.5), note=f"mean over {len(v)} CV test folds; interval = 2.5-97.5 pct of fold values (not a CI)"))
    # ---- TCGA in-sample (descriptive only)
    r_all = model_risk(m, X, "train")
    sm = surv_metrics(c.time.values, c.event.values, r_all, cfg, seed, m["cut_median"])
    metrics.append(dict(cohort=f"{COHORT_TCGA} (training, in-sample; optimistic)", endpoint="OS", model="molecular", scheme="in_sample", **sm))
    for nm, rr in [(COHORT_TCGA, r_all)]:
        sub = c.assign(risk_score=rr, split="training_in_sample")
        scores.append(sub)
    # ---- frozen external validation
    for name, (Xe, ce) in ext.items():
        ce = valid_surv(ce)
        Xe = Xe.loc[ce.index]
        miss_g = [g for g in frozen["features"] if g not in Xe.columns]
        if miss_g or len(ce) < 20:
            metrics.append(dict(cohort=name, endpoint="NA", model="molecular", scheme="external_frozen", n=len(ce), events=int(ce.event.sum()) if len(ce) else 0,
                                note=f"not evaluated: missing genes {miss_g} or n<20"))
            continue
        ep = ce.endpoint.iloc[0]
        r = model_risk(frozen_as_model(frozen), Xe, "cohort")
        sm = surv_metrics(ce.time.values, ce.event.values, r, cfg, seed, frozen["risk_median"])
        sm["note"] = sm.get("note", "") + ("; endpoint mismatch (DFS vs OS-trained) - exploratory" if ep != "OS" else "")
        metrics.append(dict(cohort=name, endpoint=ep, model="molecular", scheme="external_frozen", **sm))
        scores.append(ce.assign(risk_score=r, split="external_frozen"))
        ece = ce.assign(risk=r).dropna(subset=["age", "stage"])
        if clin_p and len(ece) >= 30:
            for mod, pp in (("clinical", clin_p), ("combined", comb_p)):
                try:
                    s2 = surv_metrics(ece.time.values, ece.event.values, lp_from(pp, ece), cfg, seed)
                    metrics.append(dict(cohort=name, endpoint=ep, model=mod, scheme="external_frozen",
                                        **{k: v for k, v in s2.items() if k not in ("logrank_p_high_vs_low", "n_high")},
                                        note=f"frozen TCGA coefficients; complete-case age+stage n={len(ece)}"))
                except Exception as ex:
                    LOG.warning("clinical model eval failed in %s: %s", name, ex)
        else:
            metrics.append(dict(cohort=name, endpoint=ep, model="clinical/combined", scheme="external_frozen", n=len(ece),
                                note="not evaluated: age and/or stage missing"))
    md = pd.DataFrame(metrics)
    md["panel_version"], md["model_version"], md["seed"] = pver, VERSION, seed
    save_csv(md, "prognosis_validation_metrics.csv")
    sc = pd.concat(scores)
    sc["risk_group"] = np.where(sc.risk_score > frozen["risk_median"], "high", "low")
    sc["risk_tertile"] = pd.cut(sc.risk_score, [-np.inf] + frozen["risk_tertiles"] + [np.inf], labels=["low", "medium", "high"]).astype(str)
    save_csv(sc[["cohort", "sample_id", "patient_id", "split", "risk_score", "risk_group", "risk_tertile", "time", "event", "endpoint", "age", "stage"]],
             "prognosis_patient_scores.csv")
    ex_ok = md[(md.scheme == "external_frozen") & (md.model == "molecular") & md.cindex.notna()]
    summ("prognosis_external_validation", "completed" if len(ex_ok) else "no cohort evaluable", n=int(ex_ok.n.sum()) if len(ex_ok) else 0,
         metric="; ".join(f"{r.cohort}: C={r.cindex:.3f} [{r.cindex_lo95:.3f},{r.cindex_hi95:.3f}]" for r in ex_ok.itertuples()),
         note=f"panel={pstatus}")


def frozen_as_model(fz):
    return dict(selected=fz["features"], coef=fz["coef"], mean=fz["train_mean"], sd=fz["train_sd"])


# ----------------------------------------------------------------------------------------
# NOVELTY 1: random-gene-set and proliferation benchmarking (identical procedure)
# ----------------------------------------------------------------------------------------
def benchmark_stage(cfg, args):
    panel, pstatus, pver = panel_gate(cfg, args)
    fp = D["proc"] / "frozen_prognosis_model.json"
    if not fp.exists():
        raise DataUnavailable("frozen model missing; run the prognosis stage first")
    fz = json.loads(fp.read_text(encoding="utf-8"))
    coh = available_cohorts(cfg)
    X, c = coh[COHORT_TCGA]
    c = valid_surv(c)
    X = X.loc[c.index]
    ext = {}
    for n, (Xe, ce) in coh.items():
        if n != COHORT_TCGA:
            ce = valid_surv(ce)
            if len(ce) >= 20:
                ext[n] = (Xe.loc[ce.index], ce)
    if not ext:
        raise DataUnavailable("no evaluable external cohort for benchmarking")
    obs = pd.read_csv(D["csv"] / "prognosis_validation_metrics.csv")
    obs = obs[(obs.scheme == "external_frozen") & (obs.model == "molecular")].set_index("cohort")
    common = set(X.columns)
    for Xe, _ in ext.values():
        common &= set(Xe.columns)
    prolif = [g for g in cfg["proliferation_markers"] if g in common]
    universe = sorted(set(X.columns[X.mean() >= 1.0]) & common - set(panel.gene_symbol) - set(prolif))
    k = len(fz["eligible_features"])
    LOG.info("Random-control universe: %d genes; set size k=%d; N=%d", len(universe), k, cfg["n_random_prognosis"])
    out_p = D["csv"] / "random_gene_set_benchmark.csv"
    done = pd.DataFrame()
    fp_payload = json.dumps(dict(
        ext={n: sorted(map(str, ce.index)) for n, (_, ce) in sorted(ext.items())}, panel=fz.get("panel_version"),
        grid=list(cfg["penalizer_grid"]), seed=cfg["seed"], k=k, p_screen=cfg["p_screen"], inner=cfg["inner_folds"],
        train_n=fz.get("n_train"), version=VERSION), sort_keys=True)
    fingerprint = hashlib.sha256(fp_payload.encode()).hexdigest()[:16]
    meta_p = D["csv"] / "random_gene_set_benchmark.meta.json"
    prev_fp = json.loads(meta_p.read_text(encoding="utf-8")).get("fingerprint") if meta_p.exists() else None
    if out_p.exists():
        old = pd.read_csv(out_p)
        if prev_fp == fingerprint:
            done = old[(old.task == "prognosis") & (old.seed == cfg["seed"]) & (old.set_size == k)]
        else:
            LOG.warning("random-control cache discarded: cohorts/panel/grid changed since it was computed (fingerprint %s -> %s)",
                        prev_fp, fingerprint)
    done_ids = set(done.set_id.unique()) if len(done) else set()
    rows = done.to_dict("records")
    rng = np.random.default_rng(cfg["seed"])
    sets = [list(rng.choice(universe, size=k, replace=False)) for _ in range(cfg["n_random_prognosis"])]

    def evaluate(genes, sid, label):
        try:
            m = fit_prog_model(X, c.time.values, c.event.values, genes, cfg, cfg["seed"])
        except ModelFitError:
            return []
        out = []
        for n, (Xe, ce) in ext.items():
            if not set(m["selected"]).issubset(Xe.columns):
                continue
            r = model_risk(m, Xe, "cohort")
            out.append(dict(task="prognosis", cohort=n, set_id=sid, set_label=label, set_size=k, n_selected=len(m["selected"]),
                            metric="cindex", value=cindex(ce.time.values, ce.event.values, r), seed=cfg["seed"]))
        return out
    for i, g in enumerate(sets):
        sid = f"random_{i:04d}"
        if sid in done_ids:
            continue
        rows += evaluate(g, sid, "random")
        if (i + 1) % 25 == 0:
            LOG.info("random controls: %d/%d", i + 1, len(sets))
            pd.DataFrame(rows).to_csv(out_p, index=False)
    prolif_rows = evaluate(prolif, "proliferation_markers", "proliferation") if len(prolif) >= 5 else []
    rows = [r for r in rows if r["set_id"] != "proliferation_markers"] + prolif_rows
    bench = pd.DataFrame(rows)
    old_cls = pd.read_csv(out_p) if out_p.exists() else pd.DataFrame()
    if len(old_cls):
        bench = pd.concat([bench, old_cls[old_cls.task != "prognosis"]]).drop_duplicates(["task", "cohort", "set_id"], keep="first")
    bench.to_csv(out_p, index=False)
    meta_p.write_text(json.dumps(dict(fingerprint=fingerprint, written=dt.datetime.now().isoformat(timespec="seconds"))), encoding="utf-8")
    comp = []
    for n in ext:
        rc = bench[(bench.task == "prognosis") & (bench.cohort == n) & (bench.set_label == "random")].value.values
        o = obs.loc[n, "cindex"] if n in obs.index else np.nan
        pr = bench[(bench.task == "prognosis") & (bench.cohort == n) & (bench.set_label == "proliferation")].value
        ep = ext[n][1].endpoint.iloc[0]
        base = dict(cohort=n, endpoint=ep, n_controls=len(rc), seed=cfg["seed"], control_set_size=k)
        comp.append(dict(base, signature="ferroptosis_panel_frozen", cindex=o, empirical_p_one_sided=empirical_p(o, rc),
                         percentile_vs_random=float(np.mean(rc <= o) * 100) if len(rc) else np.nan, note="same procedure applied to controls"))
        if len(rc):
            comp.append(dict(base, signature="random_sets_median", cindex=float(np.median(rc)), note="screening+ridge identical to panel"))
            comp.append(dict(base, signature="random_sets_95th_percentile", cindex=float(np.percentile(rc, 95))))
        if len(pr):
            comp.append(dict(base, signature="proliferation_markers (author-compiled, exploratory)", cindex=float(pr.iloc[0]),
                             note="not a published validated signature"))
        mm = pd.read_csv(D["csv"] / "prognosis_validation_metrics.csv")
        for mod in ("clinical", "combined"):
            r = mm[(mm.cohort == n) & (mm.model == mod) & (mm.scheme == "external_frozen")]
            if len(r) and pd.notna(r.cindex.iloc[0]):
                comp.append(dict(base, signature=f"{mod}_model_frozen", cindex=float(r.cindex.iloc[0]), note="TCGA-frozen age+stage" + ("+risk" if mod == "combined" else "")))
    save_csv(pd.DataFrame(comp), "signature_comparison.csv")
    sp = [r for r in comp if r["signature"] == "ferroptosis_panel_frozen"]
    summ("novelty1_specificity_prognosis", "completed", n=len(sets),
         metric="; ".join(f"{r['cohort']}: C={r['cindex']:.3f}, empirical p={r['empirical_p_one_sided']:.3f}" for r in sp if pd.notna(r["cindex"])),
         note=f"{cfg['n_random_prognosis']} random sets of size {k}, seed {cfg['seed']}; panel={pstatus}")


# ----------------------------------------------------------------------------------------
# BRANCH B: classification (dMMR/MSI) - verified labels only
# ----------------------------------------------------------------------------------------
def cls_model(seed):
    return make_pipeline(StandardScaler(), LogisticRegressionCV(Cs=[0.01, 0.1, 1, 10], cv=3, scoring="roc_auc",
                                                                 class_weight="balanced", max_iter=3000, random_state=seed))


def oof_proba(X, y, seed, reps=1, folds=5):
    p = np.zeros(len(y))
    for r in range(reps):
        for tr, te in StratifiedKFold(folds, shuffle=True, random_state=seed + r).split(X, y):
            m = cls_model(seed).fit(X[tr], y[tr])
            p[te] += m.predict_proba(X[te])[:, 1] / reps
    return p


def cls_metrics(y, p, cfg, seed, thr=0.5):
    pred = (p >= thr).astype(int)
    tn, fp, fn, tp = confusion_matrix(y, pred, labels=[0, 1]).ravel()
    pr, rc, f1, _ = precision_recall_fscore_support(y, pred, average="binary", zero_division=0)
    spec = tn / (tn + fp) if tn + fp else np.nan
    lo, hi = boot_ci(lambda a, b: roc_auc_score(a, b), len(y), cfg["n_boot"], seed, y, p)
    return dict(n=len(y), n_pos=int(y.sum()), prevalence=float(y.mean()), auc=roc_auc_score(y, p), auc_lo95=lo, auc_hi95=hi,
                precision=pr, recall=rc, f1=f1, specificity=spec, balanced_accuracy=(rc + spec) / 2, tn=tn, fp=fp, fn=fn, tp=tp, threshold=thr)


def ferro_score(E, panel, min_cov):
    drv = panel[panel.gene_role == "driver"].gene_symbol.tolist()
    sup = panel[panel.gene_role == "suppressor"].gene_symbol.tolist()
    dp, sp = [g for g in drv if g in E.columns], [g for g in sup if g in E.columns]
    info = dict(drivers_present=len(dp), drivers_total=len(drv), suppressors_present=len(sp), suppressors_total=len(sup),
                missing=";".join(sorted((set(drv) | set(sup)) - set(dp) - set(sp))))
    if not drv or not sup or len(dp) / len(drv) < min_cov or len(sp) / len(sup) < min_cov:
        raise DataUnavailable(f"FerroScore coverage below {min_cov:.0%}: {info}")
    mu, sd = zfit(E[dp + sp])
    Z = ((E[dp + sp] - mu) / sd).fillna(0)
    return Z[dp].mean(axis=1) - Z[sp].mean(axis=1), info


def classification_stage(cfg, args):
    panel, pstatus, pver = panel_gate(cfg, args)
    coh = {n: v for n, v in available_cohorts(cfg).items() if n != COHORT_TCGA}
    lab = {}
    for n, (X, c) in coh.items():
        y = c.mmr.dropna()
        if len(y) >= 40 and y.sum() >= 10 and (len(y) - y.sum()) >= 10:
            lab[n] = (X.loc[y.index], y.astype(int).values)
    if not lab:
        summ("classification_dMMR", "not run", note="no cohort with verified MMR/MSI labels (>=10 per class)")
        raise DataUnavailable("no cohort with verified MMR/MSI labels; classification not run (TCGA MSI not retrieved)")
    prim = max(lab, key=lambda k: len(lab[k][1]))
    seed = cfg["seed"]
    Xp, yp = lab[prim]
    genes = [g for g in panel.gene_symbol if g in Xp.columns]
    zc = lambda X, gs: ((X[gs] - X[gs].mean()) / X[gs].std(ddof=0).replace(0, np.nan)).fillna(0)  # cohort-wise, label-free
    LOG.info("Classification: primary cohort %s n=%d dMMR=%d (%.1f%%)", prim, len(yp), yp.sum(), 100 * yp.mean())
    rows, preds, curves = [], [], {}
    prolif = [g for g in cfg["proliferation_markers"] if g in Xp.columns]
    for label, gs, tag in [("ferroptosis_panel_LR", genes, ""), ("proliferation_markers_LR", prolif, "exploratory")]:
        if len(gs) < 3:
            continue
        p = oof_proba(zc(Xp, gs).values, yp, seed, reps=cfg["cv_repeats"])
        m = cls_metrics(yp, p, cfg, seed)
        rows.append(dict(cohort=prim, scheme="oof_cv_patient_level", model=label, label="dMMR(1) vs pMMR(0) from GEO characteristics", note=tag, **m))
        preds += [dict(cohort=prim, scheme="oof_cv", model=label, sample_id=s, y_true=int(a), prob=float(b), pred=int(b >= 0.5))
                  for s, a, b in zip(Xp.index, yp, p)]
        fpr, tpr, _ = roc_curve(yp, p)
        curves[label] = (fpr, tpr, m["auc"], m["auc_lo95"], m["auc_hi95"])
    fs, info = ferro_score(Xp, panel, cfg["min_panel_coverage"])
    rows.append(dict(cohort=prim, scheme="unsupervised_score", model="FerroScore_raw (driver-suppressor), orientation not prespecified",
                     label="dMMR vs pMMR", n=len(yp), n_pos=int(yp.sum()), prevalence=float(yp.mean()), auc=roc_auc_score(yp, fs.values),
                     note="AUC<0.5 means higher FerroScore in pMMR; descriptive only"))
    # random gene-set controls (same pipeline)
    uni = sorted(set(Xp.columns[Xp.mean() >= Xp.mean().median()]) - set(panel.gene_symbol) - set(prolif))
    rng = np.random.default_rng(seed)
    k, N = len(genes), cfg["n_random_classification"]
    ctrl = []
    for i in range(N):
        gs = list(rng.choice(uni, size=k, replace=False))
        ctrl.append(roc_auc_score(yp, oof_proba(zc(Xp, gs).values, yp, seed + 1000 + i, reps=1)))
        if (i + 1) % 25 == 0:
            LOG.info("classification random controls %d/%d", i + 1, N)
    obs_auc = rows[0]["auc"]
    ep = empirical_p(obs_auc, ctrl)
    rows[0]["empirical_p_vs_random"], rows[0]["n_random_controls"] = ep, N
    rb = pd.DataFrame([dict(task="mmr_classification", cohort=prim, set_id=f"random_{i:04d}", set_label="random", set_size=k, n_selected=k,
                            metric="auc", value=v, seed=seed) for i, v in enumerate(ctrl)])
    p0 = D["csv"] / "random_gene_set_benchmark.csv"
    base = pd.read_csv(p0) if p0.exists() else pd.DataFrame()
    if len(base):
        base = base[base.task != "mmr_classification"]
    pd.concat([base, rb]).to_csv(p0, index=False)
    # external labelled cohorts
    final = cls_model(seed).fit(zc(Xp, genes).values, yp)
    for n, (Xe, ye) in lab.items():
        if n == prim or not set(genes).issubset(Xe.columns):
            continue
        pe = final.predict_proba(zc(Xe, genes).values)[:, 1]
        rows.append(dict(cohort=n, scheme="external_frozen", model="ferroptosis_panel_LR", label="dMMR vs pMMR", note=f"trained on {prim}", **cls_metrics(ye, pe, cfg, seed)))
        preds += [dict(cohort=n, scheme="external", model="ferroptosis_panel_LR", sample_id=s, y_true=int(a), prob=float(b), pred=int(b >= 0.5)) for s, a, b in zip(Xe.index, ye, pe)]
    md = pd.DataFrame(rows)
    md["panel_version"], md["model_version"], md["seed"] = pver, VERSION, seed
    save_csv(md, "classification_metrics.csv")
    save_csv(pd.DataFrame(preds), "classification_predictions.csv")
    pd.to_pickle(curves, D["proc"] / "roc_curves.pkl")
    summ("classification_dMMR", "completed", n=len(yp), metric=f"{prim}: AUC={obs_auc:.3f} (OOF CV), empirical p vs {N} random sets={ep:.3f}",
         note=f"labels from GEO characteristics; FerroScore coverage {info}; panel={pstatus}")


# ----------------------------------------------------------------------------------------
# BRANCH C: drug response (GDSC2 primary, CTRP exploratory replication)
# ----------------------------------------------------------------------------------------
PRESPECIFIED = ("primary", "clinical")


def tier_map(cfg):
    """User drug_panel.csv wins; drugs from the code defaults that the user file lacks (e.g. new clinical tier) are added."""
    user = pd.read_csv(D["cfg"] / "drug_panel.csv")
    base = pd.read_csv(io.StringIO(DRUG_PANEL))
    dp = pd.concat([user, base], ignore_index=True)
    dp["k_norm"] = dp.drug_name.map(norm)
    dp = dp.drop_duplicates("k_norm", keep="first")
    return {r.k_norm: r.tier for r in dp.itertuples()}


def drug_stage(cfg, args):
    panel, pstatus, pver = panel_gate(cfg, args)
    ep, rp = D["proc"] / "gdsc2_expr_crc.csv.gz", D["proc"] / "gdsc2_response_crc.csv"
    if not (ep.exists() and rp.exists()):
        raise DataUnavailable("GDSC2 processed data missing; run fetch+preprocess (or place manual files)")
    E = pd.read_csv(ep, index_col=0)
    E.index = E.index.astype(str)
    R = pd.read_csv(rp, dtype={"COSMIC_ID": str})
    fs, info = ferro_score(E, panel, cfg["min_panel_coverage"])
    fsd = fs.rename("ferroscore").reset_index()
    fsd.to_csv(D["proc"] / "gdsc_ferroscores.csv", index=False)
    save_csv(E[[g for g in panel.gene_symbol if g in E.columns]], "gdsc2_expression_processed.csv", index=True)
    save_csv(R, "gdsc2_drug_response_processed.csv")
    tiers = tier_map(cfg)
    R["tier"] = R.DRUG_NAME.map(lambda x: tiers.get(norm(x), "unlisted"))
    pres = {t: sorted(R[R.tier == t].DRUG_NAME.unique()) for t in ("primary", "clinical", "exploratory")}
    LOG.info("Prespecified compounds present in GDSC2 CRC data: %s", pres)
    rows = []
    n_lines = E.shape[0]
    for (did, dn), g in R.groupby(["DRUG_ID", "DRUG_NAME"]):
        d = g.merge(fsd, on="COSMIC_ID")
        for metric in ("LN_IC50", "AUC"):
            dd = d.dropna(subset=[metric])
            base = dict(drug_id=did, drug_name=dn, tier=tiers.get(norm(dn), "unlisted"), response_metric=metric, n_cell_lines=len(dd),
                        n_crc_lines_with_expression=n_lines, n_missing_response=n_lines - len(dd), panel_version=pver)
            if len(dd) < cfg["gdsc_min_cell_lines"]:
                rows.append(dict(base, rho=np.nan, p=np.nan, note=f"n<{cfg['gdsc_min_cell_lines']}: not tested"))
                continue
            rho, p, lo, hi = spearman_ci(dd.ferroscore, dd[metric])
            rows.append(dict(base, rho=rho, rho_lo95=lo, rho_hi95=hi, p=p, note="LN_IC50 is natural log (uM) as given by GDSC" if metric == "LN_IC50" else "AUC lower = more sensitive"))
    res = pd.DataFrame(rows)
    for metric in ("LN_IC50", "AUC"):
        m = (res.response_metric == metric) & res.p.notna()
        res.loc[m, "q_BH_all_drugs"] = bh_adjust(res.loc[m, "p"].values)
        pm = m & res.tier.isin(PRESPECIFIED)
        res.loc[pm, "q_BH_prespecified_primary"] = bh_adjust(res.loc[pm, "p"].values)
    res["analysis_class"] = res.tier.map({"primary": "prespecified (primary)", "clinical": "prespecified (clinical)"}).fillna("exploratory")
    save_csv(res, "gdsc2_drug_correlations.csv")
    tested = res[res.p.notna() & (res.response_metric == "LN_IC50")]
    n_sig = int((tested.q_BH_all_drugs < 0.05).sum())
    summ("novelty3_gdsc2_ferroscore_vs_response", "completed", n=n_lines,
         metric=f"{len(tested)} drugs tested; nominal p<0.05: {int((tested.p < 0.05).sum())}; BH q<0.05: {n_sig}; primary compounds present: {pres['primary'] or 'none'}; clinical: {pres['clinical'] or 'none'}",
         note=f"FerroScore = mean z(drivers) - mean z(suppressors); coverage {info['drivers_present']}/{info['drivers_total']} drivers, {info['suppressors_present']}/{info['suppressors_total']} suppressors; panel={pstatus}")
    if not pres["primary"]:
        deviation("None of ML162/ML210/RSL3 present in the retrieved GDSC2 CRC data: canonical-ferroptosis-inducer hypothesis could not be tested in GDSC2.")
    # ---- regression models (sample size permitting; drug selection by n only, never by outcome)
    cnt = R.groupby("DRUG_NAME").COSMIC_ID.nunique()
    cand = [d for d in sum(pres.values(), []) if cnt.get(d, 0) >= cfg["model_min_cell_lines"]]
    cand += [d for d in cnt[cnt >= cfg["model_min_cell_lines"]].sort_values(ascending=False, kind="stable").index if d not in cand]
    cand = cand[:cfg["max_model_drugs"]]
    mrows, pdata = [], []
    gcols = [g for g in panel.gene_symbol if g in E.columns]
    for dn in cand:
        d = R[R.DRUG_NAME == dn].merge(fsd, on="COSMIC_ID").dropna(subset=["LN_IC50"]).drop_duplicates("COSMIC_ID")
        Xg, Xs, y = E.loc[d.COSMIC_ID, gcols].fillna(E[gcols].mean()).values, d[["ferroscore"]].values, d.LN_IC50.values
        res_by = {"baseline_mean": [], "ferroscore_linear": [], "panel_ridge": []}
        oof_keep = None
        for rep in range(3):
            oof = {k: np.zeros(len(y)) for k in res_by}
            for tr, te in KFold(5, shuffle=True, random_state=cfg["seed"] + rep).split(Xg):  # one row per cell line => cell-line-level splits
                oof["baseline_mean"][te] = y[tr].mean()
                A = np.c_[np.ones(len(tr)), Xs[tr]]
                b = np.linalg.lstsq(A, y[tr], rcond=None)[0]
                oof["ferroscore_linear"][te] = np.c_[np.ones(len(te)), Xs[te]] @ b
                mm = make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-1, 3, 9))).fit(Xg[tr], y[tr])
                oof["panel_ridge"][te] = mm.predict(Xg[te])
            for k in res_by:
                res_by[k].append((mean_absolute_error(y, oof[k]), math.sqrt(mean_squared_error(y, oof[k])), r2_score(y, oof[k])))
            if rep == 0:
                oof_keep = oof
        for k, v in res_by.items():
            a = np.array(v)
            mrows.append(dict(drug_name=dn, tier=tiers.get(norm(dn), "unlisted"), model=k, n_cell_lines=len(y), cv="5-fold x3 repeats, cell-line-level",
                              MAE=a[:, 0].mean(), RMSE=a[:, 1].mean(), R2=a[:, 2].mean(), R2_sd_over_repeats=a[:, 2].std(), target="LN_IC50 (cell-line, not a clinical dose)"))
        for k, o in oof_keep.items():
            pdata += [dict(drug_name=dn, model=k, cosmic_id=c, observed=a, predicted=b) for c, a, b in zip(d.COSMIC_ID, y, o)]
    if mrows:
        save_csv(pd.DataFrame(mrows), "drug_model_metrics.csv")
        pd.DataFrame(pdata).to_csv(D["proc"] / "drug_model_predictions.csv", index=False)
        summ("drug_response_models", "completed", n=len(cand), metric=f"{len(cand)} drugs modelled", note="exploratory baseline regressions; cell-line IC50 only")
    else:
        summ("drug_response_models", "not run", note=f"no drug with >= {cfg['model_min_cell_lines']} CRC cell lines")
    # ---- drug-sensitivity classification (prespecified drugs only; sensitive = LN_IC50 below the drug's median across CRC lines)
    crow = []
    for dn in pres["primary"] + pres["clinical"]:
        d = R[R.DRUG_NAME == dn].merge(fsd, on="COSMIC_ID").dropna(subset=["LN_IC50"]).drop_duplicates("COSMIC_ID")
        if len(d) < cfg["model_min_cell_lines"]:
            continue
        yb = (d.LN_IC50.values < np.median(d.LN_IC50.values)).astype(int)
        if yb.sum() < 5 or (len(yb) - yb.sum()) < 5:
            continue
        feats_ = {"panel_LR": E.loc[d.COSMIC_ID, gcols].fillna(E[gcols].mean()).values, "ferroscore_LR": d[["ferroscore"]].values}
        for mname, Xm in feats_.items():
            pr_ = oof_proba(Xm, yb, cfg["seed"], reps=3)
            lo_, hi_ = boot_ci(lambda a, b: roc_auc_score(a, b), len(yb), cfg["n_boot"], cfg["seed"], yb, pr_)
            crow.append(dict(drug_name=dn, tier=tiers.get(norm(dn), "unlisted"), model=mname, n_cell_lines=len(yb), n_sensitive=int(yb.sum()),
                             auc=roc_auc_score(yb, pr_), auc_lo95=lo_, auc_hi95=hi_,
                             label="sensitive = LN_IC50 below drug median (cell-line level)",
                             note="out-of-fold; median split uses all lines (label only); exploratory; not a clinical response"))
    if crow:
        save_csv(pd.DataFrame(crow), "drug_sensitivity_classification.csv")
        summ("drug_sensitivity_classification", "completed", n=len(crow),
             metric="; ".join(f"{r['drug_name']}/{r['model']}: AUC={r['auc']:.2f}" for r in crow[:6]),
             note="exploratory; prespecified drugs only; cell-line sensitivity, not dose")
    else:
        summ("drug_sensitivity_classification", "not run", note=f"no prespecified drug with >= {cfg['model_min_cell_lines']} CRC lines and >=5 per class")
    # ---- CTRP replication
    cp = D["proc"] / "ctrp_response_crc.csv"
    if not cp.exists():
        summ("ctrp_replication", "unavailable", note="CTRP data not retrieved/processed")
        return
    C = pd.read_csv(cp)
    names = R[["COSMIC_ID", "CELL_LINE_NAME"]].drop_duplicates()
    names["ccl_norm"] = names.CELL_LINE_NAME.map(norm)
    C = C.merge(names.merge(fsd, on="COSMIC_ID"), on="ccl_norm")
    save_csv(C[["ccl_name", "cpd_name", "AUC", "n_experiments", "COSMIC_ID", "ferroscore"]], "ctrp_drug_response_processed.csv")
    gd = res[res.response_metric == "AUC"].copy()
    gd["cpd_norm"] = gd.drug_name.map(norm)
    out = []
    for cn, g in C.groupby("cpd_norm"):
        if len(g) < cfg["gdsc_min_cell_lines"]:
            continue
        rho, p, lo, hi = spearman_ci(g.ferroscore, g.AUC)
        row = dict(compound=g.cpd_name.iloc[0], tier=tiers.get(cn, "unlisted"), ctrp_n_cell_lines=len(g), ctrp_rho=rho, ctrp_p=p,
                   ctrp_rho_lo95=lo, ctrp_rho_hi95=hi, ctrp_metric="AUC (mean over experiments)")
        gm = gd[gd.cpd_norm == cn]
        if len(gm) and gm.p.notna().any():
            gm = gm.sort_values("n_cell_lines", ascending=False).iloc[0]
            row.update(gdsc_drug=gm.drug_name, gdsc_n=gm.n_cell_lines, gdsc_rho=gm.rho, gdsc_p=gm.p, gdsc_metric="AUC",
                       same_compound_by_name=True, direction_concordant=bool(np.sign(gm.rho) == np.sign(rho)),
                       gdsc_nominal_discovery=bool(gm.p < 0.05),
                       replicated=bool(gm.p < 0.05 and np.sign(gm.rho) == np.sign(rho) and p < 0.05))
        else:
            row.update(same_compound_by_name=False, note="no matching GDSC2 compound with testable data")
        out.append(row)
    if out:
        o = pd.DataFrame(out)
        o["ctrp_q_BH"] = bh_adjust(o.ctrp_p.values)
        o["caveat"] = "compounds matched by normalised name; cell lines matched by normalised name; expression from GDSC; exploratory"
        o.to_csv(D["csv"] / "ctrp_replication_results.csv", index=False)
        LOG.info("saved results/csv/ctrp_replication_results.csv (%d rows)", len(o))
        nm = int(o.same_compound_by_name.sum())
        summ("ctrp_replication", "completed", n=len(o), metric=f"{len(o)} CTRP compounds tested; {nm} matched a GDSC2 compound; replicated (GDSC p<.05 & concordant & CTRP p<.05): {int(o.get('replicated', pd.Series(dtype=bool)).fillna(False).sum())}",
             note="exploratory")
    else:
        summ("ctrp_replication", "not run", note="no CTRP compound with enough matched CRC lines")


# ----------------------------------------------------------------------------------------
# FIGURES
# ----------------------------------------------------------------------------------------
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 100})


def savefig(fig, name):
    D["fig"].mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(D["fig"] / f"{name}.png", dpi=300)
    fig.savefig(D["fig"] / f"{name}.pdf")
    plt.close(fig)
    LOG.info("figure saved: %s", name)


def fig_flow(cfg):
    d = read_csv_opt("cohort_inclusion_exclusion.csv")
    if d is None:
        raise DataUnavailable("no inclusion table")
    cohorts = list(dict.fromkeys(d.cohort))
    fig, axs = plt.subplots(1, len(cohorts), figsize=(3.6 * len(cohorts), 5.5), squeeze=False)
    for ax, c in zip(axs[0], cohorts):
        s = d[d.cohort == c].reset_index(drop=True)
        ax.axis("off")
        ax.set_title(c, fontweight="bold")
        for i, r in s.iterrows():
            y = 1 - i * (1 / max(len(s), 1)) * 0.95
            ax.text(0.5, y, f"{r.step}\nn = {int(r.n_after)}" + (f"\n(excluded {int(r.n_excluded)})" if r.n_excluded > 0 else ""), ha="center", va="top",
                    bbox=dict(boxstyle="round", fc="#eef3fb", ec="#4a6fa5"), fontsize=7.5, transform=ax.transAxes)
    fig.suptitle("Cohort inclusion / exclusion flow (counts from downloaded data)")
    savefig(fig, "fig01_cohort_flowchart")


def fig_missing(cfg):
    d = read_csv_opt("clinical_data_clean.csv")
    if d is None:
        raise DataUnavailable("no clinical table")
    f = ["time", "event", "age", "stage", "mmr"]
    m = d.groupby("cohort")[f].apply(lambda x: 100 * x.isna().mean())
    ax = m.plot(kind="bar", figsize=(8, 4))
    ax.set_ylabel("% missing"); ax.set_title("Clinical-field missingness by cohort (n per cohort in table)")
    ax.set_xticklabels([f"{c}\n(n={int((d.cohort == c).sum())})" for c in m.index], rotation=0)
    savefig(ax.figure, "fig02_missing_data")


def fig_expr(cfg, panel):
    coh = available_cohorts(cfg)
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    meds = []
    for n, (X, _) in coh.items():
        g = [x for x in panel.gene_symbol if x in X.columns]
        axs[0].hist(X[g].values.ravel(), bins=60, histtype="step", density=True, label=f"{n} (n={len(X)})")
        meds.append(pd.Series(X[g].median(axis=1).values, name=n))
    axs[0].set_xlabel("log2 expression (panel genes; platform-specific scale)"); axs[0].set_ylabel("density"); axs[0].legend(fontsize=7)
    axs[1].boxplot([m.values for m in meds], labels=[m.name for m in meds]); axs[1].set_ylabel("per-sample median log2 expression")
    axs[1].tick_params(axis="x", rotation=20)
    fig.suptitle("Expression distributions (RNA-seq and microarray scales are NOT directly comparable)")
    savefig(fig, "fig03_expression_distributions")


def fig_heat(cfg, panel):
    X, c = load_cohort(COHORT_TCGA)
    g = [x for x in panel.gene_symbol if x in X.columns]
    Z = ((X[g] - X[g].mean()) / X[g].std(ddof=0)).clip(-3, 3)
    fs, _ = ferro_score(X, panel, cfg["min_panel_coverage"])
    Z = Z.loc[fs.sort_values().index]
    role = panel.set_index("gene_symbol").gene_role.reindex(g)
    order = list(role.sort_values(kind="stable").index)
    fig, ax = plt.subplots(figsize=(9, 5))
    im = ax.imshow(Z[order].T.values, aspect="auto", cmap="RdBu_r", vmin=-3, vmax=3)
    ax.set_yticks(range(len(order))); ax.set_yticklabels([f"{x} ({role[x][:3]})" for x in order], fontsize=6)
    ax.set_xlabel(f"TCGA-COAD primary tumours (n={len(Z)}), ordered by FerroScore"); fig.colorbar(im, label="z-score")
    ax.set_title("Ferroptosis gene-panel expression")
    savefig(fig, "fig04_panel_heatmap")


def fig_forest(cfg):
    u = read_csv_opt("univariate_cox_results.csv")
    if u is None:
        raise DataUnavailable("no univariate results")
    u = u.dropna(subset=["HR"]).sort_values("HR")
    fig, ax = plt.subplots(figsize=(6, 0.28 * len(u) + 1.5))
    y = np.arange(len(u))
    ax.errorbar(u.HR, y, xerr=[u.HR - u.HR_lo95, u.HR_hi95 - u.HR], fmt="o", ms=3, color="k", lw=0.8)
    ax.scatter(u.HR[u.selected_for_model], y[u.selected_for_model.values], color="crimson", zorder=3, s=14, label="selected")
    ax.axvline(1, ls="--", c="grey"); ax.set_xscale("log"); ax.set_yticks(y); ax.set_yticklabels(u.gene, fontsize=6)
    ax.set_xlabel(f"HR per SD (95% CI), TCGA-COAD OS, n={int(u.n.iloc[0])}, events={int(u.events.iloc[0])}"); ax.legend(fontsize=7)
    ax.set_title("Univariate Cox (training cohort)")
    savefig(fig, "fig05_univariate_cox_forest")


def fig_km(cfg):
    s = read_csv_opt("prognosis_patient_scores.csv")
    m = read_csv_opt("prognosis_validation_metrics.csv")
    if s is None:
        raise DataUnavailable("no scores")
    cs = list(dict.fromkeys(s.cohort))
    fig, axs = plt.subplots(1, len(cs), figsize=(4.2 * len(cs), 4), squeeze=False)
    for ax, c in zip(axs[0], cs):
        d = s[s.cohort == c].dropna(subset=["time", "event"])
        for grp, col in (("low", "#2a7fba"), ("high", "#d62728")):
            g = d[d.risk_group == grp]
            if len(g):
                KaplanMeierFitter().fit(g.time, g.event, label=f"{grp} (n={len(g)})").plot_survival_function(ax=ax, color=col, ci_show=False)
        h, l = d[d.risk_group == "high"], d[d.risk_group == "low"]
        p = logrank_test(h.time, l.time, h.event, l.event).p_value if len(h) > 4 and len(l) > 4 else np.nan
        ax.set_title(f"{c}\n{d.endpoint.iloc[0]}; log-rank p={p:.3g}" + ("\n(training, in-sample)" if d.split.iloc[0].startswith("training") else "\n(frozen external)"), fontsize=8)
        ax.set_xlabel("months")
    fig.suptitle("Kaplan-Meier: risk groups split at the frozen TCGA training median")
    savefig(fig, "fig06_kaplan_meier")


def fig_perf(cfg):
    m = read_csv_opt("prognosis_validation_metrics.csv")
    if m is None:
        raise DataUnavailable("no metrics")
    a = m[(m.model == "molecular") & m.cindex.notna()]
    fig, ax = plt.subplots(figsize=(7, 4))
    y = np.arange(len(a))
    lo = (a.cindex - a.cindex_lo95).clip(lower=0).fillna(0)
    hi = (a.cindex_hi95 - a.cindex).clip(lower=0).fillna(0)
    ax.errorbar(a.cindex, y, xerr=[lo, hi], fmt="o", color="k")
    ax.axvline(0.5, ls="--", c="grey"); ax.set_yticks(y); ax.set_yticklabels([f"{r.cohort}\n(n={int(r.n)}, ev={int(r.events)})" for r in a.itertuples()], fontsize=7)
    ax.set_xlabel("C-index (bootstrap 95% CI; CV rows show fold range)"); ax.set_title("Prognostic model performance across cohorts")
    savefig(fig, "fig07_prognosis_performance")
    b = m[m.model.isin(["clinical", "molecular", "combined"]) & m.cindex.notna() & (m.scheme != "in_sample")]
    if len(b):
        pv = b.pivot_table(index="cohort", columns="model", values="cindex")
        ax = pv.plot(kind="bar", figsize=(8, 4)); ax.axhline(0.5, ls="--", c="grey"); ax.set_ylabel("C-index")
        ax.set_title("Clinical-only vs molecular-only vs clinical+molecular"); plt.xticks(rotation=15, ha="right", fontsize=7)
        savefig(ax.figure, "fig08_clinical_vs_combined")


def fig_cls(cfg):
    cv = D["proc"] / "roc_curves.pkl"
    if not cv.exists():
        raise DataUnavailable("no classification curves")
    curves = pd.read_pickle(cv)
    fig, ax = plt.subplots(figsize=(4.6, 4.4))
    for k, (f, t, a, lo, hi) in curves.items():
        ax.plot(f, t, label=f"{k}\nAUC {a:.2f} [{lo:.2f},{hi:.2f}]")
    ax.plot([0, 1], [0, 1], "k--", lw=0.7); ax.set_xlabel("1 - specificity"); ax.set_ylabel("sensitivity"); ax.legend(fontsize=6)
    ax.set_title("dMMR vs pMMR ROC (out-of-fold, patient-level)")
    savefig(fig, "fig09_roc")
    p = read_csv_opt("classification_predictions.csv")
    p = p[(p.model == "ferroptosis_panel_LR") & (p.scheme == "oof_cv")]
    cm = confusion_matrix(p.y_true, p.pred, labels=[0, 1])
    fig, ax = plt.subplots(figsize=(3.6, 3.4))
    ax.imshow(cm, cmap="Blues")
    for i in range(2):
        for j in range(2):
            ax.text(j, i, cm[i, j], ha="center", va="center")
    ax.set_xticks([0, 1]); ax.set_yticks([0, 1]); ax.set_xticklabels(["pMMR", "dMMR"]); ax.set_yticklabels(["pMMR", "dMMR"])
    ax.set_xlabel("predicted"); ax.set_ylabel("true"); ax.set_title(f"Confusion matrix (thr 0.5, n={len(p)})")
    savefig(fig, "fig10_confusion_matrix")


def fig_random(cfg):
    b = read_csv_opt("random_gene_set_benchmark.csv")
    c = read_csv_opt("signature_comparison.csv")
    cls = read_csv_opt("classification_metrics.csv")
    if b is None:
        raise DataUnavailable("no random benchmark")
    panels = [(k, g) for k, g in b[b.set_label == "random"].groupby(["task", "cohort"])]
    fig, axs = plt.subplots(1, max(len(panels), 1), figsize=(3.8 * max(len(panels), 1), 3.6), squeeze=False)
    for ax, ((task, coh), g) in zip(axs[0], panels):
        ax.hist(g.value, bins=25, color="#bbb", edgecolor="k")
        if task == "prognosis" and c is not None:
            o = c[(c.cohort == coh) & (c.signature == "ferroptosis_panel_frozen")]
            if len(o):
                ax.axvline(o.cindex.iloc[0], c="crimson", lw=2, label=f"panel {o.cindex.iloc[0]:.3f}\nemp. p={o.empirical_p_one_sided.iloc[0]:.3f}")
        elif cls is not None:
            o = cls[(cls.cohort == coh) & (cls.model == "ferroptosis_panel_LR") & (cls.scheme == "oof_cv_patient_level")]
            if len(o):
                ax.axvline(o.auc.iloc[0], c="crimson", lw=2, label=f"panel {o.auc.iloc[0]:.3f}\nemp. p={o.empirical_p_vs_random.iloc[0]:.3f}")
        ax.set_title(f"{task}\n{coh} (N={len(g)} controls)", fontsize=8); ax.set_xlabel("C-index" if task == "prognosis" else "AUC"); ax.legend(fontsize=7)
    fig.suptitle("Random gene-set benchmark (identical procedure, fixed seed)")
    savefig(fig, "fig11_random_benchmark")


def fig_drug(cfg):
    r = read_csv_opt("gdsc2_drug_correlations.csv")
    if r is None:
        raise DataUnavailable("no drug correlations")
    t = r[(r.response_metric == "LN_IC50") & r.p.notna()].copy()
    if t.empty:
        raise DataUnavailable("no tested drugs")
    R = pd.read_csv(D["proc"] / "gdsc2_response_crc.csv", dtype={"COSMIC_ID": str})
    fsd = pd.read_csv(D["proc"] / "gdsc_ferroscores.csv", dtype={"COSMIC_ID": str})
    top = pd.concat([t[t.tier.isin(PRESPECIFIED)], t.sort_values("p")]).drop_duplicates("drug_id").head(6)
    fig, axs = plt.subplots(2, 3, figsize=(11, 6.5), squeeze=False)
    for ax, r_ in zip(axs.ravel(), top.itertuples()):
        d = R[R.DRUG_ID == r_.drug_id].merge(fsd, on="COSMIC_ID").dropna(subset=["LN_IC50"])
        ax.scatter(d.ferroscore, d.LN_IC50, s=12)
        ax.set_title(f"{r_.drug_name} [{r_.analysis_class}]\nn={r_.n_cell_lines}, rho={r_.rho:.2f}, p={r_.p:.3g}, q={r_.q_BH_all_drugs:.2g}", fontsize=7)
        ax.set_xlabel("FerroScore"); ax.set_ylabel("LN_IC50 (log uM, cell line)")
    for ax in axs.ravel()[len(top):]:
        ax.axis("off")
    fig.suptitle("FerroScore vs GDSC2 response in CRC cell lines (shown: prespecified compounds first, then smallest p; selection for display only)")
    savefig(fig, "fig12_ferroscore_vs_response")
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    for tier, col in (("unlisted", "#aaa"), ("exploratory", "#f5a623"), ("clinical", "#2a7fba"), ("primary", "crimson")):
        s = t[t.tier == tier]
        ax.scatter(s.rho, -np.log10(s.q_BH_all_drugs.clip(lower=1e-6)), s=18, c=col, label=f"{tier} (n={len(s)})")
    ax.axhline(-np.log10(0.05), ls="--", c="k", lw=0.7); ax.set_xlabel("Spearman rho (FerroScore vs LN_IC50)"); ax.set_ylabel("-log10 BH q")
    ax.legend(fontsize=7); ax.set_title(f"All {len(t)} tested GDSC2 drugs (nonsignificant results included)")
    savefig(fig, "fig13_drug_correlation_summary")


def fig_drugmodel(cfg):
    p = D["proc"] / "drug_model_predictions.csv"
    m = read_csv_opt("drug_model_metrics.csv")
    if not p.exists() or m is None:
        raise DataUnavailable("no drug models")
    d = pd.read_csv(p)
    d = d[d.model == "panel_ridge"]
    ds = list(dict.fromkeys(d.drug_name))[:6]
    fig, axs = plt.subplots(2, 3, figsize=(10, 6), squeeze=False)
    for ax, dn in zip(axs.ravel(), ds):
        g = d[d.drug_name == dn]
        r2 = m[(m.drug_name == dn) & (m.model == "panel_ridge")].R2.iloc[0]
        ax.scatter(g.observed, g.predicted, s=10); lim = [g.observed.min(), g.observed.max()]; ax.plot(lim, lim, "k--", lw=0.7)
        ax.set_title(f"{dn} (n={len(g)}, CV R2={r2:.2f})", fontsize=8); ax.set_xlabel("observed LN_IC50"); ax.set_ylabel("predicted (out-of-fold)")
    for ax in axs.ravel()[len(ds):]:
        ax.axis("off")
    savefig(fig, "fig14_drug_model_pred_vs_obs")


def figures_stage(cfg, args):
    panel, _, _ = load_panel(cfg)
    jobs = [("flowchart", lambda: fig_flow(cfg)), ("missing", lambda: fig_missing(cfg)), ("expression", lambda: fig_expr(cfg, panel)),
            ("heatmap", lambda: fig_heat(cfg, panel)), ("forest", lambda: fig_forest(cfg)), ("km", lambda: fig_km(cfg)),
            ("performance", lambda: fig_perf(cfg)), ("classification", lambda: fig_cls(cfg)), ("random", lambda: fig_random(cfg)),
            ("drug", lambda: fig_drug(cfg)), ("drugmodel", lambda: fig_drugmodel(cfg))]
    for nm, fn in jobs:
        try:
            fn()
        except Exception as e:
            plt.close("all")
            LOG.warning("figure '%s' skipped: %s", nm, str(e)[:160])


# ----------------------------------------------------------------------------------------
# REPORTING
# ----------------------------------------------------------------------------------------
def md_table(df, maxrows=40, cols=None):
    if df is None or len(df) == 0:
        return "_not available_"
    d = df.head(maxrows)
    if cols:
        d = d[[c for c in cols if c in d.columns]]
    d = d.copy()
    for c in d.columns:
        if d[c].dtype.kind == "f":
            d[c] = d[c].map(lambda v: "" if pd.isna(v) else f"{v:.4g}")
        else:
            d[c] = d[c].map(lambda v: "" if (not isinstance(v, (list, dict)) and pd.isna(v)) else v)
    h = "| " + " | ".join(map(str, d.columns)) + " |\n|" + "---|" * len(d.columns) + "\n"
    return h + "\n".join("| " + " | ".join(str(x) for x in r) + " |" for r in d.itertuples(index=False))


def methodology(cfg, panel_status, pver):
    return f"""# Methodology (generated {dt.datetime.now():%Y-%m-%d %H:%M}; pipeline v{VERSION})

**Gene panel:** {panel_status}; version hash `{pver}`. A PROVISIONAL panel is NOT the original project panel.

## Data
- TCGA-COAD: GDC API (`/files`, `/cases`), open-access *{cfg['workflow_type']}* gene-level files; value used = `tpm_unstranded`, protein-coding genes, duplicate gene names averaged, transformed log2(TPM+1). Primary tumour = barcode sample code 01 and API sample_type 'Primary Tumor'; per patient the lexicographically first sample barcode (then file name) is kept. File MD5 checked against GDC.
- OS: days_to_death if Dead; if Alive, diagnoses.days_to_last_follow_up, else the maximum follow_ups.days_to_follow_up; months = days/30.4375; valid if time>0 and event known. Stage: AJCC I-IV -> 1-4 (sub-stages merged). Age: age_at_diagnosis/365.25. TCGA DFS and MSI are not retrieved.
- GEO ({', '.join(cfg['geo']['cohorts'])}): series matrix + platform annotation. Expression scale inferred (max>100 -> log2(x+1)). Probes -> gene symbols; probes without symbol or with multiple genes dropped; per gene the probe with highest mean expression kept (outcome-free). Clinical fields mapped by regular expressions over characteristics keys (see data_dictionary.csv; override in config `geo.overrides`). Time units taken from field names, otherwise assumed months (flagged).
- GDSC2: fitted dose-response (LN_IC50 natural-log µM, AUC) + RMA basal expression; CRC lines = {cfg['gdsc_lineage_field']} == '{cfg['gdsc_lineage_value']}'; joined on verified COSMIC IDs; duplicate (line,drug) pairs averaged; missing responses never imputed.
- CTRPv2: curves/experiment/cell-line/compound tables joined by IDs; CRC = site '{cfg['ctrp_site']}'; AUC averaged over experiments; cell lines matched to GDSC expression by normalised name; exploratory.

## Harmonisation
RNA-seq and microarray values are never pooled. Each cohort is analysed separately. For external prognosis evaluation, the selected genes are z-scored **within each external cohort** (label-free, cohort-specific standardisation) and multiplied by frozen TCGA coefficients. No outcome from a validation cohort is used for preprocessing, selection, tuning or thresholds.

## Prognosis (Novelty 2)
Eligible features = panel genes present in TCGA and all loaded GEO cohorts (gene availability only). In training data: z-score (training mean/SD) -> univariate Cox per gene -> keep p < {cfg['p_screen']} (fallback: lowest-p {cfg['min_features']} genes, flagged) -> ridge Cox (lifelines `penalizer`, l1_ratio=0) with penalizer chosen from {cfg['penalizer_grid']} by {cfg['inner_folds']}-fold inner CV C-index. Risk = sum(coef x z). Model, preprocessing parameters and thresholds are frozen to `data/processed/frozen_prognosis_model.json` before any external cohort is touched. TCGA CV ({cfg['cv_folds']}-fold x {cfg['cv_repeats']}) repeats the whole screening+tuning+fit inside each training fold. Clinical model: Cox on age (years) + stage (1-4); combined: age + stage + risk. Metrics: Harrell C-index (bootstrap {cfg['n_boot']} 95% CI), HR per SD of risk, calibration slope (Cox coefficient on the linear predictor; ideal 1), log-rank for high vs low at the training median. Endpoint: OS; if a cohort has no OS, DFS/RFS is used and flagged as endpoint-mismatched/exploratory.

## Specificity benchmarking (Novelty 1)
{cfg['n_random_prognosis']} random gene sets (same size as the eligible panel, drawn without replacement from expressed genes [TCGA mean log2TPM>=1] present in all cohorts, excluding panel and proliferation genes; seed {cfg['seed']}) pass through the **identical** procedure (screening, tuning, fit in TCGA, frozen application to each external cohort). Empirical one-sided p = (1 + #controls with C >= observed)/(1 + N). A proliferation-marker set (author-compiled, not a published validated signature) is run identically as an exploratory comparator. An association with survival alone is not evidence of ferroptosis specificity.

## Classification
Target: dMMR (1) vs pMMR (0) from GEO characteristics (verified labels only). L2 logistic regression (inner CV for C) with in-fold standardisation; patient-level stratified {cfg['cv_folds']}-fold CV x {cfg['cv_repeats']} repeats (out-of-fold probabilities averaged); metrics at probability threshold 0.5 (prespecified). AUC CI by bootstrap over patients. {cfg['n_random_classification']} random gene-set controls with the same pipeline. FerroScore AUC is descriptive (orientation not prespecified).

## Drug response (Novelty 3)
FerroScore = mean z(driver genes) - mean z(suppressor genes), z-scored across the cell lines analysed; requires >= {cfg['min_panel_coverage']:.0%} coverage per role. Spearman correlation vs LN_IC50 and AUC per drug (n >= {cfg['gdsc_min_cell_lines']}), 95% CI by Fisher/Bonett-Wright, Benjamini-Hochberg q across all tested drugs per metric (and separately across prespecified primary compounds). Prespecified: ML162, ML210, RSL3; everything else exploratory (erastin explicitly exploratory). Models (n >= {cfg['model_min_cell_lines']} lines): mean baseline, FerroScore-linear, panel ridge; 5-fold x 3 CV, one row per cell line so splits are cell-line-level; drugs chosen by sample size only. Prespecified tiers: primary (ML162, ML210, RSL3) and clinical (5-fluorouracil, oxaliplatin, irinotecan); all other drugs exploratory. For prespecified drugs with >= {cfg['model_min_cell_lines']} lines a sensitivity CLASSIFIER is also fitted (sensitive = LN_IC50 below that drug's median; L2 logistic regression on the panel genes or on FerroScore alone; out-of-fold AUC, 3 repeats). The endpoint is in-vitro cell-line sensitivity, never a clinical dose or response. **Cell-line IC50 is not a clinical dose; no dosing or treatment recommendation is made.**
The dMMR vs pMMR classifier above tests whether the panel carries MSI information; it is not a classifier of a ferroptosis state.
CTRP replication: for compounds shared by normalised name, directional concordance and p in CTRP; 'replicated' = GDSC AUC p<0.05, same sign, CTRP p<0.05. Different assays/cell-line panels are not interchangeable; treat as exploratory.
"""


def reports_stage(cfg, args):
    D["rep"].mkdir(parents=True, exist_ok=True)
    panel, pstatus, pver = load_panel(cfg)
    s = read_state()
    man = pd.read_csv(D["prov"] / "download_manifest.csv") if (D["prov"] / "download_manifest.csv").exists() else pd.DataFrame()
    (D["rep"] / "methodology.md").write_text(methodology(cfg, pstatus, pver), encoding="utf-8")
    summ_df = pd.DataFrame(list(s["summary"].values()))
    if len(summ_df):
        save_csv(summ_df, "analysis_summary.csv")
    g = lambda n: read_csv_opt(n)
    ds = {}
    if len(man):
        for d, gdf in man.groupby("dataset"):
            ds[d] = dict(files=len(gdf), ok=int(gdf.retrieval_status.isin(["downloaded", "cached"]).sum()), failed=int((gdf.retrieval_status == "failed").sum()))
    L = [f"# Results summary ({dt.datetime.now():%Y-%m-%d %H:%M})", "",
         f"**Gene panel status: {pstatus}** (hash `{pver}`). " + ("All ferroptosis-specific statements below are provisional." if pstatus.startswith("PROV") else ""), "",
         "> This summary reports only quantities computed in this run. A finished script does not imply a successful analysis; see per-analysis status.", "",
         "## 1. Data retrieval", md_table(pd.DataFrame([dict(dataset=k, **v) for k, v in ds.items()])) if ds else "_no downloads recorded_", "",
         "## 2. Dataset inventory", md_table(g("dataset_inventory.csv")), "", "## 3. Cohort inclusion / exclusion", md_table(g("cohort_inclusion_exclusion.csv"), 60), "",
         "## 4. Analysis status", md_table(summ_df), ""]
    pv = g("prognosis_validation_metrics.csv")
    L += ["## 5. Novelty 2: frozen external validation (prognosis)", md_table(pv, 40, ["cohort", "endpoint", "model", "scheme", "n", "events", "cindex", "cindex_lo95", "cindex_hi95", "hr_per_sd", "hr_p", "calibration_slope", "logrank_p_high_vs_low", "note"]), ""]
    if pv is not None:
        ex = pv[(pv.scheme == "external_frozen") & (pv.model == "molecular") & pv.cindex.notna()]
        k = int((ex.cindex_lo95 > 0.5).sum())
        L += [f"External cohorts evaluated: {len(ex)}; cohorts whose bootstrap 95% CI for C-index excludes 0.5: {k}. "
              "Association with survival is not evidence of clinical utility.", ""]
    sc = g("signature_comparison.csv")
    L += ["## 6. Novelty 1: ferroptosis specificity", md_table(sc, 60, ["cohort", "endpoint", "signature", "cindex", "empirical_p_one_sided", "percentile_vs_random", "n_controls", "note"]), ""]
    if sc is not None:
        fe = sc[sc.signature == "ferroptosis_panel_frozen"]
        for r in fe.itertuples():
            L.append(f"- {r.cohort}: panel C={r.cindex:.3f}; empirical p={r.empirical_p_one_sided:.3f} vs {int(r.n_controls)} random sets "
                     + ("(not distinguishable from random at 0.05)" if r.empirical_p_one_sided >= 0.05 else "(better than most random sets; still not proof of ferroptosis-specific biology)"))
        L.append("")
    cm = g("classification_metrics.csv")
    L += ["## 7. Classification (dMMR vs pMMR)", md_table(cm, 20, ["cohort", "scheme", "model", "n", "n_pos", "auc", "auc_lo95", "auc_hi95", "balanced_accuracy", "f1", "empirical_p_vs_random", "note"]), ""]
    dr = g("gdsc2_drug_correlations.csv")
    L += ["## 8. Novelty 3: FerroScore vs drug response (GDSC2, CRC cell lines)"]
    if dr is not None:
        t = dr[dr.p.notna()]
        L += [f"Drug x metric tests: {len(t)}; nominal p<0.05: {int((t.p < 0.05).sum())}; BH q<0.05: {int((t.q_BH_all_drugs < 0.05).sum())}.", "",
              "Prespecified/exploratory compounds:", md_table(dr[dr.tier != "unlisted"], 20, ["drug_name", "tier", "response_metric", "n_cell_lines", "rho", "p", "q_BH_all_drugs", "q_BH_prespecified_primary", "note"]), "",
              "Ten smallest p-values (exploratory, multiplicity-adjusted q shown):", md_table(t.sort_values("p"), 10, ["drug_name", "response_metric", "n_cell_lines", "rho", "p", "q_BH_all_drugs"])]
    else:
        L.append("_not run_")
    L += ["", "### Drug models", md_table(g("drug_model_metrics.csv"), 40), "", "### Drug-sensitivity classification (prespecified drugs)", md_table(g("drug_sensitivity_classification.csv"), 20), "", "### CTRP replication", md_table(g("ctrp_replication_results.csv"), 30), "",
          "## 9. Deviations / limitations flagged during the run"] + [f"- {d}" for d in s["deviations"]] + [
          "", "## 10. Stage execution status", md_table(pd.DataFrame([dict(stage=k, **v) for k, v in s["stages"].items()])), "",
          "No dose prediction or treatment recommendation is made; cell-line IC50 is not a clinical dose."]
    (D["rep"] / "results_summary.md").write_text("\n".join(L), encoding="utf-8")
    lim = ["# Data limitations", "",
           "- TCGA MSI/MMR and DFS were not retrieved; MSI classification relies on MMR labels parsed from GEO characteristics (GSE39582 typically).",
           "- RNA-seq (TCGA) and microarray (GEO) differ in technology and scale; external evaluation uses cohort-specific z-scoring of selected genes, which assumes similar cohort case-mix.",
           "- GEO clinical fields are parsed heuristically (regex over characteristics keys); check `data_dictionary.csv`. Time units are assumed months when not stated.",
           "- GSE14333 typically lacks OS; if so, DFS is used and flagged (endpoint mismatch).",
           "- GEO cohorts may share patients (see train_validation_overlap_checks.csv); independence is not guaranteed.",
           "- The gene panel may be PROVISIONAL; replace it with the original project panel before drawing conclusions.",
           "- GDSC2/CTRP: few CRC lines (tens); ML162/ML210/RSL3 may be absent from GDSC2 and/or CTRP; CTRP lines/compounds matched by name.",
           "- Cell-line IC50/AUC is not a clinical dose; association is not mechanism.",
           "- TCGA-COAD contains colon tumours only, whereas GDSC lineage COREAD also includes rectal lines.",
           "- GEO samples duplicating an earlier cohort (r>0.9) are removed from analysis (data/processed/overlap_drop.csv); near-duplicates below the threshold are not detected.",
           "- The drug-sensitivity classifier uses a median split of cell-line LN_IC50 (label only); n is small, so AUCs have wide intervals.", "",
           "## Deviations recorded in this run"] + [f"- {d}" for d in s["deviations"]]
    (D["rep"] / "data_limitations.md").write_text("\n".join(lim), encoding="utf-8")
    lp = D["log"] / "run.log"
    if lp.exists():
        shutil.copy(lp, D["rep"] / "run_log.txt")
    LOG.info("reports written to %s", D["rep"])


# ----------------------------------------------------------------------------------------
# DRY RUN, FETCH / PREPROCESS STAGES, RUNNER
# ----------------------------------------------------------------------------------------
def dry_run(cfg):
    LOG.info("=== DRY RUN: validating sources (no large downloads) ===")
    res = []
    try:
        inv, total = tcga_file_inventory(cfg)
        LOG.info("GDC OK: %d files / %d cases / %d samples / %d aliquots", len(inv), inv.case_id.nunique(), inv.sample_id.nunique(), inv.aliquot_id.nunique())
        sel, _ = select_tcga_files(inv)
        LOG.info("Would download %d primary-tumour files (one per patient)", len(sel))
        res.append(("GDC TCGA-COAD", "OK", f"{len(sel)} files"))
    except Exception as e:
        res.append(("GDC TCGA-COAD", "FAIL", str(e)[:150]))
    urls = {f"GEO {a}": [f"https://ftp.ncbi.nlm.nih.gov/geo/series/{geo_stub(a)}/{a}/matrix/{a}_series_matrix.txt.gz"] for a in cfg["geo"]["cohorts"]}
    urls.update({"GDSC2 response": cfg["urls"]["gdsc2_response"], "GDSC expression": cfg["urls"]["gdsc_expression"], "CTRP": cfg["urls"]["ctrp"]})
    for k, ul in urls.items():
        ok = False
        for u in ul:
            try:
                r = SESSION.head(u, allow_redirects=True, timeout=30)
                if r.status_code == 200 and "html" not in r.headers.get("content-type", "").lower():
                    res.append((k, "OK", f"{u} ({r.headers.get('content-length', '?')} bytes)"))
                    ok = True
                    break
            except Exception as e:
                last = str(e)[:100]
        if not ok:
            res.append((k, "UNREACHABLE", f"{ul[0]} -> manual fallback (see README)"))
    for r in res:
        LOG.info("%-18s %-12s %s", *r)


def stage_fetch(cfg, args):
    ok = {}
    for label, fn in [("TCGA-COAD", lambda: fetch_tcga(cfg, args.refresh, args.dry_run))] + \
                     [(a, (lambda a=a: fetch_geo(a, cfg, args.refresh))) for a in cfg["geo"]["cohorts"]] + \
                     [("GDSC", lambda: fetch_gdsc(cfg, args.refresh)), ("CTRP", lambda: fetch_ctrp(cfg, args.refresh)),
                      ("DepMap", lambda: fetch_depmap(cfg, args.refresh))]:
        try:
            fn()
            ok[label] = "available"
        except Exception as e:
            ok[label] = "UNAVAILABLE: " + str(e)[:200]
            LOG.error("fetch %s failed: %s", label, str(e)[:300])
        finally:
            flush_manifest()
    s = read_state()
    s["datasets"].update(ok)
    write_state(s)
    if not any(v == "available" for v in ok.values()):
        raise DataUnavailable("no dataset could be retrieved")


def _processed(label):
    pr = D["proc"]
    need = {"GDSC2": ["gdsc2_expr_crc.csv.gz", "gdsc2_response_crc.csv"], "CTRP": ["ctrp_response_crc.csv"]}.get(
        label, [f"{label}_expr.csv.gz", f"{label}_clin.csv"])
    return all((pr / f).exists() for f in need)


def stage_preprocess(cfg, args):
    jobs = [(COHORT_TCGA, lambda: process_tcga(cfg))] + [(a, (lambda a=a: process_geo(a, cfg))) for a in cfg["geo"]["cohorts"]] + \
           [("GDSC2", lambda: process_gdsc(cfg)), ("CTRP", lambda: process_ctrp(cfg))]
    for label, fn in jobs:
        if not getattr(args, "reprocess", False) and _processed(label):
            LOG.info("preprocess %s: processed files exist, skipped (use --reprocess to redo)", label)
            continue
        try:
            fn()
        except Exception as e:
            LOG.error("preprocess %s failed: %s", label, str(e)[:300])
            LOG.debug(traceback.format_exc())
    preprocess_finalize(cfg, args)


# ----------------------------------------------------------------------------------------
# EXTENSION (v1.4.0): ferroptosis subtypes, cell-line GPX4-inhibitor sensitivity classifier + projection,
# overview figure, extension report
# ----------------------------------------------------------------------------------------
EXT_VERSION = "1.4.0"
COH_T = COHORT_TCGA
GPX4_TARGETS = ["ML162", "ML210", "RSL3"]


def zc(X, genes):
    """cohort-wise, label-free z-score (same rule the pipeline uses for external cohorts)."""
    A = X[genes].astype(float)
    sd = A.std(ddof=0).where(lambda s: s > 1e-8)
    return ((A - A.mean()) / sd).fillna(0)


def cox_hr(t, e, x, adj=None):
    nan = dict(logHR=np.nan, HR=np.nan, lo=np.nan, hi=np.nan, p=np.nan)
    try:
        d = pd.DataFrame({"x": np.asarray(x, float), "T": np.asarray(t, float), "E": np.asarray(e, int)})
        if adj is not None:
            for k in adj.columns:
                d[k] = np.asarray(adj[k], float)
        s = CoxPHFitter(penalizer=1e-4).fit(d, "T", "E").summary.loc["x"]
        return dict(logHR=s["coef"], HR=s["exp(coef)"], lo=s["exp(coef) lower 95%"], hi=s["exp(coef) upper 95%"], p=s["p"])
    except Exception:
        return nan


def assign_labels(Xc, genes, centroids, hi_id):
    """nearest frozen centroid in cohort-z space -> 1 = 'high ferroptosis' subtype."""
    Z = zc(Xc, genes).values
    d = ((Z[:, None, :] - centroids[None, :, :]) ** 2).sum(-1)
    return (d.argmin(1) == hi_id).astype(int)


def assoc(c, g):
    """Association of a binary group g (index = sample_id) with survival, stage, age and MMR in one cohort."""
    c = c.loc[g.index]
    o = dict(n=len(c), n_high=int(g.sum()))
    v = valid_surv(c)
    gv = g.loc[v.index]
    if len(v) >= 20 and int(v.event.sum()) >= 5 and gv.nunique() == 2:
        h, l = v[gv == 1], v[gv == 0]
        o.update(endpoint=v.endpoint.iloc[0] if "endpoint" in v else "", n_surv=len(v), events=int(v.event.sum()),
                 logrank_p=logrank_test(h.time, l.time, h.event, l.event).p_value)
        r = cox_hr(v.time, v.event, gv)
        o.update(logHR=r["logHR"], HR_high_vs_low=r["HR"], HR_lo95=r["lo"], HR_hi95=r["hi"], cox_p=r["p"])
        va = v.dropna(subset=["age", "stage"])
        if len(va) >= 30 and int(va.event.sum()) >= 5 and gv.loc[va.index].nunique() == 2:
            r2 = cox_hr(va.time, va.event, gv.loc[va.index], va[["age", "stage"]])
            o.update(HR_adj_age_stage=r2["HR"], HR_adj_lo95=r2["lo"], HR_adj_hi95=r2["hi"], p_adj=r2["p"], n_adj=len(va))
    st = c.stage.dropna()
    if len(st) >= 20:
        tab = pd.crosstab(g.loc[st.index], (st >= 3).astype(int)).reindex(index=[0, 1], columns=[0, 1], fill_value=0)
        o["pct_stage34_high"] = 100 * tab.loc[1, 1] / max(tab.loc[1].sum(), 1)
        o["pct_stage34_low"] = 100 * tab.loc[0, 1] / max(tab.loc[0].sum(), 1)
        o["stage34_fisher_p"] = stats.fisher_exact(tab.values)[1]
    ag = c.age.dropna()
    if len(ag) >= 20 and g.loc[ag.index].nunique() == 2:
        gg = g.loc[ag.index]
        o["age_median_high"], o["age_median_low"] = ag[gg == 1].median(), ag[gg == 0].median()
        o["age_mwu_p"] = stats.mannwhitneyu(ag[gg == 1], ag[gg == 0]).pvalue
    mm = c.mmr.dropna()
    if len(mm) >= 40 and mm.nunique() == 2 and g.loc[mm.index].nunique() == 2:
        gg = g.loc[mm.index]
        tab = pd.crosstab(gg, mm.astype(int)).reindex(index=[0, 1], columns=[0, 1], fill_value=0)
        o["pct_dMMR_high"], o["pct_dMMR_low"] = 100 * mm[gg == 1].mean(), 100 * mm[gg == 0].mean()
        o["mmr_fisher_p"] = stats.fisher_exact(tab.values)[1]
    return o


def assoc_light(c, g):
    """(logHR high-vs-low, dMMR fraction difference high-low) used for the random-set benchmark."""
    c = c.loc[g.index]
    v = valid_surv(c)
    gv = g.loc[v.index]
    lh = cox_hr(v.time, v.event, gv)["logHR"] if (len(v) >= 20 and int(v.event.sum()) >= 5 and gv.nunique() == 2) else np.nan
    mm, dm = c.mmr.dropna(), np.nan
    if len(mm) >= 40 and mm.nunique() == 2:
        gg = g.loc[mm.index]
        if gg.nunique() == 2:
            dm = mm[gg == 1].mean() - mm[gg == 0].mean()
    return lh, dm


def random_universe(cfg, coh, panel):
    X = coh[COH_T][0]
    common = set.intersection(*[set(v[0].columns) for v in coh.values()])
    return sorted((set(X.columns[X.mean() >= 1.0]) & common) - set(panel.gene_symbol))


# ----------------------------------------------------------------------------------------------------------------
# STAGE 1: ferroptosis subtypes
# ----------------------------------------------------------------------------------------------------------------
def subtype_stage(cfg, args):
    panel, pstatus, pver = panel_gate(cfg, args)
    coh = available_cohorts(cfg)
    if COH_T not in coh:
        raise DataUnavailable("TCGA-COAD processed data unavailable")
    genes = [g for g in panel.gene_symbol if all(g in v[0].columns for v in coh.values())]
    if len(genes) < 5:
        raise DataUnavailable("fewer than 5 panel genes shared by all cohorts")
    seed, N = cfg["seed"], cfg["n_random_classification"]
    X, c = coh[COH_T]
    Z = zc(X, genes)
    fs, _ = ferro_score(X, panel, cfg["min_panel_coverage"])
    km = KMeans(2, n_init=50, random_state=seed).fit(Z.values)
    lab = km.labels_
    hi_id = 1 if fs.values[lab == 1].mean() > fs.values[lab == 0].mean() else 0
    cent = km.cluster_centers_
    # stability (descriptive): bootstrap adjusted Rand index; silhouette for k=2..5 (k=2 is prespecified, not selected)
    rng = np.random.default_rng(seed)
    aris = []
    for b in range(100):
        i = rng.integers(0, len(Z), len(Z))
        kb = KMeans(2, n_init=10, random_state=seed + b).fit(Z.values[i])
        aris.append(adjusted_rand_score(lab, kb.predict(Z.values)))
    sil = {k: float(silhouette_score(Z.values, KMeans(k, n_init=20, random_state=seed).fit_predict(Z.values))) for k in (2, 3, 4, 5)}
    save_csv(pd.DataFrame([dict(k_prespecified=2, bootstrap_ARI_mean=np.mean(aris), bootstrap_ARI_p05=np.percentile(aris, 5),
                                   n_genes=len(genes), panel_version=pver, **{f"silhouette_k{k}": v for k, v in sil.items()})]),
                "subtype_stability.csv")
    frozen = dict(version=EXT_VERSION, genes=genes, centroids=cent.tolist(), high_cluster_id=hi_id, panel_version=pver,
                  rule="k=2 k-means on TCGA z-scores; high = cluster with higher mean (driver z - suppressor z); "
                       "other cohorts: cohort-wise z-score, nearest frozen centroid")
    (D["proc"] / "frozen_subtype_model.json").write_text(json.dumps(frozen, indent=1), encoding="utf-8")
    groups, rows, pat = {}, [], []
    for n, (Xc, cc) in coh.items():
        g = pd.Series((lab == hi_id).astype(int), index=X.index) if n == COH_T else \
            pd.Series(assign_labels(Xc, genes, cent, hi_id), index=Xc.index)
        groups[n] = g
        try:
            fsc, _ = ferro_score(Xc, panel, cfg["min_panel_coverage"])
        except DataUnavailable:
            fsc = pd.Series(np.nan, index=Xc.index)
        pat.append(pd.DataFrame(dict(cohort=n, sample_id=Xc.index, ferro_subtype=np.where(g.values == 1, "high", "low"),
                                     ferroscore=fsc.reindex(Xc.index).values)))
        rows.append(dict(cohort=n, **assoc(cc, g)))
    obs = pd.DataFrame(rows)
    obs["panel_version"] = pver
    save_csv(obs, "ferroptosis_subtype_associations.csv")
    save_csv(pd.concat(pat), "ferroptosis_subtype_patients.csv")
    # ---- random gene-set benchmark: identical clustering + projection + tests
    universe = random_universe(cfg, coh, panel)
    k = len(genes)
    obs_stat = {n: assoc_light(coh[n][1], groups[n]) for n in coh}
    rr = np.random.default_rng(seed + 7)
    ctrl = []
    for i in range(N):
        gs = list(rr.choice(universe, size=k, replace=False))
        kmr = KMeans(2, n_init=10, random_state=seed + i).fit(zc(X, gs).values)
        for n, (Xc, cc) in coh.items():
            g = pd.Series(kmr.labels_, index=Xc.index) if n == COH_T else \
                pd.Series(assign_labels(Xc, gs, kmr.cluster_centers_, 1), index=Xc.index)
            lh, dm = assoc_light(cc, g)
            ctrl.append(dict(set_id=f"random_{i:04d}", cohort=n, abs_logHR=abs(lh) if pd.notna(lh) else np.nan,
                             abs_dMMR_diff=abs(dm) if pd.notna(dm) else np.nan))
        if (i + 1) % 25 == 0:
            LOG.info("subtype random controls %d/%d", i + 1, N)
    ct = pd.DataFrame(ctrl)
    ct.to_csv(D["csv"] / "subtype_random_benchmark.csv", index=False)
    comp = []
    for n, (lh, dm) in obs_stat.items():
        cn = ct[ct.cohort == n]
        row = dict(cohort=n, n_controls=int(cn.set_id.nunique()), control_set_size=k, seed=seed)
        if pd.notna(lh):
            row.update(obs_logHR=lh, obs_HR=float(np.exp(lh)), emp_p_abs_logHR=empirical_p(abs(lh), cn.abs_logHR.values),
                       random_median_abs_logHR=float(cn.abs_logHR.median()))
        if pd.notna(dm):
            row.update(obs_dMMR_diff=dm, emp_p_abs_dMMR_diff=empirical_p(abs(dm), cn.abs_dMMR_diff.values),
                       random_median_abs_dMMR_diff=float(cn.abs_dMMR_diff.median()))
        comp.append(row)
    save_csv(pd.DataFrame(comp), "subtype_vs_random.csv")
    t = obs[obs.cohort == COH_T]
    msg = "; ".join(f"{r.cohort}: HR(high vs low)={r.HR_high_vs_low:.2f}, log-rank p={r.logrank_p:.3g}"
                    for r in obs.itertuples() if pd.notna(getattr(r, "HR_high_vs_low", np.nan)))
    summ("ext_ferroptosis_subtypes", "completed", n=len(X), metric=msg,
            note=f"k=2 prespecified; stability ARI={np.mean(aris):.2f}; {N} random sets; panel={pstatus}")


# ----------------------------------------------------------------------------------------------------------------
# cell-line data (GDSC2 response / CTRP / long file; expression from GDSC RMA, DepMap or processed GDSC)
# ----------------------------------------------------------------------------------------------------------------
def _cos(s):
    v = pd.to_numeric(pd.Series(s), errors="coerce").astype("Int64").astype(str)
    return v.where(v != "<NA>", "").values


def read_response_sources(cfg, metric):
    out = {}
    gd = D["raw"] / "gdsc"
    gd.mkdir(parents=True, exist_ok=True)
    gp = gd / "GDSC2_fitted_dose_response.xlsx"
    adopt_manual(gd, ["*GDSC2*dose*.xlsx", "*fitted*dose*.xlsx"], gp)
    if gp.exists():
        try:
            r = pd.read_excel(gp)
            fld = cfg["gdsc_lineage_field"]
            need = {"COSMIC_ID", "CELL_LINE_NAME", "DRUG_NAME", metric, fld}
            if need.issubset(r.columns):
                r = r[r[fld].astype(str) == cfg["gdsc_lineage_value"]]
                out["gdsc"] = pd.DataFrame(dict(cell=r.CELL_LINE_NAME.astype(str).values, cosmic=_cos(r.COSMIC_ID),
                                                drug=r.DRUG_NAME.astype(str).values, value=r[metric].astype(float).values))
            else:
                LOG.warning("GDSC2 file schema unexpected; missing %s", need - set(r.columns))
        except Exception as e:
            LOG.warning("could not read GDSC2 response: %s", str(e)[:150])
    cp = D["proc"] / "ctrp_response_crc.csv"
    if cp.exists():
        r = pd.read_csv(cp)
        out["ctrp"] = pd.DataFrame(dict(cell=r.ccl_name.astype(str).values, cosmic="", drug=r.cpd_name.astype(str).values,
                                        value=r.AUC.astype(float).values))
    lp = D["raw"] / "response" / "response_long.csv"
    if lp.exists():
        r = pd.read_csv(lp)
        r.columns = [x.strip().lower() for x in r.columns]
        if {"cell_line", "drug", "auc"}.issubset(r.columns):
            out["long"] = pd.DataFrame(dict(cell=r.cell_line.astype(str).values,
                                            cosmic=_cos(r["cosmic_id"]) if "cosmic_id" in r else "",
                                            drug=r.drug.astype(str).values, value=r.auc.astype(float).values))
        else:
            LOG.warning("response_long.csv needs columns cell_line, drug, auc [, cosmic_id]; found %s", list(r.columns))
    cr = read_crispr_gpx4(cfg)
    if cr is not None:
        out["crispr"] = cr   # tried last: genetic GPX4 dependency
    return out


def load_cl_expression(cfg):
    """returns (E lines x genes, meta[cosmic,name_norm], source). Order: processed GDSC RMA, DepMap, none."""
    pe = D["proc"] / "gdsc2_expr_crc.csv.gz"
    if pe.exists():
        E = pd.read_csv(pe, index_col=0)
        E.index = E.index.astype(str)
        return E.astype("float32"), pd.DataFrame(dict(cosmic=E.index, name_norm=""), index=E.index), "GDSC RMA (processed)"
    dm = D["raw"] / "depmap"
    dm.mkdir(parents=True, exist_ok=True)
    efs = sorted(p for p in dm.glob("*.csv") if re.search(r"expression", p.name, re.I))
    mfs = sorted(p for p in dm.glob("*.csv") if re.fullmatch(r"model(\.csv|_.*\.csv)", p.name, re.I) or p.name.lower() == "model.csv")
    if not efs or not mfs:
        raise DataUnavailable("no cell-line expression: GDSC RMA file missing and no DepMap files in data/raw/depmap/ "
                              "(need an *Expression*.csv and Model.csv; run --check-data)")
    efs.sort(key=lambda p: (0 if re.search("proteincoding", p.name, re.I) else 1, 0 if re.search("tpm", p.name, re.I) else 1))
    raw = pd.read_csv(efs[0], low_memory=False)
    gcols = [x for x in raw.columns if re.match(r"^.+ \(\d+\)$", str(x))]
    if not gcols:
        raise DataUnavailable(f"{efs[0].name}: no 'GENE (id)' columns found")
    if "IsDefaultEntryForModel" in raw.columns:
        raw = raw[raw["IsDefaultEntryForModel"].astype(str).str.lower().isin(["yes", "true"])]
    idc = "ModelID" if "ModelID" in raw.columns else raw.columns[0]
    E = raw[gcols].astype("float32")
    E.index = raw[idc].astype(str).values
    E.columns = [re.sub(r" \(\d+\)$", "", str(x)).upper() for x in gcols]
    E = E.loc[:, ~E.columns.duplicated()]
    E = E[~E.index.duplicated()]
    M = pd.read_csv(mfs[0], low_memory=False)
    mid = "ModelID" if "ModelID" in M.columns else M.columns[0]
    M = M.drop_duplicates(mid).set_index(mid)
    M.index = M.index.astype(str)
    nm = M["StrippedCellLineName"] if "StrippedCellLineName" in M else M.get("CellLineName", pd.Series("", index=M.index))
    meta = pd.DataFrame(dict(cosmic=_cos(M["COSMICID"]) if "COSMICID" in M else "", name_norm=nm.map(norm).values), index=M.index)
    meta = meta.reindex(E.index).fillna("")
    return E, meta, f"DepMap ({efs[0].name})"


def build_labels(cfg, args):
    E, meta, esrc = load_cl_expression(cfg)
    srcs = read_response_sources(cfg, args.metric)
    if not srcs:
        raise DataUnavailable("no drug-response data: GDSC2_fitted_dose_response.xlsx missing/unreadable, no processed CTRP, "
                              "no data/raw/response/response_long.csv (run --check-data)")
    targets = [norm(t) for t in (args.drugs or GPX4_TARGETS)]
    notes, c2k, n2k = [], {}, {}
    for k, r in meta.iterrows():
        if r.cosmic:
            c2k.setdefault(r.cosmic, k)
        if r.name_norm:
            n2k.setdefault(r.name_norm, k)
    for sname, L in srcs.items():
        if args.response_source not in ("auto", sname):
            continue
        L = L.copy()
        tg = ["gpx4crispr"] if sname == "crispr" else targets
        L["tok"] = L.drug.map(lambda d: next((t for t in tg if t in norm(d)), None))
        present = sorted(L.loc[L.tok.notna(), "drug"].unique())
        L = L[L.tok.notna()]
        if L.empty:
            notes.append(f"{sname}: none of {args.drugs or GPX4_TARGETS} present")
            continue
        L["key"] = [c2k.get(c) or n2k.get(norm(n)) for c, n in zip(L.cosmic, L.cell)]
        n_unmatched = int(L.key.isna().sum())
        L = L.dropna(subset=["key", "value"])
        P = L.groupby(["key", "tok"]).value.mean().unstack("tok")
        keep = [t for t in P.columns if P[t].notna().sum() >= args.min_lines]
        notes.append(f"{sname}: present {present}; lines per drug {P.notna().sum().to_dict()}; unmatched-to-expression rows {n_unmatched}")
        if not keep:
            continue
        pct = P[keep].rank(pct=True)          # lower AUC / IC50 = more sensitive -> low percentile
        comp = pct.mean(axis=1).dropna()
        y = (comp < comp.median()).astype(int)   # 1 = sensitive
        cons = P[keep].corr(method="spearman") if len(keep) > 1 else pd.DataFrame()
        return dict(source=sname, expr_source=esrc, E=E, y=y, comp=comp, drugs=keep, per_drug=P[keep], drug_corr=cons, notes=notes)
    raise DataUnavailable("no usable sensitivity labels. " + " | ".join(notes) +
                          " | Use --drugs to name other compounds, or place a long file (see --check-data).")


# ----------------------------------------------------------------------------------------------------------------
# STAGE 2: cell-line sensitivity classifier + projection to tumours
# ----------------------------------------------------------------------------------------------------------------
def _project(model, genes, coh):
    return {n: pd.Series(model.predict_proba(zc(Xc, genes).values)[:, 1], index=Xc.index) for n, (Xc, _) in coh.items()}


def _proj_logHR(scores, coh):
    out = {}
    for n, s in scores.items():
        v = valid_surv(coh[n][1])
        if len(v) >= 20 and int(v.event.sum()) >= 5:
            x = s.loc[v.index]
            out[n] = cox_hr(v.time, v.event, (x - x.mean()) / (x.std(ddof=0) or 1))["logHR"]
    return out


def cellline_stage(cfg, args):
    panel, pstatus, pver = panel_gate(cfg, args)
    coh = available_cohorts(cfg)
    if COH_T not in coh:
        raise DataUnavailable("TCGA-COAD processed data unavailable")
    B = build_labels(cfg, args)
    E, y = B["E"], B["y"]
    genes = [g for g in panel.gene_symbol if g in E.columns and all(g in v[0].columns for v in coh.values())]
    if len(genes) < 5:
        raise DataUnavailable(f"only {len(genes)} panel genes shared by cell-line expression and all tumour cohorts")
    lines = [k for k in y.index if k in E.index]
    y = y.loc[lines]
    npos, nneg = int(y.sum()), int((1 - y).sum())
    if min(npos, nneg) < 8:
        raise DataUnavailable(f"too few cell lines per class (sensitive={npos}, resistant={nneg}; need >=8 each)")
    seed, N, NP = cfg["seed"], cfg["n_random_classification"], args.n_perm
    Zc = zc(E.loc[lines], genes)
    yv = y.values
    LOG.info("Cell-line classifier: source=%s expression=%s drugs=%s lines=%d sensitive=%d genes=%d",
             B["source"], B["expr_source"], B["drugs"], len(lines), npos, len(genes))
    folds = int(min(5, npos, nneg))
    p = oof_proba(Zc.values, yv, seed, reps=3, folds=folds)
    auc = roc_auc_score(yv, p)
    lo, hi = boot_ci(lambda a, b: roc_auc_score(a, b), len(yv), cfg["n_boot"], seed, yv, p)
    rng = np.random.default_rng(seed + 11)
    perm = [roc_auc_score(yv, oof_proba(Zc.values, rng.permutation(yv), seed + i, reps=1, folds=folds)) for i in range(NP)]
    fs, _ = ferro_score(E.loc[lines], panel, cfg["min_panel_coverage"])
    fs_auc = roc_auc_score(yv, -fs.loc[lines].values)   # sensitive expected at HIGH FerroScore -> orientation shown, not selected
    final = cls_model(seed).fit(Zc.values, yv)
    scores = _project(final, genes, coh)
    obs_h = _proj_logHR(scores, coh)
    # ---- random controls: identical pipeline (labels fixed, gene set random)
    uni = sorted((set(E.columns[E.loc[lines].mean() >= E.loc[lines].mean().median()]) & set(random_universe(cfg, coh, panel))))
    r2 = np.random.default_rng(seed + 13)
    ctrl_auc, ctrl_h = [], []
    for i in range(N):
        gs = list(r2.choice(uni, size=len(genes), replace=False))
        Zr = zc(E.loc[lines], gs)
        ctrl_auc.append(roc_auc_score(yv, oof_proba(Zr.values, yv, seed + 1000 + i, reps=1, folds=folds)))
        mr = cls_model(seed).fit(Zr.values, yv)
        for n, h in _proj_logHR(_project(mr, gs, coh), coh).items():
            ctrl_h.append(dict(set_id=f"random_{i:04d}", cohort=n, abs_logHR=abs(h) if pd.notna(h) else np.nan))
        if (i + 1) % 25 == 0:
            LOG.info("cell-line random controls %d/%d", i + 1, N)
    pd.DataFrame(dict(set_id=[f"random_{i:04d}" for i in range(N)], auc=ctrl_auc)).to_csv(D["csv"] / "cellline_random_auc.csv", index=False)
    ch = pd.DataFrame(ctrl_h)
    ch.to_csv(D["csv"] / "cellline_projection_random_logHR.csv", index=False)
    save_csv(pd.DataFrame([dict(
        label_source=B["source"], expression_source=B["expr_source"], drugs=";".join(B["drugs"]), metric=args.metric,
        label_rule="sensitive = composite within-drug AUC/IC50 percentile below the median across CRC lines (label only)",
        n_lines=len(lines), n_sensitive=npos, n_resistant=nneg, n_genes=len(genes), cv=f"{folds}-fold x3 repeats, cell-line level",
        auc_oof=auc, auc_lo95=lo, auc_hi95=hi, permutation_p=empirical_p(auc, perm), n_perm=NP,
        random_set_median_auc=float(np.median(ctrl_auc)), empirical_p_vs_random_sets=empirical_p(auc, ctrl_auc), n_random=N,
        ferroscore_auc_high_means_sensitive=fs_auc, panel_version=pver, panel_status=pstatus,
        note="in-vitro sensitivity, not a clinical dose or response")]), "cellline_classifier_metrics.csv")
    if len(B["drug_corr"]):
        save_csv(B["drug_corr"].reset_index().rename(columns={"tok": "drug"}), "cellline_drug_label_consistency.csv")
    cf = final.named_steps["logisticregressioncv"].coef_.ravel()
    save_csv(pd.DataFrame(dict(gene=genes, coefficient=cf)).merge(panel[["gene_symbol", "gene_role"]], left_on="gene", right_on="gene_symbol", how="left")
                .drop(columns="gene_symbol").sort_values("coefficient"), "cellline_classifier_coefficients.csv")
    save_csv(pd.DataFrame(dict(cell_line_key=lines, sensitive=yv, composite_percentile=B["comp"].loc[lines].values, oof_prob=p)),
                "cellline_labels_and_oof.csv")
    # ---- projection stats in tumours
    sub = read_csv_opt("ferroptosis_subtype_patients.csv")
    rows, pat = [], []
    for n, (Xc, cc) in coh.items():
        s = scores[n]
        pat.append(pd.DataFrame(dict(cohort=n, sample_id=Xc.index, predicted_gpx4i_sensitivity=s.values)))
        row = dict(cohort=n, n=len(s))
        v = valid_surv(cc)
        if len(v) >= 20 and int(v.event.sum()) >= 5:
            x = s.loc[v.index]
            zx = (x - x.mean()) / (x.std(ddof=0) or 1)
            r = cox_hr(v.time, v.event, zx)
            row.update(endpoint=v.endpoint.iloc[0], n_surv=len(v), events=int(v.event.sum()), HR_per_SD=r["HR"], HR_lo95=r["lo"],
                       HR_hi95=r["hi"], p=r["p"])
            va = v.dropna(subset=["age", "stage"])
            if len(va) >= 30:
                ra = cox_hr(va.time, va.event, zx.loc[va.index], va[["age", "stage"]])
                row.update(HR_per_SD_adj_age_stage=ra["HR"], p_adj=ra["p"], n_adj=len(va))
            cn = ch[ch.cohort == n].abs_logHR.values
            row.update(emp_p_abs_logHR_vs_random=empirical_p(abs(r["logHR"]), cn), random_median_abs_logHR=float(np.nanmedian(cn)) if len(cn) else np.nan)
        mm = cc.mmr.dropna()
        if len(mm) >= 40 and mm.nunique() == 2:
            row.update(auc_for_dMMR=roc_auc_score(mm.astype(int), s.loc[mm.index]),
                       mmr_mwu_p=stats.mannwhitneyu(s.loc[mm.index][mm == 1], s.loc[mm.index][mm == 0]).pvalue)
        if sub is not None:
            g = sub[sub.cohort == n].set_index("sample_id").ferro_subtype
            g = g.reindex(s.index)
            if g.notna().all() and g.nunique() == 2:
                row.update(mean_score_high_subtype=float(s[g == "high"].mean()), mean_score_low_subtype=float(s[g == "low"].mean()),
                           subtype_mwu_p=stats.mannwhitneyu(s[g == "high"], s[g == "low"]).pvalue)
        rows.append(row)
    save_csv(pd.DataFrame(rows), "cellline_projection_tumour_stats.csv")
    save_csv(pd.concat(pat), "cellline_projection_tumour_scores.csv")
    summ("ext_cellline_sensitivity_classifier", "completed", n=len(lines),
            metric=f"{B['source']}/{','.join(B['drugs'])}: OOF AUC={auc:.3f} [{lo:.3f},{hi:.3f}], perm p={empirical_p(auc, perm):.3f}, "
                   f"emp p vs {N} random sets={empirical_p(auc, ctrl_auc):.3f}",
            note=f"projection onto {len(coh)} tumour cohorts; panel={pstatus}; in-vitro sensitivity only")


# ----------------------------------------------------------------------------------------------------------------
# FIGURES
# ----------------------------------------------------------------------------------------------------------------
def fig_overview():
    fig, ax = plt.subplots(figsize=(11, 6.2))
    ax.axis("off")
    ax.set_xlim(0, 11)
    ax.set_ylim(0, 6.2)

    def box(x, y, w, h, t, fc="#eef3fb", ec="#4a6fa5"):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.04", fc=fc, ec=ec, lw=1))
        ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=7.5)

    def arr(x1, y1, x2, y2):
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="->", lw=1, color="#444"))
    box(0.1, 4.6, 2.4, 1.3, "DATA\nTCGA-COAD (RNA-seq, GDC)\nGSE39582, GSE17536, GSE14333\nGDSC2 / CTRP / DepMap cell lines", "#fff4e0", "#c98a1b")
    box(3.0, 4.6, 2.2, 1.3, "PREPROCESS\nprimary tumour selection\nprobe->gene collapse\noverlap removal, audit tables")
    box(5.7, 4.6, 2.3, 1.3, "FERROPTOSIS PANEL\n(driver / suppressor genes)\nversion hash, CONFIRMED\nor PROVISIONAL label", "#e9f7ef", "#2e8b57")
    box(8.5, 4.6, 2.3, 1.3, "RANDOM GENE-SET\nBENCHMARK\nidentical procedure\nfixed seed, empirical p", "#fdecea", "#c0392b")
    box(0.1, 2.5, 3.3, 1.3, "A. PROGNOSIS\nuni-Cox screen -> ridge Cox (TCGA)\nfrozen model -> 3 external cohorts\nclinical / molecular / combined")
    box(3.7, 2.5, 3.3, 1.3, "B1. FERROPTOSIS SUBTYPES\nk-means (k=2) on TCGA, frozen centroids\nproject cohort-wise\nsurvival / stage / age / MMR")
    box(7.3, 2.5, 3.5, 1.3, "B2. CELL-LINE SENSITIVITY\nCRC lines: GPX4-inhibitor sensitive vs resistant\nL2 logistic on panel (OOF AUC)\nfrozen -> tumours")
    box(0.1, 0.5, 4.8, 1.3, "C. DRUG RESPONSE (cell lines)\nFerroScore vs LN_IC50 / AUC (BH q)\nprespecified: ML162, ML210, RSL3 + 5-FU, oxaliplatin, irinotecan\nregression models; CTRP replication", "#f4ecf7", "#7d3c98")
    box(5.3, 0.5, 5.5, 1.3, "INTEGRATION + LIMITATIONS\nprojected sensitivity vs subtype, survival, MMR\nin-vitro sensitivity != clinical dose\nno TCGA MSI/DFS; heuristic GEO parser; overlap checks", "#fff4e0", "#c98a1b")
    for x1, x2 in [(2.5, 3.0), (5.2, 5.7), (8.0, 8.5)]:
        arr(x1, 5.25, x2, 5.25)
    for x in (1.75, 5.35, 9.05):
        arr(x, 4.6, x, 3.8)
    arr(1.75, 2.5, 1.75, 1.8)
    arr(5.35, 2.5, 6.5, 1.8)
    arr(9.05, 2.5, 8.5, 1.8)
    ax.set_title("Ferroptosis classification, prognosis and cell-line sensitivity pipeline (colorectal cancer)", fontsize=10)
    savefig(fig, "fig15_pipeline_overview")


def fig_subtype_km(cfg):
    pat, sc = read_csv_opt("ferroptosis_subtype_patients.csv"), read_csv_opt("prognosis_patient_scores.csv")
    coh = available_cohorts(cfg)
    if pat is None:
        raise DataUnavailable("no subtype table")
    names = [n for n in coh if len(valid_surv(coh[n][1])) >= 20]
    fig, axs = plt.subplots(1, len(names), figsize=(4.2 * len(names), 4), squeeze=False)
    for ax, n in zip(axs[0], names):
        v = valid_surv(coh[n][1])
        g = pat[pat.cohort == n].set_index("sample_id").ferro_subtype.reindex(v.index)
        for grp, col in (("low", "#2a7fba"), ("high", "#d62728")):
            m = g == grp
            if m.sum():
                KaplanMeierFitter().fit(v.time[m], v.event[m], label=f"{grp} (n={int(m.sum())})").plot_survival_function(ax=ax, color=col, ci_show=False)
        h, l = v[g == "high"], v[g == "low"]
        p = logrank_test(h.time, l.time, h.event, l.event).p_value if len(h) > 4 and len(l) > 4 else np.nan
        ax.set_title(f"{n} ({v.endpoint.iloc[0]})\nlog-rank p={p:.3g}", fontsize=8)
        ax.set_xlabel("months")
    fig.suptitle("Ferroptosis subtypes (TCGA-frozen k=2 centroids)")
    savefig(fig, "fig16_subtype_km")


def fig_cellline():
    m = read_csv_opt("cellline_classifier_metrics.csv")
    ca = read_csv_opt("cellline_random_auc.csv")
    ch = read_csv_opt("cellline_projection_random_logHR.csv")
    st = read_csv_opt("cellline_projection_tumour_stats.csv")
    if m is None or ca is None:
        raise DataUnavailable("cell-line stage not run")
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    axs[0].hist(ca.auc, bins=25, color="#bbb", edgecolor="k")
    axs[0].axvline(m.auc_oof.iloc[0], c="crimson", lw=2, label=f"panel AUC {m.auc_oof.iloc[0]:.3f}\nemp. p={m.empirical_p_vs_random_sets.iloc[0]:.3f}")
    axs[0].set_xlabel("out-of-fold AUC (cell lines)"); axs[0].legend(fontsize=7)
    axs[0].set_title(f"GPX4-inhibitor sensitivity: panel vs {len(ca)} random sets\n({m.drugs.iloc[0]}; n={int(m.n_lines.iloc[0])} lines)", fontsize=8)
    if ch is not None and st is not None and len(ch):
        cohs = [c for c in st.cohort if c in set(ch.cohort)]
        axs[1].boxplot([ch[ch.cohort == c].abs_logHR.dropna().values for c in cohs], labels=cohs, showfliers=False)
        for i, c in enumerate(cohs, 1):
            r = st[st.cohort == c]
            if len(r) and pd.notna(r.HR_per_SD.iloc[0]):
                axs[1].scatter([i], [abs(np.log(r.HR_per_SD.iloc[0]))], c="crimson", zorder=3, s=30)
        axs[1].set_ylabel("|log HR per SD| of projected score"); axs[1].tick_params(axis="x", rotation=15)
        axs[1].set_title("Projected score vs survival (red = panel; box = random sets)", fontsize=8)
    savefig(fig, "fig17_cellline_classifier_projection")


def ext_figures_stage(cfg, args):
    for nm, fn in [("overview", fig_overview), ("subtype_km", lambda: fig_subtype_km(cfg)), ("cellline", fig_cellline)]:
        try:
            fn()
        except Exception as e:
            plt.close("all")
            LOG.warning("figure '%s' skipped: %s", nm, str(e)[:160])


# ----------------------------------------------------------------------------------------------------------------
# REPORT
# ----------------------------------------------------------------------------------------------------------------
def ext_report_stage(cfg, args):
    g = read_csv_opt
    L = [f"# Extension results (v{EXT_VERSION}; {time.strftime('%Y-%m-%d %H:%M')})", "",
         "> Cell-line sensitivity is in-vitro sensitivity, not a clinical dose or response. Associations are not mechanism.", "",
         "## 1. Ferroptosis subtypes: associations", md_table(g("ferroptosis_subtype_associations.csv"), 10), "",
         "Stability:", md_table(g("subtype_stability.csv")), "",
         "## 2. Subtypes vs random gene sets (identical procedure)", md_table(g("subtype_vs_random.csv")), "",
         "## 3. Cell-line GPX4-inhibitor sensitivity classifier", md_table(g("cellline_classifier_metrics.csv")), "",
         "Drug-label consistency (Spearman):", md_table(g("cellline_drug_label_consistency.csv")), "",
         "Coefficients:", md_table(g("cellline_classifier_coefficients.csv"), 40), "",
         "## 4. Projection of the frozen cell-line classifier onto tumours", md_table(g("cellline_projection_tumour_stats.csv")), "",
         "## 5. How to read this",
         "- An empirical p near or above 0.05 against random gene sets means the panel is not distinguishable from arbitrary genes for that task.",
         "- Projection assumes cohort-wise z-scoring makes cell-line and tumour scales comparable; tumour stroma/purity is a confounder.",
         "- Subtype and classifier definitions are frozen before any outcome of an external cohort is examined.", ""]
    s = read_state()
    L += ["## 6. Run summary"] + [f"- {k}: {v['status']} - {str(v['key_result'])[:200]}" for k, v in s["summary"].items() if k.startswith("ext_")]
    (D["rep"] / "extension_report.md").write_text("\n".join(L), encoding="utf-8")
    LOG.info("wrote %s", D["rep"] / "extension_report.md")




# ----------------------------------------------------------------------------------------
# v1.4.1: DepMap auto-download, GPX4 CRISPR-dependency label source, drug lister
# ----------------------------------------------------------------------------------------
DEPMAP_API = "https://depmap.org/portal/api/download/files"
DEPMAP_REQUIRED = ["Model.csv", "OmicsExpressionProteinCodingGenesTPMLogp1.csv"]
DEPMAP_OPTIONAL = ["CRISPRGeneEffect.csv"]


def _pick_depmap_url(listing, filename):
    """From the DepMap file listing (DataFrame, any column case) pick (url, release) of the newest release containing filename."""
    d = listing.copy()
    d.columns = [str(c).strip().lower() for c in d.columns]
    if "filename" not in d.columns or "url" not in d.columns:
        raise DataUnavailable(f"DepMap listing has unexpected columns: {list(d.columns)[:10]}")
    d = d[d["filename"].astype(str).str.strip() == filename].dropna(subset=["url"])
    if d.empty:
        return None, None
    rel = d["release"].astype(str) if "release" in d.columns else pd.Series("", index=d.index)
    d = d.assign(_rel=rel, _public=rel.str.lower().str.contains("depmap public").astype(int))
    if "release_date" in d.columns:
        d["_date"] = pd.to_datetime(d["release_date"], errors="coerce")
    else:
        d["_date"] = pd.NaT
    d = d.sort_values(["_public", "_date", "_rel"], ascending=[False, False, False])
    r = d.iloc[0]
    return str(r["url"]), str(r["_rel"])


def fetch_depmap(cfg, refresh):
    """Download Model.csv + expression (+ optional CRISPR gene effect) from DepMap into data/raw/depmap/."""
    dm = D["raw"] / "depmap"
    dm.mkdir(parents=True, exist_ok=True)
    if (D["proc"] / "gdsc2_expr_crc.csv.gz").exists():
        LOG.info("DepMap: GDSC expression already processed; DepMap expression not needed (CRISPR file still optional)")
    have_req = all((dm / f).exists() for f in DEPMAP_REQUIRED)
    have_opt = all((dm / f).exists() for f in DEPMAP_OPTIONAL)
    if have_req and have_opt and not refresh:
        for f in DEPMAP_REQUIRED + DEPMAP_OPTIONAL:
            record("DepMap", f, "cache", dm / f, "cached", "processed")
        return
    try:
        r = http("GET", DEPMAP_API, timeout=(20, 60))
        try:
            listing = pd.read_csv(io.StringIO(r.text))
            if "filename" not in [str(c).lower() for c in listing.columns]:
                raise ValueError("not a csv listing")
        except Exception:
            listing = pd.DataFrame(r.json())
    except Exception as e:
        raise DataUnavailable(
            f"DepMap file listing unavailable ({str(e)[:150]}). Download by hand from https://depmap.org/portal/data_page/ "
            f"-> files {DEPMAP_REQUIRED + DEPMAP_OPTIONAL} -> {dm} (keep the file names).")
    errs = []
    for f in DEPMAP_REQUIRED + DEPMAP_OPTIONAL:
        if (dm / f).exists() and not refresh:
            continue
        try:
            url, rel = _pick_depmap_url(listing, f)
            if not url:
                raise DataUnavailable(f"{f} not in the DepMap listing")
            LOG.info("DepMap: downloading %s from release '%s' (large files may take several minutes)", f, rel)
            download_file([url], dm / f, "DepMap", f, "processed", "text", refresh=refresh, version=rel or "")
        except Exception as e:
            errs.append(f"{f}: {str(e)[:200]}")
            LOG.warning("DepMap %s failed: %s", f, str(e)[:200])
    missing_req = [f for f in DEPMAP_REQUIRED if not (dm / f).exists()]
    if missing_req:
        raise DataUnavailable(f"DepMap required files missing {missing_req}; place them in {dm} manually. Errors: {' | '.join(errs)}")


def read_crispr_gpx4(cfg):
    """GPX4 CRISPR gene effect (DepMap) in colorectal ('Bowel') lines as a response table (lower = more GPX4-dependent)."""
    dm = D["raw"] / "depmap"
    f, mf = dm / "CRISPRGeneEffect.csv", dm / "Model.csv"
    if not (f.exists() and mf.exists()):
        return None
    try:
        hdr = list(pd.read_csv(f, nrows=0).columns)
        col = next((c for c in hdr if str(c).upper().startswith("GPX4 (")), None)
        if col is None:
            LOG.warning("CRISPRGeneEffect.csv has no GPX4 column; GPX4-dependency labels unavailable")
            return None
        g = pd.read_csv(f, usecols=[hdr[0], col])
        g.columns = ["ModelID", "value"]
        M = pd.read_csv(mf, low_memory=False)
        mid = "ModelID" if "ModelID" in M.columns else M.columns[0]
        if "OncotreeLineage" not in M.columns:
            LOG.warning("Model.csv lacks OncotreeLineage; cannot restrict GPX4 dependency to colorectal lines")
            return None
        M = M.drop_duplicates(mid).set_index(mid)
        M.index = M.index.astype(str)
        crc = M.index[M["OncotreeLineage"].astype(str).str.lower() == "bowel"]
        g["ModelID"] = g["ModelID"].astype(str)
        g = g[g.ModelID.isin(crc)].dropna()
        if g.empty:
            return None
        nm = M["StrippedCellLineName"] if "StrippedCellLineName" in M.columns else M.get("CellLineName", pd.Series("", index=M.index))
        cos = pd.Series(_cos(M["COSMICID"]) if "COSMICID" in M.columns else "", index=M.index)
        return pd.DataFrame(dict(cell=nm.reindex(g.ModelID).astype(str).values, cosmic=cos.reindex(g.ModelID).fillna("").values,
                                 drug="GPX4_CRISPR", value=g.value.astype(float).values))
    except Exception as e:
        LOG.warning("could not read GPX4 CRISPR dependency: %s", str(e)[:150])
        return None


def list_drugs(cfg):
    gp = D["raw"] / "gdsc" / "GDSC2_fitted_dose_response.xlsx"
    adopt_manual(D["raw"] / "gdsc", ["*GDSC2*dose*.xlsx", "*fitted*dose*.xlsx"], gp)
    if not gp.exists():
        print(f"\nGDSC2 file not found: {gp}")
        return
    r = pd.read_excel(gp)
    fld = cfg["gdsc_lineage_field"]
    if fld not in r.columns or "DRUG_NAME" not in r.columns:
        print(f"unexpected GDSC2 columns: {list(r.columns)[:15]}")
        return
    crc = r[r[fld].astype(str) == cfg["gdsc_lineage_value"]]
    print(f"\nGDSC2: {r.DRUG_NAME.nunique()} drugs overall; {crc.CELL_LINE_NAME.nunique()} CRC lines; {crc.DRUG_NAME.nunique()} drugs tested in CRC lines")
    toks = ["ml162", "ml210", "rsl3", "erastin", "sulfasalazine", "sorafenib", "gpx4", "ferrostatin", "fluorouracil", "oxaliplatin", "irinotecan"]
    cnt = crc.groupby("DRUG_NAME").CELL_LINE_NAME.nunique()
    hits = [(n, int(c)) for n, c in cnt.items() if any(t in norm(n) for t in toks)]
    print("Ferroptosis-related / clinical compounds present in CRC lines (name, n lines):")
    for n, c in sorted(hits):
        print(f"   {n:30} {c}")
    if not any(any(t in norm(n) for t in ("ml162", "ml210", "rsl3")) for n, _ in hits):
        print("\n-> ML162 / ML210 / RSL3 are NOT in this GDSC2 file. The cell-line stage will use GPX4 CRISPR dependency "
              "(DepMap) instead, or a long-format response file (see --check-data).")


def check_data(cfg):
    R, P = D["raw"], D["proc"]
    items = [("GDSC2 response (xlsx)", R / "gdsc" / "GDSC2_fitted_dose_response.xlsx"),
             ("GDSC RMA expression (zip)", R / "gdsc" / "Cell_line_RMA_proc_basalExp.txt.zip"),
             ("GDSC processed expression", P / "gdsc2_expr_crc.csv.gz"), ("CTRP processed", P / "ctrp_response_crc.csv"),
             ("long response file", R / "response" / "response_long.csv")]
    print("\nCell-line data present?")
    for n, p in items:
        print(f"  {'YES' if p.exists() else 'no ':3}  {n:32} {p}")
    dm = R / "depmap"
    print(f"  {'YES' if list(dm.glob('*xpression*.csv')) else 'no ':3}  DepMap expression csv           {dm}\\*Expression*.csv")
    print(f"  {'YES' if list(dm.glob('Model*.csv')) else 'no ':3}  DepMap Model.csv                {dm}\\Model.csv")
    print("""
Minimum needed for the cellline stage = (A) a drug-response source containing ML162/ML210/RSL3 AND (B) cell-line expression.
  A) GDSC2_fitted_dose_response.xlsx  (your fetch stage appears to have downloaded it already) - NOTE: check that
     the GPX4 inhibitors are really in it (the stage prints which drugs it found). CTRP contains them; if the CTRP
     server is unreachable, download from the CTRP/DepMap portals manually.
  B) expression: either the GDSC RMA zip -> data\\raw\\gdsc\\, or DepMap (free account-less download, https://depmap.org/portal/data_page/):
        OmicsExpressionProteinCodingGenesTPMLogp1.csv  and  Model.csv  -> data\\raw\\depmap\\
     DepMap lines are matched to GDSC lines by COSMIC ID, else by normalised name.
  Any other response data: data\\raw\\response\\response_long.csv with columns cell_line,drug,auc[,cosmic_id]
     (auc: lower = more sensitive; CRC lines only).
File layouts of DepMap releases change; this loader is written for the 'GENE (id)' column format and is untested against
your files - if it errors, send me the first 3 lines of the file.
""")





# ----------------------------------------------------------------------------------------
# v1.5.0: positive controls, IC50/dose prediction, scorecard
# ----------------------------------------------------------------------------------------
V15_VERSION = "1.5.0"
# Literature-based reference signatures (fixed BEFORE looking at any result; edit only with a written reason).
# Mesenchymal / EMT state -> GPX4 dependence: Viswanathan et al., Nature 2017; Hangauer et al., Nature 2017.
SIG_MESENCHYMAL = ["VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "FN1", "CDH2", "TGFBI", "CALD1", "COL1A1", "COL1A2", "SPARC", "S100A4"]
SIG_EPITHELIAL = ["CDH1", "EPCAM", "KRT8", "KRT18", "KRT19", "CLDN4", "CLDN7", "ESRP1", "GRHL2", "OCLN"]
SIG_PROLIF = ["MKI67", "PCNA", "TOP2A", "MCM2", "MCM6", "CCNB1", "CDK1", "AURKA", "BUB1", "TYMS", "CCNA2", "CDC20"]
DOSE_FERRO_TOKENS = ["ml162", "ml210", "rsl3", "erastin", "fin56", "gpx4"]
DOSE_CHEMO_TOKENS = ["fluorouracil", "5fu", "oxaliplatin", "irinotecan", "sn38", "trifluridine", "capecitabine"]


def _sig_present(E, genes):
    return [g for g in genes if g in E.columns]


def emt_score(E):
    """mean z(mesenchymal) - mean z(epithelial); z-scored within the supplied matrix (label-free)."""
    m, e = _sig_present(E, SIG_MESENCHYMAL), _sig_present(E, SIG_EPITHELIAL)
    if len(m) < 5 or len(e) < 4:
        raise DataUnavailable(f"EMT signature coverage too low (mesenchymal {len(m)}/{len(SIG_MESENCHYMAL)}, epithelial {len(e)}/{len(SIG_EPITHELIAL)})")
    return zc(E, m).mean(axis=1) - zc(E, e).mean(axis=1)


def prolif_score(E):
    g = _sig_present(E, SIG_PROLIF)
    if len(g) < 5:
        raise DataUnavailable(f"proliferation signature coverage too low ({len(g)}/{len(SIG_PROLIF)})")
    return zc(E, g).mean(axis=1)


def nested_topk_auc(X, y, k, seed, folds, reps=1):
    """data-driven ceiling: choose the k most discriminative genes INSIDE each training fold (no leakage), then logistic regression."""
    from sklearn.feature_selection import f_classif
    p = np.zeros(len(y))
    for r in range(reps):
        for tr, te in StratifiedKFold(folds, shuffle=True, random_state=seed + r).split(X, y):
            F, _ = f_classif(X[tr], y[tr])
            idx = np.argsort(-np.nan_to_num(F))[:k]
            m = cls_model(seed).fit(X[tr][:, idx], y[tr])
            p[te] += m.predict_proba(X[te][:, idx])[:, 1] / reps
    return roc_auc_score(y, p)


def poscontrol_stage(cfg, args):
    panel, pstatus, pver = panel_gate(cfg, args)
    coh = available_cohorts(cfg)
    if COH_T not in coh:
        raise DataUnavailable("TCGA-COAD processed data unavailable")
    seed, N = cfg["seed"], cfg["n_random_classification"]
    rows_cl, rows_tu = [], []
    notes = []
    # ---------------- (A) cell lines: does a known biology signal show up with this machinery?
    cl_ok = True
    try:
        B = build_labels(cfg, args)
    except DataUnavailable as e:
        cl_ok = False
        notes.append(f"cell-line positive control skipped: {str(e)[:160]}")
    if cl_ok:
        E, y = B["E"], B["y"]
        lines = [k for k in y.index if k in E.index]
        yv = y.loc[lines].values
        npos, nneg = int(yv.sum()), int((1 - yv).sum())
        if min(npos, nneg) < 8:
            raise DataUnavailable(f"too few cell lines per class (sensitive={npos}, resistant={nneg})")
        folds = int(min(5, npos, nneg))
        El = E.loc[lines]
        genes_panel = [g for g in panel.gene_symbol if g in E.columns]
        uni_all = sorted(set(El.columns[El.mean() >= El.mean().median()]) - set(panel.gene_symbol))
        sigs = {"panel (ferroptosis)": genes_panel,
                "EMT (mesenchymal+epithelial)": _sig_present(El, SIG_MESENCHYMAL + SIG_EPITHELIAL),
                "proliferation": _sig_present(El, SIG_PROLIF)}
        r2 = np.random.default_rng(seed + 31)
        for name, gl in sigs.items():
            if len(gl) < 5:
                notes.append(f"{name}: only {len(gl)} genes present in cell-line expression; skipped")
                continue
            Z = zc(El, gl).values
            auc = roc_auc_score(yv, oof_proba(Z, yv, seed, reps=3, folds=folds))
            lo, hi = boot_ci(lambda a_, b_: roc_auc_score(a_, b_), len(yv), cfg["n_boot"], seed, yv,
                             oof_proba(Z, yv, seed, reps=3, folds=folds))
            ctrl = []
            for i in range(N):
                gs = list(r2.choice(uni_all, size=len(gl), replace=False))
                ctrl.append(roc_auc_score(yv, oof_proba(zc(El, gs).values, yv, seed + 2000 + i, reps=1, folds=folds)))
            row = dict(target=f"{B['source']}/{','.join(B['drugs'])}", signature=name, n_genes=len(gl), n_lines=len(lines), auc_oof=auc,
                       auc_lo95=lo, auc_hi95=hi, random_median_auc=float(np.median(ctrl)), empirical_p_vs_random_sets=empirical_p(auc, ctrl),
                       n_random=N)
            if name == "EMT (mesenchymal+epithelial)":
                s = emt_score(El)
                row["single_score_auc_high_means_sensitive"] = roc_auc_score(yv, s.values)
            rows_cl.append(row)
            LOG.info("positive control (cell lines) %s: AUC=%.3f random median=%.3f emp p=%.3f", name, auc, np.median(ctrl), row["empirical_p_vs_random_sets"])
        k = min(len(genes_panel), 36)
        Xall = zc(El, [g for g in El.columns if El[g].std() > 1e-6]).values
        ceil = nested_topk_auc(Xall, yv, k, seed, folds, reps=3)
        rows_cl.append(dict(target=f"{B['source']}/{','.join(B['drugs'])}", signature=f"data-driven ceiling (top {k} genes chosen inside CV folds)",
                            n_genes=k, n_lines=len(lines), auc_oof=ceil))
    # ---------------- (B) tumours: do clinical stage / proliferation / EMT show up in survival?
    r3 = np.random.default_rng(seed + 37)
    uni = random_universe(cfg, coh, panel)
    for n, (Xc, cc) in coh.items():
        v = valid_surv(cc)
        if len(v) < 20 or int(v.event.sum()) < 5:
            continue
        base = dict(cohort=n, endpoint=v.endpoint.iloc[0], n_surv=len(v), events=int(v.event.sum()))
        vs = v.dropna(subset=["stage"])
        if len(vs) >= 30 and vs.stage.nunique() > 1:
            r = cox_hr(vs.time, vs.event, (vs.stage - vs.stage.mean()) / (vs.stage.std(ddof=0) or 1))
            rows_tu.append(dict(base, signature="clinical stage (per SD)", n_genes=0, HR_per_SD=r["HR"], HR_lo95=r["lo"], HR_hi95=r["hi"], p=r["p"]))
        for name, fn, gl in [("proliferation", prolif_score, SIG_PROLIF), ("EMT (mesenchymal-epithelial)", emt_score, SIG_MESENCHYMAL + SIG_EPITHELIAL)]:
            try:
                s = fn(Xc)
            except DataUnavailable as e:
                notes.append(f"{n} {name}: {str(e)[:100]}")
                continue
            x = s.loc[v.index]
            r = cox_hr(v.time, v.event, (x - x.mean()) / (x.std(ddof=0) or 1))
            ng = len(_sig_present(Xc, gl))
            ctrl = []
            for i in range(N):
                gs = list(r3.choice(uni, size=ng, replace=False))
                xr = zc(Xc, gs).mean(axis=1).loc[v.index]
                h = cox_hr(v.time, v.event, (xr - xr.mean()) / (xr.std(ddof=0) or 1))["logHR"]
                ctrl.append(abs(h) if pd.notna(h) else np.nan)
            rows_tu.append(dict(base, signature=name, n_genes=ng, HR_per_SD=r["HR"], HR_lo95=r["lo"], HR_hi95=r["hi"], p=r["p"],
                                emp_p_abs_logHR_vs_random=empirical_p(abs(r["logHR"]), ctrl), random_median_abs_logHR=float(np.nanmedian(ctrl))))
    if rows_cl:
        save_csv(pd.DataFrame(rows_cl), "positive_control_cellline.csv")
    if rows_tu:
        save_csv(pd.DataFrame(rows_tu), "positive_control_tumour.csv")
    if not rows_cl and not rows_tu:
        raise DataUnavailable("no positive-control result could be computed: " + "; ".join(notes))
    txt = []
    for r in rows_cl:
        if "empirical_p_vs_random_sets" in r and pd.notna(r.get("empirical_p_vs_random_sets", np.nan)):
            txt.append(f"{r['signature']}: AUC={r['auc_oof']:.3f} (emp p vs random={r['empirical_p_vs_random_sets']:.3f})")
        elif r["signature"].startswith("data-driven"):
            txt.append(f"ceiling AUC={r['auc_oof']:.3f}")
    st = [r for r in rows_tu if r["signature"].startswith("clinical stage")]
    if st:
        txt.append("stage HR/SD: " + ", ".join(f"{r['cohort']}={r['HR_per_SD']:.2f}" for r in st))
    summ("ext_positive_control", "completed", n=N, metric=" | ".join(txt)[:400], note="; ".join(notes)[:200])


# ---------------- IC50 / dose-response prediction
def _match_keys(L, meta):
    c2k, n2k = {}, {}
    for k, r in meta.iterrows():
        if r.cosmic:
            c2k.setdefault(r.cosmic, k)
        if r.name_norm:
            n2k.setdefault(r.name_norm, k)
    return [c2k.get(c) or n2k.get(norm(n)) for c, n in zip(L.cosmic, L.cell)]


def ridge_oof(X, y, seed, folds=5, reps=1):
    from sklearn.linear_model import RidgeCV
    p = np.zeros(len(y))
    for r in range(reps):
        for tr, te in KFold(folds, shuffle=True, random_state=seed + r).split(X):
            m = RidgeCV(alphas=np.logspace(0, 4, 9)).fit(X[tr], y[tr])
            p[te] += m.predict(X[te]) / reps
    return p


def dose_stage(cfg, args):
    panel, pstatus, pver = panel_gate(cfg, args)
    coh = available_cohorts(cfg)
    if COH_T not in coh:
        raise DataUnavailable("TCGA-COAD processed data unavailable")
    E, meta, esrc = load_cl_expression(cfg)
    srcs = read_response_sources(cfg, "LN_IC50")      # IC50 is the target here; CTRP / long files supply AUC
    if not srcs:
        raise DataUnavailable("no drug-response data (GDSC2 xlsx / CTRP / response_long.csv); run --check-data")
    tokens = [norm(t) for t in (args.dose_drugs or (DOSE_FERRO_TOKENS + DOSE_CHEMO_TOKENS))]
    genes = [g for g in panel.gene_symbol if g in E.columns and all(g in v[0].columns for v in coh.values())]
    if len(genes) < 5:
        raise DataUnavailable(f"only {len(genes)} panel genes shared by cell-line expression and all tumour cohorts")
    seed, N = cfg["seed"], cfg["n_random_classification"]
    avail, res, tum, ctrl_rows = [], [], [], []
    uni_cache = {}
    for sname, L in srcs.items():
        if sname == "crispr":
            L = L.assign(tok="gpx4crispr")
            unit = "gene effect (lower = more GPX4-dependent)"
        else:
            L = L.copy()
            L["tok"] = L.drug.map(lambda d: next((t for t in tokens if t in norm(d)), None))
            unit = "ln IC50 (uM)" if sname == "gdsc" else ("AUC" if sname in ("ctrp", "long") else "value")
        for t in tokens + (["gpx4crispr"] if sname == "crispr" else []):
            sub = L[L.tok == t] if sname != "crispr" else L
            if sname == "crispr" and t != "gpx4crispr":
                continue
            names = sorted(sub.drug.unique()) if len(sub) else []
            avail.append(dict(source=sname, token=t, compounds_found=";".join(names)[:120], n_rows=len(sub)))
        L = L[L.tok.notna()]
        if L.empty:
            continue
        L["key"] = _match_keys(L, meta)
        L = L.dropna(subset=["key", "value"])
        for t, G in L.groupby("tok"):
            P = G.groupby("key").value.mean()
            P = P[[k in E.index for k in P.index]]
            if len(P) < args.min_lines:
                avail.append(dict(source=sname, token=t, compounds_found="", n_rows=len(P), note=f"<{args.min_lines} matched lines; not modelled"))
                continue
            lines = list(P.index)
            yv = P.values.astype(float)
            Zc = zc(E.loc[lines], genes).values
            folds = int(min(5, max(3, len(lines) // 8)))
            p = ridge_oof(Zc, yv, seed, folds=folds, reps=3)
            rho, pv = stats.spearmanr(yv, p)
            lo, hi = boot_ci(lambda a_, b_: stats.spearmanr(a_, b_)[0], len(yv), cfg["n_boot"], seed, yv, p)
            El = E.loc[lines]
            if t not in uni_cache:
                uni_cache[t] = sorted(set(El.columns[El.mean() >= El.mean().median()]) - set(panel.gene_symbol))
            r2 = np.random.default_rng(seed + 41)
            ctrl = []
            for i in range(N):
                gs = list(r2.choice(uni_cache[t], size=len(genes), replace=False))
                pr = ridge_oof(zc(El, gs).values, yv, seed + 3000 + i, folds=folds, reps=1)
                ctrl.append(stats.spearmanr(yv, pr)[0])
            res.append(dict(source=sname, drug=t, compounds=";".join(sorted(G.drug.unique()))[:100], response_unit=unit, n_lines=len(lines),
                            n_genes=len(genes), spearman_oof=rho, rho_lo95=lo, rho_hi95=hi,
                            random_median_rho=float(np.nanmedian(ctrl)), empirical_p_vs_random_sets=empirical_p(rho, ctrl), n_random=N,
                            panel_version=pver, panel_status=pstatus,
                            note="in-vitro response of CRC cell lines; judge ONLY by empirical_p_vs_random_sets (out-of-fold rho is biased low for pure noise)"))
            LOG.info("dose %s/%s: n=%d rho=%.3f random median=%.3f emp p=%.3f", sname, t, len(lines), rho, np.nanmedian(ctrl), res[-1]["empirical_p_vs_random_sets"])
            fin = __import__("sklearn.linear_model", fromlist=["RidgeCV"]).RidgeCV(alphas=np.logspace(0, 4, 9)).fit(Zc, yv)
            for n, (Xc, cc) in coh.items():
                pred = fin.predict(zc(Xc, genes).values)
                d = pd.DataFrame(dict(cohort=n, sample_id=Xc.index, source=sname, drug=t, predicted_response=pred))
                if sname == "gdsc":
                    d["predicted_IC50_uM"] = np.exp(pred)
                tum.append(d)
    if avail:
        save_csv(pd.DataFrame(avail), "dose_drug_availability.csv")
    if not res:
        raise DataUnavailable("no drug had enough matched cell lines; see results/csv/dose_drug_availability.csv and use --dose-drugs / --list-drugs "
                              "(GDSC2 lacks GPX4 inhibitors; CTRP/PRISM files can be placed per --check-data)")
    R = pd.DataFrame(res)
    R["q_bh"] = bh_adjust(R.empirical_p_vs_random_sets.values)
    save_csv(R, "dose_prediction_results.csv")
    save_csv(pd.concat(tum), "dose_prediction_tumour_scores.csv")
    best = R.sort_values("empirical_p_vs_random_sets").iloc[0]
    summ("ext_dose_prediction", "completed", n=int(R.n_lines.sum()),
         metric=f"{len(R)} drug models; {int((R.empirical_p_vs_random_sets < 0.05).sum())} beat random sets (emp p<0.05); best {best.source}/{best.drug}: "
                f"rho={best.spearman_oof:.2f}, emp p={best.empirical_p_vs_random_sets:.3f}",
         note=f"ridge on {len(genes)} panel genes; panel={pstatus}; in-vitro only")


# ---------------- scorecard
def _verdict(p, lo=0.05):
    if p is None or pd.isna(p):
        return "n/a"
    return "beats random gene sets" if p < lo else "NOT distinguishable from random gene sets"


def fig_dose_poscontrol():
    d = read_csv_opt("dose_prediction_results.csv")
    pc = read_csv_opt("positive_control_cellline.csv")
    if d is None and pc is None:
        raise DataUnavailable("dose / positive-control stages not run")
    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    if d is not None and len(d):
        d = d.sort_values("spearman_oof")
        y = np.arange(len(d))
        axs[0].barh(y - 0.18, d.spearman_oof, 0.36, color="crimson", label="panel (out-of-fold rho)")
        axs[0].barh(y + 0.18, d.random_median_rho, 0.36, color="#bbb", label="random sets (median)")
        axs[0].set_yticks(y)
        axs[0].set_yticklabels([f"{r.drug} [{r.source}] p={r.empirical_p_vs_random_sets:.2f}" for r in d.itertuples()], fontsize=7)
        axs[0].axvline(0, c="k", lw=0.5)
        axs[0].set_xlabel("Spearman rho, predicted vs measured response"); axs[0].legend(fontsize=7)
        axs[0].set_title("IC50/AUC prediction from the panel (n cell lines per drug in table)", fontsize=8)
    if pc is not None and len(pc):
        q = pc.dropna(subset=["auc_oof"])
        axs[1].barh(np.arange(len(q)), q.auc_oof, color=["crimson" if s.startswith("panel") else "#4c78a8" for s in q.signature])
        axs[1].scatter(q.random_median_auc, np.arange(len(q)), c="k", marker="|", s=200, label="random-set median", zorder=3)
        axs[1].set_yticks(np.arange(len(q)))
        axs[1].set_yticklabels([s[:38] for s in q.signature], fontsize=7)
        axs[1].axvline(0.5, c="k", lw=0.5, ls=":")
        axs[1].set_xlim(0.3, 1.0); axs[1].set_xlabel("out-of-fold AUC (GPX4 dependence)"); axs[1].legend(fontsize=7)
        axs[1].set_title("Positive controls vs ferroptosis panel", fontsize=8)
    savefig(fig, "fig18_dose_and_positive_controls")


def scorecard_stage(cfg, args):
    s = read_state()
    sm = s["summary"]

    def row(keys):
        out = []
        for k in keys:
            v = sm.get(k)
            out.append(f"- `{k}`: {v['status']} - {v['key_result']}" if v else f"- `{k}`: not run")
        return out
    try:
        _, pstatus, pver = panel_gate(cfg, args)
    except Exception as e:
        pstatus, pver = f"unknown ({str(e)[:60]})", ""
    cl = read_csv_opt("cellline_classifier_metrics.csv")
    ds = read_csv_opt("dose_prediction_results.csv")
    pcc = read_csv_opt("positive_control_cellline.csv")
    pct = read_csv_opt("positive_control_tumour.csv")
    sv = read_csv_opt("subtype_vs_random.csv")
    L = [f"# Scorecard (pipeline v{VERSION}; {time.strftime('%Y-%m-%d %H:%M')})", "",
         f"**Gene panel status: {pstatus}** (version {pver}). Numbers from a PROVISIONAL panel must not be reported as the paper's results.", "",
         "Reading rule: every ferroptosis result is judged against random gene sets of the same size run through the identical procedure.", "",
         "## Output 1 - Ferroptosis classification"] + row(["classification_dMMR", "ext_ferroptosis_subtypes"])
    if sv is not None:
        L += ["- subtypes vs survival, panel vs random sets: " + "; ".join(f"{r.cohort} emp p={r.emp_p_abs_logHR:.2f}" for r in sv.itertuples()),
              "- verdict: " + _verdict(sv.emp_p_abs_logHR.min() if len(sv) else None) + " (smallest p across cohorts; not multiplicity-corrected)"]
    L += ["", "## Output 2 - Prognosis prediction"] + row(["prognosis_external_validation", "novelty1_specificity_prognosis"])
    L += ["", "## Output 3 - Drug response (sensitive vs resistant, cell lines -> tumours)"] + row(["ext_cellline_sensitivity_classifier", "drug_sensitivity_classification"])
    if cl is not None and len(cl):
        r = cl.iloc[0]
        L += [f"- verdict: {_verdict(r.empirical_p_vs_random_sets)} (panel AUC {r.auc_oof:.3f} vs random-set median {r.random_set_median_auc:.3f}; label = {r.drugs})"]
    L += ["", "## Output 4 - IC50 / dose-response prediction"] + row(["ext_dose_prediction"])
    if ds is not None and len(ds):
        L += ["", md_table(ds[["source", "drug", "response_unit", "n_lines", "spearman_oof", "rho_lo95", "rho_hi95", "random_median_rho",
                                "empirical_p_vs_random_sets", "q_bh"]], 30), "",
              f"- {int((ds.empirical_p_vs_random_sets < 0.05).sum())} of {len(ds)} drug models beat random gene sets at p<0.05 "
              f"({int((ds.q_bh < 0.05).sum())} after BH)."]
    else:
        L += ["- no dose model: see `results/csv/dose_drug_availability.csv` for which compounds exist in your response files"]
    L += ["", "## Can the pipeline detect a real signal at all? (positive controls)"] + row(["ext_positive_control"])
    if pcc is not None:
        L += ["", md_table(pcc, 10)]
    if pct is not None:
        L += ["", md_table(pct, 20)]
    L += ["", "## Still missing / checklist"]
    chk = []
    chk.append(("Original gene panel supplied (--panel my_panel.csv)", "CONFIRMED" in str(pstatus).upper() and "PROVISIONAL" not in str(pstatus).upper()))
    chk.append(("Drug-response file with a GPX4 inhibitor (ML162/ML210/RSL3)",
                bool(ds is not None and len(ds) and ds.drug.str.contains("ml162|ml210|rsl3", case=False, regex=True).any())))
    chk.append(("IC50 (not only AUC/CRISPR) models exist", bool(ds is not None and len(ds) and (ds.response_unit.str.contains("IC50")).any())))
    chk.append(("Positive control recovered a real signal (EMT or ceiling AUC clearly above random)",
                bool(pcc is not None and len(pcc) and (pcc.empirical_p_vs_random_sets < 0.05).any())))
    chk.append(("All four tumour cohorts evaluated", bool(len(s.get("inventory", {})) >= 4)))
    L += [f"- [{'x' if ok else ' '}] {t}" for t, ok in chk]
    L += ["", "In-vitro sensitivity is not a clinical dose or response. Associations are not mechanism."]
    (D["rep"] / "scorecard.md").write_text("\n".join(L), encoding="utf-8")
    LOG.info("wrote %s", D["rep"] / "scorecard.md")
    try:
        fig_dose_poscontrol()
    except Exception as e:
        plt.close("all")
        LOG.warning("figure 'dose/positive controls' skipped: %s", str(e)[:160])
    print("\n" + "\n".join(L[:60]))


STAGE_SUMMARY_KEYS = {
    "prognosis": ["prognosis_external_validation"],
    "benchmark": ["novelty1_specificity_prognosis"],
    "classification": ["classification_dMMR"],
    "drug": ["novelty3_gdsc2_ferroscore_vs_response", "drug_response_models", "drug_sensitivity_classification", "ctrp_replication"],
    "subtypes": ["ext_ferroptosis_subtypes"],
    "cellline": ["ext_cellline_sensitivity_classifier"],
    "poscontrol": ["ext_positive_control"],
    "dose": ["ext_dose_prediction"],
}


def mark_stale_summaries(stage, reason):
    """A stage that did not complete in THIS run must not keep showing an old 'completed' result."""
    s = read_state()
    for k in STAGE_SUMMARY_KEYS.get(stage, []):
        old = s["summary"].get(k)
        if old and old.get("status") == "completed":
            s["summary"][k] = dict(analysis=k, status="NOT RUN in this run (older result discarded)", n="", key_result="",
                                   note=reason[:200])
    write_state(s)


STAGE_FUNCS = [("fetch", stage_fetch), ("preprocess", stage_preprocess), ("prognosis", prognosis_stage), ("benchmark", benchmark_stage),
               ("classification", classification_stage), ("drug", drug_stage), ("subtypes", subtype_stage), ("cellline", cellline_stage), ("poscontrol", poscontrol_stage), ("dose", dose_stage), ("figures", figures_stage),
               ("ext_figures", ext_figures_stage), ("report", reports_stage), ("ext_report", ext_report_stage), ("scorecard", scorecard_stage)]


def final_prepare():
    """--final step 1+2: drop stale config keys that override code defaults; clear cached random controls."""
    cfgp = D["cfg"] / "config.yaml"
    if cfgp.exists():
        user = yaml.safe_load(cfgp.read_text(encoding="utf-8")) or {}
        if "penalizer_grid" in user:
            shutil.copy(cfgp, cfgp.with_name(f"config.bak_{dt.datetime.now():%Y%m%d_%H%M%S}.yaml"))
            user.pop("penalizer_grid")
            cfgp.write_text(yaml.safe_dump(user, sort_keys=False), encoding="utf-8")
            LOG.info("--final: removed stale penalizer_grid from config.yaml (backup saved)")
    for f in ("random_gene_set_benchmark.csv", "random_gene_set_benchmark.meta.json"):
        p = D["csv"] / f
        if p.exists():
            p.unlink()
            LOG.info("--final: cleared cached %s", f)


def final_report():
    """--final step 4: print the key result tables + sanity checks, and save them to results/reports/final_tables.txt."""
    out = []

    def emit(t=""):
        print(t)
        out.append(t)

    def show(name, cols):
        d = read_csv_opt(name)
        emit(f"\n=== {name} ===")
        if d is None:
            emit("(not found - the stage that writes it did not run)")
            return
        d = d[[c for c in cols if c in d.columns]]
        emit(d.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    show("signature_comparison.csv", ["cohort", "endpoint", "signature", "cindex", "empirical_p_one_sided", "n_controls"])
    show("prognosis_validation_metrics.csv", ["cohort", "endpoint", "model", "scheme", "n", "events", "cindex", "cindex_lo95", "cindex_hi95"])
    show("classification_metrics.csv", ["cohort", "scheme", "model", "n", "auc", "auc_lo95", "auc_hi95", "empirical_p_vs_random"])
    # sanity checks on the external cohorts
    emit("\n=== DATA SANITY ===")
    pm, inv = read_csv_opt("prognosis_validation_metrics.csv"), read_csv_opt("dataset_inventory.csv")
    problems = 0
    if pm is not None and inv is not None:
        tot = dict(zip(inv.dataset, inv.n_samples))
        ex = pm[(pm.scheme == "external_frozen") & (pm.model == "molecular") & pm.n.notna()]
        for r in ex.itertuples():
            n, ev, tn = int(r.n), int(r.events) if pd.notna(r.events) else -1, tot.get(r.cohort)
            if ev == n and n > 0:
                emit(f"WARNING {r.cohort}: all {n} evaluated patients are events (no censored patients) -> event coding is probably wrong; run --inspect-geo {r.cohort}")
                problems += 1
            elif tn and n < 0.4 * tn:
                emit(f"WARNING {r.cohort}: only {n} of {int(tn)} samples evaluated -> check missing survival fields / overlap removal")
                problems += 1
            else:
                emit(f"ok      {r.cohort}: {n} patients, {ev} events" + (f" (of {int(tn)} samples)" if tn else ""))
    else:
        emit("(metrics or inventory missing - cannot check)")
    if problems == 0:
        emit("No sanity problems detected.")
    emit("\nFull report: results/reports/results_summary.md   Figures: results/figures/")
    try:
        (D["rep"] / "final_tables.txt").write_text("\n".join(out), encoding="utf-8")
    except Exception as e:
        LOG.warning("could not save final_tables.txt: %s", e)


def run_tests():
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(PipelineTests)
    r = unittest.TextTestRunner(verbosity=2).run(suite)
    return 0 if r.wasSuccessful() else 1


def main():
    ap = argparse.ArgumentParser(description="Ferroptosis CRC pipeline")
    ap.add_argument("--stage", nargs="*", default=["all"], choices=["all"] + [s for s, _ in STAGE_FUNCS])
    ap.add_argument("--refresh", action="store_true", help="re-download files")
    ap.add_argument("--dry-run", action="store_true", help="validate data sources only")
    ap.add_argument("--accept-provisional-panel", action="store_true")
    ap.add_argument("--test", action="store_true")
    ap.add_argument("--reprocess", action="store_true", help="redo preprocessing even if processed files exist")
    ap.add_argument("--final", action="store_true", help="one command: clean config, clear cached controls, run everything, print key results")
    ap.add_argument("--n-random", type=int, default=None, help="override number of random gene sets")
    ap.add_argument("--panel", default=None, help="CSV with your original panel (gene_symbol, gene_role[, source, evidence_or_notes]); marks it CONFIRMED")
    ap.add_argument("--inspect-geo", default=None, help="print characteristics keys of a GEO cohort (e.g. GSE14333) and exit")
    ap.add_argument("--check-data", action="store_true", help="report which cell-line data files are present and exit")
    ap.add_argument("--drugs", nargs="*", default=None, help="target compounds for cell-line labels (default ML162 ML210 RSL3)")
    ap.add_argument("--response-source", default="auto", choices=["auto", "gdsc", "ctrp", "long", "crispr"])
    ap.add_argument("--metric", default="AUC", choices=["AUC", "LN_IC50"], help="GDSC2 response metric for sensitivity labels")
    ap.add_argument("--min-lines", type=int, default=20, help="min CRC cell lines with data per target drug")
    ap.add_argument("--n-perm", type=int, default=200, help="label permutations for the cell-line classifier")
    ap.add_argument("--dose-drugs", nargs="*", default=None, help="compounds for the IC50/dose stage (default: ferroptosis inducers + CRC chemotherapy)")
    ap.add_argument("--list-drugs", action="store_true", help="list ferroptosis-related compounds in the GDSC2 CRC data and exit")
    args = ap.parse_args()
    if args.test:
        sys.exit(run_tests())
    setup_logging()
    if args.inspect_geo:
        ensure_dirs()
        inspect_geo(args.inspect_geo)
        return
    if args.panel:
        ensure_project()
        import_panel(args.panel)
    if args.final:
        args.reprocess = True
        ensure_project()
        final_prepare()
        if not args.panel:
            args.accept_provisional_panel = True
    cfg = load_cfg()
    if args.list_drugs:
        list_drugs(cfg)
        return
    if args.check_data:
        check_data(cfg)
        return
    if args.n_random:
        cfg["n_random_prognosis"] = cfg["n_random_classification"] = args.n_random
    LOG.info("Ferroptosis pipeline v%s | file=%s | root=%s | seed=%s", VERSION, Path(__file__).resolve(), ROOT, cfg["seed"])
    if args.dry_run:
        dry_run(cfg)
        return
    wanted = [s for s, _ in STAGE_FUNCS] if "all" in args.stage else args.stage
    status = {}
    for name, fn in STAGE_FUNCS:
        if name not in wanted:
            continue
        LOG.info("=== STAGE: %s ===", name)
        t0 = time.time()
        try:
            fn(cfg, args)
            status[name] = dict(status="ok", seconds=round(time.time() - t0, 1), note="")
        except DataUnavailable as e:
            status[name] = dict(status="skipped/unavailable", seconds=round(time.time() - t0, 1), note=str(e)[:300])
            mark_stale_summaries(name, str(e))
            LOG.warning("%s: %s", name, e)
        except Exception as e:
            status[name] = dict(status="FAILED", seconds=round(time.time() - t0, 1), note=f"{type(e).__name__}: {str(e)[:250]}")
            mark_stale_summaries(name, f"{type(e).__name__}: {str(e)}")
            LOG.error("%s failed:\n%s", name, traceback.format_exc())
        s = read_state()
        s["stages"][name] = status[name]
        write_state(s)
    if "report" not in wanted:
        try:
            reports_stage(cfg, args)
        except Exception:
            pass
    s = read_state()
    print("\n" + "=" * 78 + "\nEXECUTION SUMMARY (a finished script is not a successful analysis)\n" + "=" * 78)
    for k, v in status.items():
        print(f"  {k:<15} {v['status']:<22} {v['seconds']:>7}s  {v['note'][:90]}")
    print("\nDatasets:")
    for k, v in s["datasets"].items():
        print(f"  {k:<12} {v[:100]}")
    print("\nRetained sample counts:")
    for k, v in s["inventory"].items():
        print(f"  {k:<12} samples={v.get('n_samples')} genes={v.get('n_genes')} events={v.get('events', '')}")
    print("\nAnalyses:")
    for k, v in s["summary"].items():
        print(f"  {k:<42} {v['status']:<20} {str(v['key_result'])[:80]}")
    print(f"\nOutputs: {D['csv']}\n         {D['fig']}\n         {D['rep']}\n         {D['prov']}")
    if args.final:
        final_report()


# ----------------------------------------------------------------------------------------
# UNIT TESTS (python ferro_pipeline.py --test). Synthetic data are used ONLY here, never in the pipeline.
# ----------------------------------------------------------------------------------------
class PipelineTests(unittest.TestCase):
    def test_v15_helpers(self):
        rng = np.random.default_rng(3)
        X = rng.normal(size=(60, 5)); y = X[:, 0] * 2 + rng.normal(0, 0.5, 60)
        p = ridge_oof(X, y, 1, folds=5)
        self.assertGreater(stats.spearmanr(y, p)[0], 0.7)
        genes = SIG_MESENCHYMAL + SIG_EPITHELIAL
        E = pd.DataFrame(rng.normal(size=(40, len(genes))), columns=genes)
        E[SIG_MESENCHYMAL] += 2.0
        self.assertEqual(len(emt_score(E)), 40)
        with self.assertRaises(DataUnavailable):
            prolif_score(E)
        yy = (X[:, 0] > 0).astype(int)
        self.assertGreater(nested_topk_auc(X, yy, 2, 1, 5), 0.9)

    def test_bh(self):
        q = bh_adjust([0.01, 0.04, 0.03, 0.005, np.nan])
        np.testing.assert_allclose(q[:4], [0.02, 0.04, 0.04, 0.02], atol=1e-9)
        self.assertTrue(np.isnan(q[4]))

    def test_stage_event_mmr_parsing(self):
        self.assertEqual(parse_stage("Stage IIA"), 2.0)
        self.assertEqual(parse_stage("Stage IV"), 4.0)
        self.assertEqual(parse_stage("3"), 3.0)
        self.assertTrue(np.isnan(parse_stage("Not Reported")))
        self.assertEqual(list(parse_event(["1", "Dead", "0", "alive", "x"]).fillna(-1)), [1, 1, 0, 0, -1])
        self.assertEqual(list(parse_event(["death", "no death", "recurrence", "no recurrence", "NA"]).fillna(-1)), [1, 0, 1, 0, -1])
        self.assertEqual(parse_mmr("dMMR"), 1.0)
        self.assertEqual(parse_mmr("pMMR"), 0.0)
        self.assertEqual(parse_mmr("MSS"), 0.0)
        self.assertTrue(np.isnan(parse_mmr("NA")))

    def test_spearman_ci(self):
        r, p, lo, hi = spearman_ci(np.arange(30), np.arange(30) + np.random.default_rng(1).normal(0, 3, 30))
        self.assertTrue(lo < r < hi and p < 0.05)

    def test_probe_collapse(self):
        e = pd.DataFrame({"s1": [1, 5, 2, 9], "s2": [1, 6, 3, 8]}, index=["p1", "p2", "p3", "p4"])
        a = pd.DataFrame({"ID": ["p1", "p2", "p3", "p4"], "SYMBOL": ["A", "A", "B///C", "---"]})
        g, info = probes_to_genes(e, a)
        self.assertEqual(list(g.index), ["A"])
        self.assertEqual(g.loc["A", "s1"], 5)

    def test_ferroscore(self):
        E = pd.DataFrame({"D1": [1, 2, 3.0], "D2": [1, 2, 3.0], "S1": [3, 2, 1.0], "S2": [3, 2, 1.0]})
        panel = pd.DataFrame({"gene_symbol": ["D1", "D2", "S1", "S2"], "gene_role": ["driver", "driver", "suppressor", "suppressor"]})
        s, _ = ferro_score(E, panel, 0.8)
        self.assertTrue(s.iloc[2] > s.iloc[0])
        with self.assertRaises(DataUnavailable):
            ferro_score(E.drop(columns=["D1", "D2"]), panel, 0.8)

    def test_series_matrix_parse(self):
        txt = ('!Series_platform_id\t"GPL570"\n!Sample_geo_accession\t"GSM1"\t"GSM2"\n!Sample_title\t"a"\t"b"\n'
               '!Sample_characteristics_ch1\t"os.event: 1"\t"os.event: 0"\n!Sample_characteristics_ch1\t"os.delay (months): 12"\t"os.delay (months): 30"\n'
               '!series_matrix_table_begin\n"ID_REF"\t"GSM1"\t"GSM2"\n"p1"\t1.5\t2.5\n!series_matrix_table_end\n')
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "x.txt.gz"
            with gzip.open(p, "wt") as f:
                f.write(txt)
            P = parse_series_matrix(p)
        self.assertEqual(P["platform"], "GPL570")
        self.assertEqual(P["expr"].shape, (1, 2))
        cl, dd = map_geo_clinical("T", P["chars"], P["titles"], DEFAULT_CFG)
        self.assertEqual(list(cl.os_event), [1.0, 0.0])
        self.assertEqual(list(cl.os_time), [12.0, 30.0])

    def test_geo_override(self):
        chars = pd.DataFrame({"rfs.time": ["10", "20"], "rfs.flag": ["1", "0"]}, index=["GSM1", "GSM2"])
        cfg = copy.deepcopy(DEFAULT_CFG)
        cfg["geo"]["overrides"] = {"T": {"dfs_time": "rfs.time", "dfs_event": "rfs.flag"}}
        cl, _ = map_geo_clinical("T", chars, {}, cfg)
        self.assertEqual(list(cl.dfs_event), [1.0, 0.0])
        self.assertEqual(list(cl.dfs_time), [10.0, 20.0])

    def test_import_panel_validation(self):
        import tempfile
        global D
        old = D
        with tempfile.TemporaryDirectory() as td:
            D = {k: Path(td) / v for k, v in dict(
                cfg="config", raw="data/raw", proc="data/processed", csv="results/csv", fig="results/figures",
                rep="results/reports", log="results/logs", prov="results/provenance", tests="tests").items()}
            try:
                ensure_dirs()
                src = Path(td) / "p.csv"
                pd.DataFrame({"Gene_Symbol": ["gpx4", "ACSL4", "GPX4"], "Gene_Role": ["Suppressor", "Driver", "suppressor"]}).to_csv(src, index=False)
                import_panel(src)
                out = pd.read_csv(D["cfg"] / "ferroptosis_genes.csv")
                self.assertEqual(sorted(out.gene_symbol), ["ACSL4", "GPX4"])
                self.assertTrue(out.confirmed.all())
            finally:
                D = old

    def test_adopt_manual(self):
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            (Path(td) / "my_CTRPv2_download.zip").write_bytes(b"PK\x03\x04x")
            dest = Path(td) / "CTRPv2_ExpandedDataset.zip"
            self.assertIsNotNone(adopt_manual(td, ["*CTRPv2*.zip"], dest))
            self.assertTrue(dest.exists())

    def test_prognosis_machinery_synthetic(self):
        rng = np.random.default_rng(0)
        n = 220
        X = pd.DataFrame(rng.normal(size=(n, 12)), columns=[f"G{i}" for i in range(12)])
        lin = 0.9 * X["G0"] - 0.6 * X["G1"]
        t = rng.exponential(30, n) * np.exp(-lin)
        cens = rng.exponential(60, n)
        time_, ev = np.minimum(t, cens), (t <= cens).astype(int)
        cfg = copy.deepcopy(DEFAULT_CFG)
        cfg["n_boot"] = 50
        m = fit_prog_model(X, time_, ev, list(X.columns), cfg, 1)
        self.assertIn("G0", m["selected"])
        r = model_risk(m, X, "cohort")
        c = cindex(time_, ev, r)
        self.assertGreater(c, 0.6)
        sm = surv_metrics(time_, ev, r, cfg, 1, m["cut_median"])
        self.assertGreater(sm["calibration_slope"], 0)
        self.assertLess(sm["logrank_p_high_vs_low"], 0.05)

    def test_tier_map_merges_defaults(self):
        import tempfile
        global D
        old = D
        with tempfile.TemporaryDirectory() as td:
            D = {k: Path(td) / v for k, v in dict(
                cfg="config", raw="data/raw", proc="data/processed", csv="results/csv", fig="results/figures",
                rep="results/reports", log="results/logs", prov="results/provenance", tests="tests").items()}
            try:
                ensure_dirs()
                (D["cfg"] / "drug_panel.csv").write_text("drug_name,tier,note\nML162,primary,x\nErastin,primary,user override\n", encoding="utf-8")
                t = tier_map({})
                self.assertEqual(t[norm("ML162")], "primary")
                self.assertEqual(t[norm("Erastin")], "primary")  # user file wins over default
                self.assertEqual(t[norm("5-Fluorouracil")], "clinical")  # missing from user file -> added from defaults
                self.assertEqual(t[norm("Irinotecan")], "clinical")
            finally:
                D = old

    def test_empirical_p(self):
        self.assertAlmostEqual(empirical_p(0.7, [0.5] * 99), 1 / 100)
        self.assertAlmostEqual(empirical_p(0.4, [0.5] * 99), 1.0)


if __name__ == "__main__":
    main()