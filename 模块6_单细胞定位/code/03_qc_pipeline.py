"""Module 6 Step 1 - scRNA QC pipeline for GSE266873 (PHE single-cell).

Defensive, cache-aware QC. Reads 10X triplets (barcodes/features/matrix) per
sample, computes per-cell QC (nCount, nFeature, percent.mt), applies filters
(nFeature 200-6000, percent.mt < 15%), then runs:
  - DRIFT detection : robust z-score (MAD) of per-sample median percent.mt and
                      median nFeature; flags outlier libraries.
  - OVERFIT/over-filtering detection : flags samples whose retention rate after
                      QC falls below MIN_RETENTION (thresholds too aggressive or
                      library degraded).

Data source resolution (defensive fallback):
  1. real GSE266873 files under data/  (preferred; auto-used when present)
  2. synthetic stand-in under data/synthetic/  (DRY-RUN; clearly tagged)
"""
import gzip
import hashlib
import json
import logging
import re
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io as sio
from scipy import sparse

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.cwd()
DATA_DIR = ROOT / "data"
SYN_DIR = DATA_DIR / "synthetic"
OUT_DIR = ROOT / "output"
LOG_DIR = ROOT / "log"
for d in (OUT_DIR, LOG_DIR):
    d.mkdir(parents=True, exist_ok=True)

LOG_FILE = LOG_DIR / "qc.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    handlers=[
        logging.FileHandler(LOG_FILE, mode="w", encoding="utf-8"),
        logging.StreamHandler(sys.stdout),
    ],
)
log = logging.getLogger("module6_qc")

MIN_NFEATURE = 200
MAX_NFEATURE = 6000
MAX_PCT_MT = 15.0
MIN_RETENTION = 0.50          # below this => over-filtering / degraded library
MAD_Z_DRIFT = 3.0             # robust z threshold for cross-sample drift
# Minimum effect-size gate: robust-z alone is hypersensitive when MAD ~ 0
# (a 1% natural spread can yield z>3). Require a meaningful absolute/relative
# deviation in addition to the z-threshold before flagging drift.
MIN_MTPCT_ABS_DEV = 0.5       # percentage points; MT% drift must exceed this
MIN_NFEATURE_REL_DEV = 0.15   # 15%; nFeature drift must exceed this fraction of median
EXPECTED_CELLS_RANGE = (500, 20000)   # per-sample sane range (immune-sorted)

SAMPLE_GROUP = {
    "GSM8255340": "G1", "GSM8255341": "G1", "GSM8255342": "G1",
    "GSM8255343": "G2", "GSM8255344": "G2", "GSM8255345": "G2",
    "GSM8255346": "G3", "GSM8255347": "G3", "GSM8255348": "G3",
}


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def discover_samples():
    """Return (mode, {gsm: {barcodes,features,matrix}})."""
    def scan(d: Path):
        found = defaultdict(dict)
        if not d.exists():
            return found
        for p in d.glob("GSM*_*"):
            name = p.name
            m = re.match(r"(GSM\d+)_.*_(barcodes\.tsv|features\.tsv|matrix\.mtx)\.gz$", name)
            if not m:
                continue
            gsm, kind = m.group(1), m.group(2).split(".")[0]
            found[gsm][kind] = p
        complete = {g: v for g, v in found.items()
                    if {"barcodes", "features", "matrix"} <= set(v)}
        return complete

    real = scan(DATA_DIR)
    if real:
        return "REAL", real
    syn = scan(SYN_DIR)
    if syn:
        return "DRY-RUN(synthetic)", syn
    return "NONE", {}


def load_10x(files: dict):
    with gzip.open(files["matrix"], "rb") as fh:
        M = sio.mmread(fh).tocsr()            # genes x cells
    with gzip.open(files["features"], "rt") as fh:
        feats = [ln.rstrip("\n").split("\t") for ln in fh]
    symbols = np.array([f[1] if len(f) > 1 else f[0] for f in feats])
    n_bc = 0
    with gzip.open(files["barcodes"], "rt") as fh:
        for _ in fh:
            n_bc += 1
    if M.shape[0] != len(symbols) and M.shape[1] == len(symbols):
        M = M.T.tocsr()
    return M, symbols, n_bc


def per_cell_qc(M, symbols):
    """M: genes x cells (csr). Returns DataFrame per cell."""
    mt_mask = np.array([s.upper().startswith("MT-") for s in symbols])
    Mc = M.tocsc()
    n_count = np.asarray(Mc.sum(axis=0)).ravel()
    n_feature = np.asarray((Mc > 0).sum(axis=0)).ravel()
    mt_counts = np.asarray(Mc[mt_mask, :].sum(axis=0)).ravel() if mt_mask.any() else np.zeros(Mc.shape[1])
    with np.errstate(divide="ignore", invalid="ignore"):
        pct_mt = np.where(n_count > 0, 100.0 * mt_counts / n_count, 0.0)
    return pd.DataFrame({
        "nCount": n_count.astype(int),
        "nFeature": n_feature.astype(int),
        "percent_mt": np.round(pct_mt, 3),
    }), int(mt_mask.sum())


def robust_z(x):
    x = np.asarray(x, dtype=float)
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    if mad == 0:
        return np.zeros_like(x), med, mad
    return 0.6745 * (x - med) / mad, med, mad


def main():
    log.info("=" * 70)
    log.info("Module 6 Step 1: GSE266873 scRNA QC")
    mode, samples = discover_samples()
    log.info(f"[DATA] source mode = {mode}; complete samples = {len(samples)}")
    if mode == "NONE":
        log.error("[DATA] no 10X triplets found in data/ or data/synthetic/ — abort")
        return 1
    if mode.startswith("DRY-RUN"):
        log.warning("[DATA] running on SYNTHETIC stand-in. Real GSE266873 download "
                    "is blocked by this host's NCBI DNS hijack. Outputs are tagged DRY-RUN. "
                    "Drop real *.mtx.gz triplets into data/ and rerun to get real results.")

    per_sample = []
    all_cell_rows = []
    for gsm in sorted(samples):
        files = samples[gsm]
        M, symbols, n_bc = load_10x(files)
        df, n_mt = per_cell_qc(M, symbols)
        n_cells = df.shape[0]
        if not (EXPECTED_CELLS_RANGE[0] <= n_cells <= EXPECTED_CELLS_RANGE[1]):
            log.warning(f"[{gsm}] cell count {n_cells} OUTSIDE sane range {EXPECTED_CELLS_RANGE}")

        keep = (
            (df["nFeature"] >= MIN_NFEATURE)
            & (df["nFeature"] <= MAX_NFEATURE)
            & (df["percent_mt"] < MAX_PCT_MT)
        )
        n_keep = int(keep.sum())
        retention = n_keep / n_cells if n_cells else 0.0

        sha_mtx = sha256_file(files["matrix"])[:12]
        log.info(
            f"[{gsm}] grp={SAMPLE_GROUP.get(gsm,'?')} genes={M.shape[0]} "
            f"cells={n_cells} (bc_file={n_bc}) MTgenes={n_mt} "
            f"medFeat={int(df['nFeature'].median())} medMT%={df['percent_mt'].median():.2f} "
            f"-> keep={n_keep} ({retention:.1%}) mtx_sha={sha_mtx}"
        )
        if M.shape[1] != n_bc and M.shape[0] != n_bc:
            log.warning(f"[{gsm}] barcode-file count {n_bc} != matrix cells {n_cells} — format check")

        per_sample.append({
            "gsm": gsm,
            "group": SAMPLE_GROUP.get(gsm, "?"),
            "n_genes": int(M.shape[0]),
            "n_cells_raw": n_cells,
            "n_mt_genes": n_mt,
            "median_nFeature": int(df["nFeature"].median()),
            "median_nCount": int(df["nCount"].median()),
            "median_pct_mt": round(float(df["percent_mt"].median()), 3),
            "n_kept": n_keep,
            "retention": round(retention, 4),
        })
        df.insert(0, "gsm", gsm)
        df.insert(1, "kept", keep.values.astype(int))
        all_cell_rows.append(df)

    ss = pd.DataFrame(per_sample)

    # ---- DRIFT detection (cross-sample, robust z on median MT% and median nFeature) ----
    z_mt, med_mt, mad_mt = robust_z(ss["median_pct_mt"])
    z_ft, med_ft, mad_ft = robust_z(ss["median_nFeature"])
    ss["z_median_pct_mt"] = np.round(z_mt, 3)
    ss["z_median_nFeature"] = np.round(z_ft, 3)
    # effect-size gates: |dev| must be meaningful, not just statistically large
    mt_abs_dev = np.abs(ss["median_pct_mt"].values - med_mt)
    ft_rel_dev = np.abs(ss["median_nFeature"].values - med_ft) / max(med_ft, 1)
    drift_mt = (np.abs(z_mt) > MAD_Z_DRIFT) & (mt_abs_dev > MIN_MTPCT_ABS_DEV)
    drift_ft = (np.abs(z_ft) > MAD_Z_DRIFT) & (ft_rel_dev > MIN_NFEATURE_REL_DEV)
    # record which samples z-tripped but were gated out (false-drift suppressed)
    z_only = ((np.abs(z_mt) > MAD_Z_DRIFT) | (np.abs(z_ft) > MAD_Z_DRIFT)) & ~(drift_mt | drift_ft)
    ss["drift_flag"] = drift_mt | drift_ft
    drift_samples = ss.loc[ss["drift_flag"], "gsm"].tolist()
    suppressed = ss.loc[z_only, "gsm"].tolist()
    log.info(f"[DRIFT] median MT%%: med={med_mt:.2f} MAD={mad_mt:.2f}; "
             f"median nFeature: med={med_ft:.0f} MAD={mad_ft:.0f}; "
             f"threshold |z|>{MAD_Z_DRIFT} AND (MT%% abs-dev>{MIN_MTPCT_ABS_DEV}pp "
             f"or nFeature rel-dev>{MIN_NFEATURE_REL_DEV:.0%})")
    if suppressed:
        log.info(f"[DRIFT] effect-size gate suppressed false-drift on {suppressed} "
                 f"(z tripped but deviation below meaningful threshold; tiny-MAD artefact).")
    if drift_samples:
        for g in drift_samples:
            r = ss.loc[ss.gsm == g].iloc[0]
            log.warning(f"[DRIFT] {g} flagged: z(MT%)={r.z_median_pct_mt} "
                        f"z(nFeat)={r.z_median_nFeature} medMT%={r.median_pct_mt}")
        log.warning(f"[DRIFT] adjustment: {drift_samples} should be down-weighted or "
                    f"sample-specific MT cutoff applied before integration in step 6.2.")
    else:
        log.info("[DRIFT] no cross-sample drift detected.")

    # ---- OVERFIT / over-filtering detection ----
    over = ss.loc[ss["retention"] < MIN_RETENTION, "gsm"].tolist()
    if over:
        for g in over:
            r = ss.loc[ss.gsm == g].iloc[0]
            log.warning(f"[OVERFIT] {g} retention {r.retention:.1%} < {MIN_RETENTION:.0%} "
                        f"— QC may be too aggressive OR library degraded (medMT%={r.median_pct_mt})")
        log.warning("[OVERFIT] adjustment: inspect whether MAX_PCT_MT/nFeature bounds are "
                    "appropriate for flagged samples; consider per-group thresholds.")
    else:
        log.info(f"[OVERFIT] all samples retain >= {MIN_RETENTION:.0%}; thresholds not over-aggressive.")

    tag = "DRYRUN" if mode.startswith("DRY-RUN") else "REAL"
    ss.to_csv(OUT_DIR / f"qc_per_sample_{tag}.csv", index=False)
    cells = pd.concat(all_cell_rows, ignore_index=True)
    cells.to_csv(OUT_DIR / f"qc_per_cell_{tag}.csv.gz", index=False, compression="gzip")
    log.info(f"[OUTPUT] qc_per_sample_{tag}.csv  qc_per_cell_{tag}.csv.gz")

    # group-level rollup (stage totals for downstream temporal analysis)
    grp = ss.groupby("group").agg(
        n_samples=("gsm", "count"),
        cells_raw=("n_cells_raw", "sum"),
        cells_kept=("n_kept", "sum"),
        med_pct_mt=("median_pct_mt", "median"),
    ).reset_index()
    grp["retention"] = (grp["cells_kept"] / grp["cells_raw"]).round(4)
    grp.to_csv(OUT_DIR / f"qc_per_group_{tag}.csv", index=False)
    log.info("[OUTPUT] qc_per_group_%s.csv\n%s", tag, grp.to_string(index=False))

    # ---- run_log.md ----
    md = [
        f"# Module 6 Step 1 - QC run log ({tag})",
        f"timestamp: {datetime.now().isoformat(timespec='seconds')}",
        f"data source mode: **{mode}**",
        f"QC thresholds: nFeature in [{MIN_NFEATURE}, {MAX_NFEATURE}], percent.mt < {MAX_PCT_MT}%",
        "",
        "## Per-sample QC",
        "| GSM | grp | cells raw | medFeat | medMT% | kept | retention | z(MT%) | drift |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for _, r in ss.iterrows():
        md.append(
            f"| {r.gsm} | {r.group} | {r.n_cells_raw} | {r.median_nFeature} | "
            f"{r.median_pct_mt} | {r.n_kept} | {r.retention:.1%} | "
            f"{r.z_median_pct_mt} | {'YES' if r.drift_flag else '-'} |"
        )
    md += ["", "## Stage rollup", "| group | samples | cells raw | cells kept | retention | med MT% |",
           "|---|---|---|---|---|---|"]
    for _, r in grp.iterrows():
        md.append(f"| {r.group} | {r.n_samples} | {r.cells_raw} | {r.cells_kept} | "
                  f"{r.retention:.1%} | {r.med_pct_mt} |")
    md += ["", "## Flags",
           f"- drift samples: {drift_samples or 'none'}",
           f"- over-filtered samples (retention<{MIN_RETENTION:.0%}): {over or 'none'}"]
    (OUT_DIR / f"run_log_{tag}.md").write_text("\n".join(md), encoding="utf-8")
    log.info(f"[OUTPUT] run_log_{tag}.md")
    log.info("DONE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
