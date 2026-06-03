"""Module 7 Step 1 - binary staging signature: LASSO feature selection.

Scope (1/3 of module 7): 7.1 binary label setup (acute vs subacute) + 7.2
LASSO-logistic with grouped CV. Nomogram (7.3), external validation (7.4),
DCA (7.5) are the remaining 2/3.

Real data: GSE296792 (module-1 preprocessed expr_mrna.csv + metadata.csv).

Defensive strategy:
  - feature space restricted to module-2 ferroptosis gene pool (avoid p>>n LASSO
    on all ~38k genes with n=50). Dual-track: full pool vs high-conf.
  - NUMERICAL STABILITY: drop ferroptosis genes detectable (nonzero) in <3 of 50
    blood samples. These near-constant columns carry no staging signal and make
    StandardScaler emit huge z-scores -> liblinear matmul overflow (observed in
    v1). Scaling is done INSIDE each CV fold via a Pipeline (no full-data leakage).
  - LEAKAGE detection: GSE296792 is longitudinal; 23/27 patients contribute BOTH
    an acute (72h) and a subacute (7day) sample. Plain StratifiedKFold leaks
    patient identity. We run BOTH plain and patient-grouped CV and report the gap;
    the grouped (honest) estimate is authoritative downstream.
  - OVERFIT detection: C is tuned HONESTLY against grouped-CV AUC (not a hardcoded
    value). Overfit gap = resubstitution AUC at selected C minus grouped-CV AUC at
    the SAME C; large gap => LASSO fits the 50 samples but does not generalise.
"""
import logging
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import (StratifiedKFold, StratifiedGroupKFold,
                                     cross_val_predict)
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.cwd()
M1 = ROOT.parent / "模块1_预检" / "output"
M2 = ROOT.parent / "模块2_铁死亡基因集" / "output"
OUT_DIR = ROOT / "output"
LOG_DIR = ROOT / "log"
for d in (OUT_DIR, LOG_DIR):
    d.mkdir(parents=True, exist_ok=True)

LOG_FILE = LOG_DIR / "lasso.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    handlers=[logging.FileHandler(LOG_FILE, mode="w", encoding="utf-8"),
              logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger("module7")

POS_CLASS = "subacute"      # y=1
NEG_CLASS = "acute"         # y=0
N_FOLDS = 10
CS_GRID = np.logspace(-3, 1, 25)
RNG = 20260530
MIN_NONZERO_SAMPLES = 3     # drop genes detected in fewer than this many samples
MIN_FEATS, MAX_FEATS = 1, 20  # admissible LASSO sparsity for C selection
LEAK_GAP_WARN = 0.05        # AUC(plain) - AUC(grouped) above this => leakage matters
OVERFIT_GAP_WARN = 0.15     # resubAUC - groupedCV above this => overfit


def load_expr(path):
    df = pd.read_csv(path, index_col=0)
    df.columns = [c.strip().strip('"') for c in df.columns]
    df.index = [str(i).strip().strip('"') for i in df.index]
    return df


def build_dataset():
    md = pd.read_csv(M1 / "metadata.csv")
    md["sample_id"] = md["sample_id"].astype(str).str.strip().str.strip('"')
    expr = load_expr(M1 / "expr_mrna.csv")          # genes x samples
    sel = md[md["stage"].isin([POS_CLASS, NEG_CLASS])].copy()
    sel = sel[sel["sample_id"].isin(expr.columns)]
    sel["patient"] = sel["sample_id"].str.split("_").str[0]
    X = expr[sel["sample_id"].tolist()].T            # samples x genes
    y = (sel["stage"].values == POS_CLASS).astype(int)
    groups = sel["patient"].values
    paired = sum(1 for _, g in sel.groupby("patient") if g["stage"].nunique() == 2)
    log.info(f"[DATA] {X.shape[0]} samples x {X.shape[1]} genes; "
             f"pos({POS_CLASS})={y.sum()} neg({NEG_CLASS})={len(y)-y.sum()}; "
             f"patients={len(set(groups))}; patients_with_both_stages={paired}")
    return X, y, groups, sel


def restrict_features(X, track):
    if track == "fullpool":
        genes = set(pd.read_csv(M2 / "ferroptosis_geneset.csv")["symbol"].astype(str))
    else:
        genes = set(pd.read_csv(M2 / "ferroptosis_geneset_high_confidence.csv")["symbol"].astype(str))
    present = [g for g in X.columns if g in genes]
    Xf = X[present]
    detect = (Xf != 0).sum(0)
    keep = detect[detect >= MIN_NONZERO_SAMPLES].index.tolist()
    dropped = len(present) - len(keep)
    log.info(f"[FEATURES/{track}] pool={len(genes)} present={len(present)} "
             f"dropped_low_detect(<{MIN_NONZERO_SAMPLES} samples)={dropped} "
             f"-> usable={len(keep)}")
    return Xf[keep]


def make_clf(C):
    return make_pipeline(
        StandardScaler(),
        LogisticRegression(penalty="l1", solver="liblinear", C=C,
                           max_iter=5000, random_state=RNG),
    )


def cv_oof_auc(Xraw, y, splitter, C, groups=None):
    clf = make_clf(C)
    kw = {"groups": groups} if groups is not None else {}
    oof = cross_val_predict(clf, Xraw, y, cv=splitter,
                            method="predict_proba", **kw)[:, 1]
    return roc_auc_score(y, oof)


def n_selected_at(Xraw, y, C):
    clf = make_clf(C).fit(Xraw, y)
    coef = clf.named_steps["logisticregression"].coef_.ravel()
    return int((coef != 0).sum()), clf


def run_track(X, y, groups, track):
    log.info("=" * 70)
    log.info(f"[TRACK] {track}")
    Xf = restrict_features(X, track)
    Xraw = Xf.values

    n_splits = min(N_FOLDS, int(np.bincount(y).min()))
    plain = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=RNG)
    n_group_splits = min(N_FOLDS, len(set(groups[y == 0])), len(set(groups[y == 1])))
    grouped = StratifiedGroupKFold(n_splits=max(2, n_group_splits))

    # honest C selection: maximise grouped-CV AUC over the grid, restricted to
    # models whose full-data fit yields an admissible sparse signature.
    grid = []
    for C in CS_GRID:
        nz, _ = n_selected_at(Xraw, y, C)
        if nz < MIN_FEATS or nz > MAX_FEATS:
            continue
        auc_g = cv_oof_auc(Xraw, y, grouped, C, groups=groups)
        grid.append((C, auc_g, nz))
    if not grid:                       # fallback: ignore sparsity constraint
        for C in CS_GRID:
            auc_g = cv_oof_auc(Xraw, y, grouped, C, groups=groups)
            nz, _ = n_selected_at(Xraw, y, C)
            grid.append((C, auc_g, nz))
    C_star, auc_group, _ = max(grid, key=lambda t: t[1])

    auc_plain = cv_oof_auc(Xraw, y, plain, C_star)
    nz, full_clf = n_selected_at(Xraw, y, C_star)
    resub_auc = roc_auc_score(y, full_clf.predict_proba(Xraw)[:, 1])
    coef = full_clf.named_steps["logisticregression"].coef_.ravel()
    mask = coef != 0
    feats = Xf.columns[mask].tolist()
    coefs = coef[mask]

    leak_gap = auc_plain - auc_group
    overfit_gap = resub_auc - auc_group

    log.info(f"[{track}] selected C={C_star:.4g}  resub(apparent)AUC={resub_auc:.3f}  "
             f"plainCV-AUC={auc_plain:.3f}  groupedCV-AUC={auc_group:.3f}  n_sel={nz}")
    log.info(f"[{track}] signature -> {feats}")
    if leak_gap > LEAK_GAP_WARN:
        log.warning(f"[LEAKAGE/{track}] plain−grouped AUC gap={leak_gap:.3f} "
                    f"(>{LEAK_GAP_WARN}). Patient identity leaks under random CV. "
                    f"ADJUST: grouped-by-patient CV is authoritative; plain CV discarded.")
    else:
        log.info(f"[LEAKAGE/{track}] plain−grouped gap={leak_gap:.3f} within tolerance "
                 f"(no optimistic inflation from paired samples at this C).")
    if overfit_gap > OVERFIT_GAP_WARN:
        log.warning(f"[OVERFIT/{track}] apparent−groupedCV gap={overfit_gap:.3f} "
                    f"(>{OVERFIT_GAP_WARN}). LASSO fits the {len(y)} training samples but "
                    f"does not generalise. ADJUST: signature reported as EXPLORATORY; "
                    f"grouped-CV AUC={auc_group:.3f} is the honest performance.")
    else:
        log.info(f"[OVERFIT/{track}] apparent−groupedCV gap={overfit_gap:.3f} acceptable.")

    return {
        "track": track,
        "n_features_usable": Xf.shape[1],
        "resub_auc": round(float(resub_auc), 4),
        "cv_auc_plain": round(float(auc_plain), 4),
        "cv_auc_grouped": round(float(auc_group), 4),
        "leak_gap": round(float(leak_gap), 4),
        "overfit_gap": round(float(overfit_gap), 4),
        "chosen_C": float(C_star),
        "n_selected": int(nz),
        "selected_features": feats,
        "selected_coefs": [round(float(c), 4) for c in coefs],
    }


def main():
    log.info("=" * 70)
    log.info("Module 7 Step 1: binary staging signature (acute vs subacute)")
    X, y, groups, meta = build_dataset()
    meta.to_csv(OUT_DIR / "binary_labels.csv", index=False)

    results = [run_track(X, y, groups, t) for t in ("fullpool", "highconf")]

    df = pd.DataFrame([{k: v for k, v in r.items()
                        if k not in ("selected_features", "selected_coefs")}
                       for r in results])
    df.to_csv(OUT_DIR / "lasso_track_summary.csv", index=False)

    sel_rows = []
    for r in results:
        for f, c in zip(r["selected_features"], r["selected_coefs"]):
            sel_rows.append({"track": r["track"], "feature": f, "coef": c})
    pd.DataFrame(sel_rows).to_csv(OUT_DIR / "lasso_selected_features.csv", index=False)

    md = [
        "# Module 7 Step 1 - LASSO staging signature (run log)",
        f"timestamp: {datetime.now().isoformat(timespec='seconds')}",
        f"binary task: {POS_CLASS}(y=1) vs {NEG_CLASS}(y=0); training=GSE296792",
        f"samples: {X.shape[0]} ; patients: {len(set(groups))} ; "
        f"paired-leakage risk: longitudinal (same patient both stages)",
        "",
        "## Track comparison (C tuned against grouped-CV AUC)",
        "| track | usable feats | apparent AUC | plain CV | grouped CV | leak gap | overfit gap | C | n_sel |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for r in results:
        md.append(f"| {r['track']} | {r['n_features_usable']} | {r['resub_auc']} | "
                  f"{r['cv_auc_plain']} | {r['cv_auc_grouped']} | {r['leak_gap']} | "
                  f"{r['overfit_gap']} | {r['chosen_C']:.4g} | {r['n_selected']} |")
    md += ["", "## Selected features"]
    for r in results:
        pairs = ", ".join(f"{f}({c:+.2f})" for f, c in
                          zip(r["selected_features"], r["selected_coefs"]))
        md.append(f"- **{r['track']}** ({r['n_selected']}): {pairs or 'none'}")
    md += ["", "## Defensive notes",
           "- Authoritative metric = **grouped-by-patient CV AUC** (plain CV would be "
           "inflated by paired leakage under a generalising model).",
           "- Feature space pre-restricted to ferroptosis pool; genes detected in "
           f"<{MIN_NONZERO_SAMPLES}/50 samples dropped (numerical stability).",
           "- Scaling performed inside each CV fold (Pipeline) — no full-data leakage.",
           "- External validation (GSE125512) deferred to 2/3; note: external labels are "
           "'24h_or_7day_unresolved' (all mapped to acute) — cannot resolve early-vs-subacute "
           "there, so external validation will need an alternative target or relabelling."]
    (OUT_DIR / "run_log.md").write_text("\n".join(md), encoding="utf-8")
    log.info("[OUTPUT] lasso_track_summary.csv, lasso_selected_features.csv, run_log.md")
    log.info("DONE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
