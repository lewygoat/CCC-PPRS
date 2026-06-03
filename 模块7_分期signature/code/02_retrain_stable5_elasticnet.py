"""Module 7 Step 2 - retrain on bootstrap-stable 5-gene subset with ElasticNet.

Motivation: step-1 highconf LASSO produced grouped-CV AUC=0.599 with 16
features and an overfit gap of 0.353. Module-5 bootstrap (333/1000 iters)
shows only 5 hub candidates carry data-supported stability:

  ALOX15  (acute sig_rate=0.856, subacute=0.952; both 95% CI strictly off zero)
  ACSL1   (acute sig_rate=0.441, subacute=0.868)
  GSTP1   (acute sig_rate=0.399, subacute=0.604)
  HMOX1   (acute sig_rate=0.423)
  GPX4    (acute sig_rate=0.294 — kept as canonical ferroptosis axis)

The other 11 highconf-LASSO features are unstable. This step retrains on the
5-gene subset with ElasticNet (l1_ratios=[0.3,0.5,0.7]) under
StratifiedGroupKFold(5). Same leakage/overfit diagnostics as step 1.
"""
import json
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
OUT_DIR = ROOT / "output"
LOG_DIR = ROOT / "log"
for d in (OUT_DIR, LOG_DIR):
    d.mkdir(parents=True, exist_ok=True)

LOG_FILE = LOG_DIR / "retrain_stable5.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    handlers=[logging.FileHandler(LOG_FILE, mode="w", encoding="utf-8"),
              logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger("module7.step2")

SELECTED_GENES = ["ALOX15", "ACSL1", "GSTP1", "HMOX1", "GPX4"]
POS_CLASS = "subacute"
NEG_CLASS = "acute"
N_FOLDS = 5
CS_GRID = np.logspace(-3, 1, 25)
L1_RATIOS = [0.3, 0.5, 0.7]
RNG = 20260530
LEAK_GAP_WARN = 0.05
OVERFIT_GAP_WARN = 0.15


def load_expr(path):
    df = pd.read_csv(path, index_col=0)
    df.columns = [c.strip().strip('"') for c in df.columns]
    df.index = [str(i).strip().strip('"') for i in df.index]
    return df


def build_dataset():
    md = pd.read_csv(M1 / "metadata.csv")
    md["sample_id"] = md["sample_id"].astype(str).str.strip().str.strip('"')
    expr = load_expr(M1 / "expr_mrna.csv")
    sel = md[md["stage"].isin([POS_CLASS, NEG_CLASS])].copy()
    sel = sel[sel["sample_id"].isin(expr.columns)]
    sel["patient"] = sel["sample_id"].str.split("_").str[0]

    missing = [g for g in SELECTED_GENES if g not in expr.index]
    if missing:
        log.warning(f"[GENE] missing from matrix: {missing}")
    present = [g for g in SELECTED_GENES if g in expr.index]
    X = expr.loc[present, sel["sample_id"].tolist()].T
    y = (sel["stage"].values == POS_CLASS).astype(int)
    groups = sel["patient"].values
    paired = sum(1 for _, g in sel.groupby("patient") if g["stage"].nunique() == 2)
    log.info(f"[DATA] samples={X.shape[0]} genes={X.shape[1]} pos({POS_CLASS})={y.sum()} "
             f"neg({NEG_CLASS})={len(y)-y.sum()} patients={len(set(groups))} "
             f"patients_with_both_stages={paired}")
    log.info(f"[DATA] gene order: {list(X.columns)}")
    return X, y, groups, sel


def make_clf(C, l1_ratio):
    return make_pipeline(
        StandardScaler(),
        LogisticRegression(penalty="elasticnet", solver="saga", C=C,
                           l1_ratio=l1_ratio, max_iter=20000,
                           random_state=RNG, tol=1e-4),
    )


def cv_oof_auc(Xraw, y, splitter, C, l1_ratio, groups=None):
    clf = make_clf(C, l1_ratio)
    kw = {"groups": groups} if groups is not None else {}
    oof = cross_val_predict(clf, Xraw, y, cv=splitter,
                            method="predict_proba", **kw)[:, 1]
    return roc_auc_score(y, oof)


def fit_full(Xraw, y, C, l1_ratio):
    clf = make_clf(C, l1_ratio).fit(Xraw, y)
    coef = clf.named_steps["logisticregression"].coef_.ravel()
    return clf, coef


def main():
    log.info("=" * 70)
    log.info(f"Module 7 step 2 — stable 5-gene ElasticNet retrain | {datetime.now()}")
    log.info(f"Selected genes (Module-5 bootstrap stable): {SELECTED_GENES}")

    X, y, groups, meta = build_dataset()
    Xraw = X.values

    n_group_splits = min(N_FOLDS, len(set(groups[y == 0])), len(set(groups[y == 1])))
    grouped = StratifiedGroupKFold(n_splits=max(2, n_group_splits))
    plain = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=RNG)
    log.info(f"[CV] grouped n_splits={grouped.get_n_splits()}; plain n_splits={N_FOLDS}")

    grid_results = []
    for l1r in L1_RATIOS:
        for C in CS_GRID:
            try:
                auc_g = cv_oof_auc(Xraw, y, grouped, C, l1r, groups=groups)
            except Exception as exc:
                log.warning(f"[GRID] C={C:.4g} l1={l1r}: {exc}")
                continue
            grid_results.append((C, l1r, auc_g))
    if not grid_results:
        log.error("Hyperparameter grid empty — abort.")
        sys.exit(2)

    C_star, l1_star, auc_group = max(grid_results, key=lambda t: t[2])
    auc_plain = cv_oof_auc(Xraw, y, plain, C_star, l1_star)
    clf, coef = fit_full(Xraw, y, C_star, l1_star)
    resub_auc = roc_auc_score(y, clf.predict_proba(Xraw)[:, 1])
    nz = int((coef != 0).sum())
    leak_gap = auc_plain - auc_group
    overfit_gap = resub_auc - auc_group

    log.info("=" * 70)
    log.info(f"[CHOSEN] C={C_star:.4g} l1_ratio={l1_star} "
             f"resubAUC={resub_auc:.4f} plainCV={auc_plain:.4f} "
             f"groupedCV={auc_group:.4f} nonzero={nz}/{len(SELECTED_GENES)}")
    feats_table = pd.DataFrame({
        "gene": X.columns,
        "coef_standardized": np.round(coef, 4),
        "selected": coef != 0,
    })
    log.info("[COEF]\n" + feats_table.to_string(index=False))
    if leak_gap > LEAK_GAP_WARN:
        log.warning(f"[LEAKAGE] plain−grouped gap={leak_gap:.4f}; grouped is authoritative.")
    else:
        log.info(f"[LEAKAGE] gap={leak_gap:.4f} within tolerance.")
    if overfit_gap > OVERFIT_GAP_WARN:
        log.warning(f"[OVERFIT] apparent−groupedCV gap={overfit_gap:.4f} (>{OVERFIT_GAP_WARN}).")
    else:
        log.info(f"[OVERFIT] gap={overfit_gap:.4f} acceptable.")

    # also compute baseline reference: each single gene's grouped CV AUC
    log.info("=" * 70)
    log.info("[SINGLE-GENE BASELINE] grouped CV AUC per gene (univariate logistic):")
    single_rows = []
    for g_name in X.columns:
        Xg = X[[g_name]].values
        try:
            auc_g = cv_oof_auc(Xg, y, grouped, C=1.0, l1_ratio=0.5, groups=groups)
        except Exception as exc:
            log.warning(f"  {g_name}: {exc}")
            continue
        single_rows.append({"gene": g_name, "grouped_cv_auc": round(float(auc_g), 4)})
        log.info(f"  {g_name}: grouped CV AUC = {auc_g:.4f}")
    pd.DataFrame(single_rows).to_csv(OUT_DIR / "stable5_single_gene_auc.csv", index=False)

    summary = {
        "n_samples": int(X.shape[0]),
        "n_genes_input": int(X.shape[1]),
        "n_genes_selected": nz,
        "chosen_C": float(C_star),
        "chosen_l1_ratio": float(l1_star),
        "resub_auc": round(float(resub_auc), 4),
        "cv_auc_plain": round(float(auc_plain), 4),
        "cv_auc_grouped": round(float(auc_group), 4),
        "leak_gap": round(float(leak_gap), 4),
        "overfit_gap": round(float(overfit_gap), 4),
        "selected_features": [g for g, c in zip(X.columns, coef) if c != 0],
        "selected_coefs": [round(float(c), 4) for c in coef if c != 0],
        "all_coefs": {g: round(float(c), 4) for g, c in zip(X.columns, coef)},
    }
    (OUT_DIR / "stable5_retrain_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    feats_table.to_csv(OUT_DIR / "stable5_coefficients.csv", index=False)
    log.info(f"[DONE] outputs -> {OUT_DIR}/stable5_*")


if __name__ == "__main__":
    main()
