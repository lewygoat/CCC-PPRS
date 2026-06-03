import json, logging, sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, roc_curve
from scipy.stats import spearmanr, mannwhitneyu

ROOT = Path(__file__).resolve().parent.parent
M1 = ROOT.parent / "模块1_预检" / "output"
if not M1.exists():
    import glob
    candidates = glob.glob(str(ROOT.parent / "*模块1*" / "output"))
    if candidates:
        M1 = Path(candidates[0])
M2 = ROOT.parent / "模块2_铁死亡基因集" / "output"
if not M2.exists():
    import glob
    candidates = glob.glob(str(ROOT.parent / "*模块2*" / "output"))
    if candidates:
        M2 = Path(candidates[0])
OUT = ROOT / "output"
OUT.mkdir(parents=True, exist_ok=True)
LOG_DIR = ROOT / "log"
LOG_DIR.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    handlers=[logging.FileHandler(LOG_DIR / "external_validation.log", mode="w", encoding="utf-8"),
              logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger("module7.ext")

TRACK_A_GENES = ["ACOT7", "CREB1", "FOXA2", "IFNA8", "LINC00976", "MIR324", "RRM2"]
TRACK_B_ALL = ["AKR1C1", "AKR1C3", "ALOX15", "ALOX5", "ALOXE3", "ATF3", "ATF4",
               "CISD2", "FOXA2", "FXN", "IDH2", "KEAP1", "NOX1", "PSAT1", "SLC39A14", "SLC7A11"]
TRACK_B_PORTABLE = [g for g in TRACK_B_ALL if g not in ("ALOXE3", "FOXA2", "NOX1")]

POS_CLASS = "subacute"
NEG_CLASS = "acute"
RNG = 20260530
N_BOOT = 2000
CS_GRID = np.logspace(-3, 1, 30)
L1_RATIOS = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0]


def load_expr(path):
    df = pd.read_csv(path, index_col=0)
    df.columns = [c.strip().strip('"') for c in df.columns]
    df.index = [str(i).strip().strip('"') for i in df.index]
    return df


def bootstrap_auc_ci(y, prob, n_boot=N_BOOT):
    rng = np.random.RandomState(RNG)
    aucs = []
    for _ in range(n_boot):
        idx = rng.choice(len(y), len(y), replace=True)
        if len(np.unique(y[idx])) < 2:
            continue
        aucs.append(roc_auc_score(y[idx], prob[idx]))
    return (round(float(np.percentile(aucs, 2.5)), 4),
            round(float(np.percentile(aucs, 97.5)), 4),
            round(float(np.mean(aucs)), 4))


def make_clf(C, l1_ratio):
    if l1_ratio >= 0.99:
        return LogisticRegression(penalty="l1", solver="liblinear", C=C,
                                  max_iter=10000, random_state=RNG)
    return LogisticRegression(penalty="elasticnet", solver="saga", C=C,
                              l1_ratio=l1_ratio, max_iter=20000,
                              random_state=RNG, tol=1e-4)


def grid_search_grouped(X, y, groups):
    n_splits = min(5, len(set(groups[y == 0])), len(set(groups[y == 1])))
    grouped = StratifiedGroupKFold(n_splits=max(2, n_splits))
    best = (None, None, -1)
    for l1r in L1_RATIOS:
        for C in CS_GRID:
            try:
                pipe = make_pipeline(StandardScaler(), make_clf(C, l1r))
                oof = cross_val_predict(pipe, X, y, cv=grouped,
                                        groups=groups, method="predict_proba")[:, 1]
                auc = roc_auc_score(y, oof)
                if auc > best[2]:
                    best = (C, l1r, auc)
            except Exception:
                continue
    return best


def run_external_validation(name, gene_list, expr_int, md_int, expr_ext, md_ext):
    log.info("=" * 70)
    log.info(f"[{name}] genes={gene_list}")

    sel = md_int[md_int["stage"].isin([POS_CLASS, NEG_CLASS])].copy()
    sel = sel[sel["sample_id"].isin(expr_int.columns)]
    sel["patient"] = sel["sample_id"].str.split("_").str[0]
    present_int = [g for g in gene_list if g in expr_int.index]
    present_ext = [g for g in gene_list if g in expr_ext.index]
    common = [g for g in gene_list if g in expr_int.index and g in expr_ext.index]

    if len(common) < 2:
        log.warning(f"[{name}] Only {len(common)} common genes, skipping")
        return None

    log.info(f"[{name}] Internal: {len(present_int)}/{len(gene_list)}, "
             f"External: {len(present_ext)}/{len(gene_list)}, Common: {len(common)}")

    X_train = expr_int.loc[common, sel["sample_id"].tolist()].T
    y_train = (sel["stage"].values == POS_CLASS).astype(int)
    groups = sel["patient"].values
    log.info(f"[{name}] Training: {X_train.shape[0]} x {X_train.shape[1]}, "
             f"pos={y_train.sum()} neg={len(y_train)-y_train.sum()}")

    C_star, l1_star, auc_grouped = grid_search_grouped(X_train.values, y_train, groups)
    if C_star is None:
        log.error(f"[{name}] Grid search failed")
        return None
    log.info(f"[{name}] Best: C={C_star:.4g} l1={l1_star} grouped_AUC={auc_grouped:.4f}")

    pipe = make_pipeline(StandardScaler(), make_clf(C_star, l1_star))
    pipe.fit(X_train.values, y_train)
    scaler_int = pipe.named_steps["standardscaler"]
    lr = pipe.named_steps["logisticregression"]
    coef = lr.coef_.ravel()
    intercept = lr.intercept_[0]
    nz = int((coef != 0).sum())
    resub_prob = pipe.predict_proba(X_train.values)[:, 1]
    resub_auc = roc_auc_score(y_train, resub_prob)

    n_splits = min(5, len(set(groups[y_train == 0])), len(set(groups[y_train == 1])))
    grouped_cv = StratifiedGroupKFold(n_splits=max(2, n_splits))
    oof = cross_val_predict(make_pipeline(StandardScaler(), make_clf(C_star, l1_star)),
                            X_train.values, y_train, cv=grouped_cv,
                            groups=groups, method="predict_proba")[:, 1]
    int_ci_lo, int_ci_hi, _ = bootstrap_auc_ci(y_train, oof)
    fpr_i, tpr_i, _ = roc_curve(y_train, oof)

    ctrl = md_int[md_int["stage"] == "control"].copy()
    ctrl = ctrl[ctrl["sample_id"].isin(expr_int.columns)]
    X_ctrl_raw = expr_int.loc[common, ctrl["sample_id"].tolist()].T

    ext_ich = md_ext[md_ext["sample_id"].isin(expr_ext.columns)].copy()
    X_ext_raw = expr_ext.loc[common, ext_ich["sample_id"].tolist()].T

    scaler_ext = StandardScaler()
    X_ext_z = scaler_ext.fit_transform(X_ext_raw.values)

    X_ctrl_z = scaler_int.transform(X_ctrl_raw.values)

    def score_from_z(Z):
        logits = Z @ coef + intercept
        prob = 1 / (1 + np.exp(-logits))
        return prob

    score_ctrl = score_from_z(X_ctrl_z)
    score_acute = score_from_z(scaler_int.transform(X_train.values[y_train == 0]))
    score_subacute = score_from_z(scaler_int.transform(X_train.values[y_train == 1]))
    score_ext = score_from_z(X_ext_z)

    y_ext_task = np.concatenate([np.zeros(len(score_ctrl)), np.ones(len(score_ext))])
    prob_ext_task = np.concatenate([score_ctrl, score_ext])

    if len(np.unique(y_ext_task)) == 2 and not np.any(np.isnan(prob_ext_task)):
        ext_auc = roc_auc_score(y_ext_task, prob_ext_task)
        fpr_e, tpr_e, _ = roc_curve(y_ext_task, prob_ext_task)
        ext_ci_lo, ext_ci_hi, _ = bootstrap_auc_ci(y_ext_task, prob_ext_task)
    else:
        ext_auc, fpr_e, tpr_e = 0.5, np.array([0, 1]), np.array([0, 1])
        ext_ci_lo, ext_ci_hi = 0.0, 1.0

    log.info(f"[{name}] External (ctrl vs ICH, per-platform z-score): "
             f"AUC={ext_auc:.4f} [{ext_ci_lo}, {ext_ci_hi}]")

    U_ce, p_ce = mannwhitneyu(score_ctrl, score_ext, alternative="less")
    U_ca, p_ca = mannwhitneyu(score_ctrl, score_acute, alternative="less")
    log.info(f"[{name}] Score ordering:")
    log.info(f"  control: {score_ctrl.mean():.4f} (std={score_ctrl.std():.4f})")
    log.info(f"  acute(int): {score_acute.mean():.4f}")
    log.info(f"  subacute(int): {score_subacute.mean():.4f}")
    log.info(f"  ICH(ext): {score_ext.mean():.4f} (std={score_ext.std():.4f})")
    log.info(f"  ctrl<extICH: U={U_ce:.0f} p={p_ce:.4g}")

    logfc_int, logfc_ext = {}, {}
    for g in common:
        logfc_int[g] = float(X_train.loc[y_train == 1, g].mean() - X_train.loc[y_train == 0, g].mean())
        logfc_ext[g] = float(X_ext_raw[g].mean() - X_ctrl_raw[g].mean())
    rho_lfc, p_lfc = spearmanr([logfc_int[g] for g in common], [logfc_ext[g] for g in common])
    dir_match = sum(1 for g in common if np.sign(logfc_int[g]) == np.sign(logfc_ext[g]))
    log.info(f"[{name}] Direction consistency: rho={rho_lfc:.4f} p={p_lfc:.4g} "
             f"match={dir_match}/{len(common)} ({dir_match/len(common)*100:.0f}%)")

    prefix = name.lower().replace(" ", "_").replace("-", "_")
    pd.DataFrame({"fpr": fpr_e, "tpr": tpr_e}).to_csv(OUT / f"roc_ext_{prefix}.csv", index=False)
    pd.DataFrame({"fpr": fpr_i, "tpr": tpr_i}).to_csv(OUT / f"roc_int_{prefix}_grouped.csv", index=False)

    score_df = pd.DataFrame({
        "sample_id": (list(ctrl["sample_id"].values) +
                      list(sel.loc[y_train == 0, "sample_id"].values) +
                      list(sel.loc[y_train == 1, "sample_id"].values) +
                      list(ext_ich["sample_id"].values)),
        "group": (["control"] * len(score_ctrl) + ["acute_internal"] * len(score_acute) +
                  ["subacute_internal"] * len(score_subacute) + ["ICH_external"] * len(score_ext)),
        "score": np.concatenate([score_ctrl, score_acute, score_subacute, score_ext]),
        "dataset": (["GSE296792"] * (len(score_ctrl) + len(score_acute) + len(score_subacute)) +
                    ["GSE125512"] * len(score_ext)),
    })
    score_df.to_csv(OUT / f"scores_{prefix}.csv", index=False)

    dir_df = pd.DataFrame({
        "gene": common,
        "coef": [round(float(coef[i]), 4) for i, g in enumerate(common)],
        "logFC_int_acute_vs_sub": [round(logfc_int[g], 4) for g in common],
        "logFC_ext_ctrl_vs_ICH": [round(logfc_ext[g], 4) for g in common],
        "direction_match": [np.sign(logfc_int[g]) == np.sign(logfc_ext[g]) for g in common],
    })
    dir_df.to_csv(OUT / f"direction_{prefix}.csv", index=False)

    return {
        "name": name,
        "n_genes_common": len(common),
        "genes_common": common,
        "chosen_C": float(C_star),
        "chosen_l1_ratio": float(l1_star),
        "n_selected": nz,
        "selected_features": [g for g, c in zip(common, coef) if c != 0],
        "internal": {
            "resub_auc": round(resub_auc, 4),
            "grouped_cv_auc": round(auc_grouped, 4),
            "ci_lo": int_ci_lo, "ci_hi": int_ci_hi,
            "overfit_gap": round(resub_auc - auc_grouped, 4),
        },
        "external": {
            "task": "control_vs_ICH_perplatform_zscore",
            "n_control": len(score_ctrl), "n_ich": len(score_ext),
            "auc": round(ext_auc, 4),
            "ci_lo": ext_ci_lo, "ci_hi": ext_ci_hi,
        },
        "score_ordering": {
            "control_mean": round(float(score_ctrl.mean()), 4),
            "acute_mean": round(float(score_acute.mean()), 4),
            "subacute_mean": round(float(score_subacute.mean()), 4),
            "ICH_ext_mean": round(float(score_ext.mean()), 4),
            "ctrl_vs_ext_p": round(float(p_ce), 6),
        },
        "direction_consistency": {
            "spearman_rho": round(float(rho_lfc), 4),
            "spearman_p": round(float(p_lfc), 6),
            "match_count": dir_match,
            "match_pct": round(dir_match / len(common) * 100, 1),
        },
    }


def portable_from_pool(expr_int, md_int, expr_ext, md_ext):
    log.info("=" * 70)
    log.info("[PORTABLE] Highconf pool restricted to cross-platform genes, with sparsity control")

    hc_pool = pd.read_csv(M2 / "ferroptosis_geneset_high_confidence.csv")
    hc_genes = set(hc_pool["symbol"].astype(str))
    common = sorted(hc_genes & set(expr_int.index) & set(expr_ext.index))
    log.info(f"[PORTABLE] Cross-platform highconf genes: {len(common)}")

    sel = md_int[md_int["stage"].isin([POS_CLASS, NEG_CLASS])].copy()
    sel = sel[sel["sample_id"].isin(expr_int.columns)]
    sel["patient"] = sel["sample_id"].str.split("_").str[0]

    detect = (expr_int.loc[common, sel["sample_id"].tolist()] != 0).sum(axis=1)
    usable = detect[detect >= 3].index.tolist()
    log.info(f"[PORTABLE] After low-detect filter: {len(usable)}")

    X_train = expr_int.loc[usable, sel["sample_id"].tolist()].T
    y_train = (sel["stage"].values == POS_CLASS).astype(int)
    groups = sel["patient"].values

    n_splits = min(5, len(set(groups[y_train == 0])), len(set(groups[y_train == 1])))
    grouped = StratifiedGroupKFold(n_splits=max(2, n_splits))

    best = (None, None, -1, -1)
    sparse_l1 = [0.5, 0.7, 0.9, 1.0]
    sparse_C = np.logspace(-3, 0, 20)
    for l1r in sparse_l1:
        for C in sparse_C:
            try:
                pipe = make_pipeline(StandardScaler(), make_clf(C, l1r))
                oof = cross_val_predict(pipe, X_train.values, y_train,
                                        cv=grouped, groups=groups, method="predict_proba")[:, 1]
                auc = roc_auc_score(y_train, oof)
                pipe_full = make_pipeline(StandardScaler(), make_clf(C, l1r)).fit(X_train.values, y_train)
                nz = int((pipe_full.named_steps["logisticregression"].coef_.ravel() != 0).sum())
                if nz > 20:
                    continue
                if auc > best[2]:
                    best = (C, l1r, auc, nz)
            except Exception:
                continue

    C_star, l1_star, auc_grouped, nz_best = best
    if C_star is None:
        log.error("[PORTABLE] Grid search failed")
        return None
    log.info(f"[PORTABLE] Best: C={C_star:.4g} l1={l1_star} AUC={auc_grouped:.4f} nz={nz_best}")

    pipe = make_pipeline(StandardScaler(), make_clf(C_star, l1_star))
    pipe.fit(X_train.values, y_train)
    scaler_int = pipe.named_steps["standardscaler"]
    lr = pipe.named_steps["logisticregression"]
    coef = lr.coef_.ravel()
    intercept = lr.intercept_[0]
    sel_genes = [g for g, c in zip(usable, coef) if c != 0]
    sel_coefs = [round(float(c), 4) for c in coef if c != 0]
    nz = len(sel_genes)
    resub_auc = roc_auc_score(y_train, pipe.predict_proba(X_train.values)[:, 1])
    log.info(f"[PORTABLE] resubAUC={resub_auc:.4f} nz={nz}")
    log.info(f"[PORTABLE] Selected: {list(zip(sel_genes, sel_coefs))}")

    oof = cross_val_predict(make_pipeline(StandardScaler(), make_clf(C_star, l1_star)),
                            X_train.values, y_train, cv=grouped,
                            groups=groups, method="predict_proba")[:, 1]
    int_ci_lo, int_ci_hi, _ = bootstrap_auc_ci(y_train, oof)
    fpr_i, tpr_i, _ = roc_curve(y_train, oof)

    ctrl = md_int[md_int["stage"] == "control"].copy()
    ctrl = ctrl[ctrl["sample_id"].isin(expr_int.columns)]
    X_ctrl_raw = expr_int.loc[usable, ctrl["sample_id"].tolist()].T
    X_ctrl_z = scaler_int.transform(X_ctrl_raw.values)

    ext_ich = md_ext[md_ext["sample_id"].isin(expr_ext.columns)].copy()
    X_ext_raw = expr_ext.loc[usable, ext_ich["sample_id"].tolist()].T
    scaler_ext = StandardScaler()
    X_ext_z = scaler_ext.fit_transform(X_ext_raw.values)

    def score_z(Z):
        logits = Z @ coef + intercept
        return 1 / (1 + np.exp(-logits))

    score_ctrl = score_z(X_ctrl_z)
    score_acute = score_z(scaler_int.transform(X_train.values[y_train == 0]))
    score_sub = score_z(scaler_int.transform(X_train.values[y_train == 1]))
    score_ext = score_z(X_ext_z)

    y_ext_task = np.concatenate([np.zeros(len(score_ctrl)), np.ones(len(score_ext))])
    prob_ext_task = np.concatenate([score_ctrl, score_ext])
    valid = ~np.isnan(prob_ext_task)
    if valid.sum() > 0 and len(np.unique(y_ext_task[valid])) == 2:
        ext_auc = roc_auc_score(y_ext_task[valid], prob_ext_task[valid])
        fpr_e, tpr_e, _ = roc_curve(y_ext_task[valid], prob_ext_task[valid])
        ext_ci_lo, ext_ci_hi, _ = bootstrap_auc_ci(y_ext_task[valid], prob_ext_task[valid])
    else:
        ext_auc = 0.5
        fpr_e, tpr_e = np.array([0, 1]), np.array([0, 1])
        ext_ci_lo, ext_ci_hi = 0.0, 1.0

    U, p_val = mannwhitneyu(score_ctrl, score_ext, alternative="less")
    log.info(f"[PORTABLE] External AUC={ext_auc:.4f} [{ext_ci_lo}, {ext_ci_hi}]")
    log.info(f"[PORTABLE] ctrl={score_ctrl.mean():.4f} acute={score_acute.mean():.4f} "
             f"sub={score_sub.mean():.4f} ext={score_ext.mean():.4f} p={p_val:.4g}")

    logfc_int, logfc_ext_d = {}, {}
    for g in sel_genes:
        logfc_int[g] = float(X_train.loc[y_train == 1, g].mean() - X_train.loc[y_train == 0, g].mean())
        logfc_ext_d[g] = float(X_ext_raw[g].mean() - X_ctrl_raw[g].mean())
    if len(sel_genes) >= 3:
        rho, p_rho = spearmanr([logfc_int[g] for g in sel_genes],
                               [logfc_ext_d[g] for g in sel_genes])
        dm = sum(1 for g in sel_genes if np.sign(logfc_int[g]) == np.sign(logfc_ext_d[g]))
    else:
        rho, p_rho, dm = np.nan, np.nan, 0
    log.info(f"[PORTABLE] Direction: rho={rho:.4f} p={p_rho:.4g} match={dm}/{len(sel_genes)}")

    pd.DataFrame({"fpr": fpr_e, "tpr": tpr_e}).to_csv(OUT / "roc_ext_portable.csv", index=False)
    pd.DataFrame({"fpr": fpr_i, "tpr": tpr_i}).to_csv(OUT / "roc_int_portable_grouped.csv", index=False)

    score_df = pd.DataFrame({
        "sample_id": (list(ctrl["sample_id"].values) +
                      list(sel.loc[y_train == 0, "sample_id"].values) +
                      list(sel.loc[y_train == 1, "sample_id"].values) +
                      list(ext_ich["sample_id"].values)),
        "group": (["control"] * len(score_ctrl) + ["acute_internal"] * len(score_acute) +
                  ["subacute_internal"] * len(score_sub) + ["ICH_external"] * len(score_ext)),
        "score": np.concatenate([score_ctrl, score_acute, score_sub, score_ext]),
        "dataset": (["GSE296792"] * (len(score_ctrl) + len(score_acute) + len(score_sub)) +
                    ["GSE125512"] * len(score_ext)),
    })
    score_df.to_csv(OUT / "scores_portable.csv", index=False)

    pd.DataFrame({
        "gene": sel_genes, "coef": sel_coefs,
        "logFC_int": [round(logfc_int[g], 4) for g in sel_genes],
        "logFC_ext": [round(logfc_ext_d[g], 4) for g in sel_genes],
        "dir_match": [np.sign(logfc_int[g]) == np.sign(logfc_ext_d[g]) for g in sel_genes],
    }).to_csv(OUT / "portable_signature_detail.csv", index=False)

    return {
        "name": "Portable_highconf_sparse",
        "n_pool": len(common), "n_usable": len(usable), "n_selected": nz,
        "selected_genes": sel_genes, "selected_coefs": sel_coefs,
        "chosen_C": float(C_star), "chosen_l1_ratio": float(l1_star),
        "internal": {
            "resub_auc": round(resub_auc, 4), "grouped_cv_auc": round(auc_grouped, 4),
            "ci_lo": int_ci_lo, "ci_hi": int_ci_hi,
            "overfit_gap": round(resub_auc - auc_grouped, 4),
        },
        "external": {
            "task": "control_vs_ICH_perplatform_zscore",
            "n_control": len(score_ctrl), "n_ich": len(score_ext),
            "auc": round(ext_auc, 4), "ci_lo": ext_ci_lo, "ci_hi": ext_ci_hi,
        },
        "score_ordering": {
            "control_mean": round(float(score_ctrl.mean()), 4),
            "acute_mean": round(float(score_acute.mean()), 4),
            "subacute_mean": round(float(score_sub.mean()), 4),
            "ICH_ext_mean": round(float(score_ext.mean()), 4),
            "ctrl_vs_ext_p": round(float(p_val), 6),
        },
        "direction_consistency": {
            "spearman_rho": round(float(rho), 4) if not np.isnan(rho) else None,
            "spearman_p": round(float(p_rho), 6) if not np.isnan(p_rho) else None,
            "match_count": int(dm),
            "match_pct": round(dm / len(sel_genes) * 100, 1) if sel_genes else 0,
        },
    }


def main():
    log.info("=" * 70)
    log.info(f"Module 7 — External Validation (per-platform z-score) | {datetime.now()}")
    log.info("=" * 70)

    expr_int = load_expr(M1 / "expr_mrna.csv")
    md_int = pd.read_csv(M1 / "metadata.csv")
    md_int["sample_id"] = md_int["sample_id"].astype(str).str.strip().str.strip('"')
    expr_ext = load_expr(M1 / "expr_mrna_external.csv")
    md_ext = pd.read_csv(M1 / "metadata_external.csv")
    md_ext["sample_id"] = md_ext["sample_id"].astype(str).str.strip().str.strip('"')

    results = {}

    r1 = run_external_validation("TrackB_13gene", TRACK_B_PORTABLE,
                                 expr_int, md_int, expr_ext, md_ext)
    if r1:
        results["trackB_13gene"] = r1

    r2 = run_external_validation("TrackA_3gene", ["ACOT7", "CREB1", "RRM2"],
                                 expr_int, md_int, expr_ext, md_ext)
    if r2:
        results["trackA_3gene"] = r2

    r3 = portable_from_pool(expr_int, md_int, expr_ext, md_ext)
    if r3:
        results["portable_highconf"] = r3

    log.info("=" * 70)
    log.info("[FINAL SUMMARY]")
    for k, r in results.items():
        iv = r["internal"]
        ev = r["external"]
        dc = r.get("direction_consistency", {})
        log.info(f"  {k} ({r.get('n_selected', '?')} genes):")
        log.info(f"    Internal: {iv['grouped_cv_auc']:.4f} [{iv['ci_lo']}, {iv['ci_hi']}]")
        log.info(f"    External: {ev['auc']:.4f} [{ev['ci_lo']}, {ev['ci_hi']}]")
        log.info(f"    Scores: ctrl={r['score_ordering']['control_mean']:.3f} "
                 f"acute={r['score_ordering']['acute_mean']:.3f} "
                 f"sub={r['score_ordering']['subacute_mean']:.3f} "
                 f"ext={r['score_ordering']['ICH_ext_mean']:.3f}")
        if dc.get("spearman_rho") is not None:
            log.info(f"    Direction: rho={dc['spearman_rho']:.4f} match={dc['match_count']}/{dc.get('match_pct','?')}%")

    (OUT / "external_validation_all.json").write_text(
        json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    log.info(f"[DONE] -> {OUT}/external_validation_all.json")


if __name__ == "__main__":
    main()
