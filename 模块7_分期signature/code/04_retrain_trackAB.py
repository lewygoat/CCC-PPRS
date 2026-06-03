import json, logging, sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import (StratifiedKFold, StratifiedGroupKFold,
                                     cross_val_predict)
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, roc_curve, brier_score_loss
from sklearn.calibration import calibration_curve

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent.parent
M1 = ROOT.parent / "模块1_预检" / "output"
if not M1.exists():
    import glob
    candidates = glob.glob(str(ROOT.parent / "*模块1*" / "output"))
    if candidates:
        M1 = Path(candidates[0])
OUT = ROOT / "output"
OUT.mkdir(parents=True, exist_ok=True)
LOG_DIR = ROOT / "log"
LOG_DIR.mkdir(parents=True, exist_ok=True)

LOG_FILE = LOG_DIR / "retrain_trackAB.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    handlers=[logging.FileHandler(LOG_FILE, mode="w", encoding="utf-8"),
              logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger("module7.trackAB")

TRACK_A_GENES = ["ACOT7", "CREB1", "FOXA2", "IFNA8", "LINC00976", "MIR324", "RRM2"]
TRACK_B_GENES = ["AKR1C1", "AKR1C3", "ALOX15", "ALOX5", "ALOXE3", "ATF3", "ATF4",
                 "CISD2", "FOXA2", "FXN", "IDH2", "KEAP1", "NOX1", "PSAT1",
                 "SLC39A14", "SLC7A11"]

POS_CLASS = "subacute"
NEG_CLASS = "acute"
CS_GRID = np.logspace(-3, 1, 30)
L1_RATIOS = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
RNG = 20260530
N_BOOT = 2000


def load_expr(path):
    df = pd.read_csv(path, index_col=0)
    df.columns = [c.strip().strip('"') for c in df.columns]
    df.index = [str(i).strip().strip('"') for i in df.index]
    return df


def build_dataset(genes):
    md = pd.read_csv(M1 / "metadata.csv")
    md["sample_id"] = md["sample_id"].astype(str).str.strip().str.strip('"')
    expr = load_expr(M1 / "expr_mrna.csv")
    sel = md[md["stage"].isin([POS_CLASS, NEG_CLASS])].copy()
    sel = sel[sel["sample_id"].isin(expr.columns)]
    sel["patient"] = sel["sample_id"].str.split("_").str[0]
    present = [g for g in genes if g in expr.index]
    missing = [g for g in genes if g not in expr.index]
    if missing:
        log.warning(f"Missing genes: {missing}")
    X = expr.loc[present, sel["sample_id"].tolist()].T
    y = (sel["stage"].values == POS_CLASS).astype(int)
    groups = sel["patient"].values
    return X, y, groups, sel, md, expr


def make_clf(C, l1_ratio):
    if l1_ratio == 1.0:
        return make_pipeline(
            StandardScaler(),
            LogisticRegression(penalty="l1", solver="liblinear", C=C,
                               max_iter=10000, random_state=RNG),
        )
    return make_pipeline(
        StandardScaler(),
        LogisticRegression(penalty="elasticnet", solver="saga", C=C,
                           l1_ratio=l1_ratio, max_iter=20000,
                           random_state=RNG, tol=1e-4),
    )


def bootstrap_auc_ci(y, prob, n_boot=N_BOOT):
    rng = np.random.RandomState(RNG)
    aucs = []
    for _ in range(n_boot):
        idx = rng.choice(len(y), len(y), replace=True)
        if len(np.unique(y[idx])) < 2:
            continue
        aucs.append(roc_auc_score(y[idx], prob[idx]))
    aucs = np.array(aucs)
    return float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5)), float(np.mean(aucs))


def run_track(track_name, genes):
    log.info("=" * 70)
    log.info(f"[TRACK {track_name}] genes={genes}")

    X, y, groups, sel, md_full, expr_full = build_dataset(genes)
    Xraw = X.values
    log.info(f"[DATA] {X.shape[0]} samples x {X.shape[1]} genes; "
             f"pos={y.sum()} neg={len(y)-y.sum()}; patients={len(set(groups))}")
    log.info(f"[DATA] gene order: {list(X.columns)}")

    n_group_splits = min(5, len(set(groups[y == 0])), len(set(groups[y == 1])))
    grouped = StratifiedGroupKFold(n_splits=max(2, n_group_splits))
    plain = StratifiedKFold(n_splits=5, shuffle=True, random_state=RNG)

    grid = []
    for l1r in L1_RATIOS:
        for C in CS_GRID:
            try:
                clf = make_clf(C, l1r)
                oof = cross_val_predict(clf, Xraw, y, cv=grouped,
                                        groups=groups, method="predict_proba")[:, 1]
                auc_g = roc_auc_score(y, oof)
                clf_full = make_clf(C, l1r).fit(Xraw, y)
                lr = clf_full.named_steps["logisticregression"]
                nz = int((lr.coef_.ravel() != 0).sum())
                grid.append((C, l1r, auc_g, nz))
            except Exception:
                continue

    if not grid:
        log.error(f"[TRACK {track_name}] Grid search produced no valid results")
        return None

    C_star, l1_star, auc_group_best, nz_best = max(grid, key=lambda t: t[2])
    log.info(f"[BEST] C={C_star:.4g} l1={l1_star} grouped_AUC={auc_group_best:.4f} nz={nz_best}")

    clf_pipe = make_clf(C_star, l1_star)
    oof_grouped = cross_val_predict(clf_pipe, Xraw, y, cv=grouped,
                                    groups=groups, method="predict_proba")[:, 1]
    auc_grouped = roc_auc_score(y, oof_grouped)

    oof_plain = cross_val_predict(make_clf(C_star, l1_star), Xraw, y,
                                  cv=plain, method="predict_proba")[:, 1]
    auc_plain = roc_auc_score(y, oof_plain)

    clf_final = make_clf(C_star, l1_star).fit(Xraw, y)
    resub_prob = clf_final.predict_proba(Xraw)[:, 1]
    auc_resub = roc_auc_score(y, resub_prob)
    coef = clf_final.named_steps["logisticregression"].coef_.ravel()
    nz = int((coef != 0).sum())

    leak_gap = auc_plain - auc_grouped
    overfit_gap = auc_resub - auc_grouped

    log.info(f"[RESULT] resubAUC={auc_resub:.4f} plainCV={auc_plain:.4f} "
             f"groupedCV={auc_grouped:.4f} nonzero={nz}/{len(genes)}")
    log.info(f"[LEAKAGE] gap={leak_gap:.4f}")
    log.info(f"[OVERFIT] gap={overfit_gap:.4f}")

    feats_table = pd.DataFrame({
        "gene": X.columns,
        "coef_standardized": np.round(coef, 4),
        "selected": coef != 0,
    })
    log.info(f"[COEF]\n{feats_table.to_string(index=False)}")

    log.info("[SINGLE-GENE BASELINE]")
    single_rows = []
    for g_name in X.columns:
        Xg = X[[g_name]].values
        try:
            sg_clf = make_pipeline(StandardScaler(),
                                   LogisticRegression(C=1.0, max_iter=5000, random_state=RNG))
            sg_oof = cross_val_predict(sg_clf, Xg, y, cv=grouped,
                                       groups=groups, method="predict_proba")[:, 1]
            sg_auc = roc_auc_score(y, sg_oof)
        except Exception:
            sg_auc = 0.5
        single_rows.append({"gene": g_name, "grouped_cv_auc": round(sg_auc, 4)})
        log.info(f"  {g_name}: {sg_auc:.4f}")

    prefix = f"track{track_name}"

    fpr_g, tpr_g, _ = roc_curve(y, oof_grouped)
    pd.DataFrame({"fpr": fpr_g, "tpr": tpr_g}).to_csv(OUT / f"roc_{prefix}_grouped.csv", index=False)

    fpr_p, tpr_p, _ = roc_curve(y, oof_plain)
    pd.DataFrame({"fpr": fpr_p, "tpr": tpr_p}).to_csv(OUT / f"roc_{prefix}_plain.csv", index=False)

    fpr_r, tpr_r, _ = roc_curve(y, resub_prob)
    pd.DataFrame({"fpr": fpr_r, "tpr": tpr_r}).to_csv(OUT / f"roc_{prefix}_apparent.csv", index=False)

    ci_lo_g, ci_hi_g, mean_g = bootstrap_auc_ci(y, oof_grouped)
    ci_lo_p, ci_hi_p, mean_p = bootstrap_auc_ci(y, oof_plain)
    ci_lo_r, ci_hi_r, mean_r = bootstrap_auc_ci(y, resub_prob)
    boot_ci = {
        "grouped_OOF": {"point_auc": round(auc_grouped, 4), "ci_lo": round(ci_lo_g, 4),
                        "ci_hi": round(ci_hi_g, 4), "boot_mean": round(mean_g, 4), "n_boot": N_BOOT},
        "plain_OOF": {"point_auc": round(auc_plain, 4), "ci_lo": round(ci_lo_p, 4),
                      "ci_hi": round(ci_hi_p, 4), "boot_mean": round(mean_p, 4), "n_boot": N_BOOT},
        "resub": {"point_auc": round(auc_resub, 4), "ci_lo": round(ci_lo_r, 4),
                  "ci_hi": round(ci_hi_r, 4), "boot_mean": round(mean_r, 4), "n_boot": N_BOOT},
    }

    try:
        frac_pos, mean_pred = calibration_curve(y, oof_grouped, n_bins=5, strategy="uniform")
        pd.DataFrame({"mean_predicted": mean_pred, "fraction_positive": frac_pos}).to_csv(
            OUT / f"calibration_{prefix}_grouped.csv", index=False)
        brier = brier_score_loss(y, oof_grouped)
        log.info(f"[CALIBRATION] Brier={brier:.4f}")
    except Exception as e:
        log.warning(f"[CALIBRATION] {e}")

    thresholds = np.arange(0.01, 0.99, 0.01)
    n = len(y)
    nb_model, nb_all = [], []
    for t in thresholds:
        tp = np.sum((oof_grouped >= t) & (y == 1))
        fp = np.sum((oof_grouped >= t) & (y == 0))
        nb_m = tp / n - fp / n * (t / (1 - t))
        nb_a = np.mean(y) - (1 - np.mean(y)) * (t / (1 - t))
        nb_model.append(nb_m)
        nb_all.append(max(nb_a, 0))
    pd.DataFrame({"threshold": thresholds, "net_benefit_model": nb_model,
                  "net_benefit_treat_all": nb_all}).to_csv(OUT / f"dca_{prefix}_grouped.csv", index=False)

    ext_expr_path = M1 / "expr_mrna_external.csv"
    ext_md_path = M1 / "metadata_external.csv"
    ext_result = None
    if ext_expr_path.exists() and ext_md_path.exists():
        ext_expr = load_expr(ext_expr_path)
        ext_md = pd.read_csv(ext_md_path)
        ext_md["sample_id"] = ext_md["sample_id"].astype(str).str.strip().str.strip('"')

        ctrl = md_full[md_full["stage"] == "control"].copy()
        ctrl = ctrl[ctrl["sample_id"].isin(expr_full.columns)]
        present_genes = [g for g in X.columns if g in expr_full.index]
        X_ctrl = expr_full.loc[present_genes, ctrl["sample_id"].tolist()].T
        y_ctrl = np.zeros(len(X_ctrl), dtype=int)

        ext_ich = ext_md[ext_md["sample_id"].isin(ext_expr.columns)].copy()
        present_ext = [g for g in X.columns if g in ext_expr.index]
        missing_ext = [g for g in X.columns if g not in ext_expr.index]
        if missing_ext:
            log.warning(f"[EXT] Missing in external: {missing_ext}")

        if len(present_ext) >= max(1, len(X.columns) - 2):
            X_ext_ich = ext_expr.loc[present_ext, ext_ich["sample_id"].tolist()].T
            y_ext_ich = np.ones(len(X_ext_ich), dtype=int)

            if missing_ext:
                for mg in missing_ext:
                    X_ctrl[mg] = 0.0 if mg not in X_ctrl.columns else X_ctrl[mg]
                    X_ext_ich[mg] = 0.0

            X_ctrl_aligned = X_ctrl[X.columns] if all(c in X_ctrl.columns for c in X.columns) else X_ctrl.reindex(columns=X.columns, fill_value=0)
            X_ext_aligned = X_ext_ich.reindex(columns=X.columns, fill_value=0)
            X_ext_all = pd.concat([X_ctrl_aligned, X_ext_aligned], axis=0)
            y_ext_all = np.concatenate([y_ctrl, y_ext_ich])

            ext_prob = clf_final.predict_proba(X_ext_all.values)[:, 1]
            fpr_e, tpr_e, _ = roc_curve(y_ext_all, ext_prob)
            auc_e = roc_auc_score(y_ext_all, ext_prob)
            pd.DataFrame({"fpr": fpr_e, "tpr": tpr_e}).to_csv(
                OUT / f"roc_{prefix}_external.csv", index=False)
            ci_lo_e, ci_hi_e, mean_e = bootstrap_auc_ci(y_ext_all, ext_prob)
            ext_result = {
                "task": "control_vs_ICH",
                "n_control": int(len(y_ctrl)),
                "n_ich": int(len(y_ext_ich)),
                "n_genes_present": len(present_ext),
                "n_genes_missing": len(missing_ext),
                "auc": round(auc_e, 4),
                "ci_lo": round(ci_lo_e, 4),
                "ci_hi": round(ci_hi_e, 4),
            }
            log.info(f"[EXTERNAL] AUC={auc_e:.4f} [{ci_lo_e:.4f}, {ci_hi_e:.4f}] "
                     f"ctrl={len(y_ctrl)} ICH={len(y_ext_ich)}")
        else:
            log.warning(f"[EXTERNAL] Too many missing genes ({len(missing_ext)}/{len(X.columns)}), skipped")

    pd.DataFrame({
        "sample_id": sel["sample_id"].values,
        "y_true": y,
        "prob_grouped_oof": oof_grouped,
        "prob_plain_oof": oof_plain,
        "prob_resub": resub_prob,
        "stage": sel["stage"].values,
        "patient": groups,
    }).to_csv(OUT / f"sample_predictions_{prefix}.csv", index=False)

    summary = {
        "track": track_name,
        "genes_input": list(X.columns),
        "n_samples": int(X.shape[0]),
        "n_genes_input": int(X.shape[1]),
        "n_genes_selected": nz,
        "chosen_C": float(C_star),
        "chosen_l1_ratio": float(l1_star),
        "resub_auc": round(float(auc_resub), 4),
        "cv_auc_plain": round(float(auc_plain), 4),
        "cv_auc_grouped": round(float(auc_grouped), 4),
        "leak_gap": round(float(leak_gap), 4),
        "overfit_gap": round(float(overfit_gap), 4),
        "selected_features": [g for g, c in zip(X.columns, coef) if c != 0],
        "selected_coefs": [round(float(c), 4) for c in coef if c != 0],
        "all_coefs": {g: round(float(c), 4) for g, c in zip(X.columns, coef)},
        "bootstrap_ci": boot_ci,
        "external_validation": ext_result,
        "single_gene_auc": single_rows,
    }
    (OUT / f"{prefix}_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    feats_table.to_csv(OUT / f"{prefix}_coefficients.csv", index=False)
    pd.DataFrame(single_rows).to_csv(OUT / f"{prefix}_single_gene_auc.csv", index=False)
    (OUT / f"{prefix}_bootstrap_ci.json").write_text(
        json.dumps(boot_ci, indent=2), encoding="utf-8")

    log.info(f"[DONE] Track {track_name} outputs -> {OUT}/{prefix}_*")
    return summary


def main():
    log.info("=" * 70)
    log.info(f"Module 7 Step 2-3 — Track A/B retrain | {datetime.now()}")
    log.info("Replacing stable5 (AUC=0.55) with original LASSO-selected features")
    log.info("=" * 70)

    results = {}
    for name, genes in [("A", TRACK_A_GENES), ("B", TRACK_B_GENES)]:
        r = run_track(name, genes)
        if r:
            results[name] = r

    log.info("=" * 70)
    log.info("[COMPARISON]")
    for name, r in results.items():
        log.info(f"  Track {name}: groupedCV={r['cv_auc_grouped']:.4f} "
                 f"[{r['bootstrap_ci']['grouped_OOF']['ci_lo']:.4f}, "
                 f"{r['bootstrap_ci']['grouped_OOF']['ci_hi']:.4f}] "
                 f"nz={r['n_genes_selected']}/{r['n_genes_input']} "
                 f"overfit_gap={r['overfit_gap']:.4f}")
        if r.get("external_validation"):
            ev = r["external_validation"]
            log.info(f"          external({ev['task']}): AUC={ev['auc']:.4f} "
                     f"[{ev['ci_lo']:.4f}, {ev['ci_hi']:.4f}]")

    (OUT / "trackAB_comparison.json").write_text(
        json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    log.info(f"[ALL DONE] comparison -> {OUT}/trackAB_comparison.json")


if __name__ == "__main__":
    main()
