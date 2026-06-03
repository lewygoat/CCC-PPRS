import json, sys
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, roc_curve, brier_score_loss
from sklearn.calibration import calibration_curve

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent
M1 = ROOT.parent / "模块1_预检" / "output"
OUT = ROOT / "output"
OUT.mkdir(parents=True, exist_ok=True)

GENES = ["ALOX15", "ACSL1", "GSTP1", "HMOX1", "GPX4"]
RNG = 20260530


def load_expr(path):
    df = pd.read_csv(path, index_col=0)
    df.columns = [c.strip().strip('"') for c in df.columns]
    df.index = [str(i).strip().strip('"') for i in df.index]
    return df


def main():
    md = pd.read_csv(M1 / "metadata.csv")
    md["sample_id"] = md["sample_id"].astype(str).str.strip().str.strip('"')
    expr = load_expr(M1 / "expr_mrna.csv")

    sel = md[md["stage"].isin(["acute", "subacute"])].copy()
    sel = sel[sel["sample_id"].isin(expr.columns)]
    sel["patient"] = sel["sample_id"].str.split("_").str[0]
    X = expr.loc[[g for g in GENES if g in expr.index], sel["sample_id"].tolist()].T
    y = (sel["stage"].values == "subacute").astype(int)
    groups = sel["patient"].values

    best = json.loads((OUT / "stable5_retrain_summary.json").read_text("utf-8"))
    C_star = best["chosen_C"]
    l1_star = best.get("chosen_l1_ratio", best.get("l1_ratio", 0.3))

    clf_pipe = make_pipeline(
        StandardScaler(),
        LogisticRegression(penalty="elasticnet", solver="saga", C=C_star,
                           l1_ratio=l1_star, max_iter=20000,
                           random_state=RNG, tol=1e-4),
    )

    grouped = StratifiedGroupKFold(n_splits=5)
    oof_grouped = cross_val_predict(clf_pipe, X.values, y,
                                    cv=grouped, groups=groups,
                                    method="predict_proba")[:, 1]
    fpr_g, tpr_g, _ = roc_curve(y, oof_grouped)
    auc_g = roc_auc_score(y, oof_grouped)
    pd.DataFrame({"fpr": fpr_g, "tpr": tpr_g}).to_csv(
        OUT / "roc_internal_grouped.csv", index=False)
    print(f"Internal grouped CV AUC = {auc_g:.4f}  ({len(fpr_g)} ROC points)")

    plain = StratifiedKFold(n_splits=5, shuffle=True, random_state=RNG)
    oof_plain = cross_val_predict(clf_pipe, X.values, y,
                                  cv=plain, method="predict_proba")[:, 1]
    fpr_p, tpr_p, _ = roc_curve(y, oof_plain)
    auc_p = roc_auc_score(y, oof_plain)
    pd.DataFrame({"fpr": fpr_p, "tpr": tpr_p}).to_csv(
        OUT / "roc_internal_plain.csv", index=False)
    print(f"Internal plain CV AUC  = {auc_p:.4f}")

    clf_pipe.fit(X.values, y)
    resub_prob = clf_pipe.predict_proba(X.values)[:, 1]
    fpr_r, tpr_r, _ = roc_curve(y, resub_prob)
    pd.DataFrame({"fpr": fpr_r, "tpr": tpr_r}).to_csv(
        OUT / "roc_internal_apparent.csv", index=False)
    print(f"Internal apparent AUC  = {roc_auc_score(y, resub_prob):.4f}")

    frac_pos, mean_pred = calibration_curve(y, oof_grouped, n_bins=5, strategy="uniform")
    pd.DataFrame({"mean_predicted": mean_pred, "fraction_positive": frac_pos}).to_csv(
        OUT / "calibration_grouped.csv", index=False)
    brier = brier_score_loss(y, oof_grouped)
    print(f"Brier score (grouped OOF) = {brier:.4f}")

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
    pd.DataFrame({
        "threshold": thresholds,
        "net_benefit_model": nb_model,
        "net_benefit_treat_all": nb_all,
    }).to_csv(OUT / "dca_grouped.csv", index=False)
    print(f"DCA exported ({len(thresholds)} thresholds)")

    ext_expr_path = M1 / "expr_mrna_external.csv"
    ext_md_path = M1 / "metadata_external.csv"
    if ext_expr_path.exists() and ext_md_path.exists():
        ext_expr = load_expr(ext_expr_path)
        ext_md = pd.read_csv(ext_md_path)
        ext_md["sample_id"] = ext_md["sample_id"].astype(str).str.strip().str.strip('"')

        ctrl = md[md["stage"] == "control"].copy()
        ctrl = ctrl[ctrl["sample_id"].isin(expr.columns)]
        X_ctrl = expr.loc[[g for g in GENES if g in expr.index], ctrl["sample_id"].tolist()].T
        y_ctrl = np.zeros(len(X_ctrl), dtype=int)

        ext_ich = ext_md.copy()
        ext_ich = ext_ich[ext_ich["sample_id"].isin(ext_expr.columns)]
        present_ext = [g for g in GENES if g in ext_expr.index]
        X_ext_ich = ext_expr.loc[present_ext, ext_ich["sample_id"].tolist()].T
        y_ext_ich = np.ones(len(X_ext_ich), dtype=int)

        X_ext = pd.concat([X_ctrl, X_ext_ich], axis=0)
        y_ext = np.concatenate([y_ctrl, y_ext_ich])

        ext_prob = clf_pipe.predict_proba(X_ext.values)[:, 1]
        fpr_e, tpr_e, _ = roc_curve(y_ext, ext_prob)
        auc_e = roc_auc_score(y_ext, ext_prob)
        pd.DataFrame({"fpr": fpr_e, "tpr": tpr_e}).to_csv(
            OUT / "roc_external.csv", index=False)
        print(f"External AUC (ctrl vs ICH) = {auc_e:.4f}  "
              f"(ctrl={len(y_ctrl)} ICH={len(y_ext_ich)})")
    else:
        print("[WARN] External data not found, skipping external ROC.")

    pd.DataFrame({
        "sample_id": sel["sample_id"].values,
        "y_true": y,
        "prob_grouped_oof": oof_grouped,
        "prob_plain_oof": oof_plain,
        "prob_resub": resub_prob,
        "stage": sel["stage"].values,
        "patient": groups,
    }).to_csv(OUT / "sample_predictions.csv", index=False)
    print(f"Sample predictions exported ({len(y)} samples)")

    print("\n=== Output files ===")
    for f in sorted(OUT.glob("roc_*.csv")) + sorted(OUT.glob("calibration_*.csv")) + \
             sorted(OUT.glob("dca_*.csv")) + sorted(OUT.glob("sample_*.csv")):
        print(f"  {f.name}")


if __name__ == "__main__":
    main()
