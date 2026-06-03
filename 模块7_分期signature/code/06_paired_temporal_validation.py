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
from scipy.stats import wilcoxon, spearmanr, mannwhitneyu

ROOT = Path(__file__).resolve().parent.parent
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

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    handlers=[logging.FileHandler(LOG_DIR / "paired_temporal.log", mode="w", encoding="utf-8"),
              logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger("module7.paired")

TRACK_B_PORTABLE = ["AKR1C1", "AKR1C3", "ALOX15", "ALOX5", "ATF3", "ATF4",
                    "CISD2", "FXN", "IDH2", "KEAP1", "PSAT1", "SLC39A14", "SLC7A11"]
TRACK_A_GENES = ["ACOT7", "CREB1", "FOXA2", "IFNA8", "LINC00976", "MIR324", "RRM2"]

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


def make_clf(C, l1_ratio):
    if l1_ratio >= 0.99:
        return LogisticRegression(penalty="l1", solver="liblinear", C=C,
                                  max_iter=10000, random_state=RNG)
    return LogisticRegression(penalty="elasticnet", solver="saga", C=C,
                              l1_ratio=l1_ratio, max_iter=20000,
                              random_state=RNG, tol=1e-4)


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


def reconstruct_ext_metadata():
    GSM_MAP = {
        "ICH1": {"gsm": "GSM3576253", "patient": "Patient1", "draw": "first", "time": "<24h"},
        "ICH2": {"gsm": "GSM3576254", "patient": "Patient1", "draw": "second", "time": "72h"},
        "ICH3": {"gsm": "GSM3576255", "patient": "Patient2", "draw": "first", "time": "<24h"},
        "ICH4": {"gsm": "GSM3576256", "patient": "Patient2", "draw": "second", "time": "72h"},
        "ICH5": {"gsm": "GSM3576257", "patient": "Patient3", "draw": "first", "time": "<24h"},
        "ICH6": {"gsm": "GSM3576258", "patient": "Patient3", "draw": "second", "time": "72h"},
        "ICH7": {"gsm": "GSM3576259", "patient": "Patient4", "draw": "first", "time": "<24h"},
        "ICH8": {"gsm": "GSM3576260", "patient": "Patient4", "draw": "second", "time": "72h"},
        "ICH9": {"gsm": "GSM3576261", "patient": "Patient5", "draw": "first", "time": "<24h"},
        "ICH10": {"gsm": "GSM3576262", "patient": "Patient5", "draw": "second", "time": "72h"},
        "ICH11": {"gsm": "GSM3576263", "patient": "Patient6", "draw": "first", "time": "<24h"},
        "ICH12": {"gsm": "GSM3576264", "patient": "Patient6", "draw": "second", "time": "72h"},
        "ICH13": {"gsm": "GSM3576265", "patient": "Patient7", "draw": "first", "time": "<24h"},
        "ICH14": {"gsm": "GSM3576266", "patient": "Patient7", "draw": "second", "time": "72h"},
        "ICH15": {"gsm": "GSM3576267", "patient": "Patient8", "draw": "first", "time": "<24h"},
        "ICH16": {"gsm": "GSM3576268", "patient": "Patient8", "draw": "second", "time": "72h"},
        "ICH17": {"gsm": "GSM3576269", "patient": "Patient9", "draw": "first", "time": "<24h"},
        "ICH18": {"gsm": "GSM3576270", "patient": "Patient9", "draw": "second", "time": "72h"},
        "ICH19": {"gsm": "GSM3576271", "patient": "Patient10", "draw": "first", "time": "<24h"},
        "ICH20": {"gsm": "GSM3576272", "patient": "Patient10", "draw": "second", "time": "72h"},
        "ICH21": {"gsm": "GSM3576273", "patient": "Patient11", "draw": "first", "time": "<24h"},
        "ICH22": {"gsm": "GSM3576274", "patient": "Patient11", "draw": "second", "time": "72h"},
    }
    rows = []
    for sid, info in GSM_MAP.items():
        rows.append({
            "sample_id": sid, "gsm": info["gsm"], "patient": info["patient"],
            "draw": info["draw"], "time_label": info["time"],
        })
    return pd.DataFrame(rows)


def run_paired_validation(name, gene_list, expr_int, md_int, expr_ext, ext_meta):
    log.info("=" * 70)
    log.info(f"[{name}] Paired temporal validation on GSE125512")

    sel = md_int[md_int["stage"].isin([POS_CLASS, NEG_CLASS])].copy()
    sel = sel[sel["sample_id"].isin(expr_int.columns)]
    sel["patient"] = sel["sample_id"].str.split("_").str[0]
    common = [g for g in gene_list if g in expr_int.index and g in expr_ext.index]
    if len(common) < 2:
        log.warning(f"[{name}] Only {len(common)} common genes, skipping")
        return None
    log.info(f"[{name}] Common genes: {len(common)}/{len(gene_list)}")

    X_train = expr_int.loc[common, sel["sample_id"].tolist()].T
    y_train = (sel["stage"].values == POS_CLASS).astype(int)
    groups = sel["patient"].values

    n_splits = min(5, len(set(groups[y_train == 0])), len(set(groups[y_train == 1])))
    grouped = StratifiedGroupKFold(n_splits=max(2, n_splits))

    best = (None, None, -1)
    for l1r in L1_RATIOS:
        for C in CS_GRID:
            try:
                pipe = make_pipeline(StandardScaler(), make_clf(C, l1r))
                oof = cross_val_predict(pipe, X_train.values, y_train, cv=grouped,
                                        groups=groups, method="predict_proba")[:, 1]
                auc = roc_auc_score(y_train, oof)
                if auc > best[2]:
                    best = (C, l1r, auc)
            except Exception:
                continue
    C_star, l1_star, auc_int = best
    if C_star is None:
        log.error(f"[{name}] Grid search failed")
        return None
    log.info(f"[{name}] Internal grouped CV AUC={auc_int:.4f}")

    pipe = make_pipeline(StandardScaler(), make_clf(C_star, l1_star))
    pipe.fit(X_train.values, y_train)
    scaler_int = pipe.named_steps["standardscaler"]
    lr = pipe.named_steps["logisticregression"]
    coef = lr.coef_.ravel()
    intercept = lr.intercept_[0]

    ext_first = ext_meta[ext_meta["draw"] == "first"].sort_values("patient")
    ext_second = ext_meta[ext_meta["draw"] == "second"].sort_values("patient")

    first_ids = ext_first["sample_id"].tolist()
    second_ids = ext_second["sample_id"].tolist()
    patients = ext_first["patient"].tolist()

    first_present = [s for s in first_ids if s in expr_ext.columns]
    second_present = [s for s in second_ids if s in expr_ext.columns]

    X_first_raw = expr_ext.loc[common, first_present].T
    X_second_raw = expr_ext.loc[common, second_present].T

    scaler_ext = StandardScaler()
    all_ext_ids = first_present + second_present
    X_all_ext_raw = expr_ext.loc[common, all_ext_ids].T
    scaler_ext.fit(X_all_ext_raw.values)

    X_first_z = scaler_ext.transform(X_first_raw.values)
    X_second_z = scaler_ext.transform(X_second_raw.values)

    def score_z(Z):
        logits = Z @ coef + intercept
        return 1 / (1 + np.exp(-np.clip(logits, -500, 500)))

    scores_first = score_z(X_first_z)
    scores_second = score_z(X_second_z)

    log.info(f"[{name}] External scores:")
    log.info(f"  First draw (<24h):  mean={scores_first.mean():.4f} std={scores_first.std():.4f}")
    log.info(f"  Second draw (72h):  mean={scores_second.mean():.4f} std={scores_second.std():.4f}")

    n_pairs = min(len(scores_first), len(scores_second))
    diff = scores_second[:n_pairs] - scores_first[:n_pairs]
    n_increase = int((diff > 0).sum())
    n_decrease = int((diff < 0).sum())
    n_tie = int((diff == 0).sum())

    if n_pairs >= 3:
        stat, p_wilcox = wilcoxon(scores_second[:n_pairs], scores_first[:n_pairs], alternative="greater")
    else:
        stat, p_wilcox = np.nan, np.nan

    log.info(f"[{name}] Paired test (n={n_pairs}):")
    log.info(f"  Increase (second > first): {n_increase}/{n_pairs}")
    log.info(f"  Wilcoxon signed-rank (one-sided, second>first): stat={stat}, p={p_wilcox:.6f}")

    y_ext = np.concatenate([np.zeros(len(scores_first)), np.ones(len(scores_second))])
    prob_ext = np.concatenate([scores_first, scores_second])
    if len(np.unique(y_ext)) == 2:
        ext_auc = roc_auc_score(y_ext, prob_ext)
        ext_ci_lo, ext_ci_hi, _ = bootstrap_auc_ci(y_ext, prob_ext)
        fpr_e, tpr_e, _ = roc_curve(y_ext, prob_ext)
    else:
        ext_auc = 0.5
        ext_ci_lo, ext_ci_hi = 0.0, 1.0
        fpr_e, tpr_e = np.array([0, 1]), np.array([0, 1])
    log.info(f"[{name}] External first-vs-second AUC={ext_auc:.4f} [{ext_ci_lo}, {ext_ci_hi}]")

    ctrl = md_int[md_int["stage"] == "control"].copy()
    ctrl = ctrl[ctrl["sample_id"].isin(expr_int.columns)]
    X_ctrl = expr_int.loc[common, ctrl["sample_id"].tolist()].T
    score_ctrl = pipe.predict_proba(X_ctrl.values)[:, 1]
    score_acute_int = pipe.predict_proba(X_train.values[y_train == 0])[:, 1]
    score_sub_int = pipe.predict_proba(X_train.values[y_train == 1])[:, 1]

    log.info(f"[{name}] Full score ordering:")
    log.info(f"  control(int):     {score_ctrl.mean():.4f}")
    log.info(f"  acute(int,72h):   {score_acute_int.mean():.4f}")
    log.info(f"  subacute(int,7d): {score_sub_int.mean():.4f}")
    log.info(f"  first(ext,<24h):  {scores_first.mean():.4f}")
    log.info(f"  second(ext,72h):  {scores_second.mean():.4f}")

    prefix = name.lower().replace(" ", "_").replace("-", "_")
    pd.DataFrame({"fpr": fpr_e, "tpr": tpr_e}).to_csv(
        OUT / f"roc_ext_paired_{prefix}.csv", index=False)

    pair_df = pd.DataFrame({
        "patient": patients[:n_pairs],
        "score_first_24h": scores_first[:n_pairs],
        "score_second_72h": scores_second[:n_pairs],
        "delta": diff,
        "direction": ["increase" if d > 0 else "decrease" for d in diff],
    })
    pair_df.to_csv(OUT / f"paired_scores_{prefix}.csv", index=False)

    score_all = pd.DataFrame({
        "sample_id": (list(ctrl["sample_id"].values) +
                      list(sel.loc[y_train == 0, "sample_id"].values) +
                      list(sel.loc[y_train == 1, "sample_id"].values) +
                      first_present + second_present),
        "group": (["control_int"] * len(score_ctrl) +
                  ["acute_int_72h"] * len(score_acute_int) +
                  ["subacute_int_7d"] * len(score_sub_int) +
                  ["first_ext_24h"] * len(scores_first) +
                  ["second_ext_72h"] * len(scores_second)),
        "score": np.concatenate([score_ctrl, score_acute_int, score_sub_int,
                                 scores_first, scores_second]),
    })
    score_all.to_csv(OUT / f"temporal_scores_{prefix}.csv", index=False)

    oof = cross_val_predict(make_pipeline(StandardScaler(), make_clf(C_star, l1_star)),
                            X_train.values, y_train, cv=grouped,
                            groups=groups, method="predict_proba")[:, 1]
    int_ci_lo, int_ci_hi, _ = bootstrap_auc_ci(y_train, oof)

    return {
        "name": name,
        "n_genes": len(common),
        "internal": {
            "grouped_cv_auc": round(auc_int, 4),
            "ci_lo": int_ci_lo, "ci_hi": int_ci_hi,
        },
        "external_paired": {
            "task": "first_draw_vs_second_draw",
            "n_patients": n_pairs,
            "n_first": len(scores_first),
            "n_second": len(scores_second),
            "score_first_mean": round(float(scores_first.mean()), 4),
            "score_second_mean": round(float(scores_second.mean()), 4),
            "n_increase": n_increase,
            "n_decrease": n_decrease,
            "wilcoxon_stat": round(float(stat), 2) if not np.isnan(stat) else None,
            "wilcoxon_p": round(float(p_wilcox), 6) if not np.isnan(p_wilcox) else None,
            "auc_first_vs_second": round(ext_auc, 4),
            "auc_ci_lo": ext_ci_lo,
            "auc_ci_hi": ext_ci_hi,
        },
        "score_ordering": {
            "control_int": round(float(score_ctrl.mean()), 4),
            "acute_int_72h": round(float(score_acute_int.mean()), 4),
            "subacute_int_7d": round(float(score_sub_int.mean()), 4),
            "first_ext_24h": round(float(scores_first.mean()), 4),
            "second_ext_72h": round(float(scores_second.mean()), 4),
        },
    }


def main():
    log.info("=" * 70)
    log.info(f"Module 7 — Paired Temporal Validation (GSE125512) | {datetime.now()}")
    log.info("=" * 70)

    expr_int = load_expr(M1 / "expr_mrna.csv")
    md_int = pd.read_csv(M1 / "metadata.csv")
    md_int["sample_id"] = md_int["sample_id"].astype(str).str.strip().str.strip('"')
    expr_ext = load_expr(M1 / "expr_mrna_external.csv")
    ext_meta = reconstruct_ext_metadata()
    log.info(f"Reconstructed GSE125512 metadata: {len(ext_meta)} samples, "
             f"{ext_meta['patient'].nunique()} patients")
    log.info(f"  First draws: {(ext_meta['draw']=='first').sum()}")
    log.info(f"  Second draws: {(ext_meta['draw']=='second').sum()}")
    ext_meta.to_csv(OUT / "GSE125512_reconstructed_metadata.csv", index=False)

    results = {}

    r1 = run_paired_validation("TrackB_13gene", TRACK_B_PORTABLE,
                               expr_int, md_int, expr_ext, ext_meta)
    if r1:
        results["trackB_13gene"] = r1

    track_a_common = [g for g in TRACK_A_GENES if g in expr_ext.index]
    if len(track_a_common) >= 2:
        r2 = run_paired_validation("TrackA_portable", track_a_common,
                                   expr_int, md_int, expr_ext, ext_meta)
        if r2:
            results["trackA_portable"] = r2

    log.info("=" * 70)
    log.info("[FINAL RESULTS]")
    for k, r in results.items():
        ep = r["external_paired"]
        so = r["score_ordering"]
        log.info(f"\n  {k} ({r['n_genes']} genes):")
        log.info(f"    Internal: AUC={r['internal']['grouped_cv_auc']:.4f} "
                 f"[{r['internal']['ci_lo']}, {r['internal']['ci_hi']}]")
        log.info(f"    External paired (first<24h vs second 72h):")
        log.info(f"      Score first={ep['score_first_mean']:.4f} second={ep['score_second_mean']:.4f}")
        log.info(f"      Increase: {ep['n_increase']}/{ep['n_patients']} patients")
        log.info(f"      Wilcoxon p={ep['wilcoxon_p']}")
        log.info(f"      AUC(first vs second)={ep['auc_first_vs_second']:.4f} "
                 f"[{ep['auc_ci_lo']}, {ep['auc_ci_hi']}]")
        log.info(f"    Score ordering: ctrl={so['control_int']:.3f} → "
                 f"first_ext={so['first_ext_24h']:.3f} → "
                 f"acute_int={so['acute_int_72h']:.3f} ≈ "
                 f"second_ext={so['second_ext_72h']:.3f} → "
                 f"sub_int={so['subacute_int_7d']:.3f}")

    (OUT / "paired_temporal_validation.json").write_text(
        json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    log.info(f"\n[DONE] -> {OUT}/paired_temporal_validation.json")


if __name__ == "__main__":
    main()
