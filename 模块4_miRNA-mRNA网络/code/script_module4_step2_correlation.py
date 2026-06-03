import json, logging, sys, re
from pathlib import Path
from itertools import product

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

M1 = Path(__file__).resolve().parent.parent.parent / "模块1_预检" / "output"
M4 = Path(__file__).resolve().parent.parent
OUT = M4 / "output"
OUT.mkdir(parents=True, exist_ok=True)

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s",
                    handlers=[logging.StreamHandler(sys.stdout)])
log = logging.getLogger("m4.step2")

RHO_THRESH = -0.3
P_THRESH = 0.05


def load_expr(path):
    df = pd.read_csv(path, index_col=0)
    df.columns = [c.strip().strip('"') for c in df.columns]
    df.index = [str(i).strip().strip('"') for i in df.index]
    return df


def normalize_sample_id(sid):
    sid = sid.strip().strip('"')
    return (sid.replace("_base", "_bas").replace("_7day", "_7da")
              .replace("_Base", "_bas").replace("_7Day", "_7da"))


def mirbase_to_probe_candidates(mirbase_id, probe_ids):
    name = mirbase_id.lower().strip()
    let_m = re.match(r"hsa-let-7([a-z]?)-([35]p)$", name)
    mir_m = re.match(r"hsa-mir-(\d+)([a-z]?)-([35]p)$", name)
    mir_na = re.match(r"hsa-mir-(\d+)([a-z]?)$", name)
    matched = []
    for pid in probe_ids:
        for part in pid.lower().replace("*", "").split("/"):
            part = part.strip()
            hit = False
            if let_m:
                letter, arm = let_m.group(1), let_m.group(2)
                pm = re.match(r"hsa-let-7-p\d+([a-z])\d*_([35]p)$", part)
                if pm and pm.group(2) == arm and letter and pm.group(1) == letter:
                    hit = True
            elif mir_m:
                num, letter, arm = mir_m.group(1), mir_m.group(2), mir_m.group(3)
                pm = re.match(rf"hsa-mir-{num}(?:-p\d+([a-z]?)\d*)?_([35]p)$", part)
                if pm and pm.group(2) == arm:
                    plet = pm.group(1) or ""
                    if (letter == "" and plet == "") or (letter != "" and plet == letter):
                        hit = True
            elif mir_na:
                num, letter = mir_na.group(1), mir_na.group(2)
                pm = re.match(rf"hsa-mir-{num}(?:-p\d+([a-z]?)\d*)?_([35]p)$", part)
                if pm:
                    plet = pm.group(1) or ""
                    if (letter == "" and plet == "") or (letter != "" and plet == letter):
                        hit = True
            if hit:
                matched.append(pid)
                break
    return sorted(set(matched))


def main():
    log.info("Module 4 Step 2: miRNA-mRNA Spearman expression correlation")

    expr_mrna = load_expr(M1 / "expr_mrna.csv")
    expr_mirna = load_expr(M1 / "expr_mirna.csv")

    mrna_norm = {normalize_sample_id(c): c for c in expr_mrna.columns}
    mirna_norm = {normalize_sample_id(c): c for c in expr_mirna.columns}
    common_norm = set(mrna_norm.keys()) & set(mirna_norm.keys())
    log.info(f"Sample ID normalization: mRNA={len(mrna_norm)} miRNA={len(mirna_norm)} paired={len(common_norm)}")

    md = pd.read_csv(M1 / "metadata.csv")
    md["sample_id"] = md["sample_id"].astype(str).str.strip().str.strip('"')
    md["sid_norm"] = md["sample_id"].apply(normalize_sample_id)
    md_paired = md[md["sid_norm"].isin(common_norm)]
    log.info(f"Paired by stage: {md_paired['stage'].value_counts().to_dict()}")

    consensus = pd.read_csv(M4 / "output" / "hub_mirna_candidates_consensus.csv")
    log.info(f"Consensus pairs: {len(consensus)}")

    double_evidence = consensus[
        (consensus["in_validated"] == 1) & (consensus["in_predicted_ge2"] == 1)
    ]
    log.info(f"Double-evidence pairs: {len(double_evidence)}")

    probe_ids = list(expr_mirna.index)
    unique_mirnas = set(consensus["mirna"].astype(str))
    log.info(f"Unique consensus miRNAs: {len(unique_mirnas)}")

    mirna_map = {}
    mapped_count = 0
    for mirbase_id in unique_mirnas:
        probes = mirbase_to_probe_candidates(mirbase_id, probe_ids)
        if probes:
            mirna_map[mirbase_id] = probes
            mapped_count += 1
    log.info(f"miRBase -> probe mapping: {mapped_count}/{len(unique_mirnas)} mapped")

    paired_list = sorted(common_norm)
    mrna_cols = [mrna_norm[s] for s in paired_list]
    mirna_cols = [mirna_norm[s] for s in paired_list]

    results = []
    n_tested = 0
    n_sig_neg = 0

    for _, row in consensus.iterrows():
        target = str(row["target"])
        mirna_id = str(row["mirna"])
        in_val = int(row["in_validated"])
        in_pred = int(row["in_predicted_ge2"])

        if target not in expr_mrna.index:
            continue
        if mirna_id not in mirna_map:
            continue

        mrna_vals = expr_mrna.loc[target, mrna_cols].values.astype(float)

        for probe in mirna_map[mirna_id]:
            mirna_vals = expr_mirna.loc[probe, mirna_cols].values.astype(float)

            valid = ~(np.isnan(mrna_vals) | np.isnan(mirna_vals))
            if valid.sum() < 10:
                continue

            rho, pval = spearmanr(mrna_vals[valid], mirna_vals[valid])
            n_tested += 1

            results.append({
                "target": target,
                "mirna_mirbase": mirna_id,
                "mirna_probe": probe,
                "rho": round(rho, 4),
                "p_value": pval,
                "n_samples": int(valid.sum()),
                "in_validated": in_val,
                "in_predicted_ge2": in_pred,
                "sig_negative": rho < RHO_THRESH and pval < P_THRESH,
            })
            if rho < RHO_THRESH and pval < P_THRESH:
                n_sig_neg += 1

    log.info(f"Tested: {n_tested} pairs")
    log.info(f"Significant negative (rho<{RHO_THRESH}, p<{P_THRESH}): {n_sig_neg}")

    df = pd.DataFrame(results)
    df.to_csv(OUT / "mirna_mrna_correlation_all.csv", index=False)

    if len(df) > 0:
        sig_neg = df[df["sig_negative"]].sort_values("rho")
        sig_neg.to_csv(OUT / "mirna_mrna_sig_negative.csv", index=False)
        log.info(f"\n=== Significant negative pairs (rho<{RHO_THRESH}, p<{P_THRESH}) ===")
        for _, r in sig_neg.iterrows():
            tag = "DOUBLE" if r["in_validated"] and r["in_predicted_ge2"] else "single"
            log.info(f"  {r['target']:8s} <- {r['mirna_mirbase']:25s} rho={r['rho']:+.3f} p={r['p_value']:.4g} n={r['n_samples']} [{tag}]")

        per_gene = df.groupby("target").agg(
            n_tested=("rho", "count"),
            n_sig_neg=("sig_negative", "sum"),
            best_rho=("rho", "min"),
        ).reset_index()
        log.info(f"\n=== Per-gene summary ===")
        for _, r in per_gene.iterrows():
            log.info(f"  {r['target']:8s}: tested={r['n_tested']:3d} sig_neg={r['n_sig_neg']:3d} best_rho={r['best_rho']:+.3f}")

        summary = {
            "n_consensus_input": len(consensus),
            "n_mirnas_mapped": mapped_count,
            "n_mirnas_total": len(unique_mirnas),
            "n_paired_samples": len(paired_list),
            "n_pairs_tested": n_tested,
            "n_sig_negative": n_sig_neg,
            "rho_threshold": RHO_THRESH,
            "p_threshold": P_THRESH,
            "per_gene": per_gene.to_dict("records"),
        }
        (OUT / "correlation_summary.json").write_text(
            json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

    log.info("[DONE]")


if __name__ == "__main__":
    main()
