import json, logging, sys, re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

M1 = Path(__file__).resolve().parent.parent.parent / "模块1_预检" / "output"
M4 = Path(__file__).resolve().parent.parent
OUT = M4 / "output"
OUT.mkdir(parents=True, exist_ok=True)

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s",
                    handlers=[logging.StreamHandler(sys.stdout)])
log = logging.getLogger("m4.step2b")

MIMAT_TO_MIRBASE = {
    "MIMAT0000068": "hsa-miR-15a-5p",
    "MIMAT0004488": "hsa-miR-15a-3p",
    "MIMAT0000417": "hsa-miR-15b-5p",
    "MIMAT0000069": "hsa-miR-16-5p",
    "MIMAT0000461": "hsa-miR-195-5p",
    "MIMAT0002820": "hsa-miR-497-5p",
    "MIMAT0000070": "hsa-miR-17-5p",
    "MIMAT0000076": "hsa-miR-21-5p",
    "MIMAT0004494": "hsa-miR-21-3p",
    "MIMAT0000080": "hsa-miR-24-3p",
    "MIMAT0000084": "hsa-miR-27a-3p",
    "MIMAT0004501": "hsa-miR-27a-5p",
    "MIMAT0000086": "hsa-miR-29a-3p",
    "MIMAT0000100": "hsa-miR-29b-3p",
    "MIMAT0000681": "hsa-miR-29c-3p",
    "MIMAT0000092": "hsa-miR-92a-3p",
    "MIMAT0000425": "hsa-miR-130a-3p",
    "MIMAT0004593": "hsa-miR-130a-5p",
    "MIMAT0000691": "hsa-miR-130b-3p",
    "MIMAT0000433": "hsa-miR-142-5p",
    "MIMAT0000434": "hsa-miR-142-3p",
    "MIMAT0000436": "hsa-miR-144-3p",
    "MIMAT0000073": "hsa-miR-19a-3p",
    "MIMAT0000074": "hsa-miR-19b-3p",
    "MIMAT0000242": "hsa-miR-129-5p",
    "MIMAT0003393": "hsa-miR-425-5p",
    "MIMAT0000256": "hsa-miR-181a-5p",
    "MIMAT0000270": "hsa-miR-181b-5p",
    "MIMAT0000258": "hsa-miR-181c-5p",
    "MIMAT0002821": "hsa-miR-181d-5p",
    "MIMAT0000271": "hsa-miR-214-3p",
    "MIMAT0000426": "hsa-miR-132-3p",
    "MIMAT0000269": "hsa-miR-212-3p",
    "MIMAT0000427": "hsa-miR-133a-3p",
    "MIMAT0000770": "hsa-miR-133b",
    "MIMAT0000439": "hsa-miR-153-3p",
    "MIMAT0003258": "hsa-miR-590-3p",
    "MIMAT0000456": "hsa-miR-186-5p",
    "MIMAT0000762": "hsa-miR-324-3p",
    "MIMAT0000761": "hsa-miR-324-5p",
    "MIMAT0000765": "hsa-miR-335-5p",
    "MIMAT0004703": "hsa-miR-335-3p",
    "MIMAT0000753": "hsa-miR-342-3p",
    "MIMAT0004694": "hsa-miR-342-5p",
    "MIMAT0001340": "hsa-miR-423-3p",
    "MIMAT0004748": "hsa-miR-423-5p",
    "MIMAT0000646": "hsa-miR-155-5p",
    "MIMAT0000451": "hsa-miR-150-5p",
    "MIMAT0000440": "hsa-miR-191-5p",
    "MIMAT0000275": "hsa-miR-218-5p",
    "MIMAT0000690": "hsa-miR-296-5p",
    "MIMAT0003220": "hsa-miR-556-3p",
    "MIMAT0000435": "hsa-miR-143-3p",
    "MIMAT0000437": "hsa-miR-145-5p",
    "MIMAT0000750": "hsa-miR-340-3p",
    "MIMAT0000066": "hsa-let-7e-5p",
    "MIMAT0000063": "hsa-let-7b-5p",
    "MIMAT0000064": "hsa-let-7c-5p",
    "MIMAT0000062": "hsa-let-7a-5p",
    "MIMAT0000067": "hsa-let-7f-5p",
    "MIMAT0000414": "hsa-let-7g-5p",
    "MIMAT0000415": "hsa-let-7i-5p",
    "MIMAT0000415": "hsa-let-7d-5p",
    "MIMAT0000104": "hsa-miR-103a-3p",
    "MIMAT0000101": "hsa-miR-103b",
    "MIMAT0000250": "hsa-miR-139-5p",
    "MIMAT0000422": "hsa-miR-124-3p",
    "MIMAT0004185": "hsa-miR-449a",
    "MIMAT0000082": "hsa-miR-26a-5p",
    "MIMAT0000083": "hsa-miR-26b-5p",
    "MIMAT0000089": "hsa-miR-31-5p",
    "MIMAT0000097": "hsa-miR-99a-5p",
    "MIMAT0000689": "hsa-miR-99b-5p",
    "MIMAT0000098": "hsa-miR-100-5p",
    "MIMAT0000099": "hsa-miR-101-3p",
    "MIMAT0000094": "hsa-miR-95-3p",
    "MIMAT0000093": "hsa-miR-106a-5p",
    "MIMAT0004951": "hsa-miR-887-3p",
    "MIMAT0000077": "hsa-miR-22-3p",
    "MIMAT0000430": "hsa-miR-138-5p",
    "MIMAT0000764": "hsa-miR-339-5p",
    "MIMAT0004692": "hsa-miR-340-5p",
    "MIMAT0000728": "hsa-miR-375-3p",
    "MIMAT0000752": "hsa-miR-328-3p",
}

TARGETED_PAIRS = [
    ("ACSL1", "hsa-miR-142-3p", "ENCORI", 6),
    ("ACSL1", "hsa-miR-556-3p", "ENCORI", 6),
    ("ACSL1", "hsa-miR-19a-3p", "ENCORI", 4),
    ("ACSL1", "hsa-miR-19b-3p", "ENCORI", 4),
    ("ACSL1", "hsa-miR-129-5p", "ENCORI", 4),
    ("ACSL1", "hsa-miR-425-5p", "ENCORI", 4),
    ("ACSL1", "hsa-miR-15a-5p", "ENCORI", 3),
    ("ACSL1", "hsa-miR-16-5p", "ENCORI", 3),
    ("ACSL1", "hsa-miR-424-5p", "ENCORI", 3),
    ("ACSL1", "hsa-miR-181a-5p", "ENCORI", 2),
    ("ACSL1", "hsa-miR-181b-5p", "ENCORI", 2),
    ("ACSL1", "hsa-miR-181c-5p", "ENCORI", 2),
    ("ACSL1", "hsa-miR-130a-3p", "ENCORI+TargetScan", 2),
    ("ACSL1", "hsa-miR-130b-3p", "ENCORI+TargetScan", 2),
    ("ACSL1", "hsa-miR-195-5p", "ENCORI", 2),
    ("ACSL1", "hsa-miR-497-5p", "ENCORI", 2),
    ("ACSL1", "hsa-miR-340-5p", "ENCORI", 2),
    ("ACSL1", "hsa-miR-339-5p", "ENCORI", 2),
    ("ACSL1", "hsa-miR-449a", "Literature(AMI ferroptosis)", 0),
    ("ACSL1", "hsa-miR-22-3p", "ENCORI", 1),
    ("ACSL1", "hsa-miR-26a-5p", "ENCORI", 1),
    ("ACSL1", "hsa-miR-26b-5p", "ENCORI", 1),
    ("ACSL1", "hsa-miR-103a-3p", "ENCORI", 1),
    ("ACSL1", "hsa-miR-186-5p", "ENCORI", 1),
    ("ACSL1", "hsa-miR-375-3p", "ENCORI", 1),
    ("GSTP1", "hsa-miR-133a-3p", "Literature(miRTarBase)", 0),
    ("GSTP1", "hsa-miR-133b", "Literature(NSCLC)", 0),
    ("GSTP1", "hsa-miR-144-3p", "Literature(prostate)", 0),
    ("GSTP1", "hsa-miR-153-3p", "Literature(prostate)", 0),
    ("GSTP1", "hsa-miR-130b-3p", "Literature(ovarian)", 0),
    ("GSTP1", "hsa-miR-590-3p", "Literature(prostate)", 0),
    ("GSTP1", "hsa-miR-186-5p", "Literature(ovarian)", 0),
    ("GPX4", "hsa-let-7b-5p", "multiMiR_validated", 0),
    ("GPX4", "hsa-let-7c-5p", "multiMiR_validated", 0),
    ("GPX4", "hsa-miR-15a-5p", "multiMiR_validated", 0),
    ("GPX4", "hsa-miR-214-3p", "multiMiR_validated", 0),
    ("GPX4", "hsa-miR-181a-5p", "multiMiR_validated", 0),
    ("GPX4", "hsa-miR-154-5p", "multiMiR_validated", 0),
    ("GPX4", "hsa-miR-342-5p", "multiMiR_validated", 0),
    ("GPX4", "hsa-miR-423-3p", "multiMiR_validated", 0),
    ("GPX4", "hsa-miR-124-3p", "multiMiR_validated", 0),
    ("HMOX1", "hsa-miR-24-3p", "multiMiR_validated", 0),
    ("HMOX1", "hsa-miR-155-5p", "multiMiR_validated", 0),
    ("FTH1", "hsa-miR-92a-3p", "multiMiR_validated", 0),
    ("FTH1", "hsa-miR-335-3p", "multiMiR_validated", 0),
    ("FTH1", "hsa-miR-101-3p", "multiMiR_validated", 0),
    ("FTH1", "hsa-miR-218-5p", "multiMiR_validated", 0),
    ("FTH1", "hsa-miR-324-3p", "multiMiR_validated", 0),
    ("FTH1", "hsa-miR-130a-3p", "multiMiR_validated", 0),
    ("SLC7A11", "hsa-miR-130a-3p", "ENCORI+multiMiR", 0),
    ("SLC7A11", "hsa-miR-29a-3p", "multiMiR_validated", 0),
    ("SLC7A11", "hsa-miR-21-5p", "multiMiR_validated", 0),
    ("SLC7A11", "hsa-miR-27a-3p", "multiMiR_validated", 0),
    ("SLC7A11", "hsa-miR-26b-5p", "multiMiR_validated", 0),
]


def load_expr(path):
    df = pd.read_csv(path, index_col=0)
    df.columns = [c.strip().strip('"') for c in df.columns]
    df.index = [str(i).strip().strip('"') for i in df.index]
    return df


def normalize_sid(s):
    return s.strip().strip('"').replace("_base", "_bas").replace("_7day", "_7da")


def build_mimat_to_probe(expr_mirna):
    probe_ids = list(expr_mirna.index)
    mirbase_to_probe = {}
    for mirbase_name in set(MIMAT_TO_MIRBASE.values()):
        name = mirbase_name.lower().strip()
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
        if matched:
            mirbase_to_probe[mirbase_name] = sorted(set(matched))
    return mirbase_to_probe


def main():
    log.info("Module 4 Step 2b: Targeted miRNA-mRNA correlation (revised hub genes)")

    expr_mrna = load_expr(M1 / "expr_mrna.csv")
    expr_mirna = load_expr(M1 / "expr_mirna.csv")

    mrna_norm = {normalize_sid(c): c for c in expr_mrna.columns}
    mirna_norm = {normalize_sid(c): c for c in expr_mirna.columns}
    common = sorted(set(mrna_norm.keys()) & set(mirna_norm.keys()))
    log.info(f"Paired samples: {len(common)}")

    mrna_cols = [mrna_norm[s] for s in common]
    mirna_cols = [mirna_norm[s] for s in common]

    mirbase_to_probe = build_mimat_to_probe(expr_mirna)
    log.info(f"miRBase -> probe mapping: {len(mirbase_to_probe)} miRNAs mapped")
    for k, v in sorted(mirbase_to_probe.items())[:10]:
        log.info(f"  {k} -> {v}")

    results = []
    for target, mirna_id, source, clip_n in TARGETED_PAIRS:
        if target not in expr_mrna.index:
            continue
        if mirna_id not in mirbase_to_probe:
            results.append({
                "target": target, "mirna": mirna_id, "source": source, "clip_n": clip_n,
                "probe": "NOT_MAPPED", "rho": np.nan, "p_value": np.nan, "n": 0,
                "sig_negative": False
            })
            continue

        mrna_vals = expr_mrna.loc[target, mrna_cols].values.astype(float)
        for probe in mirbase_to_probe[mirna_id]:
            mirna_vals = expr_mirna.loc[probe, mirna_cols].values.astype(float)
            valid = ~(np.isnan(mrna_vals) | np.isnan(mirna_vals))
            n_valid = int(valid.sum())
            if n_valid < 10:
                continue
            std_check = mirna_vals[valid].std()
            if std_check == 0:
                continue
            rho, pval = spearmanr(mrna_vals[valid], mirna_vals[valid])
            results.append({
                "target": target, "mirna": mirna_id, "source": source, "clip_n": clip_n,
                "probe": probe, "rho": round(rho, 4), "p_value": pval, "n": n_valid,
                "sig_negative": rho < -0.3 and pval < 0.05
            })

    df = pd.DataFrame(results)
    df.to_csv(OUT / "targeted_mirna_correlation.csv", index=False)

    mapped = df[df["probe"] != "NOT_MAPPED"]
    unmapped = df[df["probe"] == "NOT_MAPPED"]
    sig_neg = mapped[mapped["sig_negative"]].sort_values("rho")

    log.info(f"\n{'='*70}")
    log.info(f"Total targeted pairs: {len(TARGETED_PAIRS)}")
    log.info(f"Mapped to probes: {len(mapped)} (unmapped: {len(unmapped)})")
    log.info(f"Significant negative (rho<-0.3, p<0.05): {len(sig_neg)}")

    log.info(f"\n=== Significant negative correlations ===")
    for _, r in sig_neg.iterrows():
        log.info(f"  {r['target']:8s} <- {r['mirna']:25s} rho={r['rho']:+.4f} p={r['p_value']:.4g} n={r['n']} [{r['source']}]")

    log.info(f"\n=== Per-gene summary ===")
    for gene in ["ALOX15", "ACSL1", "GSTP1", "HMOX1", "GPX4", "FTH1", "SLC7A11"]:
        sub = mapped[mapped["target"] == gene]
        if len(sub) == 0:
            sub_all = df[df["target"] == gene]
            if len(sub_all) > 0:
                log.info(f"  {gene:8s}: 0 tested (ALOX15 lacks miRNA regulation; uses DICE/hnRNP translational control)")
            continue
        n_sig = sub["sig_negative"].sum()
        best = sub["rho"].min()
        log.info(f"  {gene:8s}: tested={len(sub):3d} sig_neg={n_sig:3d} best_rho={best:+.4f}")

    step1_sig = pd.DataFrame()
    if (OUT / "mirna_mrna_sig_negative.csv").exists():
        step1_sig = pd.read_csv(OUT / "mirna_mrna_sig_negative.csv").rename(
            columns={"mirna_mirbase": "mirna", "mirna_probe": "probe"})
        step1_sig["source"] = "M4_step1"
    all_sig = pd.concat([sig_neg, step1_sig], ignore_index=True)
    if len(all_sig) > 0:
        all_sig = all_sig.drop_duplicates(subset=["target", "mirna", "probe"]).sort_values("rho")
        all_sig.to_csv(OUT / "all_sig_negative_pairs.csv", index=False)
        log.info(f"\n=== Combined significant negative pairs (old + new): {len(all_sig)} ===")
        for _, r in all_sig.iterrows():
            src = r.get("source", "M4_step1")
            log.info(f"  {str(r['target']):8s} <- {str(r['mirna']):25s} rho={r['rho']:+.4f} [{src}]")

    summary = {
        "n_targeted_pairs": len(TARGETED_PAIRS),
        "n_mapped": len(mapped),
        "n_unmapped": len(unmapped),
        "n_sig_negative": len(sig_neg),
        "n_paired_samples": len(common),
        "ALOX15_note": "No miRNA regulation in ENCORI/miRDB/literature; regulated by DICE/hnRNP K/E1 translational control at 3'UTR",
        "per_gene": {},
    }
    for gene in ["ACSL1", "GSTP1", "HMOX1", "GPX4", "FTH1", "SLC7A11"]:
        sub = mapped[mapped["target"] == gene]
        summary["per_gene"][gene] = {
            "tested": len(sub),
            "sig_negative": int(sub["sig_negative"].sum()),
            "best_rho": round(float(sub["rho"].min()), 4) if len(sub) > 0 else None,
        }
    (OUT / "targeted_correlation_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    log.info(f"\n[DONE] -> {OUT}/targeted_*.csv, all_sig_negative_pairs.csv")


if __name__ == "__main__":
    main()
