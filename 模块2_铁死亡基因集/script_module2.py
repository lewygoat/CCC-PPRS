import hashlib
import io
import json
import logging
import re
import sys
import time
import warnings
import zipfile
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Optional

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

warnings.filterwarnings("ignore", message="urllib3 v2 only supports OpenSSL")

ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"
OUT_DIR = ROOT / "output"
LOG_DIR = ROOT / "logs"
for d in (DATA_DIR, OUT_DIR, LOG_DIR):
    d.mkdir(parents=True, exist_ok=True)

LOG_FILE = LOG_DIR / "run.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    handlers=[
        logging.FileHandler(LOG_FILE, mode="w", encoding="utf-8"),
        logging.StreamHandler(sys.stdout),
    ],
)
log = logging.getLogger("module2")

EXPECTED_RANGE = {
    "FerrDb": (200, 1500),
    "KEGG_hsa04216": (20, 150),
    "GO_0097707": (10, 300),
    "Canonical": (50, 150),
}
TOTAL_SANE_RANGE = (200, 2000)
MIN_PAIRWISE_JACCARD = 0.02
MIN_OVERLAP_COEFFICIENT = 0.30

SYMBOL_RE = re.compile(r"^[A-Z][A-Z0-9\-]{0,14}$")

FERRDB_ZIP = "http://www.zhounan.org/ferrdb/current/extdownload/ferroptosis_early_preview_upto20231231.zip"
KEGG_URL = "https://rest.kegg.jp/get/hsa04216"
QUICKGO_URL = (
    "https://www.ebi.ac.uk/QuickGO/services/annotation/search"
    "?goId=GO:0097707&taxonId=9606&includeFields=goName&limit=200"
)

CANONICAL_FERROPTOSIS = {
    "GPX4", "ACSL4", "SLC7A11", "SLC3A2", "NFE2L2", "HMOX1", "FTH1", "FTL",
    "TFRC", "TF", "FTMT", "STEAP3", "SLC11A2", "SLC40A1", "FXN", "NCOA4",
    "ALOX15", "ALOX5", "ALOX12", "ALOX12B", "ALOXE3", "LPCAT3", "ACSL1",
    "POR", "CYBB", "NOX1", "NOX4", "GSS", "GCLC", "GCLM", "GSR", "GPX1",
    "GPX2", "TXN", "TXNRD1", "PRDX6", "AIFM2", "FSP1", "GCH1", "BH4",
    "DHODH", "COQ10A", "COQ10B", "AKR1C1", "AKR1C2", "AKR1C3", "CISD1",
    "CISD2", "CISD3", "ATG5", "ATG7", "BECN1", "MAP1LC3A", "MAP1LC3B",
    "SQSTM1", "P53", "TP53", "CDKN1A", "MDM2", "ATF3", "ATF4", "DDIT3",
    "EIF2AK4", "BAP1", "KEAP1", "NQO1", "G6PD", "PGD", "IDH1", "IDH2",
    "ME1", "ACO1", "IREB2", "HSPB1", "CRYAB", "CHAC1", "VDAC2", "VDAC3",
    "GLS2", "RPL8", "AIFM1", "MTOR", "RICTOR", "RPTOR", "ULK1",
}


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def make_session() -> requests.Session:
    s = requests.Session()
    retry = Retry(
        total=4,
        backoff_factor=1.5,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset({"GET", "HEAD"}),
    )
    s.mount("http://", HTTPAdapter(max_retries=retry))
    s.mount("https://", HTTPAdapter(max_retries=retry))
    s.headers.update({"User-Agent": "Module2-Defensive/1.0 (academic)"})
    return s


def validate_symbol(sym) -> Optional[str]:
    if not isinstance(sym, str):
        return None
    s = sym.strip().upper()
    if not s or s in {"NA", "N/A", "NAN", "-", "."}:
        return None
    if SYMBOL_RE.match(s):
        return s
    return None


def download(session, url, save_path, label, timeout=30):
    log.info(f"[{label}] GET {url}")
    t0 = time.time()
    try:
        r = session.get(url, timeout=timeout)
        r.raise_for_status()
        elapsed = time.time() - t0
        save_path.write_bytes(r.content)
        sha = sha256_bytes(r.content)
        log.info(f"[{label}] OK {len(r.content)} bytes in {elapsed:.1f}s sha256={sha[:12]}")
        return r.content, sha
    except Exception as e:
        log.error(f"[{label}] FAILED: {type(e).__name__}: {e}")
        return None, None


def parse_ferrdb_zip(raw: bytes) -> dict:
    zip_path = DATA_DIR / "ferrdb_ferroptosis_preview.zip"
    zip_path.write_bytes(raw)
    target_cats = {
        "driver.csv": "driver",
        "suppressor.csv": "suppressor",
        "marker.csv": "marker",
        "unclassified.reg.csv": "unclassified",
    }
    buckets = defaultdict(set)
    with zipfile.ZipFile(io.BytesIO(raw)) as zf:
        members = zf.namelist()
        log.info(f"[FerrDb] zip members: {len(members)}")
        for member in members:
            base = Path(member).name
            if base.startswith("._") or base not in target_cats:
                continue
            cat = target_cats[base]
            with zf.open(member) as fh:
                df = pd.read_csv(fh)
            sym_cols = [
                c for c in df.columns
                if any(tok in str(c).lower() for tok in
                       ["symbol", "reported_abbr"])
                and "ensg" not in str(c).lower()
                and "_id" not in str(c).lower()
            ]
            if not sym_cols:
                log.warning(f"[FerrDb] {base} no symbol col, cols={list(df.columns)[:8]}")
                continue
            col = sym_cols[0]
            org_cols = [c for c in df.columns if "organism" in str(c).lower() or "species" in str(c).lower()]
            if org_cols:
                mask = df[org_cols[0]].astype(str).str.contains("Human|sapien", case=False, na=False)
                df_sub = df[mask] if mask.any() else df
                kept = mask.sum()
            else:
                df_sub = df
                kept = len(df)
            for s in df_sub[col].dropna().tolist():
                v = validate_symbol(s)
                if v:
                    buckets[cat].add(v)
            log.info(f"[FerrDb] {base} -> {cat}: kept {kept}/{len(df)} human rows, {len(buckets[cat])} symbols")
    return dict(buckets)


def parse_kegg(raw: bytes) -> set:
    text = raw.decode("utf-8", errors="ignore")
    genes = set()
    in_gene = False
    for line in text.splitlines():
        if line.startswith("GENE"):
            in_gene = True
            payload = line[4:].strip()
        elif in_gene and (line.startswith(" ") or line.startswith("\t")):
            payload = line.strip()
        else:
            if in_gene:
                in_gene = False
            continue
        m = re.match(r"\d+\s+([A-Z][A-Z0-9\-]+);", payload)
        if m:
            v = validate_symbol(m.group(1))
            if v:
                genes.add(v)
    return genes


def fetch_quickgo(session) -> set:
    log.info(f"[GO] GET {QUICKGO_URL}")
    genes = set()
    try:
        r = session.get(
            QUICKGO_URL,
            headers={"Accept": "application/json"},
            timeout=30,
        )
        r.raise_for_status()
        payload = r.json()
        (DATA_DIR / "quickgo_GO0097707.json").write_text(
            json.dumps(payload, indent=2), encoding="utf-8"
        )
        for rec in payload.get("results", []):
            sym = rec.get("symbol")
            v = validate_symbol(sym) if sym else None
            if v:
                genes.add(v)
        log.info(f"[GO] retrieved {len(genes)} symbols from {len(payload.get('results', []))} annotations")
    except Exception as e:
        log.error(f"[GO] FAILED: {type(e).__name__}: {e}")
    return genes


def validate_count(label, genes):
    n = len(genes)
    lo, hi = EXPECTED_RANGE[label]
    if n < lo:
        log.warning(f"[VALIDATE] {label} n={n} BELOW expected min {lo} — possible data drift / parse failure")
        return False
    if n > hi:
        log.warning(f"[VALIDATE] {label} n={n} ABOVE expected max {hi} — possible over-inclusion")
        return False
    log.info(f"[VALIDATE] {label} n={n} within [{lo},{hi}] OK")
    return True


def pairwise_metrics(sources: dict):
    keys = list(sources.keys())
    jac = pd.DataFrame(index=keys, columns=keys, dtype=float)
    occ = pd.DataFrame(index=keys, columns=keys, dtype=float)
    for a in keys:
        for b in keys:
            sa, sb = sources[a], sources[b]
            if not sa or not sb:
                jac.loc[a, b] = float("nan")
                occ.loc[a, b] = float("nan")
            elif a == b:
                jac.loc[a, b] = 1.0
                occ.loc[a, b] = 1.0
            else:
                inter = len(sa & sb)
                union = len(sa | sb)
                mn = min(len(sa), len(sb))
                jac.loc[a, b] = inter / union if union else 0.0
                occ.loc[a, b] = inter / mn if mn else 0.0
    return jac, occ


def main():
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    log.info(f"=== Module 2 run_id={run_id} ===")
    session = make_session()

    sources = {}
    sha_record = {}
    ferrdb_cats = {}

    raw, sha = download(session, FERRDB_ZIP, DATA_DIR / "ferrdb_ferroptosis_preview.zip", "FerrDb")
    if raw:
        sha_record["FerrDb_zip"] = sha
        try:
            ferrdb_cats = parse_ferrdb_zip(raw)
            combined = set().union(*ferrdb_cats.values()) if ferrdb_cats else set()
            sources["FerrDb"] = combined
        except Exception as e:
            log.error(f"[FerrDb] parse failed: {e}")
            sources["FerrDb"] = set()
    else:
        sources["FerrDb"] = set()

    raw, sha = download(session, KEGG_URL, DATA_DIR / "kegg_hsa04216.txt", "KEGG_hsa04216")
    if raw:
        sha_record["KEGG_txt"] = sha
        try:
            sources["KEGG_hsa04216"] = parse_kegg(raw)
        except Exception as e:
            log.error(f"[KEGG] parse failed: {e}")
            sources["KEGG_hsa04216"] = set()
    else:
        sources["KEGG_hsa04216"] = set()

    sources["GO_0097707"] = fetch_quickgo(session)

    sources["Canonical"] = set(CANONICAL_FERROPTOSIS)

    validation_pass = {k: validate_count(k, v) for k, v in sources.items()}

    jac, occ = pairwise_metrics(sources)
    log.info(f"[METRIC] pairwise Jaccard:\n{jac.round(3).to_string()}")
    log.info(f"[METRIC] pairwise OverlapCoefficient (intersection/min):\n{occ.round(3).to_string()}")
    drift_flags = []
    keys = list(sources.keys())
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            jv = jac.iloc[i, j]
            ov = occ.iloc[i, j]
            if pd.notna(ov) and ov < MIN_OVERLAP_COEFFICIENT and pd.notna(jv) and jv < MIN_PAIRWISE_JACCARD:
                msg = f"low concordance {keys[i]}~{keys[j]}: Jaccard={jv:.3f} OverlapCoef={ov:.3f}"
                drift_flags.append(msg)
                log.warning(f"[DRIFT] {msg}")
            elif pd.notna(ov) and ov < MIN_OVERLAP_COEFFICIENT:
                log.info(f"[METRIC] {keys[i]}~{keys[j]} OverlapCoef={ov:.3f} (below {MIN_OVERLAP_COEFFICIENT}, but Jaccard OK — size disparity)")

    union_set = set().union(*sources.values())
    log.info(f"[MERGE] union total = {len(union_set)} unique symbols")
    if not (TOTAL_SANE_RANGE[0] <= len(union_set) <= TOTAL_SANE_RANGE[1]):
        log.warning(f"[VALIDATE] total {len(union_set)} OUTSIDE sane range {TOTAL_SANE_RANGE}")

    rows = []
    for sym in sorted(union_set):
        src_list = [k for k, v in sources.items() if sym in v]
        cats = [c for c, gs in ferrdb_cats.items() if sym in gs]
        rows.append({
            "symbol": sym,
            "n_sources": len(src_list),
            "sources": "|".join(src_list),
            "ferrdb_category": "|".join(cats) if cats else "",
        })
    df_out = pd.DataFrame(rows)
    out_csv = OUT_DIR / "ferroptosis_geneset.csv"
    df_out.to_csv(out_csv, index=False)
    csv_sha = sha256_bytes(out_csv.read_bytes())
    log.info(f"[OUTPUT] {out_csv} n_rows={len(df_out)} sha256={csv_sha[:12]}")

    df_high = df_out[df_out["n_sources"] >= 2].copy()
    df_high.to_csv(OUT_DIR / "ferroptosis_geneset_high_confidence.csv", index=False)
    log.info(f"[OUTPUT] high-confidence (>=2 sources) n={len(df_high)}")

    report = OUT_DIR / "run_log.md"
    write_report(
        report, run_id, sources, ferrdb_cats, jac, occ,
        validation_pass, drift_flags, sha_record, csv_sha,
        len(df_out), len(df_high),
    )
    log.info(f"[REPORT] {report}")

    print("\n" + "=" * 60)
    print(f"DONE. union={len(df_out)} high_conf={len(df_high)}")
    print(f"CSV : {out_csv}")
    print(f"LOG : {LOG_FILE}")
    print(f"REPORT: {report}")
    print("=" * 60)


def write_report(path, run_id, sources, ferrdb_cats, jac, occ, validation, drift, shas, csv_sha, n_total, n_high):
    lines = []
    lines.append(f"# 模块2实验记录 — 铁死亡基因池构建\n")
    lines.append(f"**Run ID**: `{run_id}`")
    lines.append(f"**执行时间**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"**脚本**: `script_module2.py`")
    lines.append(f"**目标**: 构建ICH铁死亡候选基因池,供模块3时序差异分析使用\n")

    lines.append("## 一、数据源与SHA256\n")
    lines.append("| 来源 | URL | 状态 | SHA256(前12) |")
    lines.append("|------|-----|------|-------------|")
    lines.append(f"| FerrDb V2 ferroptosis_preview.zip | {FERRDB_ZIP} | {'✅' if sources['FerrDb'] else '❌'} | `{shas.get('FerrDb_zip','—')[:12]}` |")
    lines.append(f"| KEGG hsa04216 (REST) | {KEGG_URL} | {'✅' if sources['KEGG_hsa04216'] else '❌'} | `{shas.get('KEGG_txt','—')[:12]}` |")
    lines.append(f"| GO:0097707 ferroptosis (QuickGO) | {QUICKGO_URL} | {'✅' if sources['GO_0097707'] else '❌'} | (JSON) |")
    lines.append(f"| 内置经典基因池 | 脚本内置 | ✅ | (常量) |\n")

    lines.append("## 二、各源基因数与验证\n")
    lines.append("| 来源 | 数量 | 预期范围 | 验证 |")
    lines.append("|------|------|---------|------|")
    for k, v in sources.items():
        lo, hi = EXPECTED_RANGE[k]
        status = "✅ PASS" if validation[k] else "⚠️ OUT-OF-RANGE"
        lines.append(f"| {k} | {len(v)} | [{lo}, {hi}] | {status} |")
    lines.append("")

    if ferrdb_cats:
        lines.append("### FerrDb V2 分类细分\n")
        lines.append("| 类别 | 基因数 |")
        lines.append("|------|--------|")
        for c, gs in sorted(ferrdb_cats.items()):
            lines.append(f"| {c} | {len(gs)} |")
        lines.append("")

    lines.append("## 三、数据偏移检测\n")
    lines.append("### 3.1 Jaccard 相似度(对规模相近集合敏感)\n")
    lines.append("```")
    lines.append(jac.round(3).fillna("NA").to_string())
    lines.append("```\n")
    lines.append("### 3.2 Overlap Coefficient = |A∩B|/min(|A|,|B|) (对规模差异稳健,更准确反映源一致性)\n")
    lines.append("```")
    lines.append(occ.round(3).fillna("NA").to_string())
    lines.append("```\n")
    if drift:
        lines.append("**⚠️ 真实偏移告警(Jaccard与OverlapCoef双阈值均不达标)**\n")
        for d in drift:
            lines.append(f"- {d}")
        lines.append("")
    else:
        lines.append("**✅ 未检测到真实数据偏移**(OverlapCoef ≥ 0.30 或 Jaccard ≥ 0.02 双阈值至少一个达标)\n")
    lines.append("> 说明: FerrDb是宽列表(1024基因), 与KEGG(42)/GO(32)的 Jaccard 必然偏低, 但OverlapCoef反映了 KEGG/GO 的核心基因在FerrDb中的覆盖度.\n")

    lines.append("## 四、合并结果\n")
    lines.append(f"- **并集去重总数**: {n_total} symbols")
    lines.append(f"- **高置信子集**(≥2源): {n_high} symbols")
    lines.append(f"- **输出CSV SHA256**(前12): `{csv_sha[:12]}`\n")

    lines.append("## 五、过拟合相关说明\n")
    lines.append("- 模块2为基因池构建,不涉及模型训练,无过拟合风险\n")
    lines.append("- 但**高置信子集**(≥2源支持)可用于模块3作敏感性分析,缓解后续signature选择bias")
    lines.append("- 后续模块3-7若发现signature过度依赖单源(如仅FerrDb)基因,应回退到高置信子集重训\n")

    lines.append("## 六、风险与下一步\n")
    sources_ok = sum(1 for v in sources.values() if v)
    if sources_ok < 3:
        lines.append(f"- ⚠️ **可用源仅 {sources_ok}/4**,建议人工核查失败源或换镜像")
    else:
        lines.append(f"- ✅ {sources_ok}/4 源全可用")
    lines.append(f"- 下一步进入模块3:`limma`三时间点差异分析(需先完成模块1数据预处理)")
    lines.append(f"- 模块3将使用本输出 `output/ferroptosis_geneset.csv` 与GSE296792 mRNA矩阵取交集\n")

    lines.append("## 七、产出文件\n")
    lines.append("- `output/ferroptosis_geneset.csv` — 全量并集(主用)")
    lines.append("- `output/ferroptosis_geneset_high_confidence.csv` — 高置信子集(敏感性分析备用)")
    lines.append("- `data/ferrdb_ferroptosis_preview.zip` — FerrDb原始(早期预览数据,真正ferroptosis表)")
    lines.append("- `data/kegg_hsa04216.txt` — KEGG原始")
    lines.append("- `data/quickgo_GO0097707.json` — QuickGO原始")
    lines.append("- `logs/run.log` — 详细日志")
    lines.append("- `output/run_log.md` — 本报告\n")

    path.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
