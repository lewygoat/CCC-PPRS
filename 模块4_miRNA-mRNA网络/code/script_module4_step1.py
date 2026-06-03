import gzip
import hashlib
import io
import logging
import re
import sys
import time
import warnings
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Optional

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

warnings.filterwarnings("ignore", message="urllib3 v2 only supports OpenSSL")

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data"
OUT_DIR = ROOT / "output"
LOG_DIR = ROOT / "log"
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
log = logging.getLogger("module4_step1")

HUB_GENES = ["GPX4", "ACSL4", "SLC7A11", "HMOX1", "FTH1"]

EXPECTED_PER_GENE = {
    "ENCORI": (1, 600),
    "miRDB": (5, 1500),
}
MIRDB_SCORE_CUT = 80
ENCORI_MIN_CLIP = 1
ENCORI_MIN_PROGRAMS = 2

KNOWN_AXES = {
    "GPX4": ["hsa-miR-15a-5p", "hsa-miR-15b-5p", "hsa-miR-214-3p"],
    "SLC7A11": ["hsa-miR-27a-3p", "hsa-miR-26b-5p"],
    "HMOX1": ["hsa-miR-377-3p", "hsa-miR-24-3p"],
    "ACSL4": ["hsa-miR-424-5p"],
    "FTH1": ["hsa-miR-200b-3p"],
}

MIRNA_RE = re.compile(r"^hsa-(let|miR)-[0-9a-zA-Z\-]+$")

ENCORI_TMPL = (
    "https://rnasysu.com/encori/api/miRNATarget/?assembly=hg38&geneType=mRNA"
    "&miRNA=all&clipExpNum={clip}&degraExpNum=0&pancancerNum=0"
    "&programNum={prog}&program=None&target={gene}&cellType=all"
)
MIRDB_URLS = [
    "https://mirdb.org/download/miRDB_v6.0_prediction_result.txt.gz",
    "http://mirdb.org/download/miRDB_v6.0_prediction_result.txt.gz",
]


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
    s.headers.update({"User-Agent": "Module4-Defensive/1.0 (academic)"})
    return s


def download_with_fallback(session, urls, save_path: Path, label: str, timeout=120):
    last_err = None
    for url in urls:
        log.info(f"[{label}] GET {url}")
        t0 = time.time()
        try:
            r = session.get(url, timeout=timeout)
            r.raise_for_status()
            elapsed = time.time() - t0
            save_path.write_bytes(r.content)
            sha = sha256_bytes(r.content)
            log.info(
                f"[{label}] OK {len(r.content)} bytes in {elapsed:.1f}s sha256={sha[:12]}"
            )
            return r.content, sha, url
        except Exception as e:
            last_err = e
            log.warning(f"[{label}] try failed for {url}: {type(e).__name__}: {e}")
    log.error(f"[{label}] ALL fallbacks failed; last error: {last_err}")
    return None, None, None


def validate_mirna_id(name) -> Optional[str]:
    if not isinstance(name, str):
        return None
    s = name.strip()
    if not s:
        return None
    if MIRNA_RE.match(s):
        return s
    return None


def fetch_encori_for_gene(session, gene: str) -> pd.DataFrame:
    url = ENCORI_TMPL.format(
        gene=gene, clip=ENCORI_MIN_CLIP, prog=ENCORI_MIN_PROGRAMS
    )
    log.info(f"[ENCORI] GET {gene}  (clip>={ENCORI_MIN_CLIP}, prog>={ENCORI_MIN_PROGRAMS})")
    t0 = time.time()
    try:
        r = session.get(url, timeout=60)
        r.raise_for_status()
    except Exception as e:
        log.error(f"[ENCORI] {gene} FAILED: {type(e).__name__}: {e}")
        return pd.DataFrame()
    text = r.text
    cache = DATA_DIR / f"encori_{gene}.tsv"
    cache.write_text(text, encoding="utf-8")
    sha = sha256_bytes(text.encode("utf-8"))
    log.info(
        f"[ENCORI] {gene} OK {len(text)} bytes in {time.time()-t0:.1f}s sha256={sha[:12]}"
    )
    lines = [ln for ln in text.splitlines() if not ln.startswith("#")]
    if len(lines) < 2:
        log.warning(f"[ENCORI] {gene} returned empty table")
        return pd.DataFrame()
    df = pd.read_csv(io.StringIO("\n".join(lines)), sep="\t")
    if "miRNAname" not in df.columns or "geneName" not in df.columns:
        log.error(f"[ENCORI] {gene} unexpected columns: {list(df.columns)[:8]}")
        return pd.DataFrame()
    sub = df[["miRNAname", "geneName", "clipExpNum"]].copy()
    sub.columns = ["mirna", "target", "clip_n"]
    sub["mirna"] = sub["mirna"].astype(str).str.strip()
    sub["target"] = sub["target"].astype(str).str.strip().str.upper()
    sub = sub[sub["mirna"].map(validate_mirna_id).notna()]
    sub = sub.drop_duplicates(subset=["mirna", "target"])
    return sub


def collect_encori(session, hub_genes) -> tuple:
    frames = []
    shas = []
    for g in hub_genes:
        df = fetch_encori_for_gene(session, g)
        if not df.empty:
            frames.append(df)
        time.sleep(0.5)
    if not frames:
        return pd.DataFrame(), ""
    out = pd.concat(frames, ignore_index=True)
    composite = sha256_bytes(out.to_csv(index=False).encode("utf-8"))
    return out, composite


def parse_mirdb(raw: bytes) -> pd.DataFrame:
    bio = io.BytesIO(raw)
    with gzip.GzipFile(fileobj=bio) as gz:
        df = pd.read_csv(
            gz,
            sep="\t",
            header=None,
            names=["mirna", "refseq", "score"],
            dtype={"mirna": str, "refseq": str, "score": float},
        )
    log.info(f"[miRDB] raw rows={len(df)}")
    df = df[df["mirna"].str.startswith("hsa-", na=False)].copy()
    df = df[df["score"] >= MIRDB_SCORE_CUT]
    df["mirna"] = df["mirna"].str.strip()
    df = df[df["mirna"].map(validate_mirna_id).notna()]
    log.info(f"[miRDB] hsa rows score>={MIRDB_SCORE_CUT}: {len(df)} miRNAs={df['mirna'].nunique()} refseqs={df['refseq'].nunique()}")
    return df


def refseq_to_symbol(refseqs, session) -> dict:
    """NCBI E-utilities batched lookup of RefSeq -> gene symbol."""
    out = {}
    refseqs = [r for r in refseqs if isinstance(r, str) and r]
    if not refseqs:
        return out
    BATCH = 200
    for i in range(0, len(refseqs), BATCH):
        chunk = refseqs[i : i + BATCH]
        q = " OR ".join(f"{r}[Accession]" for r in chunk)
        try:
            r = session.get(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                params={"db": "nuccore", "term": q, "retmax": BATCH, "retmode": "json"},
                timeout=30,
            )
            r.raise_for_status()
            ids = r.json().get("esearchresult", {}).get("idlist", [])
            if not ids:
                continue
            r2 = session.get(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
                params={"db": "nuccore", "id": ",".join(ids), "retmode": "json"},
                timeout=30,
            )
            r2.raise_for_status()
            for k, v in r2.json().get("result", {}).items():
                if k == "uids":
                    continue
                acc = v.get("accessionversion", "").split(".")[0]
                title = v.get("title", "")
                m = re.search(r"\(([A-Z][A-Z0-9\-]{0,14})\)", title)
                if acc and m:
                    out[acc] = m.group(1)
            time.sleep(0.34)
        except Exception as e:
            log.warning(f"[NCBI] batch {i} failed: {e}")
    return out


def per_gene_query(df_pairs: pd.DataFrame, hub_genes, label: str) -> dict:
    res = {}
    for g in hub_genes:
        mirs = sorted(df_pairs.loc[df_pairs["target"] == g, "mirna"].unique().tolist())
        res[g] = mirs
        log.info(f"[{label}] {g}: {len(mirs)} miRNAs")
        lo, hi = EXPECTED_PER_GENE[label]
        if not (lo <= len(mirs) <= hi):
            log.warning(
                f"[{label}] {g} count {len(mirs)} OUTSIDE expected [{lo},{hi}] — possible drift"
            )
    return res


def jaccard(a, b):
    A, B = set(a), set(b)
    if not A and not B:
        return 1.0
    if not A or not B:
        return 0.0
    return len(A & B) / len(A | B)


def overlap_coef(a, b):
    A, B = set(a), set(b)
    if not A or not B:
        return 0.0
    return len(A & B) / min(len(A), len(B))


def biological_sanity(consensus_by_gene: dict) -> dict:
    report = {}
    for g, known_list in KNOWN_AXES.items():
        hits = [m for m in known_list if m in consensus_by_gene.get(g, [])]
        report[g] = {
            "known_axes_total": len(known_list),
            "recovered": len(hits),
            "recovered_list": hits,
        }
        if not hits:
            log.warning(f"[SANITY] {g}: NONE of {known_list} recovered — suspect")
        else:
            log.info(f"[SANITY] {g}: recovered {len(hits)}/{len(known_list)}: {hits}")
    return report


def main():
    log.info("=" * 70)
    log.info("Module 4 Step 1: multiMiR-equivalent query for hub genes")
    log.info(f"hub genes: {HUB_GENES}")

    session = make_session()
    sources = {}

    log.info("[NOTE] miRTarBase mirrors (cuhk.edu.cn / awi.cuhk.edu.cn) all timed out — "
             "switching to ENCORI CLIP-Seq experimental evidence as the experimental source.")
    df_encori, sha_encori = collect_encori(session, HUB_GENES)
    if df_encori.empty:
        log.error("[ENCORI] DEAD — cannot proceed without at least 2 sources")
        return 1
    sources["ENCORI"] = {
        "sha256": sha_encori,
        "url": "rnasysu.com/encori/api per-gene",
        "rows": len(df_encori),
        "per_gene": per_gene_query(df_encori, HUB_GENES, "ENCORI"),
    }

    mirdb_path = DATA_DIR / "miRDB_v6.0_prediction_result.txt.gz"
    if mirdb_path.exists():
        log.info(f"[miRDB] using cached {mirdb_path.name}")
        raw_mirdb = mirdb_path.read_bytes()
        sha_mirdb = sha256_bytes(raw_mirdb)
        used_url2 = "cache"
    else:
        raw_mirdb, sha_mirdb, used_url2 = download_with_fallback(
            session, MIRDB_URLS, mirdb_path, "miRDB"
        )
    if raw_mirdb is None:
        log.error("[miRDB] DEAD — cannot proceed without at least 2 sources")
        return 1
    df_mirdb_raw = parse_mirdb(raw_mirdb)

    hub_refseqs_needed = set()
    refseq_universe = df_mirdb_raw["refseq"].unique().tolist()
    log.info(f"[miRDB] mapping {len(refseq_universe)} RefSeqs to symbols via NCBI")
    rs_to_sym = refseq_to_symbol(refseq_universe, session)
    df_mirdb_raw["target"] = df_mirdb_raw["refseq"].map(rs_to_sym).fillna("")
    df_mirdb = df_mirdb_raw[df_mirdb_raw["target"].isin(HUB_GENES)].copy()
    log.info(
        f"[miRDB] after symbol mapping, rows for hub genes: {len(df_mirdb)}; "
        f"refseq->symbol resolved: {sum(1 for v in rs_to_sym.values() if v)}/{len(refseq_universe)}"
    )
    sources["miRDB"] = {
        "sha256": sha_mirdb,
        "url": used_url2,
        "rows": len(df_mirdb),
        "per_gene": per_gene_query(df_mirdb, HUB_GENES, "miRDB"),
    }

    log.info("=" * 70)
    log.info("Cross-source consensus (>=2 sources)")
    consensus = {}
    drift_report = {}
    for g in HUB_GENES:
        sets = {src: set(sources[src]["per_gene"][g]) for src in sources}
        names = list(sets.keys())
        union = set().union(*sets.values())
        cons = sorted([m for m in union if sum(m in sets[s] for s in names) >= 2])
        consensus[g] = cons
        j = jaccard(sets[names[0]], sets[names[1]])
        oc = overlap_coef(sets[names[0]], sets[names[1]])
        drift_report[g] = {"jaccard": j, "overlap_coef": oc, "consensus_n": len(cons)}
        log.info(
            f"[CONSENSUS] {g}: {names[0]}={len(sets[names[0]])} "
            f"{names[1]}={len(sets[names[1]])} J={j:.3f} OC={oc:.3f} ≥2src={len(cons)}"
        )

    sanity = biological_sanity(consensus)

    consensus_rows = []
    for g, mirs in consensus.items():
        for m in mirs:
            in_mtb = m in sources["miRTarBase"]["per_gene"][g]
            in_mirdb = m in sources["miRDB"]["per_gene"][g]
            consensus_rows.append(
                {
                    "target": g,
                    "mirna": m,
                    "in_miRTarBase": int(in_mtb),
                    "in_miRDB": int(in_mirdb),
                    "n_sources": int(in_mtb) + int(in_mirdb),
                }
            )
    df_cons = pd.DataFrame(consensus_rows)
    df_cons.to_csv(OUT_DIR / "hub_mirna_candidates_consensus.csv", index=False)
    log.info(f"[OUTPUT] consensus pairs: {len(df_cons)} -> hub_mirna_candidates_consensus.csv")

    union_rows = []
    for g in HUB_GENES:
        for src, info in sources.items():
            for m in info["per_gene"][g]:
                union_rows.append({"target": g, "mirna": m, "source": src})
    pd.DataFrame(union_rows).to_csv(OUT_DIR / "hub_mirna_candidates_union.csv", index=False)

    overfit_flag = False
    for g in HUB_GENES:
        n_cons = drift_report[g]["consensus_n"]
        n_encori = len(sources["ENCORI"]["per_gene"][g])
        n_mirdb = len(sources["miRDB"]["per_gene"][g])
        if n_cons > 0 and max(n_encori, n_mirdb) > 0:
            shrink = 1 - n_cons / max(n_encori, n_mirdb)
            if shrink < 0.30:
                log.warning(
                    f"[OVERFIT-CHECK] {g}: consensus={n_cons} vs larger source={max(n_encori, n_mirdb)} "
                    f"shrinkage={shrink:.2f} <0.30 — consensus too permissive"
                )
                overfit_flag = True
        if n_cons == 0:
            log.warning(f"[OVERFIT-CHECK] {g}: zero consensus — too strict / source disagreement")
    log.info(f"[OVERFIT-CHECK] overall flag={overfit_flag}")

    report_lines = [
        "# Module 4 Step 1 - run log",
        f"timestamp: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "## Sources",
    ]
    for src, info in sources.items():
        report_lines += [
            f"- **{src}**: url=`{info['url']}` sha256={info['sha256'][:16]} rows={info['rows']}",
        ]
    report_lines += ["", "## Per-gene candidate counts"]
    report_lines += ["| Gene | ENCORI | miRDB | Jaccard | OverlapCoef | Consensus>=2 |", "|---|---|---|---|---|---|"]
    for g in HUB_GENES:
        d = drift_report[g]
        report_lines.append(
            f"| {g} | {len(sources['ENCORI']['per_gene'][g])} | "
            f"{len(sources['miRDB']['per_gene'][g])} | {d['jaccard']:.3f} | "
            f"{d['overlap_coef']:.3f} | {d['consensus_n']} |"
        )
    report_lines += ["", "## Biological sanity (known axes recovery)"]
    report_lines += ["| Gene | known | recovered | which |", "|---|---|---|---|"]
    for g, s in sanity.items():
        report_lines.append(
            f"| {g} | {s['known_axes_total']} | {s['recovered']} | {', '.join(s['recovered_list']) or '—'} |"
        )
    (OUT_DIR / "run_log.md").write_text("\n".join(report_lines), encoding="utf-8")
    log.info(f"[OUTPUT] run_log.md written")
    log.info("DONE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
