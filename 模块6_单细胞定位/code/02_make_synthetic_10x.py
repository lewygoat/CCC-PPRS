"""Synthetic 10X generator mirroring GSE266873 structure.

ONLY used as a DRY-RUN stand-in when the real GSE266873 supplementary files
cannot be downloaded (this host's DNS resolves *.ncbi.nlm.nih.gov to a proxy
that fails TLS). The QC pipeline (03_qc_pipeline.py) is cache-aware and will
automatically prefer real *.mtx.gz files placed under data/ over this synthetic
stand-in. Outputs generated from synthetic data are clearly tagged DRY-RUN.
"""
import gzip
import os
import sys
from pathlib import Path

import numpy as np
from scipy import sparse, io as sio

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.cwd()
SYN_DIR = ROOT / "data" / "synthetic"
SYN_DIR.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(20260530)

SAMPLES = {
    "GSM8255340": ("sample1_1_ICH_2h", "G1"),
    "GSM8255341": ("sample1_2_ICH_3h", "G1"),
    "GSM8255342": ("sample1_3_ICH_2h", "G1"),
    "GSM8255343": ("sample2_1_ICH_12h", "G2"),
    "GSM8255344": ("sample2_2_ICH_12h", "G2"),
    "GSM8255345": ("sample2_3_ICH_12h", "G2"),
    "GSM8255346": ("sample3_1_ICH_26h", "G3"),
    "GSM8255347": ("sample3_2_ICH_28h", "G3"),
    "GSM8255348": ("sample3_3_ICH_32h", "G3"),
}

MT_GENES = [
    "MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2", "MT-ATP8", "MT-ATP6",
    "MT-CO3", "MT-ND3", "MT-ND4L", "MT-ND4", "MT-ND5", "MT-ND6", "MT-CYB",
]
MICROGLIA = ["P2RY12", "TMEM119", "CX3CR1", "AIF1", "CSF1R", "C1QA", "C1QB"]
NEUTROPHIL = ["FCGR3B", "CSF3R", "S100A8", "S100A9", "FUT4", "CXCR2", "MPO"]
HUB = ["GPX4", "ACSL4", "SLC7A11", "HMOX1", "FTH1"]

N_FILLER = 1968
FILLER = [f"GENE{n:05d}" for n in range(N_FILLER)]
GENES = MT_GENES + MICROGLIA + NEUTROPHIL + HUB + FILLER
N_GENES = len(GENES)
mt_idx = np.array([i for i, g in enumerate(GENES) if g.startswith("MT-")])

# per-sample cell counts (realistic for sorted immune cells); one drifted sample
BASE_CELLS = {
    "GSM8255340": 3200, "GSM8255341": 2800, "GSM8255342": 3500,
    "GSM8255343": 4100, "GSM8255344": 2600, "GSM8255345": 3900,
    "GSM8255346": 4800, "GSM8255347": 5200, "GSM8255348": 4400,
}
# GSM8255344 is the deliberately DRIFTED sample: elevated MT% (degraded library)
DRIFT_SAMPLE = "GSM8255344"


def simulate_sample(gsm):
    n_cells = BASE_CELLS[gsm]
    # mixture: 70% microglia-like, 25% neutrophil-like, 5% low-quality debris
    n_lowq = int(n_cells * 0.06)
    n_good = n_cells - n_lowq

    # baseline expression per gene (lognormal mean depth)
    base_mu = RNG.lognormal(mean=0.2, sigma=1.1, size=N_GENES)
    # cell-level depth factor
    depth = RNG.lognormal(mean=7.6, sigma=0.45, size=n_good)  # ~2000 UMIs median

    cols, rows, vals = [], [], []
    mt_boost = 4.0 if gsm == DRIFT_SAMPLE else 1.0

    for c in range(n_good):
        lam = base_mu * (depth[c] / base_mu.sum())
        lam = lam.copy()
        lam[mt_idx] *= mt_boost  # drifted sample over-expresses MT
        counts = RNG.poisson(lam)
        nz = np.nonzero(counts)[0]
        rows.extend(nz.tolist())
        cols.extend([c] * len(nz))
        vals.extend(counts[nz].tolist())

    # low-quality cells: shallow depth, very high MT fraction
    for j in range(n_lowq):
        c = n_good + j
        lam = base_mu * (RNG.lognormal(mean=5.0, sigma=0.4) / base_mu.sum())
        lam[mt_idx] *= 12.0
        counts = RNG.poisson(lam)
        nz = np.nonzero(counts)[0]
        rows.extend(nz.tolist())
        cols.extend([c] * len(nz))
        vals.extend(counts[nz].tolist())

    M = sparse.csc_matrix(
        (np.array(vals, dtype=np.int32),
         (np.array(rows), np.array(cols))),
        shape=(N_GENES, n_cells),
    )
    return M


def write_10x(gsm, tag, M):
    prefix = f"{gsm}_{tag}_"
    # matrix.mtx.gz
    mtx_path = SYN_DIR / f"{prefix}matrix.mtx.gz"
    with gzip.open(mtx_path, "wb") as fh:
        sio.mmwrite(fh, M.tocoo(), field="integer")
    # barcodes.tsv.gz
    bc_path = SYN_DIR / f"{prefix}barcodes.tsv.gz"
    with gzip.open(bc_path, "wt") as fh:
        for i in range(M.shape[1]):
            fh.write(f"{gsm}_CELL{i:06d}-1\n")
    # features.tsv.gz  (10X v3: gene_id, symbol, type)
    ft_path = SYN_DIR / f"{prefix}features.tsv.gz"
    with gzip.open(ft_path, "wt") as fh:
        for k, g in enumerate(GENES):
            fh.write(f"ENSGSYN{k:08d}\t{g}\tGene Expression\n")
    return mtx_path.stat().st_size


def main():
    print(f"[SYNTH] writing synthetic GSE266873-like 10X into {SYN_DIR}")
    print(f"[SYNTH] {N_GENES} genes ({len(mt_idx)} MT-), {len(SAMPLES)} samples, "
          f"drift sample = {DRIFT_SAMPLE} (MT-elevated)")
    for gsm, (tag, grp) in SAMPLES.items():
        M = simulate_sample(gsm)
        sz = write_10x(gsm, tag, M)
        print(f"[SYNTH] {gsm} ({tag}, {grp}): {M.shape[0]}x{M.shape[1]} "
              f"nnz={M.nnz} mtx={sz}B")
    print("[SYNTH] done")


if __name__ == "__main__":
    main()
