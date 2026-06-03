# Module 6 — Single-cell localisation (GSE266873)

## Status: NOT part of the published results

The single-cell localisation layer is **deferred to future work** in the
manuscript. None of the contents of this module contribute to any figure,
table, or numerical result reported in the paper. The manuscript explicitly
states that GSE266873 "is referenced as a reserved single-cell resource for
planned future analysis and was not analysed in the present study."

## Why a synthetic-data generator is present

`code/02_make_synthetic_10x.py` produces **synthetic** 10x matrices that only
mirror the structure (sample/barcode/feature layout) of GSE266873. It exists
solely as a **dry-run stand-in** for pipeline development, used when the real
GSE266873 supplementary files could not be downloaded in the analysis
environment (the host's DNS resolved `*.ncbi.nlm.nih.gov` to a proxy that
failed TLS).

- The synthetic matrices contain **no real biological signal**.
- Every artefact derived from them is tagged **DRYRUN** (see
  `output/*_DRYRUN.*`) and must not be interpreted as a result.
- `code/03_qc_pipeline.py` is cache-aware: if real `*.mtx.gz` files are placed
  under `data/`, it uses those in preference to the synthetic stand-in.

## Reproducing with real data

1. Download the real GSE266873 supplementary files (see `code/01_download.R`).
2. Place the `*.mtx.gz` / barcodes / features under `data/`.
3. Re-run `code/03_qc_pipeline.py`; it will ignore the synthetic stand-in.

The synthetic generator can be deleted entirely without affecting any
published result.
