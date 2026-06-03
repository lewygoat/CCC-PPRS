# CCC-PPRS

Analysis code for the manuscript *Temporal Reprogramming of the Ferroptosis
Axis in Peripheral Blood after Intracerebral Haemorrhage: A Cross-Cohort,
Cross-Platform Phase-Resolved Signature*.

The pipeline is organised as numbered modules. Each module keeps its scripts
under `code/` (module 2 keeps a single top-level script). Run modules in
numerical order; later modules read the `output/` of earlier ones as siblings.

## Module map

| Module | Folder | Paper step |
|--------|--------|-----------|
| M1 | `模块1_预检` | Preprocessing & QC; ComBat batch correction (38,349 mRNA / 912 miRNA / 77 samples) |
| M2 | `模块2_铁死亡基因集` | Four-source ferroptosis gene pool (1,059 union) |
| M3 | `模块3_时序差异_枢纽轴` | Temporal differential expression (limma); ferroptosis-restricted DE-FRGs |
| M4 | `模块4_miRNA-mRNA网络` | miRNA–mRNA query (multiMiR/ENCORI) and patient-level expression correlation |
| M5 | `模块5_bootstrap稳健性` | 333-resample bootstrap hub-stability analysis |
| M6 | `模块6_单细胞定位` | Single-cell localisation — DEFERRED to future work (see module README) |
| M7 | `模块7_分期signature` | Dual-track Elastic-Net staging signature; external + paired-temporal validation |

## Datasets (not included; download from GEO)

- GSE296792 — longitudinal training cohort (77 paired mRNA/miRNA)
- GSE125512 — cross-platform external validation (n = 22)
- GSE266873 — single-cell resource reserved for future work (not analysed)

Curated gene lists (FerrDb V2, KEGG hsa04216, GO:0097707) are downloaded or
cached by module 2.

## Dependencies

- R 4.3.x: limma 3.58, multiMiR 1.24, glmnet 4.1, pROC 1.18, rms 6.7, rmda 1.6,
  sva (ComBat), Seurat 5.0, Harmony 1.2
- Python 3: numpy, scipy, pandas, scikit-learn

## Notes on scope and reproducibility

- **Module 6 is a dry-run on synthetic data and contributes nothing to the
  published results.** The single-cell layer is deferred to future work; see
  `模块6_单细胞定位/README.md`.
- **Functional enrichment (paper Supplementary Fig. S4) has no script here.**
  The Enrichr/GSEA-KEGG results were produced with an external tool; only the
  output tables exist in the original working tree.
- **Module 4 scripts contain absolute machine-specific input paths.** They must
  be repointed to the local module output directories before they will run on
  another machine.
