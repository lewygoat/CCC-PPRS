source(file.path(Sys.getenv("PROJ_ROOT", unset = getwd()), "code", "00_utils.R"))

ensure_packages(cran = c("matrixStats"),
                bioc = c("Biobase", "sva", "edgeR", "limma"))

suppressPackageStartupMessages({
  library(Biobase)
  library(sva)
  library(edgeR)
  library(limma)
  library(matrixStats)
})

cache_path <- file.path(DATA_DIR, "_eset_list.rds")
meta_cache <- file.path(DATA_DIR, "_meta_all.rds")
if (!file.exists(cache_path)) stop("Run 01_preflight_download.R first.", call. = FALSE)
if (!file.exists(meta_cache)) stop("Run 04_metadata_align.R first.", call. = FALSE)
results  <- readRDS(cache_path)
all_meta <- readRDS(meta_cache)

OUT_DIR <- file.path(PROJ_ROOT, "output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

CPM_THRESHOLD     <- 1
SAMPLE_FRAC       <- 0.30
COMBAT_BATCH_MIN  <- 2

detect_type <- function(mat) {
  mx  <- suppressWarnings(max(mat, na.rm = TRUE))
  mn  <- suppressWarnings(min(mat, na.rm = TRUE))
  has_int <- all(mat == round(mat), na.rm = TRUE) && mn >= 0
  if (mx < 30 && mn > -5) return("log_intensity")
  if (has_int && mx > 1000) return("counts")
  if (mx > 100 && mn >= 0 && !has_int) return("intensity")
  "unknown"
}

guess_modality <- function(eset, gse_id) {
  ann  <- tolower(annotation(eset))
  feat <- rownames(exprs(eset))
  if (grepl("mirna|microrna", ann)) return("mirna")
  mirna_like <- grepl("^(hsa-)?(let-|mir-)?\\d|^let-|^hsa-mir", feat, ignore.case = TRUE) |
                grepl("^hsa[._-]", feat, ignore.case = TRUE)
  if (mean(mirna_like) > 0.5) return("mirna")
  if (length(feat) < 3000 && mean(mirna_like) > 0.2) return("mirna")
  "mrna"
}

normalize_one <- function(mat, modality, gse_id, slot_idx) {
  tp <- detect_type(mat)
  log_msg("INFO", sprintf("[%s slot %d] modality=%s, dtype=%s, dim=%d×%d, range=[%.2f, %.2f]",
                          gse_id, slot_idx, modality, tp,
                          nrow(mat), ncol(mat),
                          min(mat, na.rm = TRUE), max(mat, na.rm = TRUE)))

  if (tp == "log_intensity") {
    out <- mat
    log_msg("INFO", sprintf("[%s slot %d] already log-scale, skip log", gse_id, slot_idx))
  } else if (modality == "mirna" && tp == "counts") {
    dge <- DGEList(counts = mat)
    keep <- rowSums(cpm(dge) > CPM_THRESHOLD) >= ceiling(SAMPLE_FRAC * ncol(mat))
    log_msg("INFO", sprintf("[%s slot %d] miRNA low-expr filter: kept %d / %d features (CPM>%g in >=%.0f%% samples)",
                            gse_id, slot_idx, sum(keep), nrow(mat),
                            CPM_THRESHOLD, 100 * SAMPLE_FRAC))
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- calcNormFactors(dge, method = "TMM")
    out <- cpm(dge, log = TRUE, prior.count = 1)
  } else if (modality == "mrna" && tp == "counts") {
    dge <- DGEList(counts = mat)
    dge <- calcNormFactors(dge, method = "TMM")
    out <- cpm(dge, log = TRUE, prior.count = 1)
    log_msg("INFO", sprintf("[%s slot %d] mRNA: TMM + log2(CPM+1)", gse_id, slot_idx))
  } else if (tp == "intensity") {
    out <- log2(mat + 1)
    out <- normalizeBetweenArrays(out, method = "quantile")
    log_msg("INFO", sprintf("[%s slot %d] intensity: log2 + quantile normalize", gse_id, slot_idx))
  } else {
    out <- log2(pmax(mat, 0) + 1)
    log_msg("WARN", sprintf("[%s slot %d] dtype=unknown, fallback to log2(x+1)", gse_id, slot_idx))
  }
  out
}

run_combat <- function(mat, meta, gse_id, slot_idx) {
  m <- meta[match(colnames(mat), meta$sample_id), , drop = FALSE]
  ok <- !is.na(m$stage) & !is.na(m$batch)
  if (sum(ok) < ncol(mat))
    log_msg("WARN", sprintf("[%s slot %d] %d samples lack stage/batch, dropped from ComBat",
                             gse_id, slot_idx, ncol(mat) - sum(ok)))
  mat <- mat[, ok, drop = FALSE]
  m   <- m[ok, , drop = FALSE]
  if (length(unique(m$batch)) < COMBAT_BATCH_MIN) {
    log_msg("INFO", sprintf("[%s slot %d] only %d batch level, skip ComBat",
                             gse_id, slot_idx, length(unique(m$batch))))
    return(mat)
  }
  mod <- model.matrix(~ stage, data = m)
  out <- tryCatch(
    ComBat(dat = mat, batch = m$batch, mod = mod, par.prior = TRUE, prior.plots = FALSE),
    error = function(e) {
      log_msg("ERROR", sprintf("[%s slot %d] ComBat failed: %s", gse_id, slot_idx, conditionMessage(e)))
      mat
    })
  log_msg("INFO", sprintf("[%s slot %d] ComBat applied with %d batch levels, stage as covariate",
                           gse_id, slot_idx, length(unique(m$batch))))
  out
}

process_dataset <- function(gse_id, role) {
  obj <- results[[gse_id]]
  if (is.null(obj)) {
    log_msg("ERROR", sprintf("[%s] no data, skip", gse_id))
    return(invisible(NULL))
  }
  meta_keys <- grep(paste0("^", gse_id), names(all_meta), value = TRUE)
  meta_combined <- do.call(rbind, all_meta[meta_keys])
  if ("inferred_stage" %in% colnames(meta_combined) && !"stage" %in% colnames(meta_combined))
    colnames(meta_combined)[colnames(meta_combined) == "inferred_stage"] <- "stage"

  for (i in seq_along(obj)) {
    e <- obj[[i]]
    if (!inherits(e, "ExpressionSet")) next
    modality <- guess_modality(e, gse_id)
    mat <- exprs(e)
    mat <- mat[complete.cases(mat), , drop = FALSE]
    norm_mat <- normalize_one(mat, modality, gse_id, i)
    norm_mat <- run_combat(norm_mat, meta_combined, gse_id, i)

    if (role == "train" && modality == "mrna") {
      fn <- "expr_mrna.csv"
    } else if (role == "train" && modality == "mirna") {
      fn <- "expr_mirna.csv"
    } else if (role == "external" && modality == "mrna") {
      fn <- "expr_mrna_external.csv"
    } else if (role == "external" && modality == "mirna") {
      fn <- "expr_mirna_external.csv"
    } else {
      fn <- sprintf("expr_%s_slot%d.csv", gse_id, i)
    }
    out_path <- file.path(OUT_DIR, fn)
    write.csv(norm_mat, out_path)
    log_msg("INFO", sprintf("[%s slot %d] saved -> %s (%d × %d)",
                             gse_id, i, out_path, nrow(norm_mat), ncol(norm_mat)))

    after_var <- rowVars(norm_mat)
    nz <- after_var > 0 & is.finite(after_var)
    if (sum(nz) >= 50) {
      mat_pca <- norm_mat[nz, , drop = FALSE]
      vv <- after_var[nz]
      top <- order(vv, decreasing = TRUE)[seq_len(min(2000, length(vv)))]
      pc  <- prcomp(t(mat_pca[top, ]), scale. = TRUE)
      pct <- round(100 * (pc$sdev^2) / sum(pc$sdev^2), 1)
      log_msg("INFO", sprintf("[%s slot %d] post-norm PCA PC1=%.1f%%, PC2=%.1f%% (dropped %d zero-var)",
                                gse_id, i, pct[1], pct[2], length(after_var) - sum(nz)))
    } else {
      pct <- c(NA, NA)
      log_msg("WARN", sprintf("[%s slot %d] post-norm PCA skipped (insufficient non-zero-var rows)",
                                gse_id, i))
    }
    if (pct[1] > 60)
      log_msg("WARN", sprintf("[%s slot %d] PC1>60%% AFTER ComBat — strong residual batch suspected; check covariates",
                               gse_id, i))
  }
}

log_msg("INFO", "==== Stage 1.2-1.3+1.5 : normalize + ComBat ====")
process_dataset("GSE296792", "train")
process_dataset("GSE125512", "external")
log_msg("INFO", "Stage 1.2-1.5 complete.")
