source(file.path(Sys.getenv("PROJ_ROOT", unset = getwd()), "code", "00_utils.R"))

ensure_packages(cran = c("ggplot2", "matrixStats"))

suppressPackageStartupMessages({
  library(Biobase)
  library(ggplot2)
  library(matrixStats)
})

cache_path <- file.path(DATA_DIR, "_eset_list.rds")
if (!file.exists(cache_path)) stop("Run 01_preflight_download.R first.", call. = FALSE)
results <- readRDS(cache_path)

FIG_DIR <- file.path(PROJ_ROOT, "log", "preflight_figs")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

maybe_log2 <- function(mat) {
  q99 <- suppressWarnings(quantile(mat, 0.99, na.rm = TRUE))
  if (is.finite(q99) && q99 > 50) {
    log_msg("INFO", sprintf("Auto log2(x+1) applied (q99=%.1f)", q99))
    mat <- log2(mat + 1)
  }
  mat
}

drift_scan_one <- function(gse_id, slot_idx, eset) {
  mat <- exprs(eset)
  if (is.null(mat) || ncol(mat) < 4) {
    log_msg("WARN", sprintf("[%s][slot %d] too few samples for drift scan", gse_id, slot_idx))
    return(invisible(NULL))
  }
  mat <- mat[complete.cases(mat), , drop = FALSE]
  if (nrow(mat) < 100) {
    log_msg("WARN", sprintf("[%s][slot %d] only %d complete features, skip", gse_id, slot_idx, nrow(mat)))
    return(invisible(NULL))
  }
  mat <- maybe_log2(mat)

  vars <- rowVars(mat)
  nz   <- vars > 0 & is.finite(vars)
  if (sum(nz) < 50) {
    log_msg("WARN", sprintf("[%s][slot %d] only %d non-zero-var rows, skip PCA", gse_id, slot_idx, sum(nz)))
    return(invisible(NULL))
  }
  if (sum(nz) < length(vars))
    log_msg("INFO", sprintf("[%s][slot %d] dropped %d zero-var rows before PCA",
                            gse_id, slot_idx, length(vars) - sum(nz)))
  mat  <- mat[nz, , drop = FALSE]
  vars <- vars[nz]
  top  <- order(vars, decreasing = TRUE)[seq_len(min(2000, length(vars)))]
  pc   <- prcomp(t(mat[top, ]), scale. = TRUE)
  pct  <- round(100 * (pc$sdev^2) / sum(pc$sdev^2), 1)

  df <- data.frame(PC1 = pc$x[,1], PC2 = pc$x[,2], sample = colnames(mat))

  d  <- dist(pc$x[, 1:min(5, ncol(pc$x))])
  md <- as.matrix(d)
  mean_d <- mean(md[upper.tri(md)])
  per_sample_mean <- rowMeans(md)
  z  <- (per_sample_mean - mean(per_sample_mean)) / sd(per_sample_mean)
  outliers <- df$sample[abs(z) > 3]
  if (length(outliers) > 0)
    log_msg("WARN", sprintf("[%s][slot %d] PCA outliers (|z|>3): %s",
                            gse_id, slot_idx, paste(outliers, collapse = ", ")))

  p <- ggplot(df, aes(PC1, PC2)) +
    geom_point(size = 2, alpha = 0.7) +
    labs(title = sprintf("%s slot %d  PCA (top-var 2000)", gse_id, slot_idx),
         x = sprintf("PC1 (%.1f%%)", pct[1]),
         y = sprintf("PC2 (%.1f%%)", pct[2])) +
    theme_minimal(base_size = 11)
  ggsave(file.path(FIG_DIR, sprintf("PCA_%s_slot%d.png", gse_id, slot_idx)),
         p, width = 5, height = 4, dpi = 150)

  miss_pct <- 100 * sum(is.na(exprs(eset))) / length(exprs(eset))
  log_msg("INFO", sprintf("[%s][slot %d] missing=%.2f%%, PC1=%.1f%%, PC2=%.1f%%, meanDist=%.2f",
                          gse_id, slot_idx, miss_pct, pct[1], pct[2], mean_d))

  invisible(list(outliers = outliers, miss_pct = miss_pct, pct = pct))
}

log_msg("INFO", "==== Stage 3/3 : drift scan (PCA + outlier) ====")

for (gse in names(results)) {
  obj <- results[[gse]]
  if (is.null(obj)) next
  for (i in seq_along(obj)) {
    e <- obj[[i]]
    if (!inherits(e, "ExpressionSet")) next
    drift_scan_one(gse, i, e)
  }
}

log_msg("INFO", "Stage 3/3 complete.")
log_msg("INFO", sprintf("Figures -> %s", FIG_DIR))
log_msg("INFO", sprintf("Full log -> %s", .LOG_FILE))
