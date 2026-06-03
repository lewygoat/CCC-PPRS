suppressPackageStartupMessages({
  library(limma)
  library(matrixStats)
})

args <- commandArgs(trailingOnly = TRUE)
M5_ROOT <- Sys.getenv("M5_ROOT")
M1_OUT  <- Sys.getenv("M1_OUT")
M2_OUT  <- Sys.getenv("M2_OUT")
M3_OUT  <- Sys.getenv("M3_OUT")
stopifnot(nzchar(M5_ROOT), nzchar(M1_OUT), nzchar(M2_OUT), nzchar(M3_OUT))

N_BOOT_TOTAL <- 333L
CHECKPOINT_EVERY <- 50L
ADJ_P_CUT <- 0.05
LOG_PATH <- file.path(M5_ROOT, "log",
  sprintf("module5_boot_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S")))
CHECKPOINT_PATH <- file.path(M5_ROOT, "output", "_boot_checkpoint.rds")
RESULT_PATH <- file.path(M5_ROOT, "output", "bootstrap_stability.csv")
DRIFT_PATH  <- file.path(M5_ROOT, "output", "bootstrap_drift_monitor.csv")
HUB_PATH    <- file.path(M5_ROOT, "output", "hub_candidate_frequency.csv")

log_msg <- function(tag, msg) {
  line <- sprintf("[%s][%s] %s",
                  format(Sys.time(), "%Y-%m-%d %H:%M:%S"), tag, msg)
  cat(line, "\n", sep = "", file = LOG_PATH, append = TRUE)
  cat(line, "\n", sep = "")
}

dir.create(dirname(LOG_PATH), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(CHECKPOINT_PATH), showWarnings = FALSE, recursive = TRUE)

log_msg("INFO", "================ MODULE 5 BOOT 1/3 START ================")
log_msg("INFO", sprintf("R: %s", R.version.string))
log_msg("INFO", sprintf("N_BOOT_TOTAL=%d, CHECKPOINT_EVERY=%d, ADJ_P_CUT=%.3f",
                        N_BOOT_TOTAL, CHECKPOINT_EVERY, ADJ_P_CUT))

read_csv_safe <- function(p) {
  stopifnot(file.exists(p))
  read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
}

expr_path <- file.path(M1_OUT, "expr_mrna.csv")
meta_path <- file.path(M1_OUT, "metadata.csv")
pool_all_path <- file.path(M2_OUT, "ferroptosis_geneset.csv")
pool_hi_path  <- file.path(M2_OUT, "ferroptosis_geneset_high_confidence.csv")
hits_path     <- file.path(M3_OUT, "stage_DE_FRG_hits.csv")

log_msg("INFO", sprintf("Load expr  : %s", expr_path))
log_msg("INFO", sprintf("Load meta  : %s", meta_path))
log_msg("INFO", sprintf("Load poolA : %s", pool_all_path))
log_msg("INFO", sprintf("Load poolH : %s", pool_hi_path))
log_msg("INFO", sprintf("Load m3hit : %s", hits_path))

expr_df <- read_csv_safe(expr_path)
meta    <- read_csv_safe(meta_path)
pool_all <- read_csv_safe(pool_all_path)
pool_hi  <- read_csv_safe(pool_hi_path)
m3_hits  <- read_csv_safe(hits_path)

stopifnot(ncol(expr_df) >= 5L)
genes <- as.character(expr_df[[1L]])
expr_mat <- as.matrix(expr_df[, -1L, drop = FALSE])
rownames(expr_mat) <- genes
log_msg("INFO", sprintf("First 3 genes parsed: %s",
                        paste(head(genes, 3), collapse = ",")))
log_msg("INFO", sprintf("expr_mat: %d genes x %d samples",
                        nrow(expr_mat), ncol(expr_mat)))

stopifnot(all(c("sample_id", "stage") %in% colnames(meta)))
meta <- meta[meta$sample_id %in% colnames(expr_mat), , drop = FALSE]
expr_mat <- expr_mat[, meta$sample_id, drop = FALSE]
log_msg("INFO", sprintf("Aligned: %d samples", ncol(expr_mat)))

stage_tab <- table(meta$stage)
log_msg("INFO", sprintf("Stage distribution: %s",
  paste(names(stage_tab), stage_tab, sep = "=", collapse = "; ")))

if (!"control" %in% names(stage_tab))
  stop("control group missing in metadata; abort.")

stage_levels <- intersect(c("acute", "subacute"), names(stage_tab))
if (length(stage_levels) < 1L)
  stop("No comparable stage (acute/subacute) found; abort.")
log_msg("INFO", sprintf("Contrasts to run: %s",
  paste(sprintf("%s_vs_control", stage_levels), collapse = ", ")))

pool_gene_col <- function(df) {
  for (cand in c("symbol", "gene", "Gene", "Symbol", "GeneSymbol", "gene_symbol")) {
    if (cand %in% colnames(df)) return(cand)
  }
  colnames(df)[1]
}
pool_a_genes <- unique(toupper(pool_all[[pool_gene_col(pool_all)]]))
pool_h_genes <- unique(toupper(pool_hi[[pool_gene_col(pool_hi)]]))
log_msg("INFO", sprintf("Pool sizes: all=%d, high_conf=%d",
                        length(pool_a_genes), length(pool_h_genes)))

orig_genes <- list()
for (lv in stage_levels) {
  sub_a <- m3_hits[m3_hits$contrast == paste0(lv, "_vs_control") &
                     m3_hits$pool == "all_1059", "gene", drop = TRUE]
  sub_h <- m3_hits[m3_hits$contrast == paste0(lv, "_vs_control") &
                     m3_hits$pool == "high_conf_92", "gene", drop = TRUE]
  orig_genes[[lv]] <- list(all = unique(sub_a),
                           hi  = unique(sub_h),
                           stable = intersect(sub_a, sub_h))
  log_msg("INFO", sprintf("Original M3 hits [%s_vs_control]: all=%d, hi=%d, stable=%d",
    lv, length(orig_genes[[lv]]$all),
    length(orig_genes[[lv]]$hi),
    length(orig_genes[[lv]]$stable)))
}

design_hub_path <- file.path(M5_ROOT, "..", "..", "选题B_技术路线.md")
hub_design <- c("GPX4", "ACSL4", "SLC7A11", "HMOX1", "FTH1")
hub_m3_pool_stable <- unique(unlist(lapply(orig_genes, `[[`, "stable")))
hub_track <- unique(c(hub_design, hub_m3_pool_stable))
hub_track <- intersect(hub_track, rownames(expr_mat))
log_msg("INFO", sprintf("HUB tracking list (n=%d): %s",
  length(hub_track), paste(hub_track, collapse = ",")))

run_limma_once <- function(emat, mvec) {
  grp <- factor(mvec, levels = c("control", unique(setdiff(mvec, "control"))))
  if (length(levels(grp)) < 2L) return(NULL)
  design <- model.matrix(~0 + grp)
  colnames(design) <- levels(grp)
  fit <- tryCatch(
    suppressWarnings(lmFit(emat, design, method = "ls")),
    error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  out <- list()
  for (lv in setdiff(levels(grp), "control")) {
    ct_str <- sprintf("%s - control", lv)
    ct_mat <- tryCatch(makeContrasts(contrasts = ct_str, levels = design),
                       error = function(e) NULL)
    if (is.null(ct_mat)) next
    fit2 <- tryCatch(suppressWarnings(contrasts.fit(fit, ct_mat)),
                     error = function(e) NULL)
    if (is.null(fit2)) next
    fit2 <- tryCatch(suppressWarnings(eBayes(fit2, robust = FALSE)),
                     error = function(e) NULL)
    if (is.null(fit2)) next
    tt <- topTable(fit2, number = Inf, sort.by = "none")
    out[[lv]] <- data.frame(gene = rownames(tt),
                            logFC = tt$logFC,
                            adj.P.Val = tt$adj.P.Val,
                            stringsAsFactors = FALSE)
  }
  out
}

if (file.exists(CHECKPOINT_PATH)) {
  ckpt <- readRDS(CHECKPOINT_PATH)
  done_b <- ckpt$last_b
  log_msg("INFO", sprintf("RESUME from checkpoint: last_b=%d", done_b))
} else {
  ckpt <- list(
    last_b = 0L,
    freq_pool_a = lapply(stage_levels, function(x) integer(0L)),
    freq_pool_h = lapply(stage_levels, function(x) integer(0L)),
    logfc_pool_a = lapply(stage_levels, function(x) list()),
    overlap_a = lapply(stage_levels, function(x) numeric(0L)),
    overlap_h = lapply(stage_levels, function(x) numeric(0L)),
    hub_logfc = lapply(stage_levels, function(x) {
      m <- matrix(NA_real_, nrow = 0, ncol = length(hub_track))
      colnames(m) <- hub_track; m
    }),
    hub_sig = lapply(stage_levels, function(x) {
      m <- matrix(NA, nrow = 0, ncol = length(hub_track))
      colnames(m) <- hub_track; m
    }),
    drift_log = data.frame(),
    elapsed = numeric(0L)
  )
  names(ckpt$freq_pool_a) <- stage_levels
  names(ckpt$freq_pool_h) <- stage_levels
  names(ckpt$logfc_pool_a) <- stage_levels
  names(ckpt$overlap_a) <- stage_levels
  names(ckpt$overlap_h) <- stage_levels
  names(ckpt$hub_logfc) <- stage_levels
  names(ckpt$hub_sig) <- stage_levels
  log_msg("INFO", "Initialized fresh state.")
}

ALL_GENES <- rownames(expr_mat)
freq_a <- lapply(stage_levels, function(x) {
  if (length(ckpt$freq_pool_a[[x]]) == 0L) {
    v <- integer(length(ALL_GENES)); names(v) <- ALL_GENES; v
  } else ckpt$freq_pool_a[[x]]
})
freq_h <- lapply(stage_levels, function(x) {
  if (length(ckpt$freq_pool_h[[x]]) == 0L) {
    v <- integer(length(ALL_GENES)); names(v) <- ALL_GENES; v
  } else ckpt$freq_pool_h[[x]]
})
names(freq_a) <- stage_levels
names(freq_h) <- stage_levels

logfc_acc <- lapply(stage_levels, function(x) {
  if (length(ckpt$logfc_pool_a[[x]]) == 0L) list() else ckpt$logfc_pool_a[[x]]
})
names(logfc_acc) <- stage_levels

overlap_a <- lapply(stage_levels, function(x) ckpt$overlap_a[[x]])
overlap_h <- lapply(stage_levels, function(x) ckpt$overlap_h[[x]])
names(overlap_a) <- stage_levels
names(overlap_h) <- stage_levels

hub_logfc <- ckpt$hub_logfc
hub_sig   <- ckpt$hub_sig
elapsed   <- ckpt$elapsed
drift_log <- ckpt$drift_log

start_idx <- ckpt$last_b + 1L
log_msg("INFO", sprintf("Starting from bootstrap iter %d", start_idx))

DRIFT_STOP <- FALSE
DRIFT_OVERLAP_BATCH_MEAN_FLOOR <- 0.20
DRIFT_RUNAWAY_TIME_SEC <- 90

batch_overlap_a <- numeric(0L)
batch_overlap_h <- numeric(0L)

for (b in start_idx:N_BOOT_TOTAL) {
  set.seed(20260530L + b)
  iter_t0 <- Sys.time()

  boot_idx <- integer(0L)
  for (lv in unique(meta$stage)) {
    pool_lv <- which(meta$stage == lv)
    if (length(pool_lv) == 0L) next
    boot_idx <- c(boot_idx, sample(pool_lv, length(pool_lv), replace = TRUE))
  }
  emat_b <- expr_mat[, boot_idx, drop = FALSE]
  mvec_b <- meta$stage[boot_idx]

  row_var <- matrixStats::rowVars(emat_b)
  keep_rows <- !is.na(row_var) & row_var > 1e-8
  if (sum(keep_rows) < 100L) {
    log_msg("WARN", sprintf("Iter %d: too few variable genes (%d), skip",
                            b, sum(keep_rows)))
    next
  }
  emat_b <- emat_b[keep_rows, , drop = FALSE]

  res <- run_limma_once(emat_b, mvec_b)
  if (is.null(res)) {
    log_msg("WARN", sprintf("Iter %d: limma failed, skip", b))
    next
  }

  for (lv in stage_levels) {
    rr <- res[[lv]]
    if (is.null(rr)) next
    rr_a <- rr[rr$gene %in% pool_a_genes & rr$adj.P.Val < ADJ_P_CUT, , drop = FALSE]
    rr_h <- rr[rr$gene %in% pool_h_genes & rr$adj.P.Val < ADJ_P_CUT, , drop = FALSE]
    if (nrow(rr_a) > 0L) {
      hits_a <- rr_a$gene
      freq_a[[lv]][hits_a] <- freq_a[[lv]][hits_a] + 1L
    } else hits_a <- character(0L)
    if (nrow(rr_h) > 0L) {
      hits_h <- rr_h$gene
      freq_h[[lv]][hits_h] <- freq_h[[lv]][hits_h] + 1L
    } else hits_h <- character(0L)

    for (g in rr_a$gene) {
      if (is.null(logfc_acc[[lv]][[g]])) logfc_acc[[lv]][[g]] <- numeric(0L)
      logfc_acc[[lv]][[g]] <- c(logfc_acc[[lv]][[g]],
                                rr_a$logFC[rr_a$gene == g])
    }

    orig_a <- orig_genes[[lv]]$all
    orig_h <- orig_genes[[lv]]$hi
    oc_a <- if (length(orig_a) > 0L && length(hits_a) > 0L) {
      length(intersect(orig_a, hits_a)) /
        min(length(orig_a), length(hits_a))
    } else NA_real_
    oc_h <- if (length(orig_h) > 0L && length(hits_h) > 0L) {
      length(intersect(orig_h, hits_h)) /
        min(length(orig_h), length(hits_h))
    } else NA_real_
    overlap_a[[lv]] <- c(overlap_a[[lv]], oc_a)
    overlap_h[[lv]] <- c(overlap_h[[lv]], oc_h)

    in_hub <- intersect(hub_track, rr$gene)
    new_logfc <- setNames(rep(NA_real_, length(hub_track)), hub_track)
    new_sig   <- setNames(rep(NA, length(hub_track)), hub_track)
    if (length(in_hub) > 0L) {
      sub <- rr[match(in_hub, rr$gene), ]
      new_logfc[in_hub] <- sub$logFC
      new_sig[in_hub] <- sub$adj.P.Val < ADJ_P_CUT
    }
    hub_logfc[[lv]] <- rbind(hub_logfc[[lv]], new_logfc)
    hub_sig[[lv]]   <- rbind(hub_sig[[lv]], new_sig)
  }

  iter_t <- as.numeric(difftime(Sys.time(), iter_t0, units = "secs"))
  elapsed <- c(elapsed, iter_t)
  if (iter_t > DRIFT_RUNAWAY_TIME_SEC) {
    log_msg("DRIFT", sprintf("Iter %d: runtime %.1fs > %ds (slowdown)",
                              b, iter_t, DRIFT_RUNAWAY_TIME_SEC))
  }

  if (b %% CHECKPOINT_EVERY == 0L || b == N_BOOT_TOTAL) {
    log_msg("INFO", sprintf("===== checkpoint at iter %d / %d (mean iter time=%.2fs) =====",
                            b, N_BOOT_TOTAL, mean(elapsed, na.rm = TRUE)))
    ckpt$last_b <- b
    ckpt$freq_pool_a <- freq_a
    ckpt$freq_pool_h <- freq_h
    ckpt$logfc_pool_a <- logfc_acc
    ckpt$overlap_a <- overlap_a
    ckpt$overlap_h <- overlap_h
    ckpt$hub_logfc <- hub_logfc
    ckpt$hub_sig   <- hub_sig
    ckpt$drift_log <- drift_log
    ckpt$elapsed   <- elapsed
    saveRDS(ckpt, CHECKPOINT_PATH)

    for (lv in stage_levels) {
      recent_a <- tail(overlap_a[[lv]], CHECKPOINT_EVERY)
      mean_a <- mean(recent_a, na.rm = TRUE)
      log_msg("INFO", sprintf("Overlap [%s_vs_control] all-pool batch mean=%.3f",
                              lv, mean_a))
      drift_log <- rbind(drift_log, data.frame(
        iter = b, contrast = paste0(lv, "_vs_control"),
        pool = "all", batch_overlap_mean = mean_a))
      if (!is.na(mean_a) && mean_a < DRIFT_OVERLAP_BATCH_MEAN_FLOOR) {
        log_msg("DRIFT", sprintf("Overlap [%s_vs_control] batch mean %.3f < floor %.2f",
                                  lv, mean_a, DRIFT_OVERLAP_BATCH_MEAN_FLOOR))
      }
    }
  }

  if (DRIFT_STOP) {
    log_msg("DRIFT", "DRIFT_STOP triggered; aborting before completion.")
    break
  }
}

log_msg("INFO", "===== aggregating results =====")

stab_rows <- list()
for (lv in stage_levels) {
  for (pool_tag in c("all", "high_conf")) {
    freq <- if (pool_tag == "all") freq_a[[lv]] else freq_h[[lv]]
    sel <- freq[freq > 0L]
    sel <- sort(sel, decreasing = TRUE)
    if (length(sel) == 0L) next
    for (g in names(sel)) {
      if (pool_tag == "all" && !is.null(logfc_acc[[lv]][[g]]) &&
          length(logfc_acc[[lv]][[g]]) > 0L) {
        lfc <- logfc_acc[[lv]][[g]]
        ci_lo <- quantile(lfc, 0.025, na.rm = TRUE)
        ci_hi <- quantile(lfc, 0.975, na.rm = TRUE)
        med   <- median(lfc, na.rm = TRUE)
      } else {
        ci_lo <- NA_real_; ci_hi <- NA_real_; med <- NA_real_
      }
      stab_rows[[length(stab_rows) + 1L]] <- data.frame(
        contrast = paste0(lv, "_vs_control"),
        pool = pool_tag,
        gene = g,
        select_freq = unname(sel[g]),
        select_rate = unname(sel[g]) / N_BOOT_TOTAL,
        logfc_median = med,
        logfc_ci_lo = ci_lo,
        logfc_ci_hi = ci_hi
      )
    }
  }
}

if (length(stab_rows) > 0L) {
  stab_df <- do.call(rbind, stab_rows)
  rownames(stab_df) <- NULL
  write.csv(stab_df, RESULT_PATH, row.names = FALSE)
  log_msg("INFO", sprintf("Wrote stability table: %d rows -> %s",
                          nrow(stab_df), RESULT_PATH))
} else {
  log_msg("WARN", "No stability rows to write.")
}

drift_rows <- list()
for (lv in stage_levels) {
  oa <- overlap_a[[lv]]
  oh <- overlap_h[[lv]]
  drift_rows[[length(drift_rows) + 1L]] <- data.frame(
    contrast = paste0(lv, "_vs_control"),
    n_iter_recorded_all = length(oa),
    overlap_all_mean = mean(oa, na.rm = TRUE),
    overlap_all_sd   = sd(oa, na.rm = TRUE),
    overlap_high_mean = mean(oh, na.rm = TRUE),
    overlap_high_sd   = sd(oh, na.rm = TRUE)
  )
}
write.csv(do.call(rbind, drift_rows), DRIFT_PATH, row.names = FALSE)
log_msg("INFO", sprintf("Wrote drift monitor -> %s", DRIFT_PATH))

hub_rows <- list()
for (lv in stage_levels) {
  m_lfc <- hub_logfc[[lv]]
  m_sig <- hub_sig[[lv]]
  for (g in hub_track) {
    lfcs <- m_lfc[, g]
    sigs <- m_sig[, g]
    n_iter <- sum(!is.na(lfcs))
    n_sig  <- sum(sigs, na.rm = TRUE)
    hub_rows[[length(hub_rows) + 1L]] <- data.frame(
      contrast = paste0(lv, "_vs_control"),
      hub = g,
      n_iter = n_iter,
      sig_freq = n_sig,
      sig_rate = if (n_iter > 0L) n_sig / n_iter else NA_real_,
      logfc_median = if (n_iter > 0L) median(lfcs, na.rm = TRUE) else NA_real_,
      logfc_ci_lo = if (n_iter > 0L) quantile(lfcs, 0.025, na.rm = TRUE) else NA_real_,
      logfc_ci_hi = if (n_iter > 0L) quantile(lfcs, 0.975, na.rm = TRUE) else NA_real_
    )
  }
}
write.csv(do.call(rbind, hub_rows), HUB_PATH, row.names = FALSE)
log_msg("INFO", sprintf("Wrote hub frequency table -> %s", HUB_PATH))

log_msg("INFO", "================ MODULE 5 BOOT 1/3 END ================")
