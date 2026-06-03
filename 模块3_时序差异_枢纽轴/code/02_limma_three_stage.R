source(file.path(Sys.getenv("PROJ_ROOT", unset = getwd()), "code", "00_utils.R"))

ensure_packages(cran = c("ggplot2"),
                bioc = c("limma", "sva"))

suppressPackageStartupMessages({
  library(limma)
  library(sva)
  library(ggplot2)
})

inputs <- readRDS(file.path(OUT_DIR, "_inputs_cache.rds"))
expr   <- inputs$expr_mrna
meta   <- inputs$meta
stages <- inputs$stage_tab

LFC_CUT  <- 0.585
PADJ_CUT <- 0.05
USE_SVA  <- TRUE
USE_ARRAY_WEIGHTS <- TRUE

stage_f <- factor(meta$stage, levels = intersect(c("control","hyperacute","acute","subacute"),
                                                  unique(meta$stage)))
design <- model.matrix(~ 0 + stage_f)
colnames(design) <- levels(stage_f)
log_msg("INFO", sprintf("Design matrix rank=%d, n=%d, p=%d",
                        qr(design)$rank, nrow(design), ncol(design)))

if (USE_SVA) {
  mod0 <- model.matrix(~ 1, data = meta)
  n.sv <- tryCatch(num.sv(expr, design, method = "be"),
                   error = function(e) { log_msg("WARN", sprintf("num.sv failed: %s", conditionMessage(e))); 0 })
  log_msg("DRIFT", sprintf("SVA-detected hidden surrogate variables n.sv = %d", n.sv))
  if (n.sv > 0 && n.sv < ncol(design)) {
    sva_fit <- tryCatch(sva(expr, design, mod0, n.sv = n.sv),
                        error = function(e) { log_msg("WARN", sprintf("sva failed: %s", conditionMessage(e))); NULL })
    if (!is.null(sva_fit)) {
      design <- cbind(design, sva_fit$sv)
      colnames(design)[(ncol(design)-n.sv+1):ncol(design)] <- sprintf("SV%d", seq_len(n.sv))
      log_msg("DRIFT", sprintf("Added %d surrogate vars to design as drift defense", n.sv))
    }
  } else if (n.sv >= ncol(design)) {
    log_msg("WARN", "n.sv >= rank(design); skip SVA to avoid collinearity")
  }
}

aw <- NULL
if (USE_ARRAY_WEIGHTS) {
  aw <- tryCatch(arrayWeights(expr, design),
                 error = function(e) { log_msg("WARN", sprintf("arrayWeights failed: %s", conditionMessage(e))); NULL })
  if (!is.null(aw)) {
    z <- (aw - mean(aw)) / sd(aw)
    susp <- meta$sample_id[abs(z) > 2]
    if (length(susp) > 0)
      log_msg("DRIFT", sprintf("Low-weight samples (|z|>2): %s",
                                paste(sprintf("%s(w=%.2f)", susp, aw[abs(z) > 2]),
                                      collapse = "; ")))
  }
}

fit <- lmFit(expr, design, weights = aw, method = "robust")

contrast_pairs <- list()
for (s in setdiff(levels(stage_f), "control"))
  contrast_pairs[[sprintf("%s_vs_control", s)]] <- sprintf("%s - control", s)

cont_mat <- makeContrasts(contrasts = unlist(contrast_pairs), levels = design)
colnames(cont_mat) <- names(contrast_pairs)

fit2 <- contrasts.fit(fit, cont_mat)
fit2 <- eBayes(fit2, robust = TRUE)

de_all <- list()
for (cname in colnames(cont_mat)) {
  tt <- topTable(fit2, coef = cname, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  tt$contrast <- cname
  de_all[[cname]] <- tt

  n_up   <- sum(tt$adj.P.Val < PADJ_CUT & tt$logFC >  LFC_CUT, na.rm = TRUE)
  n_down <- sum(tt$adj.P.Val < PADJ_CUT & tt$logFC < -LFC_CUT, na.rm = TRUE)
  log_msg("INFO", sprintf("[%s] DE: up=%d, down=%d (|lfc|>%.3f, FDR<%.2g)",
                          cname, n_up, n_down, LFC_CUT, PADJ_CUT))

  pvals <- tt$P.Value[!is.na(tt$P.Value)]
  prop_lt05 <- mean(pvals < 0.05)
  if (prop_lt05 > 0.4)
    log_msg("DRIFT", sprintf("[%s] p<0.05 proportion = %.2f (>0.4 suggests inflation, check confounders)",
                              cname, prop_lt05))
  if (prop_lt05 < 0.06)
    log_msg("DRIFT", sprintf("[%s] p<0.05 proportion = %.2f (<0.06 suggests low signal/sensitivity)",
                              cname, prop_lt05))

  png(file.path(LOG_DIR, sprintf("pHist_%s.png", cname)), width = 600, height = 400, res = 100)
  hist(pvals, breaks = 50, main = sprintf("p-value histogram: %s", cname),
       xlab = "p-value", col = "grey80", border = "white")
  dev.off()
}

de_long <- do.call(rbind, de_all)
write.csv(de_long, file.path(OUT_DIR, "stage_vs_control_DE_full.csv"), row.names = FALSE)
log_msg("INFO", sprintf("Saved full DE results to %s",
                        file.path(OUT_DIR, "stage_vs_control_DE_full.csv")))

saveRDS(list(de_all = de_all, design = design, aw = aw),
        file.path(OUT_DIR, "_limma_cache.rds"))
log_msg("INFO", "Stage 3.1 limma complete.")
