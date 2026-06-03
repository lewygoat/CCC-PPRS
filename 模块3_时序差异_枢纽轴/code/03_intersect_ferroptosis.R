source(file.path(Sys.getenv("PROJ_ROOT", unset = getwd()), "code", "00_utils.R"))

inputs <- readRDS(file.path(OUT_DIR, "_inputs_cache.rds"))
fit_cache <- readRDS(file.path(OUT_DIR, "_limma_cache.rds"))
de_all <- fit_cache$de_all

LFC_CUT  <- 0.585
PADJ_CUT <- 0.05
JACCARD_WARN <- 0.50

intersect_one <- function(de_tbl, pool, pool_name) {
  sig <- subset(de_tbl,
                !is.na(adj.P.Val) & adj.P.Val < PADJ_CUT &
                  abs(logFC) > LFC_CUT)
  sig_genes <- toupper(sig$gene)
  hit <- intersect(sig_genes, pool)
  data.frame(
    contrast = unique(de_tbl$contrast),
    pool     = pool_name,
    n_sig    = length(sig_genes),
    n_pool   = length(pool),
    n_hit    = length(hit),
    enrichment = if (length(sig_genes) > 0 && length(pool) > 0)
                   length(hit) / length(sig_genes) else NA_real_,
    stringsAsFactors = FALSE
  )
}

summary_rows <- list()
hit_lists    <- list()

for (cname in names(de_all)) {
  tbl <- de_all[[cname]]

  r_all <- intersect_one(tbl, inputs$pool_all, "all_1059")
  r_hi  <- intersect_one(tbl, inputs$pool_hi,  "high_conf_92")
  summary_rows[[paste0(cname, "_all")]] <- r_all
  summary_rows[[paste0(cname, "_hi")]]  <- r_hi

  sig <- subset(tbl, !is.na(adj.P.Val) & adj.P.Val < PADJ_CUT & abs(logFC) > LFC_CUT)
  sig$gene <- toupper(sig$gene)
  hit_all <- sig[sig$gene %in% inputs$pool_all, c("gene","logFC","adj.P.Val","contrast")]
  hit_hi  <- sig[sig$gene %in% inputs$pool_hi,  c("gene","logFC","adj.P.Val","contrast")]
  hit_all$pool <- "all_1059"
  hit_hi$pool  <- "high_conf_92"
  hit_lists[[cname]] <- rbind(hit_all, hit_hi)

  set_all <- unique(hit_all$gene)
  set_hi  <- unique(hit_hi$gene)
  jac <- if (length(union(set_all, set_hi)) > 0)
           length(intersect(set_all, set_hi)) / length(union(set_all, set_hi)) else NA_real_
  log_msg("STAB", sprintf("[%s] Jaccard(all vs high_conf) = %.3f (n_all=%d, n_hi=%d)",
                           cname, jac, length(set_all), length(set_hi)))
  if (!is.na(jac) && jac < JACCARD_WARN)
    log_msg("DRIFT", sprintf("[%s] Sensitivity divergence: hub gene list strongly depends on pool choice (Jaccard %.2f<%.2f). Pool-stable genes are safer for downstream signature.",
                              cname, jac, JACCARD_WARN))
}

summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(OUT_DIR, "stage_DE_FRG_summary.csv"), row.names = FALSE)
log_msg("INFO", sprintf("Saved DE-FRG summary -> %s",
                        file.path(OUT_DIR, "stage_DE_FRG_summary.csv")))

hits_df <- do.call(rbind, hit_lists)
write.csv(hits_df, file.path(OUT_DIR, "stage_DE_FRG_hits.csv"), row.names = FALSE)
log_msg("INFO", sprintf("Saved DE-FRG hit list -> %s",
                        file.path(OUT_DIR, "stage_DE_FRG_hits.csv")))

stable_genes <- list()
for (cname in names(de_all)) {
  pool_stable <- intersect(
    hits_df$gene[hits_df$contrast == cname & hits_df$pool == "all_1059"],
    hits_df$gene[hits_df$contrast == cname & hits_df$pool == "high_conf_92"]
  )
  stable_genes[[cname]] <- pool_stable
  log_msg("STAB", sprintf("[%s] Pool-stable DE-FRG = %d", cname, length(pool_stable)))
}
saveRDS(stable_genes, file.path(OUT_DIR, "_pool_stable_genes.rds"))

log_msg("INFO", "Stage 3.2 intersection complete.")
