suppressPackageStartupMessages({
  library(multiMiR)
})

args <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(args) >= 1) args[1] else getwd()
OUT_DIR <- file.path(ROOT, "output")
LOG_DIR <- file.path(ROOT, "log")
DATA_DIR <- file.path(ROOT, "data")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(LOG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

LOG_FILE <- file.path(LOG_DIR, "run_multiMiR.log")
log_msg <- function(msg) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("%s | %s", ts, msg)
  cat(line, "\n", sep = "")
  cat(line, "\n", file = LOG_FILE, append = TRUE, sep = "")
}

file.remove(LOG_FILE)
log_msg("======================================================================")
log_msg("Module 4 Step 1 (R/multiMiR): query 5 hub genes")

HUB <- c("GPX4", "ACSL4", "SLC7A11", "HMOX1", "FTH1")
log_msg(sprintf("hub genes: %s", paste(HUB, collapse = ", ")))

KNOWN_AXES <- list(
  GPX4    = c("hsa-miR-15a-5p", "hsa-miR-15b-5p", "hsa-miR-214-3p"),
  SLC7A11 = c("hsa-miR-27a-3p", "hsa-miR-26b-5p"),
  HMOX1   = c("hsa-miR-377-3p", "hsa-miR-24-3p"),
  ACSL4   = c("hsa-miR-424-5p"),
  FTH1    = c("hsa-miR-200b-3p")
)

EXPECTED <- list(
  validated = c(0, 600),
  predicted = c(5, 3000)
)

dbinfo <- tryCatch(multimir_dbInfoVersions(),
                   error = function(e) { log_msg(paste("dbInfoVersions FAILED:", e$message)); NULL })
if (!is.null(dbinfo)) log_msg(sprintf("remote multimir DB version: %s", dbinfo$VERSION[1]))

t0 <- Sys.time()
log_msg("[QUERY] validated tables (miRTarBase + TarBase + miRecords)")
res_val <- tryCatch(
  get_multimir(
    org = "hsa",
    target = HUB,
    table = "validated",
    summary = TRUE,
    use.tibble = FALSE
  ),
  error = function(e) { log_msg(paste("validated query FAILED:", e$message)); NULL }
)
if (is.null(res_val)) {
  log_msg("validated query returned NULL — aborting")
  quit(status = 1)
}
df_val <- res_val@data
log_msg(sprintf("[validated] raw rows=%d unique miRNAs=%d targets=%d databases=%s",
                nrow(df_val), length(unique(df_val$mature_mirna_id)),
                length(unique(df_val$target_symbol)),
                paste(sort(unique(df_val$database)), collapse=",")))
write.csv(df_val, file.path(DATA_DIR, "multimir_validated_raw.csv"), row.names = FALSE)

log_msg(sprintf("validated elapsed: %.1fs", as.numeric(Sys.time() - t0, units="secs")))

t1 <- Sys.time()
log_msg("[QUERY] predicted tables (top 20% per pair, default cutoff)")
res_pred <- tryCatch(
  get_multimir(
    org = "hsa",
    target = HUB,
    table = "predicted",
    predicted.cutoff = 20,
    predicted.cutoff.type = "p",
    predicted.site = "conserved",
    summary = TRUE,
    use.tibble = FALSE
  ),
  error = function(e) { log_msg(paste("predicted query FAILED:", e$message)); NULL }
)
if (is.null(res_pred)) {
  log_msg("predicted query returned NULL — aborting")
  quit(status = 1)
}
df_pred <- res_pred@data
log_msg(sprintf("[predicted] raw rows=%d unique miRNAs=%d targets=%d databases=%s",
                nrow(df_pred), length(unique(df_pred$mature_mirna_id)),
                length(unique(df_pred$target_symbol)),
                paste(sort(unique(df_pred$database)), collapse=",")))
write.csv(df_pred, file.path(DATA_DIR, "multimir_predicted_raw.csv"), row.names = FALSE)
log_msg(sprintf("predicted elapsed: %.1fs", as.numeric(Sys.time() - t1, units="secs")))

per_gene_per_db <- function(df, label) {
  out <- list()
  for (g in HUB) {
    sub <- df[df$target_symbol == g, ]
    by_db <- split(sub$mature_mirna_id, sub$database)
    by_db <- lapply(by_db, function(x) sort(unique(x[grepl("^hsa-", x)])))
    out[[g]] <- by_db
    rng <- EXPECTED[[label]]
    total <- length(unique(sub$mature_mirna_id))
    log_msg(sprintf("[%s] %s: total=%d  per-db=[%s]", label, g, total,
                    paste(sprintf("%s:%d", names(by_db), sapply(by_db, length)), collapse=", ")))
    if (total < rng[1] || total > rng[2]) {
      log_msg(sprintf("  WARN: %s total %d OUTSIDE expected [%d,%d]", g, total, rng[1], rng[2]))
    }
  }
  out
}

log_msg("======================================================================")
val_pg <- per_gene_per_db(df_val,  "validated")
pred_pg <- per_gene_per_db(df_pred, "predicted")

jaccard <- function(a,b){ A<-unique(a); B<-unique(b); if(length(A)+length(B)==0) return(1); if(length(A)==0||length(B)==0) return(0); length(intersect(A,B))/length(union(A,B)) }
oc      <- function(a,b){ A<-unique(a); B<-unique(b); if(length(A)==0||length(B)==0) return(0); length(intersect(A,B))/min(length(A),length(B)) }

log_msg("======================================================================")
log_msg("[CONSENSUS] canonical multiMiR rule: predicted >=2 of {TargetScan, miRDB, DIANA-microT, ElMMo, PicTar, PITA}; OR validated by any of {miRTarBase, TarBase, miRecords}")

consensus_rows <- list()
drift_rows <- list()
for (g in HUB) {
  pred_sets <- pred_pg[[g]]
  pred_sources_present <- names(pred_sets)
  if (length(pred_sources_present) < 2) {
    log_msg(sprintf("[CONSENSUS] %s: only %d predicted DB returned — relaxing to ANY",
                    g, length(pred_sources_present)))
    pred_consensus_2 <- sort(unique(unlist(pred_sets)))
  } else {
    all_pred <- unlist(pred_sets, use.names = FALSE)
    tab <- table(all_pred)
    pred_consensus_2 <- sort(names(tab[tab >= 2]))
  }
  val_set <- sort(unique(unlist(val_pg[[g]], use.names = FALSE)))
  pair_set <- sort(union(pred_consensus_2, val_set))

  if (length(pred_sources_present) >= 2) {
    pair_jac <- mean(sapply(seq_along(pred_sources_present)[-length(pred_sources_present)], function(i){
      sapply((i+1):length(pred_sources_present), function(j){
        jaccard(pred_sets[[i]], pred_sets[[j]])
      })
    }), na.rm = TRUE)
    pair_oc <- mean(sapply(seq_along(pred_sources_present)[-length(pred_sources_present)], function(i){
      sapply((i+1):length(pred_sources_present), function(j){
        oc(pred_sets[[i]], pred_sets[[j]])
      })
    }), na.rm = TRUE)
  } else {
    pair_jac <- NA; pair_oc <- NA
  }
  drift_rows[[g]] <- data.frame(
    gene = g,
    n_pred_dbs = length(pred_sources_present),
    n_val_dbs  = length(val_pg[[g]]),
    pred_total = length(unique(unlist(pred_sets))),
    val_total  = length(val_set),
    pred_consensus_ge2 = length(pred_consensus_2),
    pair_jaccard_mean = round(pair_jac, 3),
    pair_overlap_coef_mean = round(pair_oc, 3),
    final_pair_n = length(pair_set),
    stringsAsFactors = FALSE
  )

  for (m in pair_set) {
    in_val <- as.integer(m %in% val_set)
    in_pred <- as.integer(m %in% pred_consensus_2)
    consensus_rows[[length(consensus_rows)+1]] <- data.frame(
      target = g, mirna = m,
      in_validated = in_val,
      in_predicted_ge2 = in_pred,
      n_evidence = in_val + in_pred,
      stringsAsFactors = FALSE
    )
  }
}
df_drift <- do.call(rbind, drift_rows)
df_cons  <- do.call(rbind, consensus_rows)
write.csv(df_drift, file.path(OUT_DIR, "drift_per_gene.csv"), row.names = FALSE)
write.csv(df_cons,  file.path(OUT_DIR, "hub_mirna_candidates_consensus.csv"), row.names = FALSE)
log_msg("[OUTPUT] drift_per_gene.csv + hub_mirna_candidates_consensus.csv written")

sanity_rows <- list()
for (g in names(KNOWN_AXES)) {
  k <- KNOWN_AXES[[g]]
  cons_g <- df_cons$mirna[df_cons$target == g]
  val_g <- sort(unique(unlist(val_pg[[g]])))
  pred_g <- sort(unique(unlist(pred_pg[[g]])))
  any_g <- union(val_g, pred_g)
  hit_cons <- intersect(k, cons_g)
  hit_any  <- intersect(k, any_g)
  sanity_rows[[g]] <- data.frame(
    gene = g,
    known = length(k),
    in_consensus = length(hit_cons),
    in_any_source = length(hit_any),
    consensus_hits = paste(hit_cons, collapse = ";"),
    any_hits       = paste(hit_any, collapse = ";"),
    stringsAsFactors = FALSE
  )
  log_msg(sprintf("[SANITY] %s: consensus_hit=%d/%d any_hit=%d/%d",
                  g, length(hit_cons), length(k), length(hit_any), length(k)))
}
df_san <- do.call(rbind, sanity_rows)
write.csv(df_san, file.path(OUT_DIR, "biological_sanity.csv"), row.names = FALSE)

overfit <- FALSE
for (i in seq_len(nrow(df_drift))) {
  g <- df_drift$gene[i]
  shrink <- 1 - df_drift$pred_consensus_ge2[i] / max(df_drift$pred_total[i], 1)
  if (df_drift$pred_consensus_ge2[i] > 0 && shrink < 0.30) {
    log_msg(sprintf("[OVERFIT] %s: consensus/total shrinkage=%.2f (<0.30) — too permissive", g, shrink))
    overfit <- TRUE
  }
  if (df_drift$final_pair_n[i] == 0) {
    log_msg(sprintf("[OVERFIT] %s: zero final pairs — too strict / source disagreement", g))
    overfit <- TRUE
  }
}
log_msg(sprintf("[OVERFIT] flag=%s", overfit))

md <- c(
  "# Module 4 Step 1 (R/multiMiR) - run log",
  sprintf("timestamp: %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S")),
  sprintf("multiMiR DB version: %s", dbinfo$VERSION[1]),
  "",
  "## Per-gene summary",
  "| Gene | pred DBs | val DBs | pred total | val total | pred>=2 | pair Jaccard | pair OC | final pairs |",
  "|---|---|---|---|---|---|---|---|---|"
)
for (i in seq_len(nrow(df_drift))) {
  r <- df_drift[i,]
  md <- c(md, sprintf("| %s | %d | %d | %d | %d | %d | %s | %s | %d |",
                      r$gene, r$n_pred_dbs, r$n_val_dbs, r$pred_total, r$val_total,
                      r$pred_consensus_ge2,
                      ifelse(is.na(r$pair_jaccard_mean), "NA", as.character(r$pair_jaccard_mean)),
                      ifelse(is.na(r$pair_overlap_coef_mean), "NA", as.character(r$pair_overlap_coef_mean)),
                      r$final_pair_n))
}
md <- c(md, "", "## Biological sanity (known ferroptosis-related axes)",
        "| Gene | known | in consensus | in any source | consensus hits |",
        "|---|---|---|---|---|")
for (i in seq_len(nrow(df_san))) {
  r <- df_san[i,]
  md <- c(md, sprintf("| %s | %d | %d | %d | %s |",
                      r$gene, r$known, r$in_consensus, r$in_any_source,
                      ifelse(nchar(r$consensus_hits)==0, "-", r$consensus_hits)))
}
writeLines(md, file.path(OUT_DIR, "run_log.md"))
log_msg("[OUTPUT] run_log.md written")
log_msg("DONE")
