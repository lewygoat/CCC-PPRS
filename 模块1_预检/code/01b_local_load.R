source(file.path(Sys.getenv("PROJ_ROOT", unset = getwd()), "code", "00_utils.R"))

suppressPackageStartupMessages({
  library(Biobase)
})

log_msg("INFO", "==== Stage 1/3-bis : local supplementary loader ====")
log_msg("INFO", "Bypassing getGEO due to DNS hijack on this host. Reading pre-downloaded suppl files.")

read_count_csv_gz <- function(path) {
  conn <- gzfile(path, "rt")
  on.exit(close(conn))
  df <- read.csv(conn, header = TRUE, check.names = FALSE,
                 stringsAsFactors = FALSE, row.names = 1)
  m <- as.matrix(df)
  storage.mode(m) <- "numeric"
  m
}

parse_pdata_296792 <- function(sample_ids, title_lookup = NULL) {
  parts <- strsplit(sample_ids, "_", fixed = TRUE)
  patient <- vapply(parts, `[`, character(1), 1)
  tp      <- vapply(parts, `[`, character(1), 2)
  stage_map <- c(base = "control", bas = "control",
                 "72h" = "acute",
                 "7day" = "subacute", "7da" = "subacute")
  stage <- stage_map[tp]
  stage[is.na(stage)] <- NA_character_
  data.frame(
    sample_id = sample_ids,
    patient   = patient,
    timepoint = tp,
    stage     = unname(stage),
    batch     = "GSE296792",
    row.names = sample_ids,
    stringsAsFactors = FALSE
  )
}

build_eset <- function(mat, pdata, gpl_label) {
  common <- intersect(colnames(mat), rownames(pdata))
  mat   <- mat[, common, drop = FALSE]
  pdata <- pdata[common, , drop = FALSE]
  adf <- AnnotatedDataFrame(data = pdata)
  eset <- ExpressionSet(assayData = mat, phenoData = adf,
                        annotation = gpl_label)
  eset
}

results <- list()

mrna_path  <- file.path(DATA_DIR, "GSE296792_Normalized_mRNA_counts_v1.csv.gz")
mirna_path <- file.path(DATA_DIR, "GSE296792_Normalized_miRNA_counts_v1.csv.gz")
ok_296 <- file.exists(mrna_path) && file.exists(mirna_path)
assert_that(ok_296, "GSE296792 supplementary mRNA+miRNA CSV.gz present in data_raw/")

if (ok_296) {
  log_msg("INFO", "[GSE296792] reading mRNA matrix...")
  mrna_mat <- read_count_csv_gz(mrna_path)
  log_msg("INFO", sprintf("[GSE296792] mRNA matrix: %d genes x %d samples", nrow(mrna_mat), ncol(mrna_mat)))

  log_msg("INFO", "[GSE296792] reading miRNA matrix...")
  mirna_mat <- read_count_csv_gz(mirna_path)
  log_msg("INFO", sprintf("[GSE296792] miRNA matrix: %d miRNAs x %d samples", nrow(mirna_mat), ncol(mirna_mat)))

  common <- intersect(colnames(mrna_mat), colnames(mirna_mat))
  log_msg("INFO", sprintf("[GSE296792] paired samples between mRNA & miRNA: %d", length(common)))

  pd_m <- parse_pdata_296792(colnames(mrna_mat))
  pd_i <- parse_pdata_296792(colnames(mirna_mat))

  eset_mrna  <- build_eset(mrna_mat,  pd_m, "GPL_GSE296792_mRNA")
  eset_mirna <- build_eset(mirna_mat, pd_i, "GPL_GSE296792_miRNA_mir")
  results[["GSE296792"]] <- list(eset_mrna, eset_mirna)

  tab_m <- table(pd_m$stage, useNA = "ifany")
  log_msg("INFO", sprintf("[GSE296792] mRNA stage distribution: %s",
                          paste(sprintf("%s=%d", names(tab_m), as.integer(tab_m)), collapse = "; ")))
}

ext_path <- file.path(DATA_DIR, "GSE125512_SecondSample_vs_FirstSample.all.gene.result.xls.gz")
if (file.exists(ext_path)) {
  log_msg("INFO", "[GSE125512] reading external dataset (xls.gz, actually TSV)...")
  conn <- gzfile(ext_path, "rt")
  raw <- read.delim(conn, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  close(conn)
  log_msg("INFO", sprintf("[GSE125512] raw shape: %d rows x %d cols", nrow(raw), ncol(raw)))

  ich_cols <- grep("^ICH[0-9]+$", colnames(raw), value = TRUE)
  log_msg("INFO", sprintf("[GSE125512] ICH expression columns detected: %d  (%s ... %s)",
                          length(ich_cols), ich_cols[1], ich_cols[length(ich_cols)]))

  symbol <- raw[["symbol"]]
  if (is.null(symbol)) symbol <- raw[[2]]
  dup <- duplicated(symbol) | is.na(symbol)
  ext_mat <- as.matrix(raw[!dup, ich_cols, drop = FALSE])
  storage.mode(ext_mat) <- "numeric"
  rownames(ext_mat) <- symbol[!dup]
  log_msg("INFO", sprintf("[GSE125512] expression matrix: %d genes x %d samples (dups dropped)",
                          nrow(ext_mat), ncol(ext_mat)))

  # External dataset: paired design "FirstSample" (24h) vs "SecondSample" (7day) within same patient
  # The provided file collapses to one row per (patient,gene) — not enough resolution for stage.
  # We treat each ICH column as one observation, all marked 'acute' (24h post-ICH) per series description.
  # Stage refinement requires per-sample series matrix raw counts (not in this file).
  pd_ext <- data.frame(
    sample_id = ich_cols,
    patient   = ich_cols,
    timepoint = "24h_or_7day_unresolved",
    stage     = "acute",
    batch     = "GSE125512",
    row.names = ich_cols,
    stringsAsFactors = FALSE
  )
  log_msg("WARN", "[GSE125512] suppl file is DE-result with collapsed columns; stage set to 'acute' as proxy. Module 7 external validation should re-derive per-sample data from raw FASTQ or contact author.")
  eset_ext <- build_eset(ext_mat, pd_ext, "GPL_GSE125512_mRNA")
  results[["GSE125512"]] <- list(eset_ext)
} else {
  log_msg("WARN", "[GSE125512] suppl file not present, skipping external dataset")
}

sm266 <- file.path(DATA_DIR, "GSE266873_series_matrix.txt.gz")
if (file.exists(sm266)) {
  log_msg("INFO", "[GSE266873] scRNA dataset — only series_matrix metadata loaded (single-cell, processed elsewhere in Module 9)")
}

saveRDS(results, file.path(DATA_DIR, "_eset_list.rds"))
log_msg("INFO", sprintf("Cached ExpressionSet list -> %s (datasets: %s)",
                        file.path(DATA_DIR, "_eset_list.rds"),
                        paste(names(results), collapse = ", ")))
log_msg("INFO", "Stage 1/3-bis complete.")
