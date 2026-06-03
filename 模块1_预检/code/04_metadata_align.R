source(file.path(Sys.getenv("PROJ_ROOT", unset = getwd()), "code", "00_utils.R"))

suppressPackageStartupMessages({
  library(Biobase)
})

cache_path <- file.path(DATA_DIR, "_eset_list.rds")
if (!file.exists(cache_path)) stop("Run 01_preflight_download.R first.", call. = FALSE)
results <- readRDS(cache_path)

STAGE_RULES <- list(
  control    = c("control", "ctrl", "healthy", "\\bhc\\b", "normal", "baseline"),
  hyperacute = c("hyperacute", "ultra.?acute", "<\\s*24\\s*h", "0.?24\\s*h",
                 "^0\\s*h", "^6\\s*h", "early", "<24", "0-24", "0_24"),
  acute      = c("\\bacute\\b", "72\\s*h", "72-96", "48-72", "3.?day", "\\b3d\\b",
                 "24-72", "24_72", "day\\s*3"),
  subacute   = c("subacute", "7\\s*day", "\\b7d\\b", "168\\s*h", "day\\s*7",
                 "7-14", "1\\s*week")
)

infer_stage <- function(label) {
  s <- tolower(as.character(label))
  for (st in names(STAGE_RULES)) {
    pats <- STAGE_RULES[[st]]
    if (any(vapply(pats, function(p) grepl(p, s, perl = TRUE), logical(1))))
      return(st)
  }
  NA_character_
}

PRE_CHECK <- file.path(PROJ_ROOT, "pre_check")
dir.create(PRE_CHECK, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(PROJ_ROOT, "output"), showWarnings = FALSE, recursive = TRUE)

extract_one <- function(gse_id, slot_idx, eset, label_cols = NULL) {
  pdata <- pData(eset)
  sids  <- sampleNames(eset)

  # Path A: local loader pre-populated stage/batch columns
  if (all(c("stage", "batch") %in% colnames(pdata))) {
    raw_lab <- if ("timepoint" %in% colnames(pdata)) as.character(pdata$timepoint)
               else rep(NA_character_, length(sids))
    return(data.frame(
      sample_id     = sids,
      gse           = gse_id,
      slot          = slot_idx,
      raw_label     = raw_lab,
      inferred_stage = as.character(pdata$stage),
      batch         = as.character(pdata$batch),
      stringsAsFactors = FALSE
    ))
  }

  # Path B: GEOquery-derived pData with ch1 fields
  if (is.null(label_cols)) {
    label_cols <- intersect(
      c("time:ch1", "timepoint:ch1", "time point:ch1", "stage:ch1", "group:ch1",
        "time_point:ch1", "characteristics_ch1", "characteristics_ch1.1",
        "title", "source_name_ch1"),
      colnames(pdata))
  }

  raw_label <- if (length(label_cols) > 0) {
    apply(pdata[, label_cols, drop = FALSE], 1,
          function(x) paste(na.omit(unique(x)), collapse = " | "))
  } else rep(NA_character_, length(sids))

  batch_col <- intersect(c("batch:ch1", "scan_protocol:ch1", "submission_date",
                            "data_processing", "platform_id"),
                         colnames(pdata))
  batch <- if (length(batch_col) > 0) as.character(pdata[, batch_col[1]]) else rep("1", length(sids))

  stage <- vapply(raw_label, infer_stage, character(1))

  data.frame(
    sample_id     = sids,
    gse           = gse_id,
    slot          = slot_idx,
    raw_label     = raw_label,
    inferred_stage = stage,
    batch         = batch,
    stringsAsFactors = FALSE
  )
}

all_meta <- list()
for (gse in names(results)) {
  obj <- results[[gse]]
  if (is.null(obj)) {
    log_msg("WARN", sprintf("[%s] no data, skip metadata extraction", gse))
    next
  }
  for (i in seq_along(obj)) {
    e <- obj[[i]]
    if (!inherits(e, "ExpressionSet")) next
    df <- extract_one(gse, i, e)
    all_meta[[sprintf("%s_slot%d", gse, i)]] <- df

    pre_csv <- file.path(PRE_CHECK, sprintf("%s_slot%d_stage_mapping.csv", gse, i))
    write.csv(df, pre_csv, row.names = FALSE)
    n_na <- sum(is.na(df$inferred_stage))
    log_msg("INFO", sprintf("[%s slot %d] stage_mapping -> %s (NA=%d/%d)",
                            gse, i, pre_csv, n_na, nrow(df)))
    if (n_na > 0)
      log_msg("WARN", sprintf("[%s slot %d] %d samples have NO inferred stage; review %s and edit overrides",
                               gse, i, n_na, pre_csv))
  }
}

override_csv <- file.path(PRE_CHECK, "stage_override.csv")
if (file.exists(override_csv)) {
  ov <- read.csv(override_csv, stringsAsFactors = FALSE)
  stopifnot(all(c("sample_id","stage") %in% colnames(ov)))
  log_msg("INFO", sprintf("Applying override: %d entries", nrow(ov)))
  for (key in names(all_meta)) {
    idx <- match(all_meta[[key]]$sample_id, ov$sample_id)
    has <- !is.na(idx)
    all_meta[[key]]$inferred_stage[has] <- ov$stage[idx[has]]
  }
}

gse296 <- grep("^GSE296792", names(all_meta), value = TRUE)
if (length(gse296) > 0) {
  meta_main <- do.call(rbind, all_meta[gse296])
  meta_main <- unique(meta_main[, c("sample_id", "inferred_stage", "batch", "raw_label")])
  colnames(meta_main)[2] <- "stage"
  out_path <- file.path(PROJ_ROOT, "output", "metadata.csv")
  write.csv(meta_main, out_path, row.names = FALSE)
  log_msg("INFO", sprintf("Wrote training metadata -> %s (n=%d)", out_path, nrow(meta_main)))
  tab <- table(meta_main$stage, useNA = "ifany")
  log_msg("INFO", sprintf("Stage distribution (train): %s",
                           paste(sprintf("%s=%d", names(tab), as.integer(tab)),
                                 collapse = "; ")))
} else {
  log_msg("ERROR", "GSE296792 missing -> cannot generate training metadata")
}

gse125 <- grep("^GSE125512", names(all_meta), value = TRUE)
if (length(gse125) > 0) {
  meta_ext <- do.call(rbind, all_meta[gse125])
  meta_ext <- unique(meta_ext[, c("sample_id", "inferred_stage", "batch", "raw_label")])
  colnames(meta_ext)[2] <- "stage"
  out_path <- file.path(PROJ_ROOT, "output", "metadata_external.csv")
  write.csv(meta_ext, out_path, row.names = FALSE)
  log_msg("INFO", sprintf("Wrote external metadata -> %s (n=%d)", out_path, nrow(meta_ext)))
}

saveRDS(all_meta, file.path(DATA_DIR, "_meta_all.rds"))
log_msg("INFO", "Stage 1.4 metadata alignment complete.")
