source(file.path(Sys.getenv("PROJ_ROOT", unset = getwd()), "code", "00_utils.R"))

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
})

cache_path <- file.path(DATA_DIR, "_eset_list.rds")
if (!file.exists(cache_path)) stop("Run 01_preflight_download.R first.", call. = FALSE)
results <- readRDS(cache_path)

EXPECT <- list(
  GSE296792 = list(
    min_samples_total = 50,
    expected_two_matrices = TRUE,
    time_field_candidates = c("time:ch1", "timepoint:ch1", "time point:ch1",
                              "characteristics_ch1", "time_point:ch1")
  ),
  GSE125512 = list(
    min_samples_total = 20,
    expected_two_matrices = FALSE,
    time_field_candidates = c("time:ch1", "timepoint:ch1", "group:ch1",
                              "characteristics_ch1")
  ),
  GSE266873 = list(
    min_samples_total = 5,
    expected_two_matrices = FALSE,
    time_field_candidates = c("time:ch1", "timepoint:ch1", "stage:ch1",
                              "characteristics_ch1", "tissue:ch1")
  )
)

find_time_field <- function(pdata, candidates) {
  hit <- intersect(tolower(candidates), tolower(colnames(pdata)))
  if (length(hit) == 0) return(NULL)
  colnames(pdata)[tolower(colnames(pdata)) %in% hit][1]
}

scan_for_time_token <- function(pdata) {
  tokens <- c("hour","hr","\\bh\\b","day","week","acute","subacute",
              "<24","24h","72h","7d","0-6","6-24","24-48","96h")
  hits <- list()
  for (col in colnames(pdata)) {
    vals <- as.character(pdata[[col]])
    if (any(grepl(paste(tokens, collapse = "|"), vals, ignore.case = TRUE))) {
      hits[[col]] <- unique(vals)[1:min(6, length(unique(vals)))]
    }
  }
  hits
}

log_msg("INFO", "==== Stage 2/3 : structure assertion ====")

report <- list()
for (gse in names(EXPECT)) {
  exp <- EXPECT[[gse]]
  obj <- results[[gse]]
  if (is.null(obj)) {
    log_msg("ERROR", sprintf("[%s] no data, skip assertion", gse))
    next
  }

  n_mat <- length(obj)
  assert_that(n_mat >= 1, sprintf("[%s] at least one ExpressionSet present", gse))
  if (exp$expected_two_matrices) {
    pass_two <- assert_that(n_mat >= 2,
      sprintf("[%s] expects >=2 matrices (mRNA + miRNA), got %d", gse, n_mat))
    if (!pass_two)
      log_msg("WARN", sprintf("[%s] miRNA matrix may need manual download from GEO suppl files",
                              gse))
  }

  total_samples <- 0
  per_slot <- list()
  for (i in seq_along(obj)) {
    e <- obj[[i]]
    if (!inherits(e, "ExpressionSet")) next
    pdata <- pData(e)
    nsmp  <- nrow(pdata)
    total_samples <- total_samples + nsmp

    gpl <- annotation(e)
    tf  <- find_time_field(pdata, exp$time_field_candidates)
    scan_hits <- scan_for_time_token(pdata)

    log_msg("INFO", sprintf("[%s][slot %d] GPL=%s, nSample=%d", gse, i, gpl, nsmp))
    if (!is.null(tf)) {
      tab <- table(pdata[[tf]], useNA = "ifany")
      log_msg("INFO", sprintf("[%s][slot %d] time field='%s', distribution: %s",
                              gse, i, tf,
                              paste(sprintf("%s=%d", names(tab), as.integer(tab)),
                                    collapse = "; ")))
    } else if (length(scan_hits) > 0) {
      log_msg("WARN", sprintf("[%s][slot %d] no explicit time field; token scan suggests: %s",
                              gse, i, paste(names(scan_hits), collapse = ", ")))
      for (cn in names(scan_hits))
        log_msg("INFO", sprintf("    %s -> %s", cn, paste(scan_hits[[cn]], collapse = " | ")))
    } else {
      log_msg("ERROR", sprintf("[%s][slot %d] NO time-related field detected, inspect pData manually",
                               gse, i))
    }

    per_slot[[i]] <- list(gpl = gpl, n = nsmp, time_field = tf)
  }

  assert_that(total_samples >= exp$min_samples_total,
              sprintf("[%s] total samples %d >= expected min %d",
                      gse, total_samples, exp$min_samples_total))

  report[[gse]] <- list(total = total_samples, slots = per_slot)
}

saveRDS(report, file.path(DATA_DIR, "_structure_report.rds"))
log_msg("INFO", "Stage 2/3 complete.")
