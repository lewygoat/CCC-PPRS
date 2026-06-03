suppressPackageStartupMessages({
  library(GEOquery)
})

args <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(args) >= 1) args[1] else getwd()
DATA_DIR <- file.path(ROOT, "data")
LOG_DIR <- file.path(ROOT, "log")
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(LOG_DIR, "download.log")
file.remove(LOG_FILE)
log_msg <- function(msg) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("%s | %s\n", ts, msg)
  cat(line)
  cat(line, file = LOG_FILE, append = TRUE)
}

options(timeout = 1200)
options(download.file.method = "libcurl")
Sys.setenv(VROOM_CONNECTION_SIZE = 5e7)

GSE <- "GSE266873"
GSMS <- sprintf("GSM%d", 8255340:8255348)
SUFFIXES <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
SAMPLE_TAG <- c(
  GSM8255340 = "sample1_1_ICH_2h",
  GSM8255341 = "sample1_2_ICH_3h",
  GSM8255342 = "sample1_3_ICH_2h",
  GSM8255343 = "sample2_1_ICH_12h",
  GSM8255344 = "sample2_2_ICH_12h",
  GSM8255345 = "sample2_3_ICH_12h",
  GSM8255346 = "sample3_1_ICH_26h",
  GSM8255347 = "sample3_2_ICH_28h",
  GSM8255348 = "sample3_3_ICH_32h"
)

URL_TEMPLATES <- list(
  https_ftp  = "https://ftp.ncbi.nlm.nih.gov/geo/samples/%s/%s/suppl/%s_%s_%s",
  http_ftp   = "http://ftp.ncbi.nlm.nih.gov/geo/samples/%s/%s/suppl/%s_%s_%s",
  ftp_proto  = "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/%s/%s/suppl/%s_%s_%s"
)

gsm_nnn <- function(gsm) sub("(GSM\\d+)\\d{3}$", "\\1nnn", gsm)

sha256_file <- function(p) {
  tryCatch(tools::md5sum(p), error = function(e) NA)
}

dl_one <- function(url, dest, label) {
  log_msg(sprintf("[%s] GET %s", label, url))
  t0 <- Sys.time()
  res <- tryCatch(
    download.file(url, destfile = dest, mode = "wb", quiet = TRUE),
    error = function(e) { log_msg(sprintf("[%s] ERROR %s", label, e$message)); -1 },
    warning = function(w) { log_msg(sprintf("[%s] WARN %s", label, w$message)); -2 }
  )
  if (is.null(res) || (is.numeric(res) && res < 0)) return(FALSE)
  ok <- file.exists(dest) && file.size(dest) > 100
  if (ok) {
    elapsed <- as.numeric(Sys.time() - t0, units = "secs")
    log_msg(sprintf("[%s] OK %d bytes %.1fs md5=%s",
                    label, file.size(dest), elapsed,
                    substr(sha256_file(dest), 1, 12)))
  } else {
    log_msg(sprintf("[%s] FAIL (file<100B or missing)", label))
    if (file.exists(dest)) file.remove(dest)
  }
  ok
}

dl_one_with_fallback <- function(gsm, suffix, tag) {
  dest <- file.path(DATA_DIR, sprintf("%s_%s_%s", gsm, tag, suffix))
  if (file.exists(dest) && file.size(dest) > 100) {
    log_msg(sprintf("[CACHE] using %s (%d bytes)", basename(dest), file.size(dest)))
    return(TRUE)
  }
  prefix <- gsm_nnn(gsm)
  for (proto_name in names(URL_TEMPLATES)) {
    url <- sprintf(URL_TEMPLATES[[proto_name]], prefix, gsm, gsm, tag, suffix)
    ok <- dl_one(url, dest, sprintf("%s/%s/%s", gsm, suffix, proto_name))
    if (ok) return(TRUE)
    Sys.sleep(1)
  }
  log_msg(sprintf("[%s/%s] ALL FALLBACKS FAILED", gsm, suffix))
  FALSE
}

log_msg("======================================================================")
log_msg(sprintf("Module 6 Step 1: download %s 10X supp files (%d samples x %d files)",
                GSE, length(GSMS), length(SUFFIXES)))

status_rows <- list()
for (gsm in GSMS) {
  tag <- SAMPLE_TAG[[gsm]]
  log_msg(sprintf("---- %s (%s) ----", gsm, tag))
  for (suf in SUFFIXES) {
    ok <- dl_one_with_fallback(gsm, suf, tag)
    status_rows[[length(status_rows)+1]] <- data.frame(
      gsm = gsm, tag = tag, file = suf, ok = ok,
      stringsAsFactors = FALSE
    )
  }
}

df <- do.call(rbind, status_rows)
write.csv(df, file.path(DATA_DIR, "_download_status.csv"), row.names = FALSE)
n_ok <- sum(df$ok); n_total <- nrow(df)
log_msg(sprintf("DOWNLOAD SUMMARY: %d/%d OK", n_ok, n_total))
if (n_ok < n_total) {
  log_msg("NETWORK BLOCKED — see download.log for which files failed.")
  log_msg("Will continue with whichever samples are complete (need all 3 files per sample).")
}
log_msg("STAGE-1 download phase complete.")
