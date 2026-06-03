suppressPackageStartupMessages({
  if (!exists(".PROJ_ROOT")) {
    .PROJ_ROOT <- normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "."), ".."), mustWork = FALSE)
  }
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

PROJ_ROOT <- Sys.getenv("PROJ_ROOT", unset = getwd())
DATA_DIR  <- file.path(PROJ_ROOT, "data_raw")
LOG_DIR   <- file.path(PROJ_ROOT, "log")
CODE_DIR  <- file.path(PROJ_ROOT, "code")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(LOG_DIR,  showWarnings = FALSE, recursive = TRUE)

.LOG_FILE <- file.path(LOG_DIR, sprintf("preflight_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S")))

log_msg <- function(level = c("INFO","WARN","ERROR","ASSERT"), msg) {
  level <- match.arg(level)
  line  <- sprintf("[%s][%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, msg)
  cat(line, "\n", sep = "")
  cat(line, "\n", sep = "", file = .LOG_FILE, append = TRUE)
  invisible(line)
}

assert_that <- function(cond, msg) {
  if (isTRUE(cond)) {
    log_msg("ASSERT", sprintf("PASS  %s", msg))
    return(invisible(TRUE))
  }
  log_msg("ASSERT", sprintf("FAIL  %s", msg))
  invisible(FALSE)
}

with_retry <- function(expr, max_try = 3, sleep_sec = 5, label = "task") {
  for (i in seq_len(max_try)) {
    res <- tryCatch(eval.parent(substitute(expr)),
                    error = function(e) structure(list(err = e), class = "retry_err"))
    if (!inherits(res, "retry_err")) {
      if (i > 1) log_msg("INFO", sprintf("[%s] succeeded on attempt %d", label, i))
      return(res)
    }
    log_msg("WARN", sprintf("[%s] attempt %d/%d failed: %s",
                            label, i, max_try, conditionMessage(res$err)))
    if (i < max_try) Sys.sleep(sleep_sec * i)
  }
  log_msg("ERROR", sprintf("[%s] all %d attempts failed", label, max_try))
  stop(sprintf("with_retry(%s) exhausted", label), call. = FALSE)
}

ensure_packages <- function(cran = character(), bioc = character()) {
  miss_cran <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
  miss_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss_cran) > 0) {
    log_msg("WARN", sprintf("Missing CRAN packages: %s", paste(miss_cran, collapse = ", ")))
    log_msg("INFO", sprintf("install.packages(c(%s))",
                             paste(sprintf('"%s"', miss_cran), collapse = ", ")))
  }
  if (length(miss_bioc) > 0) {
    log_msg("WARN", sprintf("Missing Bioconductor packages: %s", paste(miss_bioc, collapse = ", ")))
    log_msg("INFO", sprintf('BiocManager::install(c(%s))',
                             paste(sprintf('"%s"', miss_bioc), collapse = ", ")))
  }
  ok <- length(miss_cran) == 0 && length(miss_bioc) == 0
  if (!ok) stop("Required packages missing, see log for install commands.", call. = FALSE)
  invisible(TRUE)
}

cache_load <- function(gse_id, loader_fn) {
  cache_file <- file.path(DATA_DIR, sprintf("%s.rds", gse_id))
  if (file.exists(cache_file)) {
    log_msg("INFO", sprintf("[%s] cache hit -> %s", gse_id, cache_file))
    return(readRDS(cache_file))
  }
  log_msg("INFO", sprintf("[%s] cache miss, fetching...", gse_id))
  obj <- with_retry(loader_fn(), max_try = 3, sleep_sec = 10, label = gse_id)
  saveRDS(obj, cache_file)
  log_msg("INFO", sprintf("[%s] cached -> %s (%.1f MB)",
                          gse_id, cache_file, file.size(cache_file)/1024^2))
  obj
}

session_snapshot <- function() {
  log_msg("INFO", sprintf("R version: %s", R.version.string))
  log_msg("INFO", sprintf("Platform: %s", R.version$platform))
  log_msg("INFO", sprintf("PROJ_ROOT: %s", PROJ_ROOT))
  log_msg("INFO", sprintf("Log file: %s", .LOG_FILE))
}
