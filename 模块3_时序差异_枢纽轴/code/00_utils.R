`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

PROJ_ROOT <- Sys.getenv("PROJ_ROOT", unset = getwd())
M3_ROOT   <- PROJ_ROOT
OUT_DIR   <- file.path(M3_ROOT, "output")
LOG_DIR   <- file.path(M3_ROOT, "log")
M1_OUT    <- Sys.getenv("M1_OUT", unset = file.path(M3_ROOT, "..", "模块1_预检", "output"))
M2_OUT    <- Sys.getenv("M2_OUT", unset = file.path(M3_ROOT, "..", "模块2_铁死亡基因集", "output"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(LOG_DIR, showWarnings = FALSE, recursive = TRUE)

.RUN_ID   <- format(Sys.time(), "%Y%m%d_%H%M%S")
.LOG_FILE <- file.path(LOG_DIR, sprintf("module3_%s.log", .RUN_ID))

log_msg <- function(level = c("INFO","WARN","ERROR","ASSERT","DRIFT","STAB"), msg) {
  level <- match.arg(level)
  line  <- sprintf("[%s][%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, msg)
  cat(line, "\n", sep = "")
  cat(line, "\n", sep = "", file = .LOG_FILE, append = TRUE)
  invisible(line)
}

assert_that <- function(cond, msg, fatal = FALSE) {
  if (isTRUE(cond)) {
    log_msg("ASSERT", sprintf("PASS  %s", msg))
    return(invisible(TRUE))
  }
  log_msg("ASSERT", sprintf("FAIL  %s", msg))
  if (fatal) stop(sprintf("FATAL ASSERT: %s", msg), call. = FALSE)
  invisible(FALSE)
}

ensure_packages <- function(cran = character(), bioc = character()) {
  miss_c <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
  miss_b <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss_c) > 0) {
    log_msg("ERROR", sprintf("Missing CRAN: %s", paste(miss_c, collapse = ", ")))
    log_msg("INFO",  sprintf('install.packages(c(%s))',
                              paste(sprintf('"%s"', miss_c), collapse = ", ")))
  }
  if (length(miss_b) > 0) {
    log_msg("ERROR", sprintf("Missing Bioc: %s", paste(miss_b, collapse = ", ")))
    log_msg("INFO",  sprintf('BiocManager::install(c(%s))',
                              paste(sprintf('"%s"', miss_b), collapse = ", ")))
  }
  if (length(miss_c) + length(miss_b) > 0)
    stop("Install missing packages first.", call. = FALSE)
  invisible(TRUE)
}

session_snapshot <- function() {
  log_msg("INFO", sprintf("Run ID: %s", .RUN_ID))
  log_msg("INFO", sprintf("R: %s", R.version.string))
  log_msg("INFO", sprintf("M3_ROOT: %s", normalizePath(M3_ROOT, mustWork = FALSE)))
  log_msg("INFO", sprintf("M1_OUT : %s", normalizePath(M1_OUT,  mustWork = FALSE)))
  log_msg("INFO", sprintf("M2_OUT : %s", normalizePath(M2_OUT,  mustWork = FALSE)))
  log_msg("INFO", sprintf("Log    : %s", .LOG_FILE))
}
