if (Sys.getenv("PROJ_ROOT") == "") {
  Sys.setenv(PROJ_ROOT = normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "."), "..")))
}

ROOT <- Sys.getenv("PROJ_ROOT")
source(file.path(ROOT, "code", "00_utils.R"))

log_msg("INFO", "================ PREFLIGHT START ================")
session_snapshot()

steps <- c("01_preflight_download.R",
           "02_structure_assert.R",
           "03_drift_scan.R")

for (s in steps) {
  log_msg("INFO", sprintf("---- exec %s ----", s))
  res <- tryCatch(source(file.path(ROOT, "code", s), echo = FALSE),
                  error = function(e) {
                    log_msg("ERROR", sprintf("[%s] aborted: %s", s, conditionMessage(e)))
                    structure(list(), class = "step_error")
                  })
  if (inherits(res, "step_error")) {
    log_msg("ERROR", sprintf("Pipeline halted at %s. Inspect log and rerun.", s))
    break
  }
}

log_msg("INFO", "================ PREFLIGHT END ================")
