if (Sys.getenv("PROJ_ROOT") == "") {
  arg <- commandArgs(trailingOnly = FALSE)
  src <- sub("^--file=", "", arg[grep("^--file=", arg)])
  here <- if (length(src) > 0) normalizePath(dirname(src)) else getwd()
  Sys.setenv(PROJ_ROOT = normalizePath(file.path(here, "..")))
}

ROOT <- Sys.getenv("PROJ_ROOT")
source(file.path(ROOT, "code", "00_utils.R"))

log_msg("INFO", "================ MODULE 1 FULL START ================")
session_snapshot()

steps <- c("01b_local_load.R",
           "02_structure_assert.R",
           "03_drift_scan.R",
           "04_metadata_align.R",
           "05_normalize_combat.R")

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

log_msg("INFO", "================ MODULE 1 FULL END ================")
