if (Sys.getenv("PROJ_ROOT") == "") {
  arg <- commandArgs(trailingOnly = FALSE)
  src <- sub("^--file=", "", arg[grep("^--file=", arg)])
  here <- if (length(src) > 0) normalizePath(dirname(src)) else getwd()
  Sys.setenv(PROJ_ROOT = normalizePath(file.path(here, "..")))
}

ROOT  <- Sys.getenv("PROJ_ROOT")
HERE  <- file.path(ROOT, "code")
source(file.path(HERE, "00_utils.R"))

log_msg("INFO", "================ MODULE 3 (1/3) START ================")
session_snapshot()

steps <- c("01_load_inputs.R",
           "02_limma_three_stage.R",
           "03_intersect_ferroptosis.R")

for (s in steps) {
  log_msg("INFO", sprintf("---- exec %s ----", s))
  res <- tryCatch(source(file.path(HERE, s), echo = FALSE),
                  error = function(e) {
                    log_msg("ERROR", sprintf("[%s] aborted: %s", s, conditionMessage(e)))
                    structure(list(), class = "step_error")
                  })
  if (inherits(res, "step_error")) {
    log_msg("ERROR", sprintf("Pipeline halted at %s. Inspect log and rerun.", s))
    break
  }
}

log_msg("INFO", "================ MODULE 3 (1/3) END ================")
