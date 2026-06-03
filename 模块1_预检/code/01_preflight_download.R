source(file.path(Sys.getenv("PROJ_ROOT", unset = getwd()), "code", "00_utils.R"))

ensure_packages(bioc = c("GEOquery", "Biobase"))

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
})

options(timeout = 600)
Sys.setenv(VROOM_CONNECTION_SIZE = 5e7)

GSE_IDS <- c("GSE296792", "GSE125512", "GSE266873")

fetch_geo <- function(gse_id) {
  cache_load(gse_id, function() {
    getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = FALSE, destdir = DATA_DIR)
  })
}

inspect_eset_list <- function(gse_id, eset_list) {
  if (!is.list(eset_list) || length(eset_list) == 0) {
    log_msg("ERROR", sprintf("[%s] empty getGEO return", gse_id))
    return(invisible(NULL))
  }
  for (i in seq_along(eset_list)) {
    e <- eset_list[[i]]
    if (!inherits(e, "ExpressionSet")) {
      log_msg("WARN", sprintf("[%s] slot %d not ExpressionSet (class=%s)",
                              gse_id, i, paste(class(e), collapse = "/")))
      next
    }
    log_msg("INFO", sprintf("[%s] slot %d : GPL=%s, nSample=%d, nFeature=%d",
                            gse_id, i,
                            annotation(e),
                            ncol(exprs(e)),
                            nrow(exprs(e))))
  }
}

session_snapshot()
log_msg("INFO", "==== Stage 1/3 : GEO download ====")

results <- list()
for (gse in GSE_IDS) {
  log_msg("INFO", sprintf(">>> %s", gse))
  res <- tryCatch(fetch_geo(gse),
                  error = function(e) {
                    log_msg("ERROR", sprintf("[%s] giving up: %s", gse, conditionMessage(e)))
                    NULL
                  })
  results[[gse]] <- res
  if (!is.null(res)) inspect_eset_list(gse, res)
}

ok_ids <- names(results)[!vapply(results, is.null, logical(1))]
fail_ids <- setdiff(GSE_IDS, ok_ids)
log_msg("INFO", sprintf("Downloaded OK : %s", paste(ok_ids, collapse = ", ")))
if (length(fail_ids) > 0)
  log_msg("ERROR", sprintf("Download FAIL : %s -- rerun script after fixing network", paste(fail_ids, collapse = ", ")))

saveRDS(results, file.path(DATA_DIR, "_eset_list.rds"))
log_msg("INFO", "Stage 1/3 complete.")
