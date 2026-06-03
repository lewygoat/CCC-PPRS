options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 900)
options(Ncpus = max(1, parallel::detectCores() - 1))

Sys.setenv(DISABLE_AUTOBREW = "1")

cat("[FIX] Installing curl using system libcurl (DISABLE_AUTOBREW=1)\n")
install.packages("curl", type = "source",
                 configure.vars = "DISABLE_AUTOBREW=1")

cat("[FIX] Installing httr / httr2 / rentrez / rvest\n")
install.packages(c("httr", "httr2", "rvest"), type = "source")

cat("[FIX] Re-installing Bioc packages: GEOquery, sva, KEGGREST, AnnotationDbi\n")
BiocManager::install(
  c("KEGGREST", "AnnotationDbi", "annotate", "genefilter", "GEOquery", "sva"),
  ask = FALSE, update = FALSE
)

verify <- c("curl", "httr", "GEOquery", "sva", "AnnotationDbi", "genefilter")
status <- sapply(verify, function(p) {
  if (requireNamespace(p, quietly = TRUE)) "OK" else "MISSING"
})
cat("\n[FIX] Verification:\n")
for (i in seq_along(verify)) {
  cat(sprintf("  %-20s %s\n", verify[i], status[i]))
}
if (any(status == "MISSING")) {
  quit(status = 1)
}
