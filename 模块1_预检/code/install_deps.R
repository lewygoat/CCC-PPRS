options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 600)

cran_pkgs <- c("BiocManager", "ggplot2", "matrixStats")
miss_c <- cran_pkgs[!sapply(cran_pkgs, requireNamespace, quietly = TRUE)]
if (length(miss_c) > 0) {
  cat(sprintf("[INSTALL] CRAN missing: %s\n", paste(miss_c, collapse = ", ")))
  install.packages(miss_c, Ncpus = max(1, parallel::detectCores() - 1))
} else {
  cat("[INSTALL] CRAN packages already present.\n")
}

if (!requireNamespace("BiocManager", quietly = TRUE))
  stop("BiocManager install failed")

bioc_pkgs <- c("GEOquery", "Biobase", "sva", "edgeR", "limma")
miss_b <- bioc_pkgs[!sapply(bioc_pkgs, requireNamespace, quietly = TRUE)]
if (length(miss_b) > 0) {
  cat(sprintf("[INSTALL] Bioc missing: %s\n", paste(miss_b, collapse = ", ")))
  BiocManager::install(miss_b, ask = FALSE, update = FALSE,
                       Ncpus = max(1, parallel::detectCores() - 1))
} else {
  cat("[INSTALL] Bioc packages already present.\n")
}

cat("\n[INSTALL] Final verification:\n")
for (p in c(cran_pkgs, bioc_pkgs)) {
  ok <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("  %-15s : %s\n", p, ifelse(ok, "OK", "MISSING")))
}

cat("\n[INSTALL] R version: ", R.version.string, "\n", sep = "")
cat("[INSTALL] Bioc version: ", as.character(BiocManager::version()), "\n", sep = "")
