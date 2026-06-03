source(file.path(Sys.getenv("PROJ_ROOT", unset = getwd()), "code", "00_utils.R"))

EXPECTED_META_COLS  <- c("sample_id", "stage")
EXPECTED_STAGES     <- c("control", "hyperacute", "acute", "subacute")
MIN_N_PER_STAGE     <- 3
MIN_GENES_MRNA      <- 5000
MIN_FERROPTOSIS_N   <- 100

p_mrna     <- file.path(M1_OUT, "expr_mrna.csv")
p_mirna    <- file.path(M1_OUT, "expr_mirna.csv")
p_meta     <- file.path(M1_OUT, "metadata.csv")
p_pool_all <- file.path(M2_OUT, "ferroptosis_geneset.csv")
p_pool_hi  <- file.path(M2_OUT, "ferroptosis_geneset_high_confidence.csv")

assert_that(file.exists(p_mrna),  sprintf("M1 expr_mrna.csv exists  (%s)", p_mrna),  fatal = TRUE)
assert_that(file.exists(p_meta),  sprintf("M1 metadata.csv exists   (%s)", p_meta),  fatal = TRUE)
assert_that(file.exists(p_pool_all), sprintf("M2 全集 exists           (%s)", p_pool_all), fatal = TRUE)
assert_that(file.exists(p_pool_hi),  sprintf("M2 高置信集 exists      (%s)", p_pool_hi),  fatal = TRUE)
if (!file.exists(p_mirna))
  log_msg("WARN", "miRNA matrix absent; module 3 mRNA-only branch will proceed")

expr_mrna <- as.matrix(read.csv(p_mrna, row.names = 1, check.names = FALSE))
meta      <- read.csv(p_meta, stringsAsFactors = FALSE, check.names = FALSE)
pool_all  <- read.csv(p_pool_all, stringsAsFactors = FALSE)
pool_hi   <- read.csv(p_pool_hi,  stringsAsFactors = FALSE)

assert_that(all(EXPECTED_META_COLS %in% colnames(meta)),
            sprintf("metadata cols include %s", paste(EXPECTED_META_COLS, collapse = ",")),
            fatal = TRUE)
assert_that(nrow(expr_mrna) >= MIN_GENES_MRNA,
            sprintf("mRNA features >= %d (got %d)", MIN_GENES_MRNA, nrow(expr_mrna)),
            fatal = TRUE)

common <- intersect(colnames(expr_mrna), meta$sample_id)
assert_that(length(common) >= 10,
            sprintf("samples shared between expr_mrna and metadata >= 10 (got %d)", length(common)),
            fatal = TRUE)
expr_mrna <- expr_mrna[, common, drop = FALSE]
meta      <- meta[match(common, meta$sample_id), , drop = FALSE]

meta$stage <- tolower(trimws(meta$stage))
unknown <- setdiff(unique(meta$stage), EXPECTED_STAGES)
if (length(unknown) > 0)
  log_msg("WARN", sprintf("Unknown stages dropped: %s", paste(unknown, collapse = ", ")))
keep <- meta$stage %in% EXPECTED_STAGES
expr_mrna <- expr_mrna[, keep, drop = FALSE]
meta      <- meta[keep, , drop = FALSE]

stage_tab <- table(meta$stage)
log_msg("INFO", sprintf("Stage distribution: %s",
                        paste(sprintf("%s=%d", names(stage_tab), as.integer(stage_tab)),
                              collapse = "; ")))

assert_that("control" %in% names(stage_tab),
            "control samples present (required for stage-vs-control contrasts)",
            fatal = TRUE)
for (s in setdiff(EXPECTED_STAGES, "control")) {
  if (!s %in% names(stage_tab)) {
    log_msg("WARN", sprintf("Stage '%s' missing; corresponding contrast will be skipped", s))
    next
  }
  assert_that(stage_tab[[s]] >= MIN_N_PER_STAGE,
              sprintf("Stage '%s' has >= %d samples (got %d)", s, MIN_N_PER_STAGE, stage_tab[[s]]))
}
assert_that(stage_tab[["control"]] >= MIN_N_PER_STAGE,
            sprintf("Stage 'control' has >= %d samples (got %d)",
                    MIN_N_PER_STAGE, stage_tab[["control"]]),
            fatal = TRUE)

pool_all_v <- unique(toupper(trimws(pool_all$symbol %||% pool_all[[1]])))
pool_hi_v  <- unique(toupper(trimws(pool_hi$symbol  %||% pool_hi[[1]])))
assert_that(length(pool_all_v) >= MIN_FERROPTOSIS_N,
            sprintf("Ferroptosis pool size >= %d (got %d)", MIN_FERROPTOSIS_N, length(pool_all_v)),
            fatal = TRUE)
log_msg("INFO", sprintf("Ferroptosis pools loaded: all=%d, high_conf=%d",
                        length(pool_all_v), length(pool_hi_v)))

rownames(expr_mrna) <- toupper(rownames(expr_mrna))
hit_all <- length(intersect(rownames(expr_mrna), pool_all_v))
log_msg("INFO", sprintf("Ferroptosis pool ∩ mRNA matrix : %d / %d (coverage %.1f%%)",
                        hit_all, length(pool_all_v), 100 * hit_all / length(pool_all_v)))
assert_that(hit_all >= 50,
            sprintf("Ferroptosis-mRNA overlap >= 50 (got %d) — pool/expression annotation mismatch?",
                    hit_all))

inputs <- list(expr_mrna = expr_mrna, meta = meta,
               pool_all = pool_all_v, pool_hi = pool_hi_v,
               stage_tab = stage_tab)
saveRDS(inputs, file.path(OUT_DIR, "_inputs_cache.rds"))
log_msg("INFO", "Stage 3.0 input load complete.")
