# ============================================================================
# common.R  --  shared setup for every script in this repository
# ============================================================================
# Sourced at the top of each scripts/*.R. Provides:
#   ROOT, DATA, META, CURATED, RESULTS, FIGURES   (absolute paths, repo-relative)
#   need_pkgs()        load packages, stop with a clear message if missing
#   banner(), kv()     sanity-check printing helpers
#   check_expected()   compare a computed value to metadata/expected_values.csv
#                      and append the outcome to results/reproducibility_check.csv
#   write_session()    dump sessionInfo() to results/sessionInfo_<script>.txt
#   GENO_COLORS, COND_LEVELS, DATASET1_GENOTYPES ...
#
# Repo root is found from the location of this file, so the scripts can be run
# from any working directory:   Rscript scripts/01_dataset1_DESeq2.R
# ============================================================================

ROOT <- local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f)) {
    normalizePath(file.path(dirname(f), ".."))          # scripts/x.R -> repo root
  } else if (file.exists("R/common.R")) {
    normalizePath(".")                                   # interactive, cwd = root
  } else if (file.exists("../R/common.R")) {
    normalizePath("..")                                  # interactive, cwd = scripts/
  } else stop("Cannot locate repository root. Run from the repo root or via Rscript scripts/<name>.R")
})
DATA    <- file.path(ROOT, "data")
META    <- file.path(ROOT, "metadata")
CURATED <- file.path(ROOT, "curated_genes")
RESULTS <- file.path(ROOT, "results")
FIGURES <- file.path(ROOT, "figures")
for (d in c(RESULTS, FIGURES)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Sub-folders used by the pipeline
RES_D1     <- file.path(RESULTS, "dataset1")           # DESeq2 outputs, dataset 1
RES_D1_ENR <- file.path(RESULTS, "dataset1_enrichment")
RES_D1_WG  <- file.path(RESULTS, "dataset1_WGCNA")
RES_D2     <- file.path(RESULTS, "dataset2")
RES_HM     <- file.path(RESULTS, "heatmap_matrices")   # z-score / marker tables behind heatmaps
for (d in c(RES_D1, RES_D1_ENR, RES_D1_WG, RES_D2, RES_HM)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- packages --------------------------------------------------------------
# Optional: extra library trees searched AFTER the default ones (so they never
# shadow newer system packages), e.g. an annotation package installed elsewhere:
#   export PLLP_R_EXTRA_LIBS=~/R/x86_64-pc-linux-gnu-library/4.1
local({
  extra <- Sys.getenv("PLLP_R_EXTRA_LIBS")
  if (nzchar(extra)) .libPaths(c(.libPaths(), path.expand(strsplit(extra, ":")[[1]])))
})

need_pkgs <- function(...) {
  pk <- c(...)
  miss <- pk[!vapply(pk, function(p) suppressWarnings(requireNamespace(p, quietly = TRUE)), logical(1))]
  if (length(miss))
    stop("Missing R packages: ", paste(miss, collapse = ", "),
         "\nInstall with BiocManager::install(c(", paste0('"', miss, '"', collapse = ", "), "))")
  for (p in pk) suppressPackageStartupMessages(library(p, character.only = TRUE))
  invisible(pk)
}

# ---- printing helpers -----------------------------------------------------
banner <- function(txt) cat("\n", strrep("=", 76), "\n", txt, "\n", strrep("=", 76), "\n", sep = "")
kv     <- function(k, v) cat(sprintf("  %-46s %s\n", k, paste(v, collapse = ", ")))

# ---- reproducibility assertions -------------------------------------------
# metadata/expected_values.csv holds the numbers printed in the paper / saved in
# the original run. Every script calls check_expected() for the values it can
# recompute; results accumulate in results/reproducibility_check.csv.
EXPECTED <- if (file.exists(file.path(META, "expected_values.csv")))
  read.csv(file.path(META, "expected_values.csv"), stringsAsFactors = FALSE) else NULL

check_expected <- function(key, value, script = SCRIPT_NAME, tolerance = 0) {
  exp_row <- if (!is.null(EXPECTED)) EXPECTED[EXPECTED$key == key, ] else NULL
  if (is.null(exp_row) || nrow(exp_row) == 0) {
    cat(sprintf("  [no expectation] %-40s = %s\n", key, value)); return(invisible(NA))
  }
  expv <- exp_row$expected[1]
  ok <- if (is.numeric(value) && !is.na(suppressWarnings(as.numeric(expv)))) {
    abs(value - as.numeric(expv)) <= tolerance
  } else identical(as.character(value), as.character(expv))
  cat(sprintf("  [%s] %-40s got %s | published %s%s\n", if (ok) "OK  " else "DIFF", key,
              value, expv, if (!ok && nzchar(exp_row$note[1])) paste0("  (", exp_row$note[1], ")") else ""))
  f <- file.path(RESULTS, "reproducibility_check.csv")
  row <- data.frame(script = script, key = key, published = expv, reproduced = as.character(value),
                    match = ok, checked = format(Sys.time(), "%Y-%m-%d %H:%M"), stringsAsFactors = FALSE)
  if (file.exists(f)) {
    old <- read.csv(f, stringsAsFactors = FALSE)
    old <- old[!(old$script == script & old$key == key), ]
    row <- rbind(old, row)
  }
  write.csv(row, f, row.names = FALSE)
  invisible(ok)
}

write_session <- function(script = SCRIPT_NAME) {
  f <- file.path(RESULTS, paste0("sessionInfo_", script, ".txt"))
  writeLines(c(paste("Run:", format(Sys.time())), capture.output(sessionInfo())), f)
  cat("sessionInfo written to", f, "\n")
}

# Verify every expected output exists and is non-empty (per lab convention).
verify_files <- function(paths) {
  ok <- TRUE
  for (p in paths) {
    sz <- if (file.exists(p)) file.size(p) else NA
    good <- !is.na(sz) && sz > 0; ok <- ok && good
    cat(sprintf("  %-4s %-60s %s\n", if (good) "OK" else "FAIL", basename(p),
                if (is.na(sz)) "MISSING" else sprintf("%.2f MB", sz / 1e6)))
  }
  if (!ok) stop("One or more expected output files are missing or empty.")
  invisible(ok)
}

# ---- constants shared across figures --------------------------------------
GENO_COLORS <- c("EK" = "#E8A87C", "pd1116" = "#8CD493", "pd1245" = "#789BCA")
GENO_TEXT   <- c("EK" = "#C0672A", "pd1116" = "#3E9A47", "pd1245" = "#3D63A0")
DIET_COLORS <- c("HP" = "#FF6B6B", "LP" = "#9370DB")
COND_LEVELS <- c("EK_HP", "EK_LP", "pd1116_HP", "pd1116_LP", "pd1245_HP", "pd1245_LP")
ALPHA <- 0.05      # FDR threshold for DEGs
LFC   <- 1         # |log2 fold-change| threshold for DEGs

# Name of the calling script (for logs); scripts may override.
SCRIPT_NAME <- local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f)) sub("\\.R$", "", basename(f)) else "interactive"
})

cat("common.R loaded | repo root:", ROOT, "| script:", SCRIPT_NAME, "\n")
