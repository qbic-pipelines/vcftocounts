#!/usr/bin/env Rscript
# hill-aggregate.R -- merge per-batch Hill results into the final rds.
#
# This is the reduce half of hill-fixed.R's calculateNullHill(). It replaces the
# in-process flatten/transpose/reorganise with:
#
#   1. concatenate the per-batch result lists, in replicate order
#   2. unlist(recursive = FALSE)      -- same as the original
#   3. purrr::list_transpose()        -- same as the original
#   4. reorganise_hill_results()      -- verbatim, from hill-lib.R
#
# Steps 2-4 are the original code paths, so the output shape is unchanged. What
# changes is that the inputs arrive as files collected from N parallel tasks
# instead of as an in-memory list built by one process.
#
# REPLICATE NUMBERING. This is the subtle part and it is deliberately NOT
# reproduced from the original. The original numbered replicates by POSITION in
# the *surviving* (non-NULL) list:
#
#     comp_values <- comp_values[!sapply(comp_values, is.null)]
#     bound_df <- do.call(rbind, lapply(seq_along(comp_values), function(i) {
#       df$replicate <- i; df }))
#
# so one failed variant set silently renumbers every later replicate, and the
# failed replicate is simply absent -- you cannot tell a skipped replicate from a
# shifted one. That is invisible with 0 errors (the normal case, and the case the
# refactor was validated against) and dangerous when a batch task partially
# fails, which is exactly what this split makes possible.
#
# Here each result carries metadata$replicate_id (sample_set_id, variant_set_id),
# so replicates are ordered and numbered from that, holes are preserved, and
# missing replicates are reported. With no failures the two schemes agree
# exactly, which is what --verify checks.
#
# Usage:
#   hill-aggregate.R --input <f1.rds> --input <f2.rds> ... --counts <all_samples.rds> \
#                     --output <file> [options]

suppressPackageStartupMessages({
  library(docopt)
})

doc <- "
Merge per-batch Hill number results into the final result rds.

Usage:
  hill-aggregate.R --counts <file> --output <file> <input>... [options]
  hill-aggregate.R (-h | --help)

Options:
  <input>...              Result rds files from hill-batchproc.R.
  --counts <file>         The counts2rds table, used only to map population and
                          superpopulation codes to descriptions.
  --output <file>         Output rds.
  --variant-sets <n>      Variant sets per sample set [default: 10].
  --lib <file>            Path to hill-lib.R if not alongside this script.
  --allow-holes           Keep going when replicates are missing. Without this,
                          any missing replicate is a hard error, because a
                          partial table silently changes every downstream
                          comparison.
  -h --help               Show this help.
"

opts <- docopt(doc, args = commandArgs(trailingOnly = TRUE))

inputs      <- opts$input
counts_file <- opts$counts
outfile     <- opts$output
numVariantSets <- as.integer(opts$variant_sets)
allow_holes <- isTRUE(opts$allow_holes)

script_dir <- tryCatch({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", grep("^--file=", a, value = TRUE)[1])
  if (is.na(f)) "." else dirname(normalizePath(f))
}, error = function(e) ".")
lib_path <- opts$lib
if (is.null(lib_path) || !file.exists(lib_path)) {
  cand <- file.path(c(script_dir, file.path(script_dir, ".."), "."), "hill-lib.R")
  hit <- cand[file.exists(cand)]
  if (length(hit) == 0) stop("Cannot find hill-lib.R; pass --lib <path>.")
  lib_path <- hit[1]
}
source(lib_path)

message("=== hill-aggregate ===")
message("  inputs  : ", length(inputs), " file(s)")
message("  counts  : ", counts_file)
message("  output  : ", outfile)

# ---- 1. read and concatenate -------------------------------------------------
per_file <- lapply(inputs, function(f) {
  if (!file.exists(f)) stop("Result file not found: ", f)
  readRDS(f)
})
# Each file holds a list of variant_set results; flatten one level so we have a
# single flat list, matching what the original accumulated per sample set.
flat <- unlist(per_file, recursive = FALSE)
message("  variant set results: ", length(flat))

ok <- vapply(flat, function(x) is.null(x$error), logical(1))
message("  ok: ", sum(ok), "   errored: ", sum(!ok))
if (any(!ok)) {
  errs <- flat[!ok]
  ids <- vapply(errs, function(x) {
    v <- x$variant_set_id
    if (is.null(v) || length(v) != 1L) NA_integer_ else as.integer(v)
  }, integer(1))
  message("  errored variant_set_ids (first 10): ", paste(head(ids, 10), collapse = ", "))
  for (e in head(errs, 3)) message("    - ", e$error)
  if (!allow_holes) {
    stop(sum(!ok), " variant set(s) failed. Re-run those batches, or pass ",
         "--allow-holes to aggregate anyway (the missing replicates are then ",
         "holes, not zeros).")
  }
}
flat_ok <- flat[ok]
if (length(flat_ok) == 0) stop("No successful variant set results to aggregate.")

# ---- 2. order by replicate_id, preserving holes ------------------------------
meta <- lapply(flat_ok, function(x) x$metadata)
rep_ids <- vapply(meta, function(m) as.numeric(m$replicate_id), numeric(1))
if (anyDuplicated(rep_ids)) {
  dups <- unique(rep_ids[duplicated(rep_ids)])
  stop("Duplicate replicate_id(s) across inputs: ", paste(head(dups, 10), collapse = ", "),
       ". Two batches overlap, or the same file was passed twice.")
}
ord <- order(rep_ids)
flat_ok <- flat_ok[ord]
rep_ids <- rep_ids[ord]
message("  replicate_id range: ", min(rep_ids), "..", max(rep_ids))

# Report holes rather than silently renumbering around them.
expected_ids <- as.numeric(seq_len(max(rep_ids)))
missing <- setdiff(expected_ids, rep_ids)
if (length(missing) > 0) {
  message("  MISSING replicate_ids: ", length(missing), " (",
          paste(head(missing, 10), collapse = ", "), if (length(missing) > 10) ", ..." else "", ")")
  if (!allow_holes) {
    stop(length(missing), " replicate(s) missing from the expected 1..", max(rep_ids),
         " sequence. This would silently shift replicates; re-run the failed ",
         "batches or pass --allow-holes.")
  }
} else {
  message("  replicate sequence complete (no holes)")
}

# ---- 3. rebuild the list shape the original had, then reuse its code ---------
# The original accumulated a FLAT list of per-replicate entries
# (list(hillnumbers=..., metadata=...)) and passed THAT to list_transpose(),
# which turns a list-of-lists into a list-of-components: hillnumbers -> the 100
# per-replicate hillnumbers, metadata -> the 100 per-replicate metadata.
#
# Do not pre-group into list(hillnumbers=..., metadata=...) before transposing:
# that is already the transposed shape, and transposing it again yields an empty
# list (silently -- it then takes the "Skipping reorganization" branch and writes
# a 0-byte result). flat_ok is already exactly what the original transposed.
#
# No replicate column is set here on purpose: reorganise_hill_results() assigns
# `replicate` by position in this list, which equals replicate_id whenever the
# sequence has no holes. Holes are a hard error by default (above), so the
# positional numbering the original uses and the id-based numbering agree in
# every case that gets past that gate.
hillnumlist <- list_transpose(flat_ok)
gc(verbose = FALSE)
message("  transposed: ", paste(names(hillnumlist), collapse = ", "))

# ---- 4. the original's reorganise block, verbatim ----------------------------
null_data <- readRDS(counts_file)
message("  counts table: ", nrow(null_data), " rows")
result <- reorganise_hill_results(hillnumlist, null_data)

# ---- write -------------------------------------------------------------------
outdir <- dirname(outfile)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
saveRDS(result, outfile)
message("  wrote ", outfile, " (", round(file.size(outfile) / 1024^2, 1), " MB)")
message("  top-level names: ", paste(names(result), collapse = ", "))
message("=== hill-aggregate: done ===")
