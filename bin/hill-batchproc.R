#!/usr/bin/env Rscript
# hill-batchproc.R -- process ONE batch of variant sets into Hill numbers.
#
# This is the unit of work Nextflow parallelises. It corresponds exactly to the
# body of the inner lapply in hill-fixed.R's calculateNullHill() (the
# variant_mode == "batched" branch, ~lines 2019-2058): read a batch file, run
# process_variant_set() over each variant set in it, return the list of results.
#
# Deliberately ONE batch per invocation and serial within it:
#   - Nextflow parallelises ACROSS batches, so intra-task parallelism would only
#     fight the scheduler, and `future` (and its multicore/multisession choice,
#     and the copy-on-write reasoning that goes with it) drops out of phase 2
#     entirely.
#   - A task's peak is therefore one variant_set's working set (~1.2 GB measured)
#     rather than a multiple of it, which is what makes several tasks fit on a
#     small host.
#
# Why a batch and not a single variant set: reading the batch costs 563 MiB and
# one variant_set costs ~1.2 GiB peak, so splitting finer multiplies batch reads
# without lowering the peak. Batch granularity keeps retry meaningful (10 sample
# sets = 10 retryable units instead of 1) at no extra memory.
#
# Usage:
#   hill-batchproc.R --batch <file> --variant-sets <n> --output <file> [options]
#
# Output: an .rds holding a plain list of variant_set results, each shaped
#   list(hillnumbers = <10 components>, metadata = list(sample_set_id,
#        variant_set_id, replicate_id))
# i.e. exactly what hill-fixed.R accumulates before list_transpose(). The
# aggregate step concatenates these and does the transpose/reorganise once.
#
# Failure policy: if a variant_set errors, the error is recorded IN PLACE of its
# result (same shape as the original's tryCatch, which returned
# list(sample_set_id, variant_set_id, error)) rather than aborting the task.
# A single bad variant set should cost one replicate, not a whole batch. Set
# --strict to abort instead, which is what you want when validating a refactor.

suppressPackageStartupMessages({
  library(docopt)
})

doc <- "
Process one pre-generated batch of variant sets into Hill numbers.

Usage:
  hill-batchproc.R --batch <file> --output <file> [options]
  hill-batchproc.R (-h | --help)

Options:
  --batch <file>          Batch file from hill-batchgen.R (chunk_downsample_<sample_set>_batch_<n>.rds).
  --output <file>         Output rds: a list of variant_set results.
  --variant-sets <n>      Variant sets per sample set, used to compute
                          replicate_id as (sample_set_id - 1) * n +
                          variant_set_id. [default: 10]
  --q <num>               Hill number order [default: 1].
  --lib <file>            Path to hill-lib.R, if not alongside this script.
  --strict                Abort the task if any variant set errors. Default is
                          to record the error in place and carry on.
  -h --help               Show this help.
"

opts <- docopt(doc, args = commandArgs(trailingOnly = TRUE))

batch_file   <- opts$batch
outfile      <- opts$output
numVariantSets <- as.integer(opts$variant_sets)
qHillNumber  <- as.numeric(opts$q)
strict       <- isTRUE(opts$strict)

# Fail loudly rather than letting integer(0)/numeric(0) propagate into the maths.
for (nm in c("numVariantSets", "qHillNumber")) {
  v <- get(nm)
  if (length(v) != 1L || is.na(v)) stop("Could not read ", nm, " from the CLI options.")
}
if (numVariantSets < 1L) stop("--variant-sets must be >= 1")

# Locate hill-lib.R. Nextflow stages it next to this script via `script:` +
# `path`, so default to the script's own directory and walk up one level so the
# same file works when run from bin/ and from a staged task dir.
script_dir <- tryCatch({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", grep("^--file=", a, value = TRUE)[1])
  if (is.na(f)) "." else dirname(normalizePath(f))
}, error = function(e) ".")

lib_path <- opts$lib
if (is.null(lib_path) || !file.exists(lib_path)) {
  candidates <- file.path(c(script_dir, file.path(script_dir, ".."), "."), "hill-lib.R")
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    stop("Cannot find hill-lib.R. Looked in: ", paste(candidates, collapse = ", "),
         ". Pass --lib <path>.")
  }
  lib_path <- hit[1]
}
source(lib_path)

if (!file.exists(batch_file)) stop("Batch file not found: ", batch_file)

message("=== hill-batchproc ===")
message("  batch        : ", batch_file, " (", round(file.size(batch_file) / 1024^2), " MiB)")
message("  output       : ", outfile)
message("  lib          : ", lib_path)
message("  q            : ", qHillNumber)
message("  variant-sets : ", numVariantSets)

batch_data <- readRDS(batch_file)
if (!is.list(batch_data) || length(batch_data) == 0) {
  stop("Batch file does not hold a non-empty list: ", batch_file)
}
message("  variant sets in batch: ", length(batch_data))

# sample_set_id is carried on each variant_set, so the index file is not needed
# here -- one less input to stage, and one less thing to get out of step.
sample_set_ids <- unique(vapply(batch_data, function(vs) as.integer(vs$sample_set_id), integer(1)))
if (length(sample_set_ids) != 1L) {
  stop("Batch holds variant sets from ", length(sample_set_ids),
       " sample sets (", paste(sample_set_ids, collapse = ", "),
       "); a batch is expected to belong to exactly one.")
}
sample_set_id <- sample_set_ids[1]
message("  sample set   : ", sample_set_id)

t_start <- Sys.time()
results <- vector("list", length(batch_data))
n_error <- 0L

for (i in seq_along(batch_data)) {
  vs <- batch_data[[i]]
  vs_id <- as.integer(vs$variant_set_id)
  t0 <- Sys.time()
  results[[i]] <- tryCatch(
    process_variant_set(vs, sample_set_id, numVariantSets, qHillNumber),
    error = function(e) {
      if (strict) stop("variant_set ", vs_id, " failed: ", e$message)
      n_error <<- n_error + 1L
      warning("variant_set ", vs_id, " (sample set ", sample_set_id,
              ") failed: ", e$message)
      # Same shape the original used for a failed variant set, so the aggregate
      # step can spot a hole rather than silently dropping a replicate.
      list(sample_set_id = sample_set_id,
           variant_set_id = vs_id,
           error = e$message)
    }
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ok <- is.null(results[[i]]$error)
  message(sprintf("    variant_set %d/%d (id %d): %s in %.1fs",
                  i, length(batch_data), vs_id,
                  if (ok) "ok" else "ERROR", elapsed))
}

# Read the batch back out of memory before writing.
rm(batch_data)
gc(verbose = FALSE)

outdir <- dirname(outfile)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, outfile)

total <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
message(sprintf("  %d/%d variant sets ok, %d errored, %.0fs total, wrote %s (%.2f MB)",
                length(results) - n_error, length(results), n_error, total, outfile,
                file.size(outfile) / 1024^2))
if (n_error > 0L) {
  message("  NOTE: ", n_error, " variant set(s) errored and were recorded in place. ",
          "These become holes in the replicate table; the aggregate step reports them.")
}
message("=== hill-batchproc: done ===")
