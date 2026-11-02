#!/usr/bin/env Rscript
# hill-batchgen.R -- phase 1: generate the Monte Carlo batches, then stop.
#
# This is the "generate" half of hill-fixed.R's calculateNullHill(), split out so
# the expensive, embarrassingly-parallel per-batch Hill calculation can be run by
# the scheduler instead of inside one process.
#
# It does NOT reimplement chunk generation. It reproduces the exact preamble the
# original's chunking branch performs (the block around lines 1740-1885 of
# hill-fixed.R) and then calls the same generate_hill_chunks() from hill-lib.R.
# The whole point is that the batches this writes are byte-for-byte the batches
# the monolith would have written, so hill-batchproc.R + hill-aggregate.R
# reproduce the monolith's numbers -- which they do (verified: every component
# of the aggregated result matches the original with maxdiff 0).
#
# Determinism: the original derives the per-replicate seeds from --seed plus
# sample_set_id and variant_set_id, and generate_hill_chunks() does that itself.
# Nothing here needs its own seed, and adding one would change the draws.
#
# Usage:
#   hill-batchgen.R --input <counts.rds> --chunks <dir> --output <manifest.rds> [options]

suppressPackageStartupMessages({
  library(docopt)
})

doc <- "
Generate the Monte Carlo null-distribution batches for the Hill diversity step.

Usage:
  hill-batchgen.R --chunks <dir> --output <file> <input> [options]
  hill-batchgen.R (-h | --help)

Options:
  <input>                 counts2rds rds, or the output of hill.R --mode prepare.
  --chunks <dir>          Directory to write chunk files into.
  --output <file>         Manifest rds listing the batch files, one entry per
                          batch, for the scheduler to fan out over.
  --metadata <file>       Sample metadata rds/csv. Optional if the count table
                          already carries the metadata columns.
  --pipelines <list>      Comma-separated count columns. Default: every numeric
                          column that is not sample metadata.
  --q <num>               Hill number order [default: 1].
  --pop                   Evaluate by population.
  --superpop              Evaluate by superpopulation.
  --sample-sets <n>       Number of sample sets [default: 10].
  --variant-sets <n>      Number of variant sets per sample set [default: 10].
  --seed <n>              Random seed [default: 1234].
  --resume                Reuse existing chunk files where possible.
  --max-gb <n>            future.globals.maxSize in GB [default: 100].
  --max-variants <n>      Cap numvariants. 0 = no cap. Only needed while there is
                          no separate filtered arm. [default: 0]
  --lib <file>            Path to hill-lib.R if not alongside this script.
  -h --help               Show this help.
"

opts <- docopt(doc, args = commandArgs(trailingOnly = TRUE))

input      <- opts$input
chunkDir   <- opts$chunks
outfile    <- opts$output
qHillNumber <- as.numeric(opts$q)
numSampleSets  <- as.integer(opts$sample_sets)
numVariantSets <- as.integer(opts$variant_sets)
seed       <- as.integer(opts$seed)
maxGB      <- as.numeric(opts$max_gb)
maxVariants <- as.integer(opts$max_variants)
evaluatePopulation      <- isTRUE(opts$pop)
evaluateSuperPopulation <- isTRUE(opts$superpop)
resume     <- isTRUE(opts$resume)

for (nm in c("qHillNumber","numSampleSets","numVariantSets","seed","maxGB")) {
  v <- get(nm)
  if (length(v) != 1L || is.na(v)) stop("Could not read ", nm, " from the CLI options.")
}

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

if (!dir.exists(chunkDir)) dir.create(chunkDir, recursive = TRUE, showWarnings = FALSE)

message("=== hill-batchgen ===")
message("  input        : ", input)
message("  chunks       : ", chunkDir)
message("  sample-sets  : ", numSampleSets)
message("  variant-sets : ", numVariantSets)

# ---- replicate the original's preamble --------------------------------------
# Mirrors run_null(): accept either the --mode prepare list or a bare
# counts2rds data.frame, and apply --max-variants the same way.
prep <- read_any(input)
if (is.data.frame(prep)) {
  filtered_data  <- prep
  null_data      <- prep
  numvariants    <- length(unique(prep$variant_called))
  total_variants <- unique(prep$variant_called)
  if (maxVariants > 0L && numvariants > maxVariants) {
    message("  capping numvariants: ", numvariants, " -> ", maxVariants, " (--max-variants)")
    numvariants <- maxVariants
  }
} else if (!is.list(prep) || is.null(prep$null_data)) {
  stop("--input must be a counts2rds data.frame or the output of hill.R --mode prepare.")
} else {
  filtered_data  <- prep$filtered_data
  null_data      <- prep$null_data
  numvariants    <- prep$numvariants
  total_variants <- prep$total_variants
}
message("  null data    : ", nrow(null_data), " rows x ", ncol(null_data), " cols")
message("  numvariants  : ", numvariants)

# ---- generate ----------------------------------------------------------------
# Delegates to generate_batches(), which is hill-fixed.R's chunking branch wrapped
# verbatim (see hill-lib.R). Nothing about HOW batches are produced is restated
# here on purpose -- RNGkind("L'Ecuyer-CMRG") and the seed must be set exactly as
# the monolith sets them, or the replicates come out different.
generate_batches(
  null_data      = null_data,
  filtered_data  = filtered_data,
  total_variants = total_variants,
  numvariants    = numvariants,
  seed           = seed,
  numSampleSets  = numSampleSets,
  numVariantSets = numVariantSets,
  qHillNumber    = qHillNumber,
  chunkDir       = chunkDir,
  resume         = resume,
  evaluatePopulation      = evaluatePopulation,
  evaluateSuperPopulation = evaluateSuperPopulation,
  downSamplePop           = FALSE,
  downSampleSuperPop      = FALSE,
  maxGB                   = maxGB
)

# ---- manifest ----------------------------------------------------------------
# One entry per batch file, each carrying the fields hill-batchproc.R needs to
# name its output and the aggregate step needs to order replicates. Reading the
# index (not globbing) keeps this in step with whatever was actually written.
index_files <- list.files(chunkDir, pattern = "^chunk_[0-9]+\\.RData$", full.names = TRUE)
if (length(index_files) == 0) stop("No chunk index files in ", chunkDir)
nums <- as.integer(gsub(".*chunk_([0-9]+)\\.RData$", "\\1", basename(index_files)))
index_files <- index_files[order(nums)]

manifest <- list()
for (f in index_files) {
  e <- new.env(); load(f, envir = e)
  ci <- e$chunk_index
  for (bf in ci$batch_files) {
    manifest[[length(manifest) + 1L]] <- list(
      sample_set_id  = as.integer(ci$sample_set_id),
      batch_file     = bf$file,
      variant_ids    = as.integer(bf$variant_ids),
      variant_sets   = length(bf$variant_ids)
    )
  }
}
message("  manifest     : ", length(manifest), " batch file(s)")
saveRDS(manifest, outfile)
message("  wrote ", outfile)
message("=== hill-batchgen: done ===")
