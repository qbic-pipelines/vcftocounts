#!/usr/bin/env Rscript
# =============================================================================
# hill.R -- Hill number diversity calculation
#
# Scope: the Hill number part only.
#   prepare   count table + metadata -> prepared.rds
#   null      calculateNullHill()      (Monte Carlo null replicates)
#   filtered  calculate_filtered_hill()
#
# Variance decomposition, ANOVA-style pipeline differences, filtered-vs-null
# patterns and plots are NOT here. They live in their own scripts so each can be
# a pipeline step of its own.
#
# The reference implementation (../R/*.R of the R-to-Nextflow repo) is vendored
# verbatim below the marker, so this file is self-contained: it can sit in a
# Nextflow bin/ directory and be invoked by bare name from any task work dir.
# =============================================================================

doc <- "
Hill number diversity calculation.

Usage:
  hill.R --mode <mode> --output <file> <input> [options]
  hill.R (-h | --help)

Modes:
  prepare      Count table and sample metadata -> prepared.rds.
  null         Null distribution Hill numbers (chunked, resumable).
  filtered     Filtered Hill numbers.

Options:
  <input>                 Input: for null/filtered the output of --mode prepare,
                          otherwise an rds (preferred) or csv count table.
  --mode <mode>           prepare | null | filtered.
  --output <file>         Output rds.
  --metadata <file>       Sample metadata, used by --mode prepare. Optional if
                          the count table already carries the metadata columns.
  --filtered <file>       Separate count table for the filtered arm, used by the
                          prepare mode. Omit while one table stands in for both;
                          when given, numvariants becomes the size of the
                          filtered set and the null arm resamples that many
                          variants from the full table.
  --pipelines <list>      Comma-separated count columns. Default: every numeric
                          column that is not sample metadata.
  --q <num>               Hill number order [default: 1].
  --pop                   Evaluate by population.
  --superpop              Evaluate by superpopulation.
  --downsample-pop        Downsample by population (creates replicates).
  --downsample-superpop   Downsample by superpopulation.
  --chunks <dir>          Directory for null-distribution chunk files
                          [default: chunks].
  --sample-sets <n>       Number of sample sets [default: 10].
  --variant-sets <n>      Number of variant sets per sample set [default: 10].
  --cores <n>             Number of cores [default: 1].
  --parallel <mode>       Parallel backend: auto | multicore | multisession.
                          auto = multicore where available. multicore forks, so
                          workers share the count table copy-on-write and only
                          pay for what they modify; multisession gives every
                          worker a private copy of it, which can OOM a small
                          machine. [default: auto]
  --seed <n>              Random seed [default: 1234].
  --resume                Reuse existing chunk files where possible.
  --max-gb <n>            future.globals.maxSize in GB [default: 100].
  --max-variants <n>      Cap numvariants (the number of filtered variants
                          sampled per replicate). 0 = no cap. Only needed while
                          there is no separate filtered arm: with one table
                          standing in for both, numvariants would otherwise be
                          the whole table and each replicate would resample all
                          of it. [default: 0]
  -h --help               Show this help.
"

suppressPackageStartupMessages({
  library(docopt)
  library(data.table)
  library(dplyr)
})

opts <- docopt(doc, args = commandArgs(trailingOnly = TRUE))

mode    <- opts$mode
input   <- opts$input
outfile <- opts$output
filteredPath <- opts$filtered

# NB: docopt turns hyphens in option names into underscores (--sample-sets ->
# opts$sample_sets), not into dots. Reading opts$sample.sets silently gives NULL
# and the [default: ...] in the doc never applies.
qHillNumber             <- as.numeric(opts$q)
numSampleSets           <- as.integer(opts$sample_sets)
numVariantSets          <- as.integer(opts$variant_sets)
numCores                <- as.integer(opts$cores)
parallelBackend         <- opts$parallel
seed                    <- as.integer(opts$seed)
chunkDir                <- opts$chunks
maxGB                   <- as.numeric(opts$max_gb)
evaluatePopulation      <- isTRUE(opts$pop)
evaluateSuperPopulation <- isTRUE(opts$superpop)
downSampleByPop         <- isTRUE(opts$downsample_pop)
downSampleBySuperPop    <- isTRUE(opts$downsample_superpop)
resume                  <- isTRUE(opts$resume)
maxVariants             <- as.integer(opts$max_variants)

# Fail loudly rather than propagating integer(0)/numeric(0) into the reference code.
for (nm in c("numSampleSets", "numVariantSets", "numCores", "seed")) {
  if (length(get(nm)) != 1L || is.na(get(nm))) stop("Could not read ", nm, " from the CLI options.")
}
if (length(maxGB) != 1L || is.na(maxGB)) stop("Could not read maxGB from the CLI options.")
if (!parallelBackend %in% c("auto", "multicore", "multisession")) {
  stop("Unknown --parallel '", parallelBackend,
       "'. Expected one of: auto, multicore, multisession.")
}

# Parallel backend selection, shared by every mode that goes parallel.
# multisession (PSOCK) gives each worker its own copy of the count table: at 4
# workers and a 3.5M-row table that multiplies the resident set past a 6 GB
# container limit and the task is OOM-killed (SIGKILL, exit 137). multicore
# forks, so workers inherit the parent's pages copy-on-write and only pay for
# their own writes. Fork is unavailable on Windows, and some container runtimes
# have been seen hanging on it, hence the explicit --parallel override.
select_plan <- function(ncores, backend) {
  if (identical(backend, "auto")) {
    forkable <- tryCatch(.Platform$OS.type == "unix" && future::supportsMulticore(),
                         error = function(e) FALSE)
    backend <- if (isTRUE(forkable)) "multicore" else "multisession"
  }
  if (identical(backend, "multicore")) {
    future::plan(future::multicore, workers = ncores)
  } else {
    future::plan(future::multisession, workers = ncores)
  }
  message("Setting up parallel processing with ", ncores, " workers using ", backend)
  message("Current plan: ", class(future::plan())[1])
  message("Number of workers: ", future::nbrOfWorkers())
  invisible(backend)
}

valid_modes <- c("prepare", "null", "filtered")
if (!mode %in% valid_modes) {
  stop("Unknown --mode '", mode, "'. Expected one of: ",
       paste(valid_modes, collapse = ", "))
}
if (is.null(input)) stop("<input> is required.")
if (is.null(outfile)) stop("--output is required.")

outdir <- dirname(outfile)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Input reader: rds or csv
# -----------------------------------------------------------------------------
# A prepared.rds holds a list of five objects and a null/filtered result holds a
# list of Hill components, so this only converts csv input and otherwise hands
# the object back untouched. Sub-selecting happens per mode.
read_any <- function(path) {
  if (!file.exists(path)) stop("Input not found: ", path)
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    readRDS(path)
  } else if (ext %in% c("csv", "tsv", "txt")) {
    as.data.frame(data.table::fread(path, showProgress = FALSE))
  } else {
    stop("Unsupported extension '", ext, "' for ", path, " (expected .rds or .csv)")
  }
}

# -----------------------------------------------------------------------------
# Sample metadata
# -----------------------------------------------------------------------------
# Canonical names, i.e. what hill_calculations.R get_metadata_cols() expects.
ANALYSIS_METADATA_COLS <- c("sample_name", "variant_called", "population_code",
                            "super_population_code", "population_description",
                            "super_population_description")
DEMO_COLS <- c("population_code", "super_population_code",
               "population_description", "super_population_description")

# Aliases so this accepts either the analysis-shaped metadata counts2rds.py
# writes or a raw metadata table (1000G sequence index style).
METADATA_ALIASES <- list(
  sample_name                  = c("sample_name", "Run", "SAMPLE_NAME", "sample", "Sample", "ID", "id"),
  population_code              = c("population_code", "population", "POPULATION", "pop"),
  super_population_code        = c("super_population_code", "Super_Population_Code", "superpopulation", "super_population"),
  population_description       = c("population_description", "Population_Description"),
  super_population_description = c("super_population_description", "Super_Population_Description")
)

normalise_metadata <- function(md) {
  for (target in names(METADATA_ALIASES)) {
    if (target %in% names(md)) next
    hit <- intersect(METADATA_ALIASES[[target]], names(md))
    if (length(hit) > 0) {
      md[[target]] <- md[[hit[1]]]
      message("  metadata: using column '", hit[1], "' as '", target, "'")
    } else {
      md[[target]] <- NA_character_
    }
  }
  md <- md[, names(METADATA_ALIASES), drop = FALSE]
  md <- md[!is.na(md$sample_name), , drop = FALSE]
  md <- md[!duplicated(md$sample_name), , drop = FALSE]
  for (cl in names(METADATA_ALIASES)) md[[cl]] <- as.character(md[[cl]])
  if (all(is.na(md$population_code))) {
    warning("No population column found in the metadata.")
  }
  md
}

# Which columns carry counts: explicit --pipelines, else every numeric column.
resolve_pipelines <- function(tbl, pipelines_arg) {
  if (!is.null(pipelines_arg) && nzchar(pipelines_arg)) {
    wanted <- trimws(strsplit(pipelines_arg, ",")[[1]])
    missing <- setdiff(wanted, names(tbl))
    if (length(missing) > 0) {
      stop("--pipelines names not found in input: ", paste(missing, collapse = ", "))
    }
    return(wanted)
  }
  candidates <- setdiff(names(tbl), ANALYSIS_METADATA_COLS)
  numericish <- candidates[vapply(tbl[candidates], is.numeric, logical(1))]
  if (length(numericish) == 0) {
    stop("No numeric count columns found. Pass --pipelines explicitly.")
  }
  numericish
}

# =============================================================================
# BEGIN vendored reference code -- ../R/{utils,hill_calculations}.R
# Copied verbatim from the reference implementation. Do not edit here; fix
# upstream and re-vendor.
# =============================================================================
# ---- begin utils.R -----------------------------------------------------
# utils.R
# Utility functions for diversity analysis
# These functions are used throughout the analysis pipeline

#' Sanitize count matrices for Hill number calculations
#' 
#' Ensure count matrices are numeric, non-NA, and non-negative before passing
#' to hillR functions which require non-negative community data.
#' 
#' @param mat A matrix or data frame of count data
#' @return A sanitized numeric matrix with NA replaced by 0 and negatives clamped to 0
sanitize_counts <- function(mat) {
  if(is.null(mat)) return(mat)
  m <- as.matrix(mat)
  # Coerce to numeric where possible
  storage.mode(m) <- "numeric"
  # Replace NA with zero
  m[(is.na(m)) | (m < 0)] <- 0
  return(m)
}

#' Create indices to downsample samples to match length of smallest group
#' 
#' Can downsample by population or superpopulation
#' 
#' @param demoInfo Data frame with demographic information including population_code and super_population_code
#' @param popOrSuperpop Either "population" or "Superpop"
#' @return List with pop_groups (indices per group) and min_size (minimum group size)
downsamplesamples<-function(demoInfo, popOrSuperpop="Superpop"){
  # Normalize the input string to handle both "superpop" and "superpopulation"
  level_lower <- tolower(popOrSuperpop)
  is_superpop <- grepl("^superpop", level_lower)
  
  if(is_superpop){
    # Group by superpopulation
    pop_groups <- split(seq_len(nrow(demoInfo)), demoInfo$super_population_code)
  } else {
    # Group by population
    pop_groups <- split(seq_len(nrow(demoInfo)), demoInfo$population_code)
  }
  # Find minimum group size
  min_size <- min(sapply(pop_groups, length))
  message(paste0("Minimum ", popOrSuperpop, " size: ", min_size))
  
  return(list(pop_groups=pop_groups, min_size=min_size))
}

#' Average Hill number results across multiple replicates
#' 
#' Handles dataframes from hillR by averaging numeric columns
#' 
#' @param hill_list List of Hill number results (typically data frames)
#' @return Averaged result with same structure as input elements
average_hill_results <- function(hill_list) {
  if(length(hill_list) == 0) return(NULL)
  # Get the structure from the first element
  first_result <- hill_list[[1]]
  if(is.data.frame(first_result)) {
    # Average each numeric column across all replicates
    avg_result <- first_result
    for(col in colnames(first_result)) {
      if(is.numeric(first_result[[col]])) {
        col_values <- sapply(hill_list, function(x) x[[col]])
        # Check if col_values is a vector (1D) or matrix (2D)
        if(is.null(dim(col_values))) {
          # It's a vector (single row dataframe case), use mean
          avg_result[[col]] <- mean(col_values, na.rm = TRUE)
        } else {
          # It's a matrix (multiple rows), use rowMeans
          avg_result[[col]] <- rowMeans(col_values, na.rm = TRUE)
        }
      }
    }
    return(avg_result)
  } else {
    # If not a dataframe, just average the values
    return(mean(unlist(hill_list), na.rm = TRUE))
  }
}
# ---- end utils.R ----

# ---- begin hill_calculations.R -----------------------------------------
# hill_calculations.R
suppressPackageStartupMessages(library(hillR))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(future.apply))
suppressPackageStartupMessages(library(progressr))
suppressPackageStartupMessages(library(data.table))

# =============================================================================
# Helper Functions
# =============================================================================

#' Create code-to-description mapping from data
#' @param data Data frame containing code and description columns
#' @param code_col Name of code column
#' @param desc_col Name of description column
#' @return Named vector mapping codes to descriptions
create_code_mapping <- function(data, code_col, desc_col) {
  mapping <- setNames(data[[desc_col]], data[[code_col]])
  mapping[!duplicated(names(mapping))]
}

# DEVIATION (added helper): rename a character vector by lookup, keeping the
# original name wherever the lookup yields NA. The reference does this inline in
# rename_list_components() but not where population/superpopulation list elements
# are named, so when population_description is all-NA (as it is whenever no
# population lookup is supplied) the elements get named NA and later lookups like
# data_split[[group_name]] return NULL, which surfaces as
# "no applicable method for 'group_by' applied to an object of class NULL".
rename_codes_safely <- function(values, code_to_desc) {
  out <- code_to_desc[values]
  out[is.na(out)] <- values[is.na(out)]
  unname(out)
}

#' Get metadata columns to exclude from analysis
#' @return Character vector of column names to exclude
get_metadata_cols <- function() {
  c("sample_name", "variant_called", "population_code", "super_population_code", 
    "population_description", "super_population_description")
}

#' Process Hill numbers for a group with replicates
#' @param group_data Data for a specific group (population or superpopulation)
#' @param group_name Name of the group
#' @param chunk_sample_sets List of sample sets for replicates
#' @param qHillNumber Hill number order
#' @param pairwise Whether to calculate pairwise
#' @return Data frame with Hill numbers for all replicates
process_group_with_replicates <- function(group_data, group_name, chunk_sample_sets, qHillNumber, pairwise = FALSE) {
  # Get sample sets for this group across all replicates
  group_sample_sets <- lapply(chunk_sample_sets, function(samples) {
    sort(intersect(samples, group_data$sample_name))
  })
  
  # Check if all replicates are identical
  all_identical <- length(group_sample_sets) > 0 && 
                   all(sapply(group_sample_sets, function(x) identical(x, group_sample_sets[[1]])))
  
  if(all_identical) {
    message(paste0("  ", group_name, ": All replicates identical (point estimate)"))
  } else {
    message(paste0("  ", group_name, ": ", length(group_sample_sets), " replicates"))
  }

  # Calculate Hill numbers for each replicate
  hill_results_list <- lapply(group_sample_sets, function(sample_set) {
    # DEVIATION (bug fix): subset_data was referenced below but never defined.
    # Matches how the surviving replicates are built elsewhere in this file
    # (e.g. group_data[group_data$sample_name %in% sampled_samples, ]).
    subset_data <- group_data[group_data$sample_name %in% sample_set, , drop = FALSE]
    if(nrow(subset_data) == 0) {
      if(!pairwise) warning(paste0("No samples found for ", group_name))
      return(NULL)
    }
    keep_cols <- setdiff(colnames(subset_data), get_metadata_cols())
    xmat <- as.matrix(subset_data[, keep_cols, drop=FALSE])
    xmat <- sanitize_counts(xmat)
    safe_hill(xmat, q=qHillNumber, pairwise=pairwise, stop_on_error=FALSE, 
              progress=if(pairwise) FALSE else NULL, return_on_error=NULL)
  })
  
  # Combine into data frame
  valid_results <- hill_results_list[!sapply(hill_results_list, is.null)]
  if(length(valid_results) > 0) {
    do.call(rbind, lapply(valid_results, as.data.frame))
  } else {
    NULL
  }
}

#' Calculate Hill numbers by group with or without replicates
#' @param filtered_data Filtered variant data
#' @param group_col Column name for grouping (e.g., "population_code")
#' @param code_to_desc Named vector mapping codes to descriptions
#' @param chunk_sample_sets Optional list of sample sets for replicates
#' @param qHillNumber Hill number order
#' @param pairwise Whether to calculate pairwise
#' @param use_replicates Whether to use downsampling replicates
#' @return Named list of Hill numbers by group
calculate_hill_by_group <- function(filtered_data, group_col, code_to_desc, 
                                    chunk_sample_sets = NULL, qHillNumber, 
                                    pairwise = FALSE, use_replicates = FALSE) {
  
  if(use_replicates && !is.null(chunk_sample_sets) && length(chunk_sample_sets) > 0) {
    # Process with replicates
    group_names <- unique(filtered_data[[group_col]])
    
    result <- lapply(group_names, function(group_name) {
      group_data <- filtered_data[filtered_data[[group_col]] == group_name, ]
      process_group_with_replicates(group_data, group_name, chunk_sample_sets, qHillNumber, pairwise)
    })
    names(result) <- rename_codes_safely(group_names, code_to_desc)
    
  } else {
    # Process without replicates
    data_split <- split(filtered_data, filtered_data[[group_col]])
    names(data_split) <- rename_codes_safely(names(data_split), code_to_desc)
    
    result <- future_lapply(data_split, function(x) {
      keep_cols <- setdiff(colnames(x), get_metadata_cols())
      x <- as.matrix(x[, keep_cols, drop=FALSE])
      x <- sanitize_counts(x)
      context_msg <- if(pairwise) "hill_taxa_parti_pairwise failed" else "hill_taxa_parti failed"
      safe_hill(x, q=qHillNumber, pairwise=pairwise, stop_on_error=TRUE, 
                context_msg=context_msg, progress=if(pairwise) FALSE else NULL)
    }, future.globals=FALSE, future.packages=c("hillR"), future.stdout=FALSE, future.seed=TRUE)
  }
  
  return(result)
}

#' Rename list components from codes to descriptions
#' @param result_list List with Hill number results
#' @param component_names Names of components to rename
#' @param code_to_desc Mapping from codes to descriptions
#' @return Modified result list with renamed components
rename_list_components <- function(result_list, component_names, code_to_desc) {
  for(comp_name in component_names) {
    if(!is.null(result_list[[comp_name]]) && is.list(result_list[[comp_name]])) {
      orig_names <- names(result_list[[comp_name]])
      if(!is.null(orig_names) && length(orig_names) > 0) {
        new_names <- code_to_desc[orig_names]
        new_names[is.na(new_names)] <- orig_names[is.na(new_names)]
        names(result_list[[comp_name]]) <- new_names
        message(paste0("  Renamed ", comp_name, " from codes to descriptions"))
      }
    }
  }
  return(result_list)
}

#' Apply a function to a list or vector with optional parallelism and error handling
#' 
#' @param X List or vector to iterate over
#' @param FUN Function to apply
#' @param parallel Use future_lapply if TRUE, else lapply
#' @param future_pkgs Character vector of packages to load in parallel
#' @param ... Additional arguments passed to FUN
#' @return List of results (errors returned as NA)
safe_apply <- function(X, FUN, parallel=FALSE, future_pkgs=NULL, ...) {
  safe_FUN <- function(...) {
    tryCatch(FUN(...), error=function(e) { warning(e$message); return(NA) })
  }
  if (parallel) {
    if (!requireNamespace("future.apply", quietly = TRUE)) stop("Please install the 'future.apply' package.")
    future.apply::future_lapply(X, safe_FUN, future.globals=FALSE, future.packages=future_pkgs, ...)
  } else {
    lapply(X, safe_FUN, ...)
  }
}

#' Compute by-variant Hill numbers for a grouping column (population or superpopulation)
#' @param filtered_data data.frame of filtered variants
#' @param group_col string name of grouping column ('population_code' or 'super_population_code')
#' @param group_desc_map named vector mapping codes -> descriptions (optional)
#' @param chunk_sample_sets list of sample sets for replicates (optional)
#' 
#' @param qHillNumber numeric
#' @param pairwise boolean
#' @return named list by group of hill results (data.frames or NA)
compute_filtered_by_variant <- function(filtered_data, group_col, group_desc_map=NULL, chunk_sample_sets=NULL, qHillNumber, pairwise=FALSE) {
  # If replicates (chunk sample sets) provided, average across replicates per group
  if(!is.null(chunk_sample_sets) && length(chunk_sample_sets) > 0) {
    filtered_split <- split(filtered_data, filtered_data[[group_col]])
    if(!is.null(group_desc_map)) names(filtered_split) <- rename_codes_safely(names(filtered_split), group_desc_map)
    res <- future_lapply(names(filtered_split), function(group_name) {
      group_data <- filtered_split[[group_name]]
      hill_results_list <- lapply(seq_along(chunk_sample_sets), function(rep) {
        subset_data <- group_data[group_data$sample_name %in% chunk_sample_sets[[rep]], ]
        subset_by_var <- subset_data %>% group_by(variant_called) %>% summarise(across(where(is.numeric), sum), .groups='drop')
        subset_by_var <- data.frame(subset_by_var)
        keep_cols <- setdiff(colnames(subset_by_var), c("sample_name", "variant_called", "population_code", "super_population_code", "population_description", "super_population_description"))
        xmat <- as.matrix(subset_by_var[, keep_cols, drop=FALSE])
        if(is.null(dim(xmat)) || nrow(xmat) == 0 || ncol(xmat) == 0) return(NA)
        xmat <- sanitize_counts(xmat)
        if(pairwise) {
          safe_hill(xmat, q=qHillNumber, pairwise=TRUE, stop_on_error=FALSE, progress=FALSE, return_on_error=NA)
        } else {
          safe_hill(xmat, q=qHillNumber, pairwise=FALSE, stop_on_error=FALSE, return_on_error=NA)
        }
      })
      average_hill_results(hill_results_list)
    }, future.globals=FALSE, future.packages=c("hillR", "dplyr"), future.stdout=FALSE, future.seed=TRUE)
    names(res) <- names(filtered_split)
    return(res)
  } else {
    # No replicates: aggregate once per group
    other_col <- setdiff(c("population_code", "super_population_code"), group_col)
    filtered_split_by_group <- split(filtered_data[, setdiff(colnames(filtered_data), c("sample_name", other_col))], filtered_data[[group_col]])
    filtered_split_by_group <- lapply(filtered_split_by_group, function(x) {
      x <- x %>% group_by(variant_called) %>% summarise(across(where(is.numeric), sum), .groups='drop')
      x <- data.frame(x)
      rownames(x) <- x$variant_called
      x
    })
    res <- lapply(filtered_split_by_group, function(x) {
      keep_cols <- setdiff(colnames(x), c("sample_name", "variant_called", "population_code", "super_population_code", "population_description", "super_population_description"))
      xmat <- as.matrix(x[, keep_cols, drop=FALSE])
      if(is.null(dim(xmat)) || nrow(xmat) == 0 || ncol(xmat) == 0) return(NA)
      xmat <- sanitize_counts(xmat)
      if(pairwise) {
        safe_hill(xmat, q=qHillNumber, pairwise=TRUE, stop_on_error=TRUE, context_msg=paste0("hill_taxa_parti_pairwise failed for ", group_col), progress=FALSE)
      } else {
        safe_hill(xmat, q=qHillNumber, pairwise=FALSE, stop_on_error=TRUE, context_msg=paste0("hill_taxa_parti failed for ", group_col))
      }
    })
    return(res)
  }
}

#' Safe wrapper for hillR calls (handles sanitization + errors)
#' @param x matrix-like object or data.frame
#' @param q Hill order
#' @param pairwise If TRUE, call hill_taxa_parti_pairwise
#' @param stop_on_error If TRUE, rethrow error; else return `return_on_error` and warn
#' @param context_msg Optional context to prepend to warnings/errors
#' @param progress Passed to pairwise function `.progress`
#' @param return_on_error Value to return on error (default NULL)
#' @return Result of hillR call or `return_on_error` on failure
safe_hill <- function(x, q, pairwise=FALSE, stop_on_error=FALSE, context_msg=NULL, progress=FALSE, return_on_error=NULL) {
  if(is.null(x)) return(return_on_error)
  xmat <- tryCatch(sanitize_counts(as.matrix(x)), error=function(e){ 
    if(stop_on_error) stop(paste0(if(!is.null(context_msg)) paste0(context_msg, ": ") else "", "failed to sanitize input: ", e$message))
    warning(paste0(if(!is.null(context_msg)) paste0(context_msg, ": ") else "", "failed to sanitize input: ", e$message))
    return(return_on_error)
  })
  if(pairwise) {
    tryCatch({
      result <- hill_taxa_parti_pairwise(t(xmat), q = q, .progress = progress)
      # Fix: hillR returns beta=0 for identical samples (diagonal), but in multiplicative
      # partitioning, beta = gamma/alpha, so identical samples should have beta = 1
      if(!is.null(result) && is.data.frame(result)) {
        # Identify beta diversity columns (various possible names from hillR)
        beta_cols <- grep("beta|local_similarity|turnover|TD_beta", colnames(result), 
                         ignore.case = TRUE)
        if(length(beta_cols) > 0) {
          for(col_idx in beta_cols) {
            col_name <- colnames(result)[col_idx]
            # For self-comparisons (site1 == site2), set beta to 1
            if("site1" %in% colnames(result) && "site2" %in% colnames(result)) {
              self_comparison <- result$site1 == result$site2
              if(any(self_comparison, na.rm = TRUE)) {
                result[self_comparison, col_idx] <- 1
                if(sum(self_comparison, na.rm = TRUE) > 0) {
                  message(sprintf("Corrected %d self-comparison beta values from 0 to 1 in column '%s'", 
                                sum(self_comparison, na.rm = TRUE), col_name))
                }
              }
            }
          }
        }
      }
      return(result)
    }, error = function(e) {
      full_msg <- paste0(if(!is.null(context_msg)) paste0(context_msg, ": ") else "", e$message)
      if(stop_on_error) stop(full_msg) else { warning(full_msg); return(return_on_error) }
    })
  } else {
    tryCatch({
      result <- hill_taxa_parti(t(xmat), q = q)
      return(result)
    }, error = function(e) {
      full_msg <- paste0(if(!is.null(context_msg)) paste0(context_msg, ": ") else "", e$message)
      if(stop_on_error) stop(full_msg) else { warning(full_msg); return(return_on_error) }
    })
  }
}

#' 
#' Computes Hill numbers (diversity metrics) for total data and splits by population/superpopulation
#' 
#' @param sampled_variants_mat Matrix of sampled variants (rows = observations, cols = variant sites)
#' @param sampled_variants_mat_pop_by_var List of population-specific variant matrices
#' @param sampled_variants_mat_superpop_by_var List of superpopulation-specific variant matrices

#' @param rows_keep Row indices to keep in analysis
#' @param variants_keep Variant IDs to keep (not currently used but kept for compatibility)
#' @param local_idx_pop List of population-specific row indices
#' @param local_idx_superpop List of superpopulation-specific row indices
#' @param qHillNumber Hill number order (q parameter, typically 0, 1, or 2)
#' @return List with Hill number results for different data splits
hillnum <- function(sampled_variants_mat, sampled_variants_mat_pop_by_var=NULL, 
                   sampled_variants_mat_superpop_by_var=NULL, rows_keep, variants_keep, 
                   local_idx_pop, local_idx_superpop, qHillNumber) {
  sanitized_sampled_variants_mat <- sanitize_counts(sampled_variants_mat)
  samplenulldatahill <- safe_hill(sanitized_sampled_variants_mat, q=qHillNumber, pairwise=FALSE, stop_on_error=FALSE, return_on_error=NA, context_msg="hill_taxa_parti failed on TotalData")
  samplenulldatahillpw <- safe_hill(sanitized_sampled_variants_mat, q=qHillNumber, pairwise=TRUE, stop_on_error=FALSE, return_on_error=NA, context_msg="hill_taxa_parti_pairwise failed on TotalData", progress=FALSE)
  # Split by population (only if local_idx_pop is provided)
  samplenullsplitbypop <- NULL
  if(!is.null(local_idx_pop) && length(local_idx_pop) > 0) {
    samplenullsplitbypop <- safe_apply(local_idx_pop, function(x) {
      if(length(x) == 0) {
        return(list(pop=NA, pw=NA))
      }
      submat <- sampled_variants_mat[which(rows_keep %in% x), , drop=FALSE]
      if(nrow(submat) == 0 || ncol(submat) == 0) {
        return(list(pop=NA, pw=NA))
      }
      submat <- sanitize_counts(as.matrix(submat))
      returnlist <- list(
        pop = safe_hill(submat, q=qHillNumber, pairwise=FALSE, stop_on_error=FALSE, return_on_error=NA),
        pw = safe_hill(submat, q=qHillNumber, pairwise=TRUE, stop_on_error=FALSE, progress=FALSE, return_on_error=NA)
      )
      return(returnlist)
    }, parallel=FALSE)
    samplenullsplitbypop <- list_transpose(samplenullsplitbypop)
  }
  
  # Split by superpopulation (only if local_idx_superpop is provided)
  samplenullsplitbysuperpop <- NULL
  if(!is.null(local_idx_superpop) && length(local_idx_superpop) > 0) {
    samplenullsplitbysuperpop <- safe_apply(local_idx_superpop, function(x) {
      if(length(x) == 0) {
        return(list(pop=NA, pw=NA))
      }
      submat <- sampled_variants_mat[which(rows_keep %in% x), , drop=FALSE]
      if(nrow(submat) == 0 || ncol(submat) == 0) {
        return(list(pop=NA, pw=NA))
      }
      submat <- sanitize_counts(as.matrix(submat))
      returnlist <- list(
        pop = safe_hill(submat, q=qHillNumber, pairwise=FALSE, stop_on_error=FALSE, return_on_error=NA),
        pw = safe_hill(submat, q=qHillNumber, pairwise=TRUE, stop_on_error=FALSE, progress=FALSE, return_on_error=NA)
      )
      return(returnlist)
    }, parallel=FALSE)
    samplenullsplitbysuperpop <- list_transpose(samplenullsplitbysuperpop)
  }
  
  # Population by variant
  samplenulldatasplitbypopbyvariant <- NULL
  if(!is.null(sampled_variants_mat_pop_by_var)){
    samplenulldatasplitbypopbyvariant <- lapply(sampled_variants_mat_pop_by_var, function(x) {
      keep_cols <- setdiff(colnames(x), c("sample_name", "variant_called", "population_code", "super_population_code", "population_description", "super_population_description"))
      xmat <- as.matrix(x[, keep_cols, drop=FALSE])
      xmat <- sanitize_counts(xmat)
      returnlist <- list(
        pop = safe_hill(xmat, q=qHillNumber, pairwise=FALSE, stop_on_error=TRUE, context_msg="hill_taxa_parti failed for population-by-variant matrix"),
        pw = safe_hill(xmat, q=qHillNumber, pairwise=TRUE, stop_on_error=TRUE, context_msg="hill_taxa_parti_pairwise failed for population-by-variant matrix", progress=FALSE)
      )
      return(returnlist)
    })
    samplenulldatasplitbypopbyvariant <- list_transpose(samplenulldatasplitbypopbyvariant)
  }
  
  # Superpopulation by variant
  samplenulldatasplitbysuperpopbyvariant <- NULL
  if(!is.null(sampled_variants_mat_superpop_by_var)){
    samplenulldatasplitbysuperpopbyvariant <- lapply(sampled_variants_mat_superpop_by_var, function(x) {
      keep_cols <- setdiff(colnames(x), c("sample_name", "variant_called", "population_code", "super_population_code", "population_description", "super_population_description"))
      xmat <- as.matrix(x[, keep_cols, drop=FALSE])
      xmat <- sanitize_counts(xmat)
      returnlist <- list(
        pop = safe_hill(xmat, q=qHillNumber, pairwise=FALSE, stop_on_error=TRUE, context_msg="hill_taxa_parti failed for superpopulation-by-variant matrix"),
        pw = safe_hill(xmat, q=qHillNumber, pairwise=TRUE, stop_on_error=TRUE, context_msg="hill_taxa_parti_pairwise failed for superpopulation-by-variant matrix", progress=FALSE)
      )
      return(returnlist)
    })
    samplenulldatasplitbysuperpopbyvariant <- list_transpose(samplenulldatasplitbysuperpopbyvariant)
  }
  
  hilldiversitylist <- list(
    TotalData = samplenulldatahill,
    SplitByPopulation = if(!is.null(samplenullsplitbypop)) samplenullsplitbypop$pop else NULL,
    SplitBySuperPopulation = if(!is.null(samplenullsplitbysuperpop)) samplenullsplitbysuperpop$pop else NULL,
    SplitByPopulationByVariant = if(!is.null(samplenulldatasplitbypopbyvariant)) samplenulldatasplitbypopbyvariant$pop else NULL,
    SplitBySuperPopulationByVariant = if(!is.null(samplenulldatasplitbysuperpopbyvariant)) samplenulldatasplitbysuperpopbyvariant$pop else NULL,
    PairwiseTotalData = samplenulldatahillpw,
    PairwiseSplitByPopulation = if(!is.null(samplenullsplitbypop)) samplenullsplitbypop$pw else NULL,
    PairwiseSplitBySuperPopulation = if(!is.null(samplenullsplitbysuperpop)) samplenullsplitbysuperpop$pw else NULL,
    PairwiseSplitByPopulationByVariant = if(!is.null(samplenulldatasplitbypopbyvariant)) samplenulldatasplitbypopbyvariant$pw else NULL,
    PairwiseSplitBySuperPopulationByVariant = if(!is.null(samplenulldatasplitbysuperpopbyvariant)) samplenulldatasplitbysuperpopbyvariant$pw else NULL
  )
  
  return(hilldiversitylist)
}

#' Downsample filtered data by population/superpopulation and calculate averaged Hill numbers
#' 
#' If prespecified_samples is provided, uses those exact sample sets from null distribution chunks
#' If return_all_replicates=TRUE, returns list with individual replicates instead of averaging
#' 
#' @param data_split List of data frames split by population or superpopulation
#' @param grouping_col Name of grouping column for logging
#' @param num_reps Number of downsampling replicates
#' @param qHillNumber Hill number order
#' @param seed Random seed for reproducibility
#' @param pairwise If TRUE, use pairwise Hill number calculation
#' @param prespecified_samples Optional list of pre-specified sample sets to use
#' @param parallel If TRUE, use parallel processing
#' @param return_all_replicates If TRUE, return all replicates instead of averaging
#' @return Either averaged Hill numbers or list of all replicates (depending on return_all_replicates)
calculate_downsampled_hill <- function(data_split, grouping_col, num_reps, qHillNumber, seed, 
                                      pairwise=FALSE, prespecified_samples=NULL, parallel=FALSE, 
                                      return_all_replicates=FALSE) {
  # Get sample counts per group
  sample_counts <- lapply(data_split, function(x) {
    length(unique(x$sample_name))
  })
  min_samples <- min(unlist(sample_counts))
  
  if(is.null(prespecified_samples)) {
    message(paste0("Downsampling ", grouping_col, " to ", min_samples, " samples per group over ", num_reps, " replicates"))
    set.seed(seed)
  } else {
    message(paste0("Using prespecified samples from null distribution chunks (", length(prespecified_samples), " replicates)"))
    num_reps <- length(prespecified_samples)
  }
  
  # For each group, perform downsampling and calculate Hill numbers
  result_by_group <- lapply(names(data_split), function(group_name) {
    group_data <- data_split[[group_name]]
    unique_samples <- unique(group_data$sample_name)
    
    # If this group has exactly min_samples or fewer and no prespecified samples, no need to downsample
    if(is.null(prespecified_samples) && length(unique_samples) <= min_samples) {
      keep_cols <- setdiff(colnames(group_data), c("sample_name", "variant_called", "population_code", "super_population_code", "population_description", "super_population_description"))
      mat <- as.matrix(group_data[, keep_cols, drop=FALSE])
      mat <- sanitize_counts(mat)
      if(pairwise) {
        return(safe_hill(mat, q=qHillNumber, pairwise=TRUE, stop_on_error=FALSE, context_msg=paste0("Hill pairwise failed for ", group_name), return_on_error=NA, progress=FALSE))
      } else {
        return(safe_hill(mat, q=qHillNumber, pairwise=FALSE, stop_on_error=FALSE, context_msg=paste0("Hill failed for ", group_name), return_on_error=NA))
      }
    }
    
    # Perform multiple downsampling replicates
    if(parallel) {
      hill_results <- future_lapply(1:num_reps, function(rep) {
        # Use prespecified samples if provided, otherwise sample randomly
        if(!is.null(prespecified_samples)) {
          sampled_samples <- prespecified_samples[[rep]]
        } else {
          sampled_samples <- sample(unique_samples, min_samples, replace=FALSE)
        }
        
        subset_data <- group_data[group_data$sample_name %in% sampled_samples, ]
        
        keep_cols <- setdiff(colnames(subset_data), c("sample_name", "variant_called", "population_code", "super_population_code", "population_description", "super_population_description"))
        mat <- as.matrix(subset_data[, keep_cols, drop=FALSE])
        mat <- sanitize_counts(mat)
        
        if(pairwise) {
          return(safe_hill(mat, q=qHillNumber, pairwise=TRUE, stop_on_error=FALSE, context_msg=paste0("Hill pairwise failed for ", group_name, " rep ", rep), return_on_error=NA, progress=FALSE))
        } else {
          return(safe_hill(mat, q=qHillNumber, pairwise=FALSE, stop_on_error=FALSE, context_msg=paste0("Hill failed for ", group_name, " rep ", rep), return_on_error=NA))
        }
      }, future.seed=TRUE, future.globals=FALSE, future.packages=c("hillR"), future.stdout=FALSE)
    } else {
      hill_results <- lapply(1:num_reps, function(rep) {
        # Use prespecified samples if provided, otherwise sample randomly
        if(!is.null(prespecified_samples)) {
          sampled_samples <- prespecified_samples[[rep]]
        } else {
          sampled_samples <- sample(unique_samples, min_samples, replace=FALSE)
        }
        
        subset_data <- group_data[group_data$sample_name %in% sampled_samples, ]
        
        keep_cols <- setdiff(colnames(subset_data), c("sample_name", "variant_called", "population_code", "super_population_code", "population_description", "super_population_description"))
        mat <- as.matrix(subset_data[, keep_cols, drop=FALSE])
        mat <- sanitize_counts(mat)
        
        if(pairwise) {
          return(safe_hill(mat, q=qHillNumber, pairwise=TRUE, stop_on_error=FALSE, context_msg=paste0("Hill pairwise failed for ", group_name, " rep ", rep), return_on_error=NA, progress=FALSE))
        } else {
          return(safe_hill(mat, q=qHillNumber, pairwise=FALSE, stop_on_error=FALSE, context_msg=paste0("Hill failed for ", group_name, " rep ", rep), return_on_error=NA))
        }
      })
    }
    
    # Return all replicates or average
    if(return_all_replicates) {
      # Add metadata to each replicate
      hill_results_with_meta <- lapply(seq_along(hill_results), function(i) {
        list(
          value = hill_results[[i]],
          metadata = list(
            group_name = group_name,
            replicate_id = i,
            sample_set_id = i  # Matches null distribution sample sets
          )
        )
      })
      return(hill_results_with_meta)
    } else {
      # Average across replicates
      return(average_hill_results(hill_results))
    }
  })
  
  names(result_by_group) <- names(data_split)
  
  if(return_all_replicates) {
    # Restructure to be replicate-first instead of group-first
    all_replicates <- list()
    for(rep_id in 1:num_reps) {
      rep_data <- list()
      for(group_name in names(result_by_group)) {
        rep_data[[group_name]] <- result_by_group[[group_name]][[rep_id]]$value
      }
      all_replicates[[rep_id]] <- list(
        data = rep_data,
        metadata = list(
          sample_set_id = rep_id,
          replicate_id = rep_id
        )
      )
    }
    return(all_replicates)
  }
  
  return(result_by_group)
}

#' Generate Hill number analysis chunks
#' 
#' Internal helper function for calculateNullHill. Generates chunk files containing
#' pre-sampled variant sets for each sample set. This dramatically reduces memory usage
#' by storing only sampled variants instead of the full dataset.
#' 
#' @param repsvec Vector of sample set IDs to generate
#' @param null_data Full null/unfiltered variant data
#' @param filtered_data Filtered variant data (to count variants per sample set)
#' @param null_mat_num Sanitized numeric matrix of null data
#' @param null_data_split_by_pop_by_variant Population-specific variant matrices (or NULL)
#' @param null_data_split_by_super_pop_by_variant Superpopulation-specific variant matrices (or NULL)
#' @param pop_idx Row indices split by population
#' @param super_idx Row indices split by superpopulation
#' @param sample_groups Sample groupings for downsampling (or NULL)
#' @param min_sample_size Minimum samples per group (or NULL)
#' @param downsample_samples Whether to downsample samples
#' @param evaluatePopulation Whether to create population indices
#' @param evaluateSuperPopulation Whether to create superpopulation indices
#' @param downSamplePop Whether to include population variant matrices
#' @param downSampleSuperPop Whether to include superpopulation variant matrices
#' @param numSampleSets Total number of sample sets
#' @param numVariantSets Number of variant sets per sample
#' @param numvariants Number of variants to sample per replicate
#' @param total_variants Vector of all variant IDs
#' @param seed Random seed
#' @param numCores Number of cores for parallel processing
#' @param chunkDir Directory to save chunks
#' @param maxGB Maximum global object size in GB
#' @return NULL (writes chunk files to disk)
generate_hill_chunks <- function(repsvec, filtered_data, counts = NULL, counts_type = c("matrix","rds_dir"),
                                sample_info = NULL, sample_name_by_row = NULL, variant_ids_by_row = NULL,
                                null_data_split_by_pop_by_variant = NULL,
                                null_data_split_by_super_pop_by_variant = NULL,
                                pop_idx = NULL, super_idx = NULL, sample_groups = NULL, min_sample_size = NULL,
                                downsample_samples = FALSE, evaluatePopulation = FALSE, evaluateSuperPopulation = FALSE,
                                downSamplePop = FALSE, downSampleSuperPop = FALSE,
                                numSampleSets = NULL, numVariantSets = NULL, numvariants = NULL, total_variants = NULL,
                                seed = NULL, numCores = 1, chunkDir = NULL, maxGB = 100, keep_cols = NULL,...) {
  
  message("Generating chunks for replicates...")
  message("Note: Pre-sampling variants during chunk generation dramatically reduces memory")
  
  # NOTE: data.table subsetting was removed because setDT() corrupts
  # shared memory across multicore forks (FutureInterruptError on workers).
  # The data.frame [i, j] subsetting is slower but correct. Optimization
  # should target the hillR calls themselves, not the data reshaping.

  # Calculate data structure details
  total_unique_variants <- length(total_variants)
  # Caller must provide `counts` (matrix or file-backed object) and sample metadata (`sample_info`, `sample_name_by_row`, `variant_ids_by_row`).
  # This avoids serializing large data.frames into worker globals.
  # total_observations is number of rows (observations) in the counts-backed table
  total_observations <- if(!is.null(sample_name_by_row)) length(sample_name_by_row) else if(!is.null(counts)) tryCatch(NROW(counts), error=function(e) stop("Cannot determine total observations")) else stop("Need sample_name_by_row or counts to determine total_observations")
  # Number of unique samples (used for downsampling diagnostics)
  unique_samples <- if(!is.null(sample_info)) nrow(sample_info) else if(!is.null(sample_name_by_row)) length(unique(sample_name_by_row)) else NA_integer_

  # numvariants = raw number of variants to sample (from downsampled filtered data)
  # We sample this exact number (or fewer if not available after downsampling)

  # Average observations per variant (how many samples have each variant on average)
  avg_obs_per_variant <- total_observations / total_unique_variants

  # Calculate expected metrics after downsampling
  if(downsample_samples) {
    # Calculate total samples after downsampling (min_sample_size per group × number of groups)
    num_groups <- length(sample_groups)
    total_samples_after_downsampling <- min_sample_size * num_groups
    sample_ratio <- total_samples_after_downsampling / unique_samples

    # Expected unique variants after sample downsampling (fewer samples → fewer variants observed)
    expected_unique_variants_after_sampling <- round(total_unique_variants * sample_ratio)
    # Sample the raw number from numvariants (or all available if fewer)
    expected_variants_sampled <- min(numvariants, expected_unique_variants_after_sampling)
    # Expected observations
    expected_obs_base <- round(avg_obs_per_variant * expected_variants_sampled * sample_ratio)
  } else {
    # No sample downsampling: sample numvariants from all available
    expected_variants_sampled <- min(numvariants, total_unique_variants)
    expected_obs_base <- round(avg_obs_per_variant * expected_variants_sampled)
    num_groups <- 0
    total_samples_after_downsampling <- unique_samples
  }

  # Print diagnostic information
  message(sprintf("  Unique samples: %d",  unique_samples))
  message(sprintf("  Average observations per variant: %.1f", avg_obs_per_variant))
  if(downsample_samples) {
    message(sprintf("  Downsampling to %d samples per group × %d groups = %d total samples (%.1f%% of %d)", 
                    min_sample_size, num_groups, total_samples_after_downsampling,
                    100 * sample_ratio, unique_samples))
    message("  Variant counts calculated per sample set (see individual sample set logs)")
  } else {
    message(sprintf("  Variants to sample: %d (from filtered data)", numvariants))
  }
  message(sprintf("  Expected observations per variant set: ~%d", expected_obs_base))
  message(sprintf("  Variant sets per sample set: %d", numVariantSets))

  # Calculate storage multiplier based on what's being stored
  # Base matrix: expected_obs_base rows (sampled variants × samples)
  # By-population matrices: expected_variants_sampled rows × num_populations (aggregated by variant)
  # By-superpopulation matrices: expected_variants_sampled rows × num_superpopulations

  base_matrix_size <- expected_obs_base
  additional_matrix_size <- 0

  if(!is.null(downSamplePop) && downSamplePop && !is.null(pop_idx)) {
    num_populations <- length(pop_idx)
    # Each population gets one row per variant (aggregated across samples)
    pop_matrix_size <- expected_variants_sampled * num_populations
    additional_matrix_size <- additional_matrix_size + pop_matrix_size
    message(sprintf("  Including population-by-variant matrices (%d populations × %d variants = %d rows)", 
                    num_populations, expected_variants_sampled, pop_matrix_size))
  }

  if(!is.null(downSampleSuperPop) && downSampleSuperPop && !is.null(super_idx)) {
    num_superpopulations <- length(super_idx)
    superpop_matrix_size <- expected_variants_sampled * num_superpopulations
    additional_matrix_size <- additional_matrix_size + superpop_matrix_size
    message(sprintf("  Including superpopulation-by-variant matrices (%d superpopulations × %d variants = %d rows)", 
                    num_superpopulations, expected_variants_sampled, superpop_matrix_size))
  }

  # Total effective rows across all matrices
  total_effective_rows <- base_matrix_size + additional_matrix_size

  if(base_matrix_size > 0) {
    multiplier <- total_effective_rows / base_matrix_size
  } else {
    multiplier <- 1.0
    if(additional_matrix_size > 0) {
      message(sprintf("  WARNING: base_matrix_size is 0, but additional_matrix_size is %d", 
                     additional_matrix_size))
    }
  }

  # Get number of numeric columns (pipelines) for accurate memory calculation
  numeric_cols <- tryCatch({
    if(!is.null(counts)) NCOL(counts) else sum(sapply(sample_info, is.numeric))
  }, error=function(e) {
    warning("Could not determine numeric column count from counts; defaulting to 0")
    0
  })

  
  # Calculate memory sizes (8 bytes per double + R object overhead ~5x)
  # Chunk: numVariantSets replicates × total_effective_rows per replicate
  # R stores matrices with overhead (row/col names, metadata, list structures, RData headers)
  r_overhead_factor <- 5.0
  chunk_size_mb <- (total_effective_rows * numeric_cols * 8 * numVariantSets * r_overhead_factor) / 1024^2

  
  # Full size: storing all observations once (no replication)
  # Base: total observations, Additional: all variants × all populations
  full_base_rows <- total_observations
  full_additional_rows <- 0
  
  if(!is.null(downSamplePop) && downSamplePop && !is.null(pop_idx)) {
    num_populations <- length(pop_idx)
    full_additional_rows <- full_additional_rows + (total_unique_variants * num_populations)
  }
  
  if(!is.null(downSampleSuperPop) && downSampleSuperPop && !is.null(super_idx)) {
    num_superpopulations <- length(super_idx)
    full_additional_rows <- full_additional_rows + (total_unique_variants * num_superpopulations)
  }
  
  if(downsample_samples) {
    # sample_ratio already calculated correctly earlier (includes num_groups)
    full_base_rows <- round(total_observations * sample_ratio)
  }
  
  full_total_rows <- full_base_rows + full_additional_rows
  full_size_mb <- (full_total_rows * numeric_cols * 8 * r_overhead_factor) / 1024^2

  
  # Set up parallel processing for VARIANT SETS (not sample sets)
  max_size_bytes <- 1024^3 * maxGB
  
  # Use multisession (PSOCK sockets) in containers where multicore (fork) hangs.
  # Detect container by checking for .dockerenv or /proc/1/cgroup containing "docker"/"lxc".
  select_plan(numCores, parallelBackend)
  
  chunk_results <- tryCatch({
    local({
      # Options set here are automatically restored when block exits
      options(future.globals.maxSize = max_size_bytes)
      options(future.wait.timeout = 31536000)
      
      # Process sample sets SEQUENTIALLY
      lapply(repsvec, function(sample_set_id){
      tryCatch({
        if(sample_set_id/numSampleSets*10 == floor(sample_set_id/numSampleSets*10)){
          message(paste0("  Generating sample set ", sample_set_id, " of ", numSampleSets))
        }
        
        # Determine which samples to keep for this sample set
        samples_to_keep <- NULL
        if(downsample_samples){
          set.seed(seed + sample_set_id)
          sampled_samples_idx <- unlist(lapply(sample_groups, function(grp){
            sample(grp, min_sample_size, replace = FALSE)
          }), use.names = FALSE)
          unique_samples <- sample_info
          samples_to_keep <- unique_samples$sample_name[sampled_samples_idx]
        }
        
        # Filter null data by samples
        if(downsample_samples){
          rows_keep <- which(sample_name_by_row %in% samples_to_keep)
        } else {
          rows_keep <- seq_len(total_observations)
        }
        
        # Build indices by population/superpopulation (only if flags are set)
        local_idx_pop <- NULL
        local_idx_superpop <- NULL
        
        # Always create population indices for pairwise calculations (needed for pipeline analysis)
        local_idx_pop <- lapply(pop_idx, function(map) {
          subset <- map[map %in% rows_keep]
          subset <- subset[!is.na(subset)]
          return(subset)
        })
        
        if(!is.null(evaluateSuperPopulation) && evaluateSuperPopulation) {
          local_idx_superpop <- lapply(super_idx, function(map) {
            subset <- map[map %in% rows_keep]
            subset <- subset[!is.na(subset)]
            return(subset)
          })
        }
        
        # Get full data matrices for this sample set (indexed into the provided counts backend)
        # Use data.table subsetting (no copy) instead of data.frame [i,j] which copies 3.5M rows.
        null_mat_sample_set <- counts[rows_keep, keep_cols, drop=FALSE]
        null_mat_sample_set <- sanitize_counts(null_mat_sample_set)
        # Build pop/superpop variant matrices (only if downSample flags are set)
        sampled_variants_mat_pop_by_var <- NULL
        sampled_variants_mat_superpop_by_var <- NULL
        
        if(!is.null(downSamplePop) && downSamplePop && !is.null(null_data_split_by_pop_by_variant)){
          # Use index-based approach to avoid storing full data copies
          if(is.character(null_data_split_by_pop_by_variant) && null_data_split_by_pop_by_variant == "USE_INDICES"){
            message("  Building population matrices from indices (memory-efficient)...")
            sampled_variants_mat_pop_by_var <- lapply(names(pop_idx), function(pop_code) {
              # Get rows for this population that are in our sample set
              pop_rows <- pop_idx[[pop_code]]
              pop_rows_in_sample <- intersect(pop_rows, rows_keep)
              
              if(length(pop_rows_in_sample) == 0) return(NULL)
              
              # Get the data for this population from main dataset
              # Use row indices to subset: sample_name, variant, and numeric columns
              pop_sample_names <- sample_name_by_row[pop_rows_in_sample]
              pop_variant_ids <- variant_ids_by_row[pop_rows_in_sample]
              pop_counts <- counts[pop_rows_in_sample, , drop=FALSE]
              # DEVIATION (bug fix): restrict to keep_cols. 'counts' is the full
              # null_data frame, so the default above also carries the five
              # character metadata columns into the sum below, and data.table
              # aborts with "Type 'character' is not supported by GForce sum".
              # The legacy branch below already excludes them via setdiff().
              pop_counts <- pop_counts[, intersect(keep_cols, colnames(pop_counts)), drop = FALSE]
              # DEVIATION (bug fix): restrict to keep_cols. 'counts' is the full
              # null_data frame, so the default above also carries the five
              # character metadata columns into the sum below, and data.table
              # aborts with "Type 'character' is not supported by GForce sum".
              # The legacy branch below already excludes them via setdiff().
              pop_counts <- pop_counts[, intersect(keep_cols, colnames(pop_counts)), drop = FALSE]
              
              # Create temporary data.table for aggregation
              dt <- data.table(
                variant_called = pop_variant_ids
              )
              # Add count columns
              for(col_idx in seq_len(ncol(pop_counts))) {
                dt[[colnames(pop_counts)[col_idx]]] <- pop_counts[, col_idx]
              }
              
              # Aggregate by variant (sum across samples within this population)
              sample_cols <- colnames(pop_counts)
              x_agg <- dt[, lapply(.SD, sum), by = .(variant_called), .SDcols = sample_cols]
              x_agg <- as.data.frame(x_agg)
              rownames(x_agg) <- x_agg$variant_called
              return(x_agg)
            })
            names(sampled_variants_mat_pop_by_var) <- names(pop_idx)
          } else {
            # Legacy approach: use pre-split data (high memory)
            sampled_variants_mat_pop_by_var <- lapply(null_data_split_by_pop_by_variant, function(x){
              # Downsample to selected samples if needed
              if(downsample_samples){
                x <- x[x$sample_name %in% samples_to_keep, ]
              }
              # NOW aggregate: sum across samples for each variant within this population
              if(nrow(x) > 0) {
                # Use data.table for efficient aggregation
                dt <- as.data.table(x)
                sample_cols <- setdiff(colnames(dt), c("sample_name", "variant_called", "population_code", "super_population_code", "population_description", "super_population_description"))
                x_agg <- dt[, lapply(.SD, sum), by = .(variant_called), .SDcols = sample_cols]
                x_agg <- as.data.frame(x_agg)
                rownames(x_agg) <- x_agg$variant_called
                return(x_agg)
              } else {
                return(NULL)
              }
            })
          }
        }
        
        if(!is.null(downSampleSuperPop) && downSampleSuperPop && !is.null(null_data_split_by_super_pop_by_variant)){
          # Use index-based approach to avoid storing full data copies
          if(is.character(null_data_split_by_super_pop_by_variant) && null_data_split_by_super_pop_by_variant == "USE_INDICES"){
            message("  Building superpopulation matrices from indices (memory-efficient)...")
            sampled_variants_mat_superpop_by_var <- lapply(names(super_idx), function(superpop_code) {
              # Get rows for this superpopulation that are in our sample set
              superpop_rows <- super_idx[[superpop_code]]
              superpop_rows_in_sample <- intersect(superpop_rows, rows_keep)
              
              if(length(superpop_rows_in_sample) == 0) return(NULL)
              
              # Get the data for this superpopulation from main dataset
              superpop_sample_names <- sample_name_by_row[superpop_rows_in_sample]
              superpop_variant_ids <- variant_ids_by_row[superpop_rows_in_sample]
              superpop_counts <- counts[superpop_rows_in_sample, , drop=FALSE]
              # DEVIATION (bug fix): same as the population path above.
              superpop_counts <- superpop_counts[, intersect(keep_cols, colnames(superpop_counts)), drop = FALSE]
              # DEVIATION (bug fix): same as the population path above.
              superpop_counts <- superpop_counts[, intersect(keep_cols, colnames(superpop_counts)), drop = FALSE]
              
              # Create temporary data.table for aggregation
              dt <- data.table(
                variant_called = superpop_variant_ids
              )
              # Add count columns
              for(col_idx in seq_len(ncol(superpop_counts))) {
                dt[[colnames(superpop_counts)[col_idx]]] <- superpop_counts[, col_idx]
              }
              
              # Aggregate by variant (sum across samples within this superpopulation)
              sample_cols <- colnames(superpop_counts)
              x_agg <- dt[, lapply(.SD, sum), by = .(variant_called), .SDcols = sample_cols]
              x_agg <- as.data.frame(x_agg)
              rownames(x_agg) <- x_agg$variant_called
              return(x_agg)
            })
            names(sampled_variants_mat_superpop_by_var) <- names(super_idx)
          } else {
            # Legacy approach: use pre-split data (high memory)
            sampled_variants_mat_superpop_by_var <- lapply(null_data_split_by_super_pop_by_variant, function(x){
              # Downsample to selected samples if needed
              if(downsample_samples){
                x <- x[x$sample_name %in% samples_to_keep, ]
              }
              # NOW aggregate: sum across samples for each variant within this superpopulation
              if(nrow(x) > 0) {
                # Use data.table for efficient aggregation
                dt <- as.data.table(x)
                sample_cols <- setdiff(colnames(dt), c("sample_name", "variant_called", "population_code", "super_population_code", "population_description", "super_population_description"))
                x_agg <- dt[, lapply(.SD, sum), by = .(variant_called), .SDcols = sample_cols]
                x_agg <- as.data.frame(x_agg)
                rownames(x_agg) <- x_agg$variant_called
                return(x_agg)
              } else {
                return(NULL)
              }
            })
          }
        }
        

        # --- BEGIN: DIAGNOSTIC AND ENFORCEMENT FOR FILTERED VARIANT AVAILABILITY ---
        # Instead of using null data, check the filtered variants available for the selected individuals
        filtered_variants_for_samples <- NULL
        if (!is.null(filtered_data) && !is.null(samples_to_keep)) {
          filtered_variants_for_samples <- unique(filtered_data$variant_called[filtered_data$sample_name %in% samples_to_keep])
        } else if (!is.null(filtered_data)) {
          filtered_variants_for_samples <- unique(filtered_data$variant_called)
        } else {
          filtered_variants_for_samples <- character(0)
        }

        n_variants_in_filtered <- length(filtered_variants_for_samples)
        if (n_variants_in_filtered == 0) {
          stop("No filtered variants available for selected individuals in this sample set. Cannot proceed.")
        }

        # Diagnostic: warn if fewer variants than requested
        if (n_variants_in_filtered < numvariants) {
          message(sprintf("Only %d filtered variants available for selected individuals (requested %d). Will sample all available.", n_variants_in_filtered, numvariants))
        }

        # Now, sample up to the number of filtered variants available for these individuals
        n_to_sample <- min(numvariants, n_variants_in_filtered)
        n_to_sample <- max(n_to_sample, 1)  # Always sample at least 1

        # Get variant IDs for this sample set (from precomputed per-row vector)
        variant_ids_in_chunk <- variant_ids_by_row[rows_keep]
        unique_variants_available <- unique(variant_ids_in_chunk)
        
        # Stream variant sets to disk one at a time to minimize RAM usage
        if(sample_set_id/numSampleSets*10 == floor(sample_set_id/numSampleSets*10)){
          message(paste0("    Processing sample set ", sample_set_id, " of ", numSampleSets, " with ", numVariantSets, " variant sets in parallel"))
          message(sprintf("    Filtered data: %d variants | Null data: %d variants available | Sampling: %d variants", 
                         n_variants_in_filtered, length(unique_variants_available), n_to_sample))
        }
        
        # Parallelize variant sets using future_lapply
        # Batch variant sets to reduce file I/O: process in groups of 10
        batch_size <- 10
        num_batches <- ceiling(numVariantSets / batch_size)
        
        # Save large matrices to temporary file for workers to access
        # This prevents copying large objects to each worker
        shared_data_file <- tempfile(pattern = paste0("shared_data_", sample_set_id, "_"), tmpdir = chunkDir, fileext = ".RData")
        shared_data <- list(
          null_mat_sample_set = null_mat_sample_set,
          sampled_variants_mat_pop_by_var = sampled_variants_mat_pop_by_var,
          sampled_variants_mat_superpop_by_var = sampled_variants_mat_superpop_by_var,
          local_idx_pop = local_idx_pop,
          rows_keep = rows_keep,
          variant_ids_in_chunk = variant_ids_in_chunk
        )
        save(shared_data, file = shared_data_file, compress = FALSE)
        rm(shared_data)  # Free memory
        gc(verbose = FALSE)
        
        # Only export small objects to workers
        variant_files <- future_lapply(seq_len(num_batches), function(batch_id) {
          # Load shared data in worker (each worker gets its own copy only when needed)
          e <- new.env()
          load(shared_data_file, envir = e)
          
          # Extract from environment - data is wrapped in shared_data list
          null_mat_sample_set <- e$shared_data$null_mat_sample_set
          sampled_variants_mat_pop_by_var <- e$shared_data$sampled_variants_mat_pop_by_var
          sampled_variants_mat_superpop_by_var <- e$shared_data$sampled_variants_mat_superpop_by_var
          local_idx_pop <- e$shared_data$local_idx_pop
          rows_keep <- e$shared_data$rows_keep
          variant_ids_in_chunk <- e$shared_data$variant_ids_in_chunk
          rm(e)
          gc(verbose = FALSE)
          
          # Determine which variant sets are in this batch
          start_variant <- (batch_id - 1) * batch_size + 1
          end_variant <- min(batch_id * batch_size, numVariantSets)
          variant_ids_in_batch <- start_variant:end_variant
          
          # Process all variant sets in this batch
          batch_data <- lapply(variant_ids_in_batch, function(variant_set_id) {
            set.seed(seed + sample_set_id * 10000 + variant_set_id)
            # Sample variants: always sample without replacement (empirical null)
            sampled_variants <- sample(unique_variants_available, n_to_sample, replace = FALSE)
            # Select all rows for each chosen unique variant
            rows_keep_variant <- which(variant_ids_in_chunk %in% sampled_variants)
            nullmat_variant <- null_mat_sample_set[rows_keep_variant, , drop=FALSE]

            pop_by_var_variant <- NULL
            superpop_by_var_variant <- NULL
            if(!is.null(sampled_variants_mat_pop_by_var)){
              pop_by_var_variant <- lapply(sampled_variants_mat_pop_by_var, function(x){
                x[which(rownames(x) %in% sampled_variants), , drop=FALSE]
              })
            }
            if(!is.null(sampled_variants_mat_superpop_by_var)){
              superpop_by_var_variant <- lapply(sampled_variants_mat_superpop_by_var, function(x){
                x[which(rownames(x) %in% sampled_variants), , drop=FALSE]
              })
            }

            # Map original (global) indices for downstream local index creation.
            # Positions are resolved with one vectorised %in% scan. The previous
            # form scanned original_indices once per matched index
            # (which(original_indices == m) inside lapply), which is O(n*m): at
            # full-variant scale that cost ~11 min per population per variant set
            # instead of the ~0.6 s below. Row order and multiplicity are
            # unchanged, so downstream results are identical.
            original_indices <- rows_keep[rows_keep_variant]
            local_idx_pop_variant <- lapply(local_idx_pop, function(idx){
              if(is.null(idx) || length(idx) == 0) return(integer(0))
              matching <- original_indices %in% idx
              if(!any(matching)) return(integer(0))
              as.integer(which(matching))
            })
            local_idx_superpop_variant <- lapply(local_idx_superpop, function(idx){
              if(is.null(idx) || length(idx) == 0) return(integer(0))
              matching <- original_indices %in% idx
              if(!any(matching)) return(integer(0))
              as.integer(which(matching))
            })
            # Create variant_set object for this variant set
            variant_set <- list(
              sample_set_id = sample_set_id,
              variant_set_id = variant_set_id,
              nullmat = nullmat_variant,
              pop_by_var = pop_by_var_variant,
              superpop_by_var = superpop_by_var_variant,
              local_idx_pop = local_idx_pop_variant,
              local_idx_superpop = local_idx_superpop_variant
            )
            return(variant_set)
          })
          
          # Save entire batch as a single file using qs for speed and compactness
          batch_file <- file.path(chunkDir, sprintf("chunk_downsample_%d_batch_%d.rds", sample_set_id, batch_id))
          tmp_batch_file <- tempfile(pattern = paste0("chunk_downsample_", sample_set_id, "_batch_"), tmpdir = chunkDir, fileext = ".rds")
          saveRDS(batch_data, tmp_batch_file)
          if(file.exists(batch_file)) file.remove(batch_file)
          if(!file.rename(tmp_batch_file, batch_file)) {
            warning(paste0("Failed to atomically move batch file to ", batch_file, "; writing directly."))
            saveRDS(batch_data, batch_file)
          }
          # Return the batch file path and variant set IDs in this batch
          rm(batch_data, null_mat_sample_set, sampled_variants_mat_pop_by_var, 
             sampled_variants_mat_superpop_by_var, local_idx_pop, rows_keep, 
             variant_ids_in_chunk)
          gc(verbose = FALSE)
          return(list(file = batch_file, variant_ids = variant_ids_in_batch))
        }, future.seed=TRUE, future.globals=list(
          shared_data_file = shared_data_file,
          batch_size = batch_size,
          numVariantSets = numVariantSets,
          chunkDir = chunkDir,
          sample_set_id = sample_set_id,
          seed = seed,
          unique_variants_available = unique_variants_available,
          n_to_sample = n_to_sample
        ), future.packages=c("data.table"), future.stdout=FALSE)
        
        # Clean up shared data file
        if(file.exists(shared_data_file)) {
          file.remove(shared_data_file)
        }
        
        # variant_files is now a list of batch metadata from future_lapply
        # Each element has $file (path) and $variant_ids (vector of variant set IDs)
        batch_files <- sapply(variant_files, function(x) x$file)
        
        # Verify batch files were created
        if(length(batch_files) == 0) {
          stop(sprintf("No batch files generated for sample set %d", sample_set_id))
        }
        if(!all(file.exists(batch_files))) {
          missing <- batch_files[!file.exists(batch_files)]
          stop(sprintf("Missing batch files for sample set %d: %s", sample_set_id, paste(head(missing, 5), collapse=", ")))
        }
        message(sprintf("  Sample set %d: Generated %d batch files (%d variant sets)", 
                       sample_set_id, length(batch_files), numVariantSets))
        
        # Save chunk index with batch information
        chunk_index <- list(
          sample_set_id = sample_set_id, 
          batch_files = variant_files,  # Store full list with file paths and variant_ids
          samples_to_keep = if(downsample_samples) samples_to_keep else NULL
        )
        chunk_file_tmp <- tempfile(pattern = paste0("chunk_", sample_set_id, "_"), tmpdir = chunkDir)
        chunk_file_final <- file.path(chunkDir, sprintf("chunk_%d.RData", sample_set_id))
        save(chunk_index, file = chunk_file_tmp, compress = FALSE)
        if(file.exists(chunk_file_final)) file.remove(chunk_file_final)
        if(!file.rename(chunk_file_tmp, chunk_file_final)) {
          warning(paste0("Failed to atomically move chunk index to ", chunk_file_final, "; writing directly."))
          save(chunk_index, file = chunk_file_final, compress = FALSE)
        }
        gc(verbose = FALSE)
        return(list(success=TRUE, sample_set_id=sample_set_id))
      }, error=function(e) {
        message(paste0("ERROR in sample set ", sample_set_id, ": ", e$message))
        return(list(success=FALSE, sample_set_id=sample_set_id, error=e$message))
      })
    # End of lapply over sample sets (sequential)
    })
    })
  }, FutureInterruptError = function(e) {
    message("\n=== INTERRUPTED ===")
    message("Chunk generation was interrupted. Existing chunks have been saved.")
    message(paste0("Resume from this point using --resume flag. Chunks saved in: ", chunkDir))
    return(NULL)
  }, interrupt = function(e) {
    message("\n=== INTERRUPTED ===")
    message("Chunk generation was interrupted. Existing chunks have been saved.")
    message(paste0("Resume from this point using --resume flag. Chunks saved in: ", chunkDir))
    return(NULL)
  }, error = function(e) {
    message(paste0("\nERROR during chunk generation: ", e$message))
    stop(e)
  })
  
  # Check if interrupted
  if(is.null(chunk_results)) {
    message("Exiting due to interruption...")
    quit(save="no", status=1)
  }
  
  # Check for failures
  failures <- Filter(function(x) !x$success, chunk_results)
  if(length(failures) > 0) {
    msg <- paste0("Failed to generate ", length(failures), " chunks: ", 
                  paste(sapply(failures, function(x) x$sample_set_id), collapse=", "))
    stop(msg)
  }
  
  message("Chunk generation complete.")
  gc(verbose = FALSE)
  
  return(invisible(NULL))
}

#' Calculate Hill numbers for filtered data across populations
#' 
#' @param filtered_data Filtered variant data
#' @param qHillNumber Hill number order
#' @param evaluatePopulation Whether to evaluate by population
#' @param evaluateSuperPopulation Whether to evaluate by superpopulation
#' @param downSampleByPop Whether downsampling was done by population
#' @param downSampleBySuperPop Whether downsampling was done by superpopulation
#' @param chunk_sample_sets List of sample sets from chunks (if downsampling)
#' @return List with filtered data Hill numbers
calculate_filtered_hill <- function(filtered_data, qHillNumber, evaluatePopulation, 
                                   evaluateSuperPopulation, downSampleByPop, downSampleBySuperPop, 
                                   chunk_sample_sets=NULL,
                                   numVariantSets=NULL, seed=NULL, variant_sample_fraction=NULL) {
  
  # Create code-to-description mappings
  pop_code_to_desc <- create_code_mapping(filtered_data, "population_code", "population_description")
  superpop_code_to_desc <- create_code_mapping(filtered_data, "super_population_code", "super_population_description")
  keep_cols <- setdiff(colnames(filtered_data), get_metadata_cols())
  
  # Calculate overall Hill numbers
  filtered_data_hill <- NULL
  filtered_data_hill_pw <- NULL
  
  if(!is.null(chunk_sample_sets) && length(chunk_sample_sets) > 0) {
    message("Calculating overall Hill numbers across downsampled replicates...")
    
    sample_sets <- lapply(chunk_sample_sets, function(x) sort(x))
    all_identical <- length(sample_sets) > 0 && all(sapply(sample_sets, function(x) identical(x, sample_sets[[1]])))

    if(all_identical) {
      message("  All sample sets are identical (point estimate for overall data)")
    } else {
      message(paste0("  ", length(sample_sets), " replicates for overall data"))
    }

    # Calculate Hill numbers for each replicate
    hill_results_list <- lapply(sample_sets, function(sample_set) {
      subset_data <- filtered_data[filtered_data$sample_name %in% sample_set, ]
      if(nrow(subset_data) == 0) {
        warning("No samples found in overall data for a replicate")
        return(NULL)
      }
      keep_cols_subset <- setdiff(colnames(subset_data), get_metadata_cols())
      message(paste0("    Subset has ", nrow(subset_data), " rows, ", length(keep_cols_subset), " pipeline columns"))
      xmat <- as.matrix(subset_data[, keep_cols_subset, drop=FALSE])
      message(paste0("    Matrix dims: ", nrow(xmat), " x ", ncol(xmat)))
      xmat <- sanitize_counts(xmat)
      message(paste0("    After sanitize: ", nrow(xmat), " x ", ncol(xmat)))
      result <- safe_hill(xmat, q=qHillNumber, pairwise=FALSE, stop_on_error=FALSE, return_on_error=NULL)
      message(paste0("    safe_hill result is NULL: ", is.null(result)))
      if(!is.null(result)) {
        message(paste0("    safe_hill result class: ", class(result)[1]))
      }
      return(result)
    })
    
    valid_results <- hill_results_list[!sapply(hill_results_list, is.null)]
    if(length(valid_results) > 0) {
      filtered_data_hill <- do.call(rbind, lapply(valid_results, as.data.frame))
      message(paste0("  Created ", nrow(filtered_data_hill), " overall diversity replicate(s)"))
    }
    
    # Calculate pairwise
    hill_pw_results_list <- lapply(sample_sets, function(sample_set) {
      subset_data <- filtered_data[filtered_data$sample_name %in% sample_set, ]
      if(nrow(subset_data) == 0) return(NULL)
      keep_cols_subset <- setdiff(colnames(subset_data), get_metadata_cols())
      xmat <- as.matrix(subset_data[, keep_cols_subset, drop=FALSE])
      xmat <- sanitize_counts(xmat)
      safe_hill(xmat, q=qHillNumber, pairwise=TRUE, stop_on_error=FALSE, progress=FALSE, return_on_error=NULL)
    })
    
    valid_pw_results <- hill_pw_results_list[!sapply(hill_pw_results_list, is.null)]
    if(length(valid_pw_results) > 0) {
      filtered_data_hill_pw <- do.call(rbind, lapply(valid_pw_results, as.data.frame))
    }
    
  } else {
    # No downsampling - calculate once from all samples
    message("Calculating overall Hill numbers from all samples (no downsampling)...")
    filtered_mat <- sanitize_counts(as.matrix(filtered_data[, keep_cols, drop=FALSE]))
    filtered_data_hill <- as.data.frame(safe_hill(filtered_mat, q=qHillNumber, pairwise=FALSE, 
                                                   stop_on_error=TRUE, context_msg="hill_taxa_parti failed on filtered data"))
    filtered_data_hill_pw <- as.data.frame(safe_hill(filtered_mat, q=qHillNumber, pairwise=TRUE, 
                                                      stop_on_error=TRUE, context_msg="hill_taxa_parti_pairwise failed on filtered data", progress=FALSE))
  }
  
  message("Calculated filtered data hill numbers")
  
  # Calculate by population (sample-level) - always calculate pairwise for pipeline analysis
  filtered_data_hill_by_pop <- NULL
  filtered_data_hill_by_pop_pw <- NULL
  
  # Always calculate pairwise by population for pipeline comparison analysis
  message("Calculating pairwise Hill numbers by population (for pipeline analysis)...")
  use_reps <- !is.null(downSampleByPop) && downSampleByPop
  
  if(use_reps) {
    message(paste0("  Creating ", length(chunk_sample_sets), " replicates for pairwise data by population..."))
  } else {
    message("  Processing all samples (no replicates)...")
  }
  
  filtered_data_hill_by_pop_pw <- calculate_hill_by_group(
    filtered_data, "population_code", pop_code_to_desc, chunk_sample_sets, 
    qHillNumber, pairwise=TRUE, use_replicates=use_reps
  )
  
  # Calculate non-pairwise by population only if evaluatePopulation is set
  if(!is.null(evaluatePopulation) && evaluatePopulation) {
    message("Calculating Hill numbers by population (sample-level)...")
    
    filtered_data_hill_by_pop <- calculate_hill_by_group(
      filtered_data, "population_code", pop_code_to_desc, chunk_sample_sets, 
      qHillNumber, pairwise=FALSE, use_replicates=use_reps
    )

  }
  
  # Calculate by superpopulation (sample-level) if requested
  filtered_data_hill_by_super_pop <- NULL
  filtered_data_hill_by_super_pop_pw <- NULL
  
  if(!is.null(evaluateSuperPopulation) && evaluateSuperPopulation) {
    message("Calculating Hill numbers by superpopulation (sample-level)...")
    use_reps <- !is.null(downSampleBySuperPop) && downSampleBySuperPop
    
    if(use_reps) {
      message(paste0("  Creating ", length(chunk_sample_sets), " replicates for filtered data by superpopulation..."))
    } else {
      message("  Processing all samples (no replicates)...")
    }
    
    filtered_data_hill_by_super_pop <- calculate_hill_by_group(
      filtered_data, "super_population_code", superpop_code_to_desc, chunk_sample_sets, 
      qHillNumber, pairwise=FALSE, use_replicates=use_reps
    )
    
    filtered_data_hill_by_super_pop_pw <- calculate_hill_by_group(
      filtered_data, "super_population_code", superpop_code_to_desc, chunk_sample_sets, 
      qHillNumber, pairwise=TRUE, use_replicates=use_reps
    )
  }
  
  # Calculate by-variant Hill numbers if requested
  filtered_data_hill_by_pop_by_variant <- NULL
  filtered_data_hill_by_pop_by_variant_pw <- NULL
  
  if(!is.null(downSampleByPop) && downSampleByPop) {
    message("Calculating Hill numbers by population and variant (variant-level)...")
    chunk_sets <- if(!is.null(chunk_sample_sets) && length(chunk_sample_sets) > 0) chunk_sample_sets else NULL
    
    # Capture warnings during calculation
    filtered_data_hill_by_pop_by_variant <- withCallingHandlers(
      compute_filtered_by_variant(
        filtered_data, "population_code", pop_code_to_desc, chunk_sets, qHillNumber, pairwise=FALSE
      ),
      warning = function(w) {
        message(sprintf("  [WARNING in pop-by-variant]: %s", conditionMessage(w)))
        invokeRestart("muffleWarning")
      }
    )
    
    filtered_data_hill_by_pop_by_variant_pw <- withCallingHandlers(
      compute_filtered_by_variant(
        filtered_data, "population_code", pop_code_to_desc, chunk_sets, qHillNumber, pairwise=TRUE
      ),
      warning = function(w) {
        message(sprintf("  [WARNING in pop-by-variant pairwise]: %s", conditionMessage(w)))
        invokeRestart("muffleWarning")
      }
    )
  }
  
  filtered_data_hill_by_super_pop_by_variant <- NULL
  filtered_data_hill_by_super_pop_by_variant_pw <- NULL
  
  if(!is.null(downSampleBySuperPop) && downSampleBySuperPop) {
    message("Calculating Hill numbers by superpopulation and variant (variant-level)...")
    chunk_sets <- if(!is.null(chunk_sample_sets) && length(chunk_sample_sets) > 0) chunk_sample_sets else NULL
    
    # Capture warnings during calculation
    filtered_data_hill_by_super_pop_by_variant <- withCallingHandlers(
      compute_filtered_by_variant(
        filtered_data, "super_population_code", superpop_code_to_desc, chunk_sets, qHillNumber, pairwise=FALSE
      ),
      warning = function(w) {
        message(sprintf("  [WARNING in superpop-by-variant]: %s", conditionMessage(w)))
        invokeRestart("muffleWarning")
      }
    )
    
    filtered_data_hill_by_super_pop_by_variant_pw <- withCallingHandlers(
      compute_filtered_by_variant(
        filtered_data, "super_population_code", superpop_code_to_desc, chunk_sets, qHillNumber, pairwise=TRUE
      ),
      warning = function(w) {
        message(sprintf("  [WARNING in superpop-by-variant pairwise]: %s", conditionMessage(w)))
        invokeRestart("muffleWarning")
      }
    )
  }

  return(list(
    TotalData = filtered_data_hill,
    PairwiseTotalData = filtered_data_hill_pw,
    SplitByPopulation = filtered_data_hill_by_pop,
    SplitBySuperPopulation = filtered_data_hill_by_super_pop,
    SplitByPopulationByVariant = filtered_data_hill_by_pop_by_variant,
    SplitBySuperPopulationByVariant = filtered_data_hill_by_super_pop_by_variant,
    PairwiseSplitByPopulation = filtered_data_hill_by_pop_pw,
    PairwiseSplitBySuperPopulation = filtered_data_hill_by_super_pop_pw,
    PairwiseSplitByPopulationByVariant = filtered_data_hill_by_pop_by_variant_pw,
    PairwiseSplitBySuperPopulationByVariant = filtered_data_hill_by_super_pop_by_variant_pw
  ))
}

#' Calculate null distribution Hill numbers using chunk-based approach
#' 
#' Generates null distribution by creating sample sets and variant sets,
#' storing intermediate results in chunk files for memory efficiency.
#' Supports resume functionality and parallel processing.
#' 
#' @param null_data Null/unfiltered variant data
#' @param filtered_data Filtered variant data (used to count variants per sample set)
#' @param total_variants Vector of all variant IDs
#' @param numvariants Number of variants to sample per replicate (only used for diagnostics)
#' @param numCores Number of CPU cores for parallel processing
#' @param seed Random seed for reproducibility
#' @param numSampleSets Number of sample sets (outer loop)
#' @param numVariantSets Number of variant sets per sample (inner loop)
#' @param qHillNumber Hill number order
#' @param chunkDir Directory to store chunk files
#' @param show_progress Whether to show progress messages
#' @param evaluateSuperPopulation Whether to evaluate by superpopulation
#' @param evaluatePopulation Whether to evaluate by population
#' @param resume Whether to resume from existing chunks
#' @param downSamplePop Whether to downsample by population
#' @param downSampleSuperPop Whether to downsample by superpopulation
#' @param maxGB Maximum global object size for future in GB
#' @return List with Hill number results and metadata
calculateNullHill <- function(null_data, filtered_data, total_variants, numvariants, numCores, seed, 
                             numSampleSets, numVariantSets, qHillNumber, chunkDir,
                             show_progress=TRUE, evaluateSuperPopulation=FALSE, 
                             evaluatePopulation=FALSE, resume=FALSE, downSamplePop=FALSE, 
                             downSampleSuperPop=FALSE, maxGB=100) {
  
  if(downSamplePop && downSampleSuperPop){
    stop("Cannot downsample by both population and superpopulation simultaneously")
  }
  
  # Diagnostic output
  message("\n=== calculateNullHill Input Summary ===")
  message(sprintf("  Null data dimensions: %d rows × %d cols", nrow(null_data), ncol(null_data)))
  message(sprintf("  Total variants (length): %d", length(total_variants)))
  message(sprintf("  numvariants parameter: %d", numvariants))
  message(sprintf("  Unique variants in null_data: %d", length(unique(null_data$variant_called))))
  message(sprintf("  maxGB parameter (future.globals.maxSize): %d GB", maxGB))
  message("=======================================\n")
  
  # Set global size limit IMMEDIATELY before any parallel operations
  max_size_bytes <- 1024^3 * maxGB
  options(future.globals.maxSize = max_size_bytes)
  # Determine if we're downsampling by samples
  downsample_samples <- downSamplePop || downSampleSuperPop
  sample_groups <- NULL
  min_sample_size <- NULL
  
  if(downsample_samples){
    # Determine which level to downsample at
    if(downSamplePop && downSampleSuperPop) {
      stop("Cannot downsample at both population and superpopulation levels simultaneously. Choose one.")
    }
    
    popOrSuperpop <- if(downSamplePop) "population" else "superpopulation"
    # Get unique samples and their groupings (see the note below: request only the
    # metadata columns that actually exist, so a missing superpopulation column is
    # not fatal)
    meta_cols_for_samples <- intersect(
      c("sample_name", "population_code", "super_population_code",
        "population_description", "super_population_description"),
      colnames(null_data))
    unique_samples <- unique(null_data[, meta_cols_for_samples, drop = FALSE])
    downsample_info <- downsamplesamples(unique_samples, popOrSuperpop)
    sample_groups <- downsample_info$pop_groups
    min_sample_size <- downsample_info$min_size
    message(paste0("Downsampling by ", popOrSuperpop, " to ", min_sample_size, " samples per group"))
  }
  
  # Generate chunks for sample sets only (variant sampling happens during Hill calculation)
  repsvec <- seq_len(numSampleSets)
  variant_ids <- null_data$variant_called
  pop_vec <- null_data$population_code
  super_vec <- null_data$super_population_code
  pop_idx <- split(seq_len(nrow(null_data)), pop_vec)

  # split() on an all-NA (or absent) grouping column yields an empty split and
  # aborts with "group length is 0 but data length > 0", taking the whole run down
  # over a metadata problem. Downstream everything already treats a NULL index as
  # "no such level" (see the !is.null(super_idx) guards), so degrade with a clear
  # message instead. The reference has no such guard; this is a deliberate
  # difference, and it is what makes --superpop safe to enable by default.
  if (is.null(super_vec) || all(is.na(super_vec))) {
    message("WARNING: no superpopulation column in the input, so the ")
    message("         by-superpopulation level will be empty. Populate ")
    message("         super_population_code/super_population_description in the ")
    message("         metadata to enable it.")
    super_idx <- NULL
  } else if (any(is.na(super_vec))) {
    message("WARNING: ", sum(is.na(super_vec)), " row(s) have no ")
    message("         super_population_code and are excluded from the ")
    message("         by-superpopulation level.")
    super_idx <- split(seq_len(nrow(null_data)), super_vec)
    super_idx <- super_idx[!is.na(names(super_idx))]
  } else {
    super_idx <- split(seq_len(nrow(null_data)), super_vec)
  }
  keep_cols <- setdiff(colnames(null_data), c("sample_name", "variant_called", "population_code", "super_population_code", "population_description", "super_population_description"))

  # Instead of storing full data copies per population/superpopulation,
  # we'll pass indices and reconstruct in chunks (saves massive memory)
  null_data_split_by_pop_by_variant <- NULL
  null_data_split_by_super_pop_by_variant <- NULL
  
  if(downSamplePop){
    message("Preparing population variant matrices...")
    message("  NOTE: Using index-based approach to minimize memory usage")
    message("  Data will be subset and aggregated during chunk generation")
    # Calculate memory saved
    n_pops <- length(unique(null_data$population_code))
    full_copy_mb <- (nrow(null_data) * n_pops * ncol(null_data) * 8) / (1024^2)
    message(sprintf("  Memory saved: ~%.0f MB (avoided %dx data duplication)", full_copy_mb, n_pops))
    # Store just the index information - aggregation happens in chunks
    null_data_split_by_pop_by_variant <- "USE_INDICES"  # Flag to indicate index-based approach
  }
  
  if(downSampleSuperPop){
    message("Preparing superpopulation variant matrices...")
    message("  NOTE: Using index-based approach to minimize memory usage")
    message("  Data will be subset and aggregated during chunk generation")
    # Calculate memory saved
    n_superpops <- length(unique(null_data$super_population_code))
    full_copy_mb <- (nrow(null_data) * n_superpops * ncol(null_data) * 8) / (1024^2)
    message(sprintf("  Memory saved: ~%.0f MB (avoided %dx data duplication)", full_copy_mb, n_superpops))
    # Store just the index information - aggregation happens in chunks
    null_data_split_by_super_pop_by_variant <- "USE_INDICES"  # Flag to indicate index-based approach
  }
  
  message("Sanitizing null data count matrix...")
  # Create a numeric counts backend (mori support removed)
  counts<- null_data
  
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  
  if(!dir.exists(chunkDir)){
    dir.create(chunkDir)
  }
  
  chunkofiles <- list.files(chunkDir, pattern="^chunk_[0-9]+\\.RData$", full.names=FALSE)
  skipchunk <- FALSE
  
  if(resume && length(chunkofiles) >= numSampleSets){
    message("Found all chunk files, skipping chunk generation")
    skipchunk <- TRUE
  } else if(resume && (length(chunkofiles) > 0)){
    # Extract chunk numbers from existing files
    existing_chunk_nums <- as.integer(gsub("chunk_([0-9]+)\\.RData", "\\1", chunkofiles))
    existing_chunk_nums <- existing_chunk_nums[!is.na(existing_chunk_nums)]
    # Only generate missing chunks
    repsvec <- setdiff(seq_len(numSampleSets), existing_chunk_nums)
    message(paste0("Found ", length(existing_chunk_nums), " existing chunk files"))
    message(paste0("Will generate ", length(repsvec), " missing chunks: ", 
                   if(length(repsvec) <= 20) paste(repsvec, collapse=", ") else paste(c(head(repsvec, 10), "...", tail(repsvec, 10)), collapse=", ")))
  } else if(resume && length(chunkofiles) == 0){
    message("No chunk files found, starting chunk generation from beginning")
  } else if(length(chunkofiles) > numSampleSets){
    stop(paste0("Found ", length(chunkofiles), " chunk files but expected ", numSampleSets, "; use a different folder or check for duplicates"))
  }
  
  # Generate chunks if needed
  if(!skipchunk){
    # Prepare sample-level metadata to pass to chunk generator (small objects)
    # Metadata columns are requested by name, so a metadata file lacking
    # superpopulation columns aborts the whole run with 'undefined columns
    # selected'. That matters because --pop/--superpop are on by default, so a
    # missing column must not be fatal. Take the columns that exist.
    meta_cols_for_samples <- intersect(
      c("sample_name", "population_code", "super_population_code",
        "population_description", "super_population_description"),
      colnames(null_data))
    unique_samples_df <- unique(null_data[, meta_cols_for_samples, drop = FALSE])
    sample_name_by_row <- null_data$sample_name
    variant_ids_by_row <- variant_ids

    generate_hill_chunks(
      repsvec = repsvec,
      filtered_data = filtered_data,
      counts = counts,
      counts_type = "matrix",
      sample_info = unique_samples_df,
      sample_name_by_row = sample_name_by_row,
      variant_ids_by_row = variant_ids_by_row,
      null_data_split_by_pop_by_variant = null_data_split_by_pop_by_variant,
      null_data_split_by_super_pop_by_variant = null_data_split_by_super_pop_by_variant,
      pop_idx = pop_idx,
      super_idx = super_idx,
      sample_groups = sample_groups,
      min_sample_size = min_sample_size,
      downsample_samples = downsample_samples,
      evaluatePopulation = evaluatePopulation,
      evaluateSuperPopulation = evaluateSuperPopulation,
      downSamplePop = downSamplePop,
      downSampleSuperPop = downSampleSuperPop,
      numSampleSets = numSampleSets,
      numVariantSets = numVariantSets,
      numvariants = numvariants,
      total_variants = total_variants,
      seed = seed,
      numCores = numCores,
      chunkDir = chunkDir,
      # variantEmpiricalNull removed
      maxGB = maxGB,
      keep_cols = keep_cols,
      qHillNumber = qHillNumber,
      show_progress = show_progress
    )
  }

  # Calculate Hill numbers from pre-generated chunk RData files
  filelist <- list.files(chunkDir, pattern="^chunk_[0-9]+\\.RData$", full.names = TRUE)
  # Sort chunk files by numeric chunk ID to ensure deterministic processing order
  if(length(filelist) > 0) {
    chunk_nums <- as.integer(gsub(".*chunk_([0-9]+)\\.RData$", "\\1", basename(filelist)))
    filelist <- filelist[order(chunk_nums)]
  }
  message(paste0("Calculating Hill numbers from ", length(filelist), " sample sets with ", 
                numVariantSets, " pre-generated variant sets each = ", 
                length(filelist) * numVariantSets, " total replicates"))
  message(paste0("Using ", numCores, " workers"))
  
  max_size_bytes <- 1024^3 * maxGB
  options(future.progress = FALSE)
  
  select_plan(numCores, parallelBackend)
  message(paste0("Max global size: ", round(max_size_bytes / 1024^3, 1), " GB"))
  message("Parallelization strategy: Sequential sample sets, parallel variant sets")

  # Warn about peak phase-2 memory before the workers fork, and name the
  # degenerate two-arm case. Each parallel task loads a whole batch file, i.e.
  # num_batches_per_set variant sets, and each variant set holds a nullmat with
  # one row per observation of the sampled variants.
  #
  # Measured on the full 3.48M-row table: a batch file is 1.30 GiB on disk and
  # 3.24 GiB once loaded, so deserialised objects run ~2.5x the file size, and
  # hillR's own working set adds roughly as much again. Workers are SIGKILLed
  # when that overflows the container, and the parent then reports only a bare
  # FutureInterruptError with "parallel job did not deliver a result", which
  # gives no hint that memory was the cause. Naming it here turns a silent
  # overnight death into a readable warning.
  IN_MEM_FACTOR <- 2.5
  est_gb <- tryCatch({
    batches <- list.files(chunkDir, pattern = "_batch_[0-9]+\\.RData$", full.names = TRUE)
    per_set <- if (length(filelist) > 0 && length(batches) > 0) {
      length(batches) / length(filelist)
    } else 0
    if (length(batches) > 0 && per_set > 0) {
      sum(file.size(batches)) / length(batches) * per_set * min(numCores, per_set) *
        IN_MEM_FACTOR / 1024^3
    } else 0
  }, error = function(e) 0)
  if (est_gb > 0) {
    message(sprintf("Estimated peak phase-2 memory: ~%.1f GB across %d workers (in-memory factor %.1fx)",
                    est_gb, numCores, IN_MEM_FACTOR))
    if (est_gb > 6) {
      warning(sprintf(
        paste0("Estimated phase-2 memory (~%.1f GB) is likely to exceed the ",
               "container limit. Each worker loads a whole batch file, so workers ",
               "may be killed and this phase will then fail with a bare ",
               "FutureInterruptError. Supply a real --filtered arm (which shrinks ",
               "the batch files in proportion to the filtered variant count), ",
               "lower --variant-sets, or raise the container memory."), est_gb))
    }
  }
  
  hillnumlist <- tryCatch({
    local({
      # Options set here are automatically restored when block exits
      options(future.globals.maxSize = max_size_bytes)
      options(future.wait.timeout = 31536000)
      
      # Process sample sets SEQUENTIALLY
      lapply(filelist, function(chunk_file) {
      e <- new.env()
      load(chunk_file, envir=e)
      chunk_index <- e$chunk_index
      
      # Debug: Check chunk structure
      if(is.null(chunk_index)) {
        stop(sprintf("chunk_index loaded from %s is NULL", basename(chunk_file)))
      }
      chunk_names <- names(chunk_index)
      if(is.null(chunk_names) || length(chunk_names) == 0) {
        stop(sprintf("chunk_index from %s has no named elements. Structure: %s", basename(chunk_file), paste(capture.output(str(chunk_index)), collapse="\n")))
      }
      message(sprintf("Processing %s: chunk_index has fields [%s]", basename(chunk_file), paste(chunk_names, collapse=", ")))

      # Process each pre-generated variant set. Support three chunk formats:
      # - legacy: chunk_index$variant_sets is an in-memory list of variant_set objects
      # - streaming: chunk_index$variant_files is a character vector of filenames (one per variant set)
      # - batched: chunk_index$batch_files is a list of batch metadata (multiple variant sets per file)
      variant_results <- list()
      if(!is.null(chunk_index$variant_sets)) {
        variant_iter <- chunk_index$variant_sets
        variant_mode <- "legacy"
        message(sprintf("  Using LEGACY format: in-memory variant sets"))
      } else if(!is.null(chunk_index$batch_files)) {
        variant_iter <- chunk_index$batch_files
        variant_mode <- "batched"
        message(sprintf("  Using BATCHED format: %d batch files", length(variant_iter)))
      } else if(!is.null(chunk_index$variant_files)) {
        variant_iter <- chunk_index$variant_files
        variant_mode <- "streaming"
        message(sprintf("  Using STREAMING format: %d individual variant files", length(variant_iter)))
      } else {
        stop("chunk_index does not contain variant_sets, variant_files, or batch_files")
      }
      # Helper to process a single variant_set object
      process_variant_set <- function(variant_set, chunk_sample_set_id, numVariantSets, qHillNumber) {
        sampled_mat <- variant_set$nullmat
        sampled_pop_by_var <- variant_set$pop_by_var
        sampled_superpop_by_var <- variant_set$superpop_by_var
        local_idx_pop_variant <- variant_set$local_idx_pop
        local_idx_superpop_variant <- variant_set$local_idx_superpop
        
        # Calculate Hill numbers
        vallist <- list()
        vallist$hillnumbers <- hillnum(
          sampled_variants_mat = sampled_mat,
          sampled_variants_mat_pop_by_var = sampled_pop_by_var,
          sampled_variants_mat_superpop_by_var = sampled_superpop_by_var,
          local_idx_pop = local_idx_pop_variant,
          local_idx_superpop = local_idx_superpop_variant,
          qHillNumber = qHillNumber,
          rows_keep = seq_len(nrow(sampled_mat)),
          variants_keep = NULL  # Added for compatibility
        )
        
        # Attach metadata
        vallist$metadata <- list(
          sample_set_id = chunk_sample_set_id,
          variant_set_id = variant_set$variant_set_id,
          replicate_id = (chunk_sample_set_id - 1) * numVariantSets + variant_set$variant_set_id
        )
        
        return(vallist)
      }
      
      # Parallelize variant set processing using future_lapply
      chunk_sample_set_id <- chunk_index$sample_set_id
      
      if(variant_mode == "batched") {
        # New batched format: each batch file contains multiple variant sets
        # Process each batch in parallel, then flatten results
        batch_results <- future_lapply(variant_iter, function(batch_info) {
          batch_file <- batch_info$file
          if(!file.exists(batch_file)) {
            warning(paste0("Batch file not found: ", batch_file))
            return(list())
          }
          
          # Load batch_data from disk
          batch_data <- readRDS(batch_file)
          # Process each variant set in this batch with error isolation
          # so a single failed variant set does not kill the whole batch worker
          batch_results <- lapply(batch_data, function(vs) {
            tryCatch(
              process_variant_set(vs, chunk_sample_set_id, numVariantSets, qHillNumber),
              error = function(e) {
                warning(paste0(
                  "Failed variant_set ", vs$variant_set_id,
                  " (sample set ", chunk_sample_set_id,
                  "): ", e$message
                ))
                return(list(
                  sample_set_id = vs$sample_set_id,
                  variant_set_id = vs$variant_set_id,
                  error = e$message
                ))
              }
            )
          })

          # Cleanup: batch_data is a local binding, not an environment, so drop
          # it directly. (The streaming branch below loads into an env, hence
          # the different cleanup there.)
          rm(batch_data)
          gc(verbose = FALSE)
          return(batch_results)
        }, future.seed=TRUE, future.globals=structure(TRUE, add=c("process_variant_set", "chunk_sample_set_id", "numVariantSets", "qHillNumber", "hillnum", "sanitize_counts", "safe_hill", "list_transpose")), future.packages=c("hillR", "purrr"), future.stdout=FALSE)
        
        # Flatten: batch_results is a list of lists, convert to single list
        variant_results <- unlist(batch_results, recursive = FALSE)
        
      } else if(variant_mode == "streaming") {
        # Streaming format: each file contains one variant_set
        if(!all(file.exists(variant_iter))) {
          missing <- variant_iter[!file.exists(variant_iter)]
          stop(paste0("Missing variant files: ", paste(head(missing, 5), collapse=", ")))
        }
        
        variant_results <- future_lapply(variant_iter, function(vf) {
          ve <- new.env()
          load(vf, envir=ve)
          if(!exists("variant_set", envir=ve)) return(NULL)
          vs <- get("variant_set", envir=ve)
          res <- process_variant_set(vs, chunk_sample_set_id, numVariantSets, qHillNumber)
          rm(list=ls(envir=ve), envir=ve)
          gc(verbose = FALSE)
          return(res)
        }, future.seed=TRUE, future.globals=structure(TRUE, add=c("process_variant_set", "chunk_sample_set_id", "numVariantSets", "qHillNumber", "hillnum", "sanitize_counts", "safe_hill", "list_transpose")), future.packages=c("hillR", "purrr"), future.stdout=FALSE)
        
      } else if(variant_mode == "legacy") {
        # Legacy format: variant_sets already in memory
        variant_results <- future_lapply(variant_iter, function(vs) {
          process_variant_set(vs, chunk_sample_set_id, numVariantSets, qHillNumber)
        }, future.seed=TRUE, future.globals=structure(TRUE, add=c("process_variant_set", "chunk_sample_set_id", "numVariantSets", "qHillNumber", "hillnum", "sanitize_counts", "safe_hill", "list_transpose")), future.packages=c("hillR", "purrr"), future.stdout=FALSE)
        
      } else {
        stop("Unsupported variant_mode")
      }
      
      # Clear chunk data from memory
      rm(chunk_index, e)
      gc(verbose = FALSE)
      
      return(variant_results)
    })
    # End of lapply over sample sets (sequential)
    })
  }, FutureInterruptError = function(e) {
    message("\n=== INTERRUPTED ===")
    message("Hill number calculation was interrupted.")
    message("Chunk files are saved. You can resume by running the script again with --resume flag.")
    return(NULL)
  }, interrupt = function(e) {
    message("\n=== INTERRUPTED ===")
    message("Hill number calculation was interrupted.")
    message("Chunk files are saved. You can resume by running the script again with --resume flag.")
    return(NULL)
  }, error = function(e) {
    message(paste0("\nERROR during Hill number calculation: ", e$message))
    stop(e)
  })
  
  # Check if interrupted
  if(is.null(hillnumlist)) {
    message("Exiting due to interruption...")
    quit(save="no", status=1)
  }
  
  message("=== Starting reorganization of Hill number results ===")
  message(paste0("hillnumlist class: ", class(hillnumlist)))
  message(paste0("hillnumlist length: ", length(hillnumlist)))
  if(length(hillnumlist) > 0) {
    message(paste0("First element class: ", class(hillnumlist[[1]])))
    message(paste0("First element length: ", length(hillnumlist[[1]])))
  }
  
  # Flatten the nested list structure
  message("Flattening nested list with unlist(recursive=FALSE)...")
  hillnumlist <- unlist(hillnumlist, recursive = FALSE)
  gc(verbose = FALSE)  # Clean up after flattening
  
  hillnumlist <- list_transpose(hillnumlist)
  gc(verbose = FALSE)  # Clean up after transpose
  
  message("After list_transpose:")
  message(paste0("  Names: ", paste(names(hillnumlist), collapse=", ")))
  message(paste0("  Has hillnumbers: ", !is.null(hillnumlist$hillnumbers)))
  if(!is.null(hillnumlist$hillnumbers)) {
    message(paste0("  hillnumbers length: ", length(hillnumlist$hillnumbers)))
    if(length(hillnumlist$hillnumbers) > 0) {
      message(paste0("  First hillnumber names: ", paste(names(hillnumlist$hillnumbers[[1]]), collapse=", ")))
    }
  }
  
  # Reorganize structure for easier analysis
  if(!is.null(hillnumlist$hillnumbers) && length(hillnumlist$hillnumbers) > 0) {
    message("Reorganizing Hill number structure...")
    component_names <- names(hillnumlist$hillnumbers[[1]])
    message(paste0("  Component names: ", paste(component_names, collapse=", ")))
    
    result <- list()
    
    for(comp_name in component_names) {
      message(paste0("  Processing ", comp_name, "..."))
      comp_values <- lapply(hillnumlist$hillnumbers, function(x) x[[comp_name]])
      
      non_null_indices <- which(!sapply(comp_values, is.null))
      message(paste0("    Non-NULL values: ", length(non_null_indices), " / ", length(comp_values)))
      
      if(length(non_null_indices) == 0) {
        message("    All NULL - storing as-is")
        result[[comp_name]] <- comp_values
      } else {
        first_non_null <- comp_values[[non_null_indices[1]]]
        message(paste0("    First non-NULL class: ", class(first_non_null)[1]))
        message(paste0("    Is list: ", is.list(first_non_null)))
        message(paste0("    Is data.frame: ", is.data.frame(first_non_null)))
        if(is.list(first_non_null)) {
          message(paste0("    Length: ", length(first_non_null)))
          if(length(first_non_null) > 0 && length(names(first_non_null)) > 0) {
            message(paste0("    Has names: TRUE"))
            message(paste0("    Names: ", paste(head(names(first_non_null), 3), collapse=", ")))
          }
        }
        if(is.list(first_non_null) && !is.data.frame(first_non_null)) {
          message(paste0("    Has names: ", length(names(first_non_null)) > 0))
          if(length(names(first_non_null)) > 0) {
            message(paste0("    Names: ", paste(head(names(first_non_null), 3), collapse=", ")))
          }
        }
        
        if(is.list(first_non_null) && !is.data.frame(first_non_null) && length(names(first_non_null)) > 0) {
          # This is a by-population or by-superpopulation component
          all_pop_names <- unique(unlist(lapply(comp_values, function(x) if(!is.null(x)) names(x) else NULL)))
          message(paste0("    Found ", length(all_pop_names), " populations/groups"))
          
          if(length(all_pop_names) > 0) {
            pop_transposed <- list()
            for(pop_name in all_pop_names) {
              # Extract all replicates for this population
              pop_replicates <- lapply(comp_values, function(x) {
                if(!is.null(x) && pop_name %in% names(x)) {
                  x[[pop_name]]
                } else {
                  NULL
                }
              })
              
              # Remove NULLs and bind into data frame
              pop_replicates <- pop_replicates[!sapply(pop_replicates, is.null)]
              if(length(pop_replicates) > 0) {
                # Bind rows - each replicate should be a data frame with 1 row
                # Add replicate_id column to track which replicate each row came from
                pop_transposed[[pop_name]] <- tryCatch({
                  bound_df <- do.call(rbind, lapply(seq_along(pop_replicates), function(i) {
                    df <- as.data.frame(pop_replicates[[i]])
                    df$replicate <- i
                    df
                  }))
                  rownames(bound_df) <- NULL  # Reset rownames to clean sequential numbers
                  message(paste0("      ", pop_name, ": bound ", nrow(bound_df), " rows"))
                  bound_df
                }, error = function(e) {
                  message(paste0("      ", pop_name, ": error binding - ", e$message))
                  pop_replicates
                })
              } else {
                pop_transposed[[pop_name]] <- NULL
              }
            }
            result[[comp_name]] <- pop_transposed
          } else {
            message("    No population names found")
            result[[comp_name]] <- comp_values
          }
        } else if(is.data.frame(first_non_null) && !is.null(rownames(first_non_null)) && 
                  length(rownames(first_non_null)) > 0 && rownames(first_non_null)[1] != "1") {
          # This is a data frame with population names in rownames (from binding hill_taxa_parti results)
          message("    Data frame with rownames - binding all replicates first...")
          comp_values <- comp_values[!sapply(comp_values, is.null)]
          bound_df <- do.call(rbind, lapply(seq_along(comp_values), function(i) {
            df <- as.data.frame(comp_values[[i]])
            df$replicate <- i
            df
          }))
          
          # Extract population names from rownames (strip trailing digits)
          pop_codes <- gsub("[0-9]+$", "", rownames(bound_df))
          unique_pops <- unique(pop_codes)
          message(paste0("    Detected ", length(unique_pops), " populations from rownames: ", 
                        paste(head(unique_pops, 3), collapse=", ")))
          
          # Split into list by population
          pop_list <- list()
          for(pop in unique_pops) {
            pop_rows <- bound_df[pop_codes == pop, , drop=FALSE]
            rownames(pop_rows) <- NULL
            pop_list[[pop]] <- pop_rows
            message(paste0("      ", pop, ": ", nrow(pop_rows), " rows"))
          }
          result[[comp_name]] <- pop_list
        } else {
          # This is TotalData or similar - bind all replicates
          message("    Binding all replicates as single data frame...")
          comp_values <- comp_values[!sapply(comp_values, is.null)]
          if(length(comp_values) > 0 && is.data.frame(comp_values[[1]])) {
            result[[comp_name]] <- tryCatch({
              bound_df <- do.call(rbind, lapply(seq_along(comp_values), function(i) {
                df <- as.data.frame(comp_values[[i]])
                df$replicate <- i
                df
              }))
              rownames(bound_df) <- NULL  # Reset rownames to clean sequential numbers
              message(paste0("    Bound ", nrow(bound_df), " rows"))
              bound_df
            }, error = function(e) {
              message(paste0("    Error binding: ", e$message))
              comp_values
            })
          } else {
            message("    Not a data frame or empty - storing as list")
            result[[comp_name]] <- comp_values
          }
        }
      }
    }
    
    result$metadata <- hillnumlist$metadata
    message("Reorganization complete. Result structure:")
    message(paste0("  Names: ", paste(names(result), collapse=", ")))
  } else {
    message("Skipping reorganization - hillnumbers is NULL or empty")
    result <- hillnumlist
  }
  
  message("Returning result with top-level names: ", paste(names(result), collapse=", "))
  
  # Convert population/superpopulation codes to descriptions for plot labeling
  pop_code_to_desc <- create_code_mapping(null_data, "population_code", "population_description")
  superpop_code_to_desc <- create_code_mapping(null_data, "super_population_code", "super_population_description")
  
  # Apply mappings to rename list elements
  pop_components <- c("SplitByPopulation", "SplitByPopulationByVariant", 
                      "PairwiseSplitByPopulation", "PairwiseSplitByPopulationByVariant")
  superpop_components <- c("SplitBySuperPopulation", "SplitBySuperPopulationByVariant",
                           "PairwiseSplitBySuperPopulation", "PairwiseSplitBySuperPopulationByVariant")
  
  result <- rename_list_components(result, pop_components, pop_code_to_desc)
  result <- rename_list_components(result, superpop_components, superpop_code_to_desc)
  
  gc(verbose = FALSE)
  return(result)
}
# ---- end hill_calculations.R ----

# ---- begin variance_analysis.R::extract_sample_sets_from_chunks (this function only) ----
extract_sample_sets_from_chunks <- function(chunkDir) {
  filelist <- list.files(chunkDir, pattern="^chunk_[0-9]+\\.RData$", full.names = TRUE)
  
  # Sort by chunk number to ensure consistent ordering
  if(length(filelist) > 0) {
    chunk_nums <- as.integer(gsub(".*chunk_([0-9]+)\\.RData$", "\\1", basename(filelist)))
    filelist <- filelist[order(chunk_nums)]
  }
  
  chunk_sample_sets <- lapply(filelist, function(x) {
    e <- new.env()
    load(x, envir=e)
    # Fixed: was e$chunk$samples_to_keep, should be e$chunk_index$samples_to_keep
    return(e$chunk_index$samples_to_keep)
  })
  return(chunk_sample_sets)
}
# ---- end variance_analysis.R::extract_sample_sets_from_chunks ----

# =============================================================================
# END vendored reference code
# =============================================================================

# -----------------------------------------------------------------------------
# Mode: prepare
# -----------------------------------------------------------------------------
run_prepare <- function(input, metadata_path, pipelines_arg, max_variants = 0L,
                        filtered_path = NULL) {
  counts <- read_any(input)
  if (!is.data.frame(counts)) {
    stop("--mode prepare expects a count table, got ", class(counts)[1], ".")
  }
  if (nrow(counts) == 0) stop("Count table has no rows: ", input)
  if (!"sample_name" %in% names(counts)) stop("Count table has no sample_name column.")
  if (!"variant_called" %in% names(counts)) stop("Count table has no variant_called column.")

  pipelines <- resolve_pipelines(counts, pipelines_arg)
  message("  pipelines (", length(pipelines), "): ", paste(pipelines, collapse = ", "))

  md <- if (!is.null(metadata_path)) normalise_metadata(read_any(metadata_path)) else NULL

  # Both arms go through identical treatment, so that they differ only in which
  # variants they hold. Each arm gets its own metadata merge, because the two
  # tables come from separate VCF2COUNTS/COUNTS2RDS runs and may carry different
  # sample sets.
  normalise_arm <- function(df, arm) {
    for (cl in DEMO_COLS) if (!cl %in% names(df)) df[[cl]] <- NA_character_

    if (!is.null(md)) {
      df <- merge(df, md, by = "sample_name", all.x = TRUE, sort = FALSE,
                  suffixes = c("", ".md"))
      for (cl in DEMO_COLS) {
        mdcol <- paste0(cl, ".md")
        if (mdcol %in% names(df)) {
          df[[cl]] <- ifelse(is.na(df[[cl]]) | df[[cl]] == "",
                             df[[mdcol]], df[[cl]])
          df[[mdcol]] <- NULL
        }
      }
      unmatched <- setdiff(unique(df$sample_name), md$sample_name)
      if (length(unmatched) > 0) {
        warning("No metadata for ", length(unmatched), " ", arm, " sample(s): ",
                paste(head(unmatched, 10), collapse = ", "))
      }
    }

    df <- df[, intersect(c(ANALYSIS_METADATA_COLS, pipelines), names(df)),
             drop = FALSE]
    for (cl in pipelines) df[[cl]] <- as.numeric(df[[cl]])

    # In place of the reference data_io.R dedupe, which exists only to make
    # paste0(sample_name, "_", variant_called) unique row names downstream.
    dup <- duplicated(df[, c("sample_name", "variant_called")])
    if (any(dup)) {
      message("  collapsing ", sum(dup), " duplicate (sample_name, variant_called) row(s) in the ", arm, " arm")
      num_cols <- intersect(pipelines, names(df))
      DT <- data.table::as.data.table(df)
      df <- as.data.frame(
        DT[, lapply(.SD, max), by = c("sample_name", "variant_called"),
           .SDcols = c(num_cols, setdiff(ANALYSIS_METADATA_COLS[3:6], c("sample_name", "variant_called")))])
    }
    df
  }

  null_data <- normalise_arm(counts, "null")

  # Without --filtered both arms are the same table, so numvariants equals the
  # whole variant set, every replicate resamples the entire pool without
  # replacement, and the null is a permutation of the observed data.
  if (!is.null(filtered_path)) {
    filtered_data <- normalise_arm(read_any(filtered_path), "filtered")
  } else {
    filtered_data <- null_data
  }

  demoRunandEthnicity <- unique(counts[, intersect(
    c("sample_name", DEMO_COLS), names(counts)), drop = FALSE])

  numvariants    <- length(unique(filtered_data$variant_called))
  total_variants <- unique(null_data$variant_called)

  # numvariants drives how many variants each Monte Carlo replicate resamples,
  # so it is the main cost knob of --mode null. By design it is the size of the
  # FILTERED set; --max-variants stays as an override for a run whose filtered
  # arm is still the full table.
  if (max_variants > 0L && numvariants > max_variants) {
    message("  capping numvariants: ", numvariants, " -> ", max_variants,
            " (--max-variants)")
    numvariants <- max_variants
  }

  message("  samples: ", length(unique(null_data$sample_name)),
          "  null rows: ", nrow(null_data),
          "  filtered rows: ", nrow(filtered_data),
          "  numvariants: ", numvariants,
          "  total_variants: ", length(total_variants))
  if (is.null(filtered_path)) {
    message("  NOTE: no --filtered table given, so the null arm is a permutation ")
    message("        of the observed data and numvariants is the whole table.")
  }

  list(data_list = null_data,
       filtered_data = filtered_data,
       null_data = null_data,
       demoRunandEthnicity = demoRunandEthnicity,
       numvariants = numvariants,
       total_variants = total_variants)
}

# -----------------------------------------------------------------------------
# Mode: null
# -----------------------------------------------------------------------------
# numvariants_cap is the --max-variants value. It has to be applied here as well
# as in run_prepare, because the Nextflow module feeds --mode null a bare
# counts2rds data.frame rather than the output of --mode prepare. Without the
# cap this branch derives numvariants = the whole table (~3M variants), and
# every Monte Carlo replicate resamples all of it, so null mode never finishes.
run_null <- function(input, numvariants_cap = 0L) {
  prep <- read_any(input)
  # Accept both the list from --mode prepare and a bare data.frame
  # (counts2rds output) — extract the needed pieces in the latter case.
  if (is.data.frame(prep)) {
    # counts2rds output: a plain data.frame
    filtered_data  <- prep
    null_data      <- prep
    numvariants    <- length(unique(prep$variant_called))
    total_variants <- unique(prep$variant_called)
    if (numvariants_cap > 0L && numvariants > numvariants_cap) {
      message("  capping numvariants: ", numvariants, " -> ", numvariants_cap,
              " (--max-variants)")
      numvariants <- numvariants_cap
    }
    message("  [null] derived numvariants=", numvariants,
            " total_variants=", length(total_variants))
  } else if (!is.list(prep) || is.null(prep$null_data)) {
    stop("--mode null needs a data.frame (counts2rds output) or the ",
         "output of --mode prepare (a list holding null_data, filtered_data, ",
         "numvariants and total_variants).")
  } else {
    # --mode prepare output: the prepared list
    filtered_data  <- prep$filtered_data
    null_data      <- prep$null_data
    numvariants    <- prep$numvariants
    total_variants <- prep$total_variants
  }
  if (!dir.exists(chunkDir)) {
    dir.create(chunkDir, recursive = TRUE, showWarnings = FALSE)
  }

  calculateNullHill(
    null_data      = null_data,
    filtered_data  = filtered_data,
    total_variants = total_variants,
    numvariants    = numvariants,
    numCores       = numCores,
    seed           = seed,
    numSampleSets  = numSampleSets,
    numVariantSets = numVariantSets,
    qHillNumber    = qHillNumber,
    chunkDir       = chunkDir,
    show_progress  = TRUE,
    evaluateSuperPopulation = evaluateSuperPopulation,
    evaluatePopulation      = evaluatePopulation,
    resume                  = resume,
    downSamplePop           = downSampleByPop,
    downSampleSuperPop      = downSampleBySuperPop,
    maxGB                   = maxGB)
}

# -----------------------------------------------------------------------------
# Mode: filtered
# -----------------------------------------------------------------------------
run_filtered <- function(input) {
  prep <- read_any(input)
  filtered_data <- if (is.list(prep) && !is.data.frame(prep) && !is.null(prep$filtered_data)) {
    prep$filtered_data
  } else if (is.data.frame(prep)) {
    prep
  } else {
    stop("--mode filtered needs the output of --mode prepare, or a count table.")
  }

  # Reuse the null stage's checkpointed sample sets, so filtered and null are
  # averaged over the same replicates. main.R does the same.
  chunk_sample_sets <- NULL
  if ((downSampleByPop || downSampleBySuperPop) && dir.exists(chunkDir)) {
    chunk_sample_sets <- extract_sample_sets_from_chunks(chunkDir)
    message("  reusing ", length(chunk_sample_sets), " sample set(s) from ", chunkDir)
  }

  res <- calculate_filtered_hill(
    filtered_data, qHillNumber, evaluatePopulation, evaluateSuperPopulation,
    downSampleByPop, downSampleBySuperPop, chunk_sample_sets,
    numVariantSets = numVariantSets, seed = seed, variant_sample_fraction = 0.8)

  res$has_replicates <- !is.null(chunk_sample_sets)
  res$num_replicates <- if (is.null(chunk_sample_sets)) 1L else length(chunk_sample_sets)
  res
}

# -----------------------------------------------------------------------------
# Dispatch
# -----------------------------------------------------------------------------
message("=== hill.R ", mode, " ===")
message("  input : ", input)
message("  output: ", outfile)

result <- switch(mode,
  prepare  = run_prepare(input, opts$metadata, opts$pipelines, maxVariants, filteredPath),
  null     = run_null(input, maxVariants),
  filtered = run_filtered(input))

saveRDS(result, outfile)
message("  wrote ", outfile, " (", round(file.size(outfile) / 1024^2, 1), " MB)")
message("=== hill.R ", mode, ": done ===")
