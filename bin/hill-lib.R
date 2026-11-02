# hill-lib.R -- shared diversity functions for the split hill pipeline.
#
# EXTRACTED from bin/hill-fixed.R, not rewritten. Lineage and exact source
# ranges are recorded per block below so the extract can be re-verified against
# the original after any change there.
#
# Deliberately contains NO commandArgs()/docopt and NO top-level side effects:
# sourcing this file defines functions and nothing else, so it is safe to source
# from a Nextflow task, an interactive session, or a test harness.
#
# Requires: hillR, purrr (>= 1.0.0, for list_transpose), data.table, dplyr.

suppressPackageStartupMessages({
  library(hillR)
  library(purrr)
  library(data.table)
  library(dplyr)
  # future + future.apply: generate_hill_chunks() calls future_lapply() and
  # select_plan() sets the plan. Both were attached by the original script.
  library(future)
  library(future.apply)
})

# ---- begin sanitize_counts (hill-fixed.R 243-251) ----
sanitize_counts <- function(mat) {
  if(is.null(mat)) return(mat)
  m <- as.matrix(mat)
  # Coerce to numeric where possible
  storage.mode(m) <- "numeric"
  # Replace NA with zero
  m[(is.na(m)) | (m < 0)] <- 0
  return(m)
}
# ---- end sanitize_counts ----

# ---- begin safe_apply (hill-fixed.R 474-484) ----
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
# ---- end safe_apply ----

# ---- begin safe_hill (hill-fixed.R 554-601) ----
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
# ---- end safe_hill ----

# ---- begin hillnum (hill-fixed.R 616-710) ----
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
# ---- end hillnum ----

# ---- begin process_variant_set (hill-fixed.R 1986-2014) ----
# LIFTED TO TOP LEVEL. In hill-fixed.R this is a closure defined inside
# calculateNullHill(), which is why the original cannot process a batch in a
# fresh process. Only the 6-space nesting indent was removed; the body is
# byte-identical. numVariantSets/qHillNumber stay explicit arguments, so no
# enclosing-scope state is required.
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
# ---- end process_variant_set ----
# ---- begin create_code_mapping (hill-fixed.R 331-334) ----
create_code_mapping <- function(data, code_col, desc_col) {
  mapping <- setNames(data[[desc_col]], data[[code_col]])
  mapping[!duplicated(names(mapping))]
}
# ---- end create_code_mapping ----

# ---- begin rename_list_components (hill-fixed.R 451-464) ----
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
# ---- end rename_list_components ----

# ---- begin extract_sample_sets_from_chunks (hill-fixed.R 2307-2323) ----
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
# ---- end extract_sample_sets_from_chunks ----

# ---- begin reorganise_hill_results (hill-fixed.R 2144-2299) ----
# WRAPPED VERBATIM. In hill-fixed.R this is inline inside calculateNullHill()
# and reads two things from the enclosing scope: `hillnumlist` (the
# transposed results) and `null_data` (only for the code->description
# mapping at the end). Both are now explicit arguments, and the two \\<-
# assignments that fed them are hoisted out as named locals. The body is
# otherwise unchanged, deliberately: this block is the one piece whose
# behaviour is hard to re-derive (three different binding strategies chosen
# by inspecting the first non-NULL value), so it is copied, not rewritten.
reorganise_hill_results <- function(hillnumlist, null_data) {

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
}
# ---- end reorganise_hill_results ----

# ---- begin generate_batches (hill-fixed.R 1721-1885) ----
# WRAPPED VERBATIM from calculateNullHill()'s chunking branch. This block is
# NOT re-derived: it builds the population/superpopulation row indices, sets
# RNGkind("L'Ecuyer-CMRG") and the seed, computes keep_cols, and drives
# generate_hill_chunks(). All of that is load-bearing for reproducibility --
# a reimplementation that got RNGkind or the index construction subtly wrong
# would still run and would still produce plausible Hill numbers, just not
# the same ones. So it is copied.
#
# Changes made in wrapping, all mechanical:
#   - the enclosing scope's values became arguments
#   - max_size_bytes/options() stays: it must run before any worker forks
#   - numCores is not used by this branch
#   - the trailing 'return(invisible(NULL))' became the return value
generate_batches <- function(null_data, filtered_data, total_variants, numvariants,
                             seed, numSampleSets, numVariantSets, qHillNumber,
                             chunkDir, resume = FALSE,
                             evaluatePopulation = FALSE,
                             evaluateSuperPopulation = FALSE,
                             downSamplePop = FALSE, downSampleSuperPop = FALSE,
                             maxGB = 100) {
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
    # Fixed at 1, not inherited: batch GENERATION is no longer the parallel step
    # (the scheduler fans out over the batches it writes), and the original's
    # numCores was the CLI global for phase 2. Passing 1 keeps generation
    # single-threaded and deterministic, and avoids a second, competing
    # parallelism inside a task that is already one scheduler slot.
    numCores = 1L,
    chunkDir = chunkDir,
    # variantEmpiricalNull removed
    maxGB = maxGB,
    keep_cols = keep_cols,
    qHillNumber = qHillNumber,
    show_progress = show_progress
  )
}
}
# ---- end generate_batches ----

# ---- begin read_any (hill-fixed.R 155-165) ----
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
# ---- end read_any ----

# ---- begin generate_hill_chunks (hill-fixed.R 880-1455) ----
# EXTRACTED VERBATIM (576 lines). This is the Monte Carlo sampler: the seed
# arithmetic (seed + sample_set_id, and seed + sample_set_id*10000 +
# variant_set_id), the per-batch sampling, the atomic temp-file-then-rename
# writes, and the index construction all live here. It calls nothing
# outside itself and base/hillR/data.table, so it lifts out cleanly.
#
# This replaced a hand-written reimplementation attempt, which is roughly
# what one should expect: the preamble that feeds this function builds
# population row indices and calls RNGkind("L'Ecuyer-CMRG"), and getting any
# of that subtly wrong yields plausible Hill numbers that are not the same
# ones. Copy it; do not rewrite it.
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
  # Edited for the extract: the original reads `parallelBackend` from the script's
  # CLI-derived global. Lifted out of that scope, default to "auto" so the
  # documented behaviour (multicore where the platform supports fork) is kept.
  # Batch generation is not the parallel step any more -- the scheduler is -- and
  # this is called with numCores = 1 by hill-batchgen.R, so this only ever sets a
  # sequential plan.
  select_plan(numCores, get0("parallelBackend", ifnotfound = "auto"))
  
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
# ---- end generate_hill_chunks ----

# ---- begin select_plan (hill-fixed.R 121-136) ----
# Parallel backend selection, called by generate_hill_chunks(). parallelBackend
# is an explicit argument here rather than the script's CLI-derived global.
select_plan <- function(ncores, backend = "auto") {
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
# ---- end select_plan ----

