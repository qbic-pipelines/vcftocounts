#!/usr/bin/env Rscript
# =============================================================================
# pipecomp.R -- do the pipeline differences agree between the filtered and the
#               null Hill analysis?
#
# Scope: the filtered-vs-null pattern comparison stage only.
#   Input   TWO rds files, both written by hill.R:
#             --filtered  the value of hill.R --mode filtered (a point estimate)
#             --null      the value of hill.R --mode null (Monte Carlo replicates)
#   Output  one CSV set per level analyzed (total / by population / by
#           superpopulation).
#
# This is the stage main.R runs after analyze_pipeline_differences:
#   filtered_vs_null_results$total_pattern <- compare_pipeline_patterns(
#     filtered_pairwise = filt_total, null_pairwise = null_total, ...)
#
# What it answers. analyze_pipeline_differences asks "which pipeline pairs
# diverge, on the null data?" This stage asks whether that answer survives
# filtering. For each pipeline pair it takes the filtered beta as one value and
# the null betas as a distribution, correlates the filtered vector against EACH
# null replicate, and reports the distribution of those correlations:
#   Spearman rho > 0.7   the pattern is PRESERVED -- filtering moves all pairs
#                        the same way
#   Spearman rho < -0.7  the pattern is REVERSED -- filtering flips which pair
#                        looks most divergent
#   |rho| < 0.3          the pattern is WEAK/CHANGED -- filtering hits some
#                        pairs harder than others
# A low CV of that correlation distribution means the answer is stable across
# replicates; the empirical proportions make the call without a p-value.
#
# It also reports per-pair effect sizes ((filtered - null mean) / null SD), a
# permutation test for whether those effect sizes are heterogeneous across
# pairs, and -- when two or more groups carry a correlation distribution -- a
# cross-population test of whether pattern preservation itself differs between
# populations.
#
# Unlike analyze_pipeline_differences this stage does NOT require replicates in
# the filtered arm. calculate_filtered_hill performs no resampling, so the
# filtered arm is a single point estimate per pair; that is exactly the intended
# input here. The replicates come from the null arm. If the null arm carries no
# replicate column either, the script falls back to the reference's permutation
# path (see the log line "No replicate information"), which still yields a
# p-value but no correlation distribution.
#
# The reference implementation is vendored verbatim under the marker, so the
# file is self-contained and can be invoked by bare name from any task work
# directory.
# =============================================================================

doc <- "
Do the pipeline differences agree between the filtered and the null analysis?

Usage:
  pipecomp.R --filtered <file> --null <file> --output <prefix> [options]
  pipecomp.R (-h | --help)

Options:
  --filtered <file>       Filtered-arm rds from hill.R --mode filtered. One
                          point estimate per pipeline pair, which is expected.
  --null <file>           Null-arm rds from hill.R --mode null. Supplies the
                          Monte Carlo replicates the comparison runs on.
  --output <prefix>       Base name for the output files. The level and then the
                          table name are appended, e.g.
                          <prefix>_filtered_vs_null_total_pattern_summary.csv
  --outdir <dir>          Directory for the output files [default: .].
  --n-permutations <n>    Permutations for the fallback tests, which run only
                          when the null arm carries no replicates [default: 9999].
  --seed <n>              Random seed, set before each level so a rerun gives
                          the same p-values. main.R never seeds this stage, so
                          without this the permutation p-values change on every
                          run. It affects only the permutation draws.
                          [default: 1234]
  --pop                   Also compare the per-population split. Only runs when
                          BOTH arms carry PairwiseSplitByPopulation.
  --superpop              Also compare the per-superpopulation split. Only runs
                          when BOTH arms carry PairwiseSplitBySuperPopulation.
  --q <num>               Hill number order, reported in the log only. It does
                          not change the tests. [default: 1]
  -h --help               Show this help.
"

suppressPackageStartupMessages({
  library(docopt)
})

opts <- docopt(doc, args = commandArgs(trailingOnly = TRUE))

filteredFile   <- opts$filtered
nullFile       <- opts$null
outPrefix      <- opts$output
outdir         <- opts$outdir
nPermutations  <- as.integer(opts$n_permutations)
seed           <- as.integer(opts$seed)
wantPop        <- isTRUE(opts$pop)
wantSuperpop   <- isTRUE(opts$superpop)
qHillNumber    <- as.numeric(opts$q)

# Fail loudly rather than propagating integer(0)/NA into the reference code.
for (nm in c("nPermutations", "seed", "qHillNumber")) {
  if (length(get(nm)) != 1L || is.na(get(nm))) stop("Could not read ", nm, " from the CLI options.")
}
if (nPermutations < 1) stop("--n-permutations must be at least 1.")
assert_single <- function(x, what) {
  if (length(x) != 1L || is.na(x) || !nzchar(x)) stop("Could not read ", what, " from the CLI options.")
  x
}
outPrefix <- assert_single(outPrefix, "--output")
outdir    <- assert_single(outdir, "--outdir")
for (f in c(filteredFile, nullFile)) {
  if (length(f) != 1L || is.na(f) || !nzchar(f)) stop("Missing input file argument.")
  if (!file.exists(f)) stop("Input file not found: ", f)
}
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Reading the arms
# -----------------------------------------------------------------------------
# Both arms are plain hill.R output, so readRDS is enough; a wrong file should
# be named as such rather than failing later inside the vendored code.
read_arm <- function(path, label) {
  arm <- tryCatch(readRDS(path),
                  error = function(e) stop("Could not read the ", label,
                                           " arm '", path, "': ", conditionMessage(e)))
  if (!is.list(arm) || is.data.frame(arm)) {
    stop("The ", label, " arm '", path, "' is a ", class(arm)[1],
         ", not a hill.R result list. Expected the rds written by hill.R.")
  }
  if (is.null(arm$PairwiseTotalData)) {
    stop("The ", label, " arm '", path, "' has no PairwiseTotalData. ",
         "Was hill.R run with pairwise output, i.e. --mode null or --mode filtered? ",
         "Present names: ", paste(names(arm), collapse = ", "))
  }
  arm
}

# Does this arm actually hold replicates? Decides which reference path runs, and
# is the single most useful thing to print when a result looks surprising.
#
# This mirrors the reference condition exactly: compare_pipeline_patterns takes
# the correlation-distribution path only when a `replicate` column is present AND
# there are more rows than unique pairs. Repeated rows WITHOUT that column (as in
# a hand-built frame) take the permutation fallback instead, so detecting on row
# counts alone would report replicates that the reference will not use.
arm_has_replicates <- function(arm, field) {
  d <- arm[[field]]
  if (is.null(d)) return(NA)
  df <- convert_pairwise_to_df(d)
  if (is.null(df) || !is.data.frame(df)) return(NA)
  if (!"replicate" %in% colnames(df)) return(FALSE)
  key <- paste(pmin(df$site1, df$site2), pmax(df$site1, df$site2), sep = "-")
  nrow(df) > length(unique(key))
}

# One line per group, so a wrong shape is visible before anything is computed.
report_level <- function(label, filtered_pairwise, null_pairwise) {
  groups <- intersect(names(filtered_pairwise), names(null_pairwise))
  message("  ", label, ": ", length(groups), " shared group(s)")
  for (g in groups) {
    fd <- convert_pairwise_to_df(filtered_pairwise[[g]])
    nd <- convert_pairwise_to_df(null_pairwise[[g]])
    f_pairs <- if (is.data.frame(fd)) length(unique(paste(fd$site1, fd$site2, sep = "-"))) else 0L
    n_pairs <- if (is.data.frame(nd)) length(unique(paste(nd$site1, nd$site2, sep = "-"))) else 0L
    n_reps  <- if (is.data.frame(nd) && n_pairs > 0) round(nrow(nd) / n_pairs) else 0L
    message("    ", g, ": pairs=", n_pairs,
            " filtered_obs/pair=", if (is.data.frame(fd) && f_pairs > 0) round(nrow(fd) / f_pairs) else 0L,
            " null_reps/pair=", n_reps)
  }
  groups
}

# Compare one level and write its CSVs. analyse = FALSE means the level is
# skipped, and this reports which side was missing rather than silently doing
# nothing.
run_level <- function(level_label, filtered_pairwise, null_pairwise, file_stem, n_perm) {
  message("\n=== Level: ", level_label, " ===")
  if (is.null(filtered_pairwise) || is.null(null_pairwise)) return(NULL)

  groups <- report_level(level_label, filtered_pairwise, null_pairwise)
  if (length(groups) == 0) {
    message("  skipped: the two arms share no group names at this level")
    return(NULL)
  }
  if (length(groups) == 1 && level_label != "total") {
    message("  note: a cross-population test needs 2+ groups, so it will not run")
  }

  # Seed per level so each level's p-values are reproducible, and so adding or
  # removing a level does not shift the permutations of the others.
  set.seed(seed + match(level_label, c("total", "by_population", "by_superpopulation")))

  results <- compare_pipeline_patterns(
    filtered_pairwise = filtered_pairwise,
    null_pairwise     = null_pairwise,
    save_to_file      = file.path(outdir, paste0(file_stem, ".csv")),
    n_permutations    = n_perm
  )

  if (is.null(results)) {
    message("  ", level_label, ": no comparison produced (see the messages above).")
  } else if (!is.null(results$summary)) {
    message("  ", level_label, ": ", nrow(results$summary), " group(s) summarised.")
  }
  results
}

# The banner is printed in the footer, after the vendored functions are defined,
# because arm_has_replicates() calls convert_pairwise_to_df() and R evaluates
# message() arguments in the order the statements run.

#' Convert hillR pairwise format to data frame for analysis
#' Handles both single replicate and multiple replicate (null distribution) formats
#' @param pairwise_data List with TD_beta, site1, site2 components OR list of such lists (replicates)
#' @return Data frame with TD_beta column and pipeline pair rownames
convert_pairwise_to_df <- function(pairwise_data) {
  if(is.null(pairwise_data)) return(NULL)
  
  # If already a data frame, return as-is
  if(is.data.frame(pairwise_data)) {
    return(pairwise_data)
  }
  
  # If it's a single hillR result with TD_beta, site1, site2 components
  if(is.list(pairwise_data) && all(c("TD_beta", "site1", "site2") %in% names(pairwise_data))) {
    df <- data.frame(
      TD_beta = as.numeric(pairwise_data$TD_beta),
      site1 = as.character(pairwise_data$site1),
      site2 = as.character(pairwise_data$site2),
      stringsAsFactors = FALSE
    )
    # Add pair_name column for consistent access
    df$pair_name <- paste(df$site1, df$site2, sep = "-")
    return(df)
  }
  
  # If it's a list of replicates (null distribution format)
  # Each replicate has TD_beta, site1, site2
  if(is.list(pairwise_data) && length(pairwise_data) > 0) {
    # Check if first non-NULL element has the hillR structure
    first_element <- pairwise_data[[which(!sapply(pairwise_data, is.null))[1]]]
    
    if(is.list(first_element) && all(c("TD_beta", "site1", "site2") %in% names(first_element))) {
      # This is a list of hillR results - aggregate across replicates
      all_beta <- c()
      all_site1 <- c()
      all_site2 <- c()
      all_replicate <- c()  # Track replicate ID
      
      for(i in seq_along(pairwise_data)) {
        replicate <- pairwise_data[[i]]
        if(!is.null(replicate) && all(c("TD_beta", "site1", "site2") %in% names(replicate))) {
          all_beta <- c(all_beta, as.numeric(replicate$TD_beta))
          all_site1 <- c(all_site1, as.character(replicate$site1))
          all_site2 <- c(all_site2, as.character(replicate$site2))
          all_replicate <- c(all_replicate, rep(i, length(replicate$TD_beta)))  # Add replicate ID
        }
      }
      
      if(length(all_beta) > 0) {
        df <- data.frame(
          TD_beta = all_beta,
          site1 = all_site1,
          site2 = all_site2,
          replicate = all_replicate,  # Include replicate column
          stringsAsFactors = FALSE
        )
        # Add pair_name column for easy identification (don't use rownames with replicates)
        df$pair_name <- paste(df$site1, df$site2, sep = "-")
        return(df)
      }
    }
  }
  
  # Otherwise return NULL
  return(NULL)
}

#' Test for heterogeneity in correlation patterns across populations/groups
#' 
#' Tests whether the filtered-vs-null pattern preservation differs across populations.
#' Uses the summary statistics from compare_pipeline_patterns() to test if populations
#' show different levels of pattern preservation.
#' 
#' @param summary_df Summary data frame from compare_pipeline_patterns with correlation metrics
#' @return List with heterogeneity test results and summary statistics
test_cross_population_heterogeneity <- function(summary_df) {
  if(is.null(summary_df) || nrow(summary_df) < 2) {
    return(list(
      test_performed = FALSE,
      note = "Need at least 2 groups for heterogeneity test"
    ))
  }
  
  # Only test if we have correlation distributions (not permutation-based single values)
  groups_with_dist <- summary_df[summary_df$Has_Correlation_Distribution == TRUE, ]
  
  if(nrow(groups_with_dist) < 2) {
    return(list(
      test_performed = FALSE,
      note = "Need at least 2 groups with correlation distributions"
    ))
  }
  
  message(sprintf("\n=== Testing Cross-Population Heterogeneity ==="))
  message(sprintf("Testing if pattern preservation differs across %d populations...", nrow(groups_with_dist)))
  
  # Extract correlation metrics
  spearman_means <- groups_with_dist$Spearman_Mean
  pearson_means <- groups_with_dist$Pearson_Mean
  
  # Calculate variance and range
  spearman_var <- var(spearman_means, na.rm = TRUE)
  spearman_range <- range(spearman_means, na.rm = TRUE)
  spearman_sd <- sd(spearman_means, na.rm = TRUE)
  
  pearson_var <- var(pearson_means, na.rm = TRUE)
  pearson_range <- range(pearson_means, na.rm = TRUE)
  pearson_sd <- sd(pearson_means, na.rm = TRUE)
  
  # Test for differences using Kruskal-Wallis (non-parametric)
  # Only if there's variance to test
  if(spearman_var > 0) {
    # Can't do Kruskal-Wallis with single values per group
    # Instead report descriptive statistics
    message(sprintf("  Spearman correlation across populations:"))
    message(sprintf("    Mean: %.4f", mean(spearman_means, na.rm = TRUE)))
    message(sprintf("    SD: %.4f", spearman_sd))
    message(sprintf("    Range: [%.4f, %.4f]", spearman_range[1], spearman_range[2]))
    message(sprintf("    Variance: %.6f", spearman_var))
    
    spearman_heterogeneity <- "Variance detected"
  } else {
    message(sprintf("  Spearman correlation across populations:"))
    message(sprintf("    All populations identical: %.4f", spearman_means[1]))
    message(sprintf("    Variance: 0.000000 (no heterogeneity)"))
    
    spearman_heterogeneity <- "No variance - perfectly consistent"
  }
  
  if(pearson_var > 0) {
    message(sprintf("  Pearson correlation across populations:"))
    message(sprintf("    Mean: %.4f", mean(pearson_means, na.rm = TRUE)))
    message(sprintf("    SD: %.4f", pearson_sd))
    message(sprintf("    Range: [%.4f, %.4f]", pearson_range[1], pearson_range[2]))
    message(sprintf("    Variance: %.6f", pearson_var))
    
    pearson_heterogeneity <- "Variance detected"
  } else {
    message(sprintf("  Pearson correlation across populations:"))
    message(sprintf("    All populations identical: %.4f", pearson_means[1]))
    message(sprintf("    Variance: 0.000000 (no heterogeneity)"))
    
    pearson_heterogeneity <- "No variance - perfectly consistent"
  }
  
  # Interpretation
  message(sprintf("  "))
  message(sprintf("  INTERPRETATION:"))
  if(spearman_var == 0 && pearson_var == 0) {
    message(sprintf("    ZERO heterogeneity → Pattern preservation is IDENTICAL across all populations"))
    message(sprintf("    This indicates filtering affects pipeline comparisons consistently regardless of ancestry"))
    interpretation <- "Perfect consistency across populations"
  } else if(spearman_var < 0.01 && pearson_var < 0.01) {
    message(sprintf("    MINIMAL heterogeneity → Pattern preservation is highly consistent across populations"))
    message(sprintf("    Minor differences are likely due to sampling variation"))
    interpretation <- "Highly consistent across populations"
  } else {
    message(sprintf("    DETECTABLE heterogeneity → Pattern preservation varies across populations"))
    message(sprintf("    Some populations show stronger/weaker pattern preservation than others"))
    interpretation <- "Populations differ in pattern preservation"
  }
  
  return(list(
    test_performed = TRUE,
    n_groups = nrow(groups_with_dist),
    spearman = list(
      mean = mean(spearman_means, na.rm = TRUE),
      sd = spearman_sd,
      variance = spearman_var,
      range = spearman_range,
      heterogeneity = spearman_heterogeneity
    ),
    pearson = list(
      mean = mean(pearson_means, na.rm = TRUE),
      sd = pearson_sd,
      variance = pearson_var,
      range = pearson_range,
      heterogeneity = pearson_heterogeneity
    ),
    interpretation = interpretation
  ))
}

compare_pipeline_patterns <- function(filtered_pairwise, null_pairwise, 
                                     save_to_file = NULL, n_permutations = 9999) {
  if(is.null(filtered_pairwise) || is.null(null_pairwise)) {
    message("Cannot compare patterns: missing data")
    return(NULL)
  }
  
  # NEW APPROACH: Calculate correlation for EACH null replicate
  # Example: With 3 pipeline pairs and 100 null replicates:
  # Filtered:      [1.5, 2.0, 1.8] - one value per pair (length 3)
  # Null rep 1:    [1.3, 1.8, 1.6] - one value per pair (length 3) → cor1
  # Null rep 2:    [1.4, 1.9, 1.7] - one value per pair (length 3) → cor2
  # ...
  # Null rep 100:  [1.6, 2.1, 1.9] - one value per pair (length 3) → cor100
  # Result: Distribution of 100 correlations [cor1, cor2, ..., cor100]
  # If CV of correlations is low → pattern is stable
  
  # Get common groups
  common_groups <- intersect(names(filtered_pairwise), names(null_pairwise))
  
  if(length(common_groups) == 0) {
    warning("No common groups between filtered and null data")
    return(NULL)
  }
  
  message(sprintf("\n=== Comparing Pipeline Difference Patterns (Filtered vs Null) ==="))
  message(sprintf("Testing if filtering differentially affects pipeline comparisons..."))
  message(sprintf("Analyzing %d group(s)\n", length(common_groups)))
  
  all_results <- list()
  
  for(group_name in common_groups) {
    filt_data <- convert_pairwise_to_df(filtered_pairwise[[group_name]])
    null_data <- convert_pairwise_to_df(null_pairwise[[group_name]])
    
    if(is.null(filt_data) || is.null(null_data)) {
      message(sprintf("  %s: Skipping (data conversion failed)", group_name))
      next
    }
    
    # Get unique pipeline pairs
    if("pair_name" %in% colnames(filt_data)) {
      pair_names <- unique(filt_data$pair_name)
    } else {
      # Fallback: construct from site1/site2
      pair_names <- unique(paste(filt_data$site1, filt_data$site2, sep = "-"))
    }
    
    if(length(pair_names) < 2) {
      message(sprintf("  %s: Skipping (need at least 2 pipeline pairs)", group_name))
      next
    }
    
    # Check data structure
    # Count replicates per pair
    if("pair_name" %in% colnames(filt_data)) {
      filt_counts <- table(filt_data$pair_name)
      n_filt_per_pair <- as.numeric(filt_counts[1])
    } else {
      filt_counts <- table(paste(filt_data$site1, filt_data$site2, sep = "-"))
      n_filt_per_pair <- as.numeric(filt_counts[1])
    }
    
    # More robust: check if we have multiple replicates by using pair_name or site1/site2 columns
    if("pair_name" %in% colnames(null_data)) {
      n_unique_null_pairs <- length(unique(null_data$pair_name))
    } else if("site1" %in% colnames(null_data) && "site2" %in% colnames(null_data)) {
      null_pair_combos <- paste(pmin(null_data$site1, null_data$site2), 
                                 pmax(null_data$site1, null_data$site2), sep="-")
      n_unique_null_pairs <- length(unique(null_pair_combos))
    } else {
      # Last fallback to rownames (shouldn't happen with updated convert_pairwise_to_df)
      n_unique_null_pairs <- length(unique(rownames(null_data)))
    }
    
    n_total_null_rows <- nrow(null_data)
    has_null_replicates <- n_total_null_rows > n_unique_null_pairs
    n_null_per_pair <- if(has_null_replicates) round(n_total_null_rows / n_unique_null_pairs) else 1
    
    message(sprintf("  %s: %d pairs, filtered=%d obs/pair, null=%d obs/pair", 
                   group_name, length(pair_names), n_filt_per_pair, n_null_per_pair))
    message(sprintf("    Pair names (first 3): %s", paste(head(pair_names, 3), collapse=", ")))
    
    # Diagnostic: Check for replicate column
    message(sprintf("    Null data columns: %s", paste(colnames(null_data), collapse=", ")))
    message(sprintf("    Has 'replicate' column: %s", "replicate" %in% colnames(null_data)))
    if("replicate" %in% colnames(null_data)) {
      message(sprintf("    Unique replicates: %d", length(unique(null_data$replicate))))
    }
    message(sprintf("    Total null rows: %d, Unique pairs: %d, Replicates detected: %s", 
                   n_total_null_rows, n_unique_null_pairs, has_null_replicates))
    
    # Extract filtered beta values (point estimates for each pair)
    # Creates vector: filt_betas = [beta_pair1, beta_pair2, beta_pair3, ...]
    filt_betas <- sapply(pair_names, function(pair) {
      if("pair_name" %in% colnames(filt_data)) {
        vals <- as.numeric(filt_data[filt_data$pair_name == pair, "TD_beta",drop=TRUE])
      } else {
        matched <- paste(filt_data$site1, filt_data$site2, sep = "-") == pair
        vals <- as.numeric(filt_data[matched, "TD_beta", drop=TRUE])
      }
      if(length(vals) == 0) return(NA_real_)
      mean(vals, na.rm = TRUE)
    }, USE.NAMES = TRUE)
    
    # Extract null distribution summary statistics
    # FIX: Add drop=TRUE
    null_betas_mean <- sapply(pair_names, function(pair) {
      if("pair_name" %in% colnames(null_data)) {
        vals <- null_data[null_data$pair_name == pair, "TD_beta", drop=TRUE]
      } else {
        vals <- null_data[paste(null_data$site1, null_data$site2, sep = "-") == pair, "TD_beta", drop=TRUE]
      }
      if(length(vals) == 0) return(NA_real_)
      mean(vals, na.rm = TRUE)
    })
    
    # FIX: Add drop=TRUE
    null_betas_sd <- sapply(pair_names, function(pair) {
      if("pair_name" %in% colnames(null_data)) {
        vals <- null_data[null_data$pair_name == pair, "TD_beta", drop=TRUE]
      } else {
        vals <- null_data[paste(null_data$site1, null_data$site2, sep = "-") == pair, "TD_beta", drop=TRUE]
      }
      if(length(vals) == 0) return(NA_real_)
      sd(vals, na.rm = TRUE)
    })
    
    # Better approach: Calculate correlation for EACH null replicate
    if("replicate" %in% colnames(null_data) && has_null_replicates) {
      replicate_ids <- unique(null_data$replicate)
      n_replicates <- length(replicate_ids)
      
      message(sprintf("    Computing correlation distribution across %d null replicates...", n_replicates))
      
      # FIX: Add drop=TRUE in replicate loops
      null_cors_pearson <- sapply(replicate_ids, function(rep_id) {
        null_rep_betas <- sapply(pair_names, function(pair) {
          if("pair_name" %in% colnames(null_data)) {
            vals <- null_data[null_data$pair_name == pair & null_data$replicate == rep_id, "TD_beta", drop=TRUE]
          } else {
            pair_match <- paste(null_data$site1, null_data$site2, sep = "-") == pair
            vals <- null_data[pair_match & null_data$replicate == rep_id, "TD_beta", drop=TRUE]
          }
          if(length(vals) == 0) return(NA_real_)
          vals[1]  # Take first value (should only be one per replicate per pair)
        })
        if(all(is.na(null_rep_betas)) || all(is.na(filt_betas))) return(NA_real_)
        cor(filt_betas, null_rep_betas, method = "pearson", use = "complete.obs")
      })
      null_cors_pearson <- null_cors_pearson[!is.na(null_cors_pearson)]
      
      # FIX: Add drop=TRUE
      null_cors_spearman <- sapply(replicate_ids, function(rep_id) {
        null_rep_betas <- sapply(pair_names, function(pair) {
          if("pair_name" %in% colnames(null_data)) {
            vals <- null_data[null_data$pair_name == pair & null_data$replicate == rep_id, "TD_beta", drop=TRUE]
          } else {
            pair_match <- paste(null_data$site1, null_data$site2, sep = "-") == pair
            vals <- null_data[pair_match & null_data$replicate == rep_id, "TD_beta", drop=TRUE]
          }
          if(length(vals) == 0) return(NA_real_)
          vals[1]
        })
        if(all(is.na(null_rep_betas)) || all(is.na(filt_betas))) return(NA_real_)
        cor(filt_betas, null_rep_betas, method = "spearman", use = "complete.obs")
      })
      null_cors_spearman <- null_cors_spearman[!is.na(null_cors_spearman)]
      
      # Summary of correlation distributions
      obs_cor_pearson <- mean(null_cors_pearson, na.rm = TRUE)
      obs_cor_pearson_median <- median(null_cors_pearson, na.rm = TRUE)
      obs_cor_pearson_sd <- sd(null_cors_pearson, na.rm = TRUE)
      
      obs_cor_spearman <- mean(null_cors_spearman, na.rm = TRUE)
      obs_cor_spearman_median <- median(null_cors_spearman, na.rm = TRUE)
      obs_cor_spearman_sd <- sd(null_cors_spearman, na.rm = TRUE)
      
      # Test: Is variability in correlations low (pattern stable) or high (pattern unstable)?
      # Use CV (coefficient of variation) as measure of stability
      cv_pearson <- obs_cor_pearson_sd / abs(obs_cor_pearson)
      cv_spearman <- obs_cor_spearman_sd / abs(obs_cor_spearman)
      
      # Empirical proportions - what fraction of replicates show each pattern?
      # These ARE the empirical p-values - no additional test needed!
      prop_high_positive <- sum(null_cors_spearman >= 0.7) / length(null_cors_spearman)
      prop_high_negative <- sum(null_cors_spearman <= -0.7) / length(null_cors_spearman)
      prop_weak <- sum(abs(null_cors_spearman) <= 0.3) / length(null_cors_spearman)
      
      # Store for output (no binomial tests needed)
      p_value_pearson <- NA
      p_value_spearman <- NA
      
      message(sprintf("    Pearson: mean=%.3f, median=%.3f, SD=%.3f, CV=%.3f", 
                     obs_cor_pearson, obs_cor_pearson_median, obs_cor_pearson_sd, cv_pearson))
      message(sprintf("    Spearman: mean=%.3f, median=%.3f, SD=%.3f, CV=%.3f", 
                     obs_cor_spearman, obs_cor_spearman_median, obs_cor_spearman_sd, cv_spearman))
      
      # Store correlation distributions for output
      has_cor_distribution <- TRUE
      
    } else {
      # Fallback: If no replicate column, use old method (filtered vs null mean)
      message(sprintf("    No replicate information - using mean null values"))
      
      obs_cor_pearson <- cor(filt_betas, null_betas_mean, method = "pearson")
      obs_cor_pearson_median <- obs_cor_pearson
      obs_cor_pearson_sd <- NA
      
      obs_cor_spearman <- cor(filt_betas, null_betas_mean, method = "spearman")
      obs_cor_spearman_median <- obs_cor_spearman
      obs_cor_spearman_sd <- NA
      
      cv_pearson <- NA
      cv_spearman <- NA
      
      # Permutation test
      perm_cors_pearson <- replicate(n_permutations, {
        cor(sample(filt_betas), null_betas_mean, method = "pearson")
      })
      p_value_pearson <- (sum(abs(perm_cors_pearson) >= abs(obs_cor_pearson)) + 1) / (n_permutations + 1)
      
      perm_cors_spearman <- replicate(n_permutations, {
        cor(sample(filt_betas), null_betas_mean, method = "spearman")
      })
      p_value_spearman <- (sum(abs(perm_cors_spearman) >= abs(obs_cor_spearman)) + 1) / (n_permutations + 1)
      
      null_cors_pearson <- NULL
      null_cors_spearman <- NULL
      has_cor_distribution <- FALSE
      
      # Set NA for empirical proportions (not available without distribution)
      prop_high_positive <- NA
      prop_high_negative <- NA
      prop_weak <- NA
    }
    
    # Test 3: Heterogeneity of effect sizes ACROSS pipeline pairs
    # Calculate effect size (Cohen's d) for each pair: (filtered - null_mean) / null_sd
    # This creates a vector of effect sizes (one per pair)
    # Test: Do some pairs show larger deviations from null than others?
    effect_sizes <- (filt_betas - null_betas_mean) / null_betas_sd
    effect_sizes[is.infinite(effect_sizes)] <- NA
    
    # Test if effect sizes vary significantly across pairs (F-test)
    if(length(effect_sizes[!is.na(effect_sizes)]) >= 2) {
      # Use variance of effect sizes as test statistic
      obs_var <- var(effect_sizes, na.rm = TRUE)
      
      # Permutation: shuffle pair labels to test if variance is significant
      perm_vars <- replicate(n_permutations, {
        var(sample(effect_sizes), na.rm = TRUE)
      })
      p_value_heterogeneity <- (sum(perm_vars >= obs_var) + 1) / (n_permutations + 1)
    } else {
      obs_var <- NA
      p_value_heterogeneity <- NA
    }
    
    # Per-pair effect sizes and significance
    pair_details <- data.frame(
      Group = group_name,
      Pipeline_Pair = pair_names,
      Filtered_Beta = filt_betas,
      Null_Mean_Beta = null_betas_mean,
      Null_SD_Beta = null_betas_sd,
      Difference = filt_betas - null_betas_mean,
      Effect_Size = effect_sizes,
      stringsAsFactors = FALSE
    )
    
    # Summary for this group
    group_summary <- data.frame(
      Group = group_name,
      N_Pairs = length(pair_names),
      N_Null_Replicates = if(has_cor_distribution) length(null_cors_pearson) else NA,
      Pearson_Mean = obs_cor_pearson,
      Pearson_Median = obs_cor_pearson_median,
      Pearson_SD = obs_cor_pearson_sd,
      Pearson_CV = cv_pearson,
      Pearson_P = p_value_pearson,
      Spearman_Mean = obs_cor_spearman,
      Spearman_Median = obs_cor_spearman_median,
      Spearman_SD = obs_cor_spearman_sd,
      Spearman_CV = cv_spearman,
      Spearman_P = p_value_spearman,
      # Empirical proportions from correlation distribution
      Prop_High_Positive = if(has_cor_distribution) prop_high_positive else NA,
      Prop_High_Negative = if(has_cor_distribution) prop_high_negative else NA,
      Prop_Weak = if(has_cor_distribution) prop_weak else NA,
      # Other tests
      Effect_Var = obs_var,
      Heterogeneity_P = p_value_heterogeneity,
      Mean_Effect_Size = mean(effect_sizes, na.rm = TRUE),
      SD_Effect_Size = sd(effect_sizes, na.rm = TRUE),
      Has_Correlation_Distribution = has_cor_distribution,
      # Interpretation flags
      Pattern_Preserved = obs_cor_spearman > 0.7,  # High positive = pattern preserved
      Pattern_Opposite = obs_cor_spearman < -0.7,  # High negative = opposite pattern
      Pattern_Weak = abs(obs_cor_spearman) < 0.3,   # Weak = pattern changed/none
      Pattern_Stable_Low_CV = has_cor_distribution && !is.na(cv_spearman) && cv_spearman < 0.2,
      Differential_Effect = !is.na(p_value_heterogeneity) && p_value_heterogeneity < 0.05,
      stringsAsFactors = FALSE
    )
    
    # Interpret results based on whether we have correlation distribution
    if(has_cor_distribution) {
      message(sprintf("    Correlation distribution analysis (across %d replicates):", length(null_cors_pearson)))
      message(sprintf("      Pearson: mean=%.3f, median=%.3f, SD=%.3f, CV=%.3f", 
                     obs_cor_pearson, obs_cor_pearson_median, obs_cor_pearson_sd, cv_pearson))
      message(sprintf("      Spearman: mean=%.3f, median=%.3f, SD=%.3f, CV=%.3f", 
                     obs_cor_spearman, obs_cor_spearman_median, obs_cor_spearman_sd, cv_spearman))
      message(sprintf("    "))
      message(sprintf("    EMPIRICAL PROPORTIONS:"))
      message(sprintf("      High positive (≥0.7): %.1f%%", prop_high_positive * 100))
      message(sprintf("      High negative (≤-0.7): %.1f%%", prop_high_negative * 100))
      message(sprintf("      Weak (|r|≤0.3): %.1f%%", prop_weak * 100))
      message(sprintf("    "))
      message(sprintf("    INTERPRETATION:"))
      if(obs_cor_spearman > 0.7) {
        message(sprintf("      HIGH POSITIVE correlation (%.3f) → Pattern PRESERVED", obs_cor_spearman))
        message(sprintf("      Filtering affects all pipeline pairs similarly"))
        message(sprintf("      [%.1f%% of replicates show this pattern]", prop_high_positive * 100))
      } else if(obs_cor_spearman < -0.7) {
        message(sprintf("      HIGH NEGATIVE correlation (%.3f) → OPPOSITE Pattern", obs_cor_spearman))
        message(sprintf("      Filtering reverses the pattern (most divergent becomes least)"))
        message(sprintf("      [%.1f%% of replicates show this pattern]", prop_high_negative * 100))
      } else if(abs(obs_cor_spearman) < 0.3) {
        message(sprintf("      WEAK correlation (%.3f) → Pattern CHANGED/NONE", obs_cor_spearman))
        message(sprintf("      Filtering differentially affects pipeline pairs"))
        message(sprintf("      [%.1f%% of replicates show this pattern]", prop_weak * 100))
      } else {
        message(sprintf("      MODERATE correlation (%.3f) → Mixed pattern", obs_cor_spearman))
      }
      if(cv_spearman < 0.2) {
        message(sprintf("      Low CV (%.3f) → Stable correlation across replicates", cv_spearman))
      } else {
        message(sprintf("      High CV (%.3f) → Variable correlation across replicates", cv_spearman))
      }
    } else {
      message(sprintf("    Pearson r = %.3f (p = %.4f)", obs_cor_pearson, p_value_pearson))
      message(sprintf("    Spearman rho = %.3f (p = %.4f)", obs_cor_spearman, p_value_spearman))
      message(sprintf("    "))
      message(sprintf("    INTERPRETATION:"))
      if(obs_cor_spearman > 0.7) {
        message(sprintf("      HIGH POSITIVE → Pattern preserved"))
      } else if(obs_cor_spearman < -0.7) {
        message(sprintf("      HIGH NEGATIVE → Opposite pattern"))
      } else if(abs(obs_cor_spearman) < 0.3) {
        message(sprintf("      WEAK → Pattern changed/none"))
      } else {
        message(sprintf("      MODERATE → Mixed effects"))
      }
    }
    message(sprintf("    Effect heterogeneity: Var = %.3f (p = %s)", 
                   obs_var, 
                   if(is.na(p_value_heterogeneity)) "NA" else sprintf("%.4f", p_value_heterogeneity)))
    if(!is.na(p_value_heterogeneity) && p_value_heterogeneity < 0.05) {
      message(sprintf("      Significant heterogeneity → Some pairs show larger deviations from null"))
    }
    message(sprintf("    "))
    all_results[[group_name]] <- list(
      summary = group_summary,
      pair_details = pair_details,
      correlation_distributions = if(has_cor_distribution) {
        list(pearson = null_cors_pearson, spearman = null_cors_spearman)
      } else {
        NULL
      }
    )
  }
  
  if(length(all_results) == 0) {
    message("No valid comparisons performed")
    return(NULL)
  }
  
  # Combine summaries across all groups
  summary_df <- do.call(rbind, lapply(all_results, function(x) x$summary))
  pair_details_df <- do.call(rbind, lapply(all_results, function(x) x$pair_details))
  
  # Apply multiple testing correction ONLY if we have p-values (when no correlation distribution)
  if(!is.null(summary_df) && nrow(summary_df) > 0 && any(!is.na(summary_df$Pearson_P))) {
    # Correct Pearson p-values
    summary_df$Pearson_P_Bonferroni <- pmin(summary_df$Pearson_P * nrow(summary_df), 1.0)
    summary_df$Pearson_P_FDR <- p.adjust(summary_df$Pearson_P, method = "fdr")
    
    # Correct Spearman p-values
    summary_df$Spearman_P_Bonferroni <- pmin(summary_df$Spearman_P * nrow(summary_df), 1.0)
    summary_df$Spearman_P_FDR <- p.adjust(summary_df$Spearman_P, method = "fdr")
  }
  
  # Correct heterogeneity p-values
  if(!is.null(summary_df) && nrow(summary_df) > 0 && any(!is.na(summary_df$Heterogeneity_P))) {
    summary_df$Heterogeneity_P_Bonferroni <- pmin(summary_df$Heterogeneity_P * nrow(summary_df), 1.0)
    summary_df$Heterogeneity_P_FDR <- p.adjust(summary_df$Heterogeneity_P, method = "fdr")
  }
  
  message(sprintf("\n=== Pattern Comparison Summary ==="))
  message(sprintf("  Total groups analyzed: %d", nrow(summary_df)))
  message(sprintf("  "))
  
  # Count groups with correlation distributions
  n_with_dist <- sum(summary_df$Has_Correlation_Distribution, na.rm = TRUE)
  if(n_with_dist > 0) {
    message(sprintf("  Correlation distribution analysis:"))
    message(sprintf("    Groups analyzed: %d / %d", n_with_dist, nrow(summary_df)))
    message(sprintf("    Pattern PRESERVED (Spearman > 0.7): %d / %d", 
                   sum(summary_df$Pattern_Preserved, na.rm = TRUE), nrow(summary_df)))
    message(sprintf("    Pattern OPPOSITE (Spearman < -0.7): %d / %d", 
                   sum(summary_df$Pattern_Opposite, na.rm = TRUE), nrow(summary_df)))
    message(sprintf("    Pattern WEAK/CHANGED (|Spearman| < 0.3): %d / %d", 
                   sum(summary_df$Pattern_Weak, na.rm = TRUE), nrow(summary_df)))
    message(sprintf("    Stable correlation (CV < 0.2): %d / %d", 
                   sum(summary_df$Pattern_Stable_Low_CV, na.rm = TRUE), n_with_dist))
    message(sprintf("  "))
  }
  
  message(sprintf("  Effect sizes:"))
  message(sprintf("    Differential effects (heterogeneity p<0.05): %d / %d", 
                 sum(summary_df$Differential_Effect, na.rm = TRUE), nrow(summary_df)))
  message(sprintf("  "))
  message(sprintf("  INTERPRETATION GUIDE:"))
  message(sprintf("    HIGH POSITIVE correlation (>0.7) → Filtering affects all pipeline pairs similarly"))
  message(sprintf("    HIGH NEGATIVE correlation (<-0.7) → Filtering produces OPPOSITE pattern"))
  message(sprintf("    WEAK correlation (|r|<0.3) → Filtering differentially affects certain pairs"))
  message(sprintf("    Low CV → Correlation is stable/robust across null replicates"))
  message(sprintf("    Heterogeneity → Some pairs deviate more from null than others"))
  
  # Save results if requested
  if(!is.null(save_to_file)) {
    summary_file <- sub("\\.csv$", "_pattern_summary.csv", save_to_file)
    detail_file <- sub("\\.csv$", "_pattern_details.csv", save_to_file)
    
    write.csv(summary_df, summary_file, row.names = FALSE)
    write.csv(pair_details_df, detail_file, row.names = FALSE)
    
    message(sprintf("  Summary saved to: %s", summary_file))
    message(sprintf("  Detailed pair comparisons saved to: %s", detail_file))
    
    # Save correlation distributions if available
    n_with_dist <- sum(summary_df$Has_Correlation_Distribution, na.rm = TRUE)
    if(n_with_dist > 0) {
      cor_dist_file <- sub("\\.csv$", "_correlation_distributions.csv", save_to_file)
      cor_dist_list <- lapply(names(all_results), function(group_name) {
        cor_dists <- all_results[[group_name]]$correlation_distributions
        if(!is.null(cor_dists)) {
          data.frame(
            Group = group_name,
            Replicate = seq_along(cor_dists$pearson),
            Pearson_r = cor_dists$pearson,
            Spearman_rho = cor_dists$spearman,
            stringsAsFactors = FALSE
          )
        } else {
          NULL
        }
      })
      cor_dist_df <- do.call(rbind, cor_dist_list[!sapply(cor_dist_list, is.null)])
      if(!is.null(cor_dist_df) && nrow(cor_dist_df) > 0) {
        write.csv(cor_dist_df, cor_dist_file, row.names = FALSE)
        message(sprintf("  Correlation distributions saved to: %s", cor_dist_file))
      }
    }
  }
  
  # Test for cross-population heterogeneity in correlation patterns
  cross_pop_test <- test_cross_population_heterogeneity(summary_df)
  
  return(list(
    summary = summary_df,
    pair_details = pair_details_df,
    detailed_results = all_results,
    cross_population_test = cross_pop_test
  ))
}

# -----------------------------------------------------------------------------
# Levels and dispatch
# -----------------------------------------------------------------------------
# main.R wraps each level's data as list(AllData = ...) for the total level and
# passes the split lists through as-is. compare_pipeline_patterns then intersects
# the two arms' names, so the wrapping has to match on both sides for the group
# to line up. Both arms come from hill.R, which uses the same names, so a level
# is compared whenever both arms carry it.
run_started <- Sys.time()

filteredArm <- read_arm(filteredFile, "filtered")
nullArm     <- read_arm(nullFile, "null")

# Printed here rather than at the top of the script: arm_has_replicates() calls
# convert_pairwise_to_df(), which the vendored block above defines, and R
# evaluates these arguments when the statement runs.
message("=== pipecomp.R ===")
message("  filtered: ", filteredFile)
message("  null    : ", nullFile)
message("  output  : ", file.path(outdir, paste0(outPrefix, "_filtered_vs_null_*")))
message("  permutations (fallback path only): ", nPermutations)
message("  seed    : ", seed, " (permutation draws only)")
message("  q       : ", qHillNumber)
message("  replicates: filtered=", arm_has_replicates(filteredArm, "PairwiseTotalData"),
        " null=", arm_has_replicates(nullArm, "PairwiseTotalData"))
if (identical(arm_has_replicates(nullArm, "PairwiseTotalData"), FALSE)) {
  message("  NOTE: the null arm carries no replicates, so the reference falls back ")
  message("        to its permutation path: a p-value but no correlation distribution.")
}

results <- list()

# Guard the degenerate case: if both arms are the same file, the comparison is a
# tautology (every replicate correlates perfectly with itself) and would report
# "pattern PRESERVED" that means nothing. That is what happens when one count
# table stands in for both arms, so say it plainly instead of reporting it as a
# finding.
if (identical(normalizePath(filteredFile, mustWork = FALSE),
              normalizePath(nullFile, mustWork = FALSE))) {
  message("\n!!! --filtered and --null are the SAME file.")
  message("    The two arms are identical, so every correlation will be ~1 and the")
  message("    result is a tautology, not a finding. Pass a real filtered arm, i.e.")
  message("    the output of hill.R --mode filtered on a genuinely filtered count table.")
}

# Total data is always compared when both arms have it.
if (!is.null(filteredArm$PairwiseTotalData) && !is.null(nullArm$PairwiseTotalData)) {
  results$total_pattern <- run_level(
    "total",
    list(AllData = filteredArm$PairwiseTotalData),
    list(AllData = nullArm$PairwiseTotalData),
    paste0(outPrefix, "_filtered_vs_null_total"),
    nPermutations)
} else {
  message("\n=== Level: total ===")
  message("  skipped: ", if (is.null(filteredArm$PairwiseTotalData)) "filtered arm has no PairwiseTotalData" else "null arm has no PairwiseTotalData")
}

# Splits run only when asked for AND present in BOTH arms. Compare the hill.R
# arm fields directly, since the reference splits on those names.
#
# Presence is NOT the same as usability. hill.R builds both split lists
# unconditionally, so an unusable level still arrives as a NAMED list whose
# entries are empty shells (and when the metadata has no superpopulation column
# those names are blank). is.null() therefore cannot tell an absent level from an
# empty one, so count the groups that actually carry rows.
usable_groups <- function(x) {
  if (is.null(x)) return(0L)
  if (is.data.frame(x)) return(if (nrow(x) > 0) 1L else 0L)
  sum(vapply(x, function(e) is.data.frame(e) && nrow(e) > 0, logical(1)))
}

compare_split <- function(level_label, field, file_suffix, enabled, flag) {
  message("\n=== Level: ", level_label, " ===")
  f <- filteredArm[[field]]
  n <- nullArm[[field]]
  if (is.null(f) || is.null(n)) {
    message("  skipped: ",
            if (is.null(f) && is.null(n)) "neither arm carries it"
            else if (is.null(f)) "filtered arm has no " else "null arm has no ", field)
    return(NULL)
  }
  # Present but empty: say which arm and how to fix it, rather than running the
  # reference on empty groups and reporting a meaningless result.
  nf <- usable_groups(f)
  nn <- usable_groups(n)
  if (nf == 0L || nn == 0L) {
    message("  skipped: ", field, " is present but empty (filtered: ", nf,
            " group(s) with rows, null: ", nn, "). Populate population_code / ",
            "super_population_code in the metadata to enable it.")
    return(NULL)
  }
  if (!enabled) {
    message("  skipped: ", flag, " not set")
    return(NULL)
  }
  run_level(level_label, f, n, paste0(outPrefix, "_filtered_vs_null_", file_suffix), nPermutations)
}

results$by_population_pattern <- compare_split(
  "by_population", "PairwiseSplitByPopulation", "by_population", wantPop, "--pop")
results$by_superpopulation_pattern <- compare_split(
  "by_superpopulation", "PairwiseSplitBySuperPopulation", "by_superpopulation", wantSuperpop, "--superpop")

# Drop the skipped levels so the metadata reports only what ran.
results <- results[!vapply(results, is.null, logical(1))]

# -----------------------------------------------------------------------------
# Metadata
# -----------------------------------------------------------------------------
# Mirrors hill.R's use of a metadata entry in the output, and records the seed
# so a p-value that looks odd can be reproduced.
CompMetadata <- list(
  script                = "pipecomp.R",
  filtered_file         = filteredFile,
  null_file             = nullFile,
  qHillNumber           = qHillNumber,
  nPermutations         = nPermutations,
  seed                  = seed,
  levelsAnalyzed        = names(results),
  evaluatePopulation    = wantPop,
  evaluateSuperPopulation = wantSuperpop,
  filteredHasReplicates = arm_has_replicates(filteredArm, "PairwiseTotalData"),
  nullHasReplicates     = arm_has_replicates(nullArm, "PairwiseTotalData"),
  runStarted            = format(run_started),
  runFinished           = format(Sys.time())
)
saveRDS(CompMetadata, file.path(outdir, paste0(outPrefix, "_filtered_vs_null_metadata.rds")))

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
message("\n=== pipecomp.R: done ===")
if (length(results) == 0) {
  message("  no level could be compared.")
  message("  Both arms must come from hill.R and carry the same level. The total ")
  message("  level needs PairwiseTotalData on both sides.")
} else {
  message("  levels compared: ", paste(names(results), collapse = ", "))
}

# Collect what this run wrote: same prefix AND written after we started. Matching
# on the prefix alone would also list a previous run's results, since reruns
# reuse the same file names.
saved <- list.files(outdir, pattern = paste0("^", outPrefix, "_filtered_vs_null.*\\.csv$"),
                    full.names = FALSE)
saved <- saved[file.mtime(file.path(outdir, saved)) >= run_started]
message("  CSV files written: ", length(saved))
for (f in sort(saved)) message("    ", f)
message("  metadata: ", file.path(outdir, paste0(outPrefix, "_filtered_vs_null_metadata.rds")))
