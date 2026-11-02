#!/usr/bin/env Rscript
# =============================================================================
# plots.R -- plot the Hill diversity distributions
# =============================================================================
#
# Scope
# -----
# This is the plotting stage of the diversity analysis, downstream of hill.R.
# It consumes the null arm (and optionally the filtered arm) produced by
# DIVANALYSIS_AGGREGATE and produces:
#
#  1. Density plots of the null Hill diversity distributions, overall, per
#     population, per superpopulation, for q = TD_gamma, TD_alpha, TD_beta;
#  2. Overlay plots: null density curves with the filtered arm drawn as a
#     vertical line (filtered is a point estimate from calculate_filtered_hill);
#  3. Pairwise beta-diversity heatmaps: filtered, null-mean, and difference
#     (filtered - null mean) for population and superpopulation splits;
#  4. Kruskal-Wallis tests of the null distributions across populations with
#     pairwise permutation follow-ups, written as CSVs.
#
# Statistical tests comparing the filtered arm to the null arm are NOT here.
# They live in pipecomp.R (correlation distributions, Spearman rho, effect
# sizes).  This script is pure visualization + null-only tests.
#
# Vendored reference code
# -----------------------
# The functions below are vendored byte-faithfully from ../R/plotting.R:
#   create_density_plots_internal
#   create_null_vs_filtered_comparison
#   create_histogram_plot
#   create_plot_annotations
#   calculate_two_sided_pvalue
#   has_valid_split
#   plot_pairwise_beta_heatmaps
#   compare_null_distributions_by_population
#
# See the marker below for the vendored block boundary.
#
# MultiQC compatibility
# ---------------------
# MultiQC discovers images by FILENAME SUFFIX: a file must end in '_mqc.png'
# to be picked up. CSVs are written as '<prefix>_nulltest_<level>_<metric>.csv'
# and '<prefix>_comparison.csv' for the filtered-vs-null p-value table.
#
# Usage
# -----
#   plots.R --null <file> --output <prefix> [options]
#   plots.R --null <file> --filtered <file> --output <prefix> [options]
#
# Options
# -------
#   --null=<file>             null-arm rds written by hill.R --mode null
#   --filtered=<file>         filtered-arm rds (optional; if omitted only
#                             null-only plots are produced)
#   --output=<prefix>         filename prefix for the outputs [default: hill_plots]
#   --outdir=<dir>            directory for the outputs [default: .]
#   --pop                     plot / test the by-population split
#   --superpop                plot / test the by-superpopulation split
#   --n-permutations=<int>    permutations for the pairwise follow-ups [default: 10000]
#   --seed=<int>              RNG seed for those permutations [default: 1234]
#   --q=<int>                 Hill order, carried through for reporting [default: 1]

suppressPackageStartupMessages({
  library(docopt)
  library(ggplot2)
  library(dplyr)
  library(pals)
})

"Usage:
  plots.R --null <file> --output <prefix> [options]
  plots.R --null <file> --filtered <file> --output <prefix> [options]
  plots.R --help

Options:
  -h --help                  Show this screen.
  --null <file>              Null-arm rds from hill.R --mode null.
  --filtered <file>          Filtered-arm rds (optional).
  --output <prefix>          Filename prefix for outputs [default: hill_plots].
  --outdir <dir>             Output directory [default: .].
  --pop                      Plot and test the by-population split.
  --superpop                 Plot and test the by-superpopulation split.
  --n-permutations <int>     Permutations for pairwise tests [default: 10000].
  --seed <int>               RNG seed for permutation tests [default: 1234].
  --q <int>                  Hill order, for reporting only [default: 1].
" -> doc

opts <- docopt(doc)

null_file   <- opts$null
filtered_file <- opts[["--filtered"]]
out_prefix  <- opts$output
outdir      <- opts$outdir
want_pop    <- isTRUE(opts$pop)
want_superpop <- isTRUE(opts$superpop)
n_permutations <- suppressWarnings(as.integer(opts$n_permutations))
seed <- suppressWarnings(as.integer(opts$seed))
q_hill <- suppressWarnings(as.integer(opts$q))

if (is.na(n_permutations) || n_permutations < 1) {
  stop("--n-permutations must be a positive integer")
}
if (is.na(seed)) {
  stop("--seed must be an integer")
}
if (!file.exists(null_file)) {
  stop("Null input file not found: ", null_file)
}
if (!grepl("\\.rds$", null_file, ignore.case = TRUE)) {
  stop("Null input must be an .rds file, got: ", null_file)
}
has_filtered <- isTRUE(!is.null(filtered_file) && file.exists(filtered_file) &&
                        grepl("\\.rds$", filtered_file, ignore.case = TRUE))
if (has_filtered && !grepl("\\.rds$", filtered_file, ignore.case = TRUE)) {
  stop("Filtered input must be an .rds file, got: ", filtered_file)
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Vendored reference code -- ../R/plotting.R
#
# Vendored byte-faithfully from the reference implementation.  See the file
# header for the provenance of each function.
# =============================================================================

# Reference: plotting.R line 17
calculate_two_sided_pvalue <- function(filtered_value, null_values, exclude_index = NULL) {
  if (is.null(exclude_index) || length(exclude_index) == 0) {
    p_value <- (sum(null_values >= filtered_value) + sum(null_values <= filtered_value) + 1) / (length(null_values) + 1)
  } else {
    null_values_excl <- null_values[-exclude_index]
    p_value <- (sum(null_values_excl >= filtered_value) + sum(null_values_excl <= filtered_value) + 1) / (length(null_values_excl) + 1)
  }
  return(p_value)
}

# Reference: plotting.R line 86
create_plot_annotations <- function(test_this_metric, p_value, p_value_bonferroni,
                                     pvalue_results, n_tests) {
  p_label <- ""
  caption_parts <- character(0)
  if (test_this_metric) {
    p_label <- sprintf("p = %.4f", p_value)
    if (!is.na(p_value_bonferroni)) {
      caption_parts <- c(caption_parts, sprintf("Bonferroni-corrected p = %.4f", p_value_bonferroni))
    }
  }
  return(list(p_label = p_label, caption = paste(caption_parts, collapse = " | ")))
}

# Reference: plotting.R line 125
create_histogram_plot <- function(null_values, filtered_df, metric, pop) {
  null_vals <- na.omit(null_values)
  if (length(null_vals) == 0) return(NULL)

  # Scale y-axis to match histogram density
  max_freq <- max(hist(null_vals, breaks = 30, plot = FALSE)$counts)

  p <- ggplot() +
    geom_histogram(data = data.frame(x = null_vals), aes(x = x, y = ..density..),
                   bins = 30, fill = "#0072B2", alpha = 0.5, color = "black") +
    geom_vline(data = filtered_df, aes(xintercept = !!sym(metric)),
               color = "#D55E00", linewidth = 1.2, linetype = "dashed") +
    scale_y_continuous(
      name = "Density",
      sec.axis = dup_axis(labels = paste(seq(0, max_freq, length.out = 5), collapse = ", "),
                          name = NULL)
    )

  return(p)
}

# Reference: plotting.R line 211
has_valid_split <- function(x) {
  if (is.null(x) || !is.list(x)) return(FALSE)
  if (length(x) == 0) return(FALSE)
  # Must have at least one element that is a data.frame with rows
  any(vapply(x, function(e) is.data.frame(e) && nrow(e) > 0, logical(1)))
}

# Reference: plotting.R line 488
create_density_plots_internal <- function(data_list, title_prefix, save_plot) {
  if (is.null(data_list) || length(data_list) == 0) {
    return(NULL)
  }

  # Remove metadata if present
  if (!is.null(data_list$metadata)) {
    data_list <- data_list[names(data_list) != "metadata"]
  }

  # Determine grouping column based on title_prefix
  if (grepl("SuperPopulation", title_prefix, ignore.case = TRUE)) {
    group_col <- "SuperPopulation"
    legend_label <- "SuperPopulation"
  } else {
    group_col <- "Population"
    legend_label <- "Population"
  }

  # Bind all data with correct .id
  combined_df <- bind_rows(data_list, .id = group_col)

  metrics <- c("TD_beta", "TD_alpha", "TD_gamma")
  metric_labels <- c("Beta Diversity", "Alpha Diversity", "Gamma Diversity")

  # Ensure grouping column is a factor with alphabetical levels for consistent ordering
  combined_df[[group_col]] <- factor(combined_df[[group_col]], levels = sort(unique(combined_df[[group_col]])))

  for (i in seq_along(metrics)) {
    metric <- metrics[i]
    if (!metric %in% colnames(combined_df)) next

    n_groups <- length(unique(combined_df[[group_col]]))
    colors <- unname(pals::glasbey(n_groups))
    names(colors) <- levels(combined_df[[group_col]])

    # Determine per-group replicate and constancy status for this metric
    stats_per_group <- combined_df %>%
      group_by(.data[[group_col]]) %>%
      summarize(n = n(), n_unique = length(unique(na.omit(!!sym(metric)))), .groups = "drop")

    groups_with_single <- stats_per_group[[group_col]][stats_per_group$n_unique <= 1]
    groups_with_replicates <- stats_per_group[[group_col]][stats_per_group$n_unique > 1]

    has_replicates <- length(groups_with_replicates) > 0
    has_single_values <- length(groups_with_single) > 0

    if (has_replicates && has_single_values) {
      # Mixed case - some groups have replicates, others are point-like
      df_replicates <- combined_df[combined_df[[group_col]] %in% groups_with_replicates, ]
      df_single <- combined_df[combined_df[[group_col]] %in% groups_with_single, ]

      all_groups <- levels(combined_df[[group_col]])
      plotted_groups <- unique(c(as.character(df_replicates[[group_col]]),
                                  as.character(df_single[[group_col]])))
      plotted_colors <- colors[plotted_groups]

      p <- ggplot() +
        geom_density(data = df_replicates,
                     aes(x = !!sym(metric), fill = .data[[group_col]], color = .data[[group_col]]),
                     alpha = 0.3) +
        geom_vline(data = df_single,
                  aes(xintercept = !!sym(metric), color = .data[[group_col]]),
                  linetype = "dashed", linewidth = 1, show.legend = TRUE) +
        scale_color_manual(
          values = plotted_colors,
          breaks = all_groups[all_groups %in% names(plotted_colors)],
          guide = guide_legend(
            override.aes = list(
              linetype = ifelse(all_groups[all_groups %in% names(plotted_colors)] %in% groups_with_single,
                                "dashed", "solid"),
              fill = ifelse(all_groups[all_groups %in% names(plotted_colors)] %in% groups_with_single,
                            NA, plotted_colors[all_groups[all_groups %in% names(plotted_colors)]]),
              alpha = ifelse(all_groups[all_groups %in% names(plotted_colors)] %in% groups_with_single,
                             1, 0.5)
            )
          )
        ) +
        scale_fill_manual(
          values = plotted_colors,
          breaks = all_groups[all_groups %in% names(plotted_colors)],
          guide = "none"
        ) +
        labs(title = paste0(title_prefix, ": ", metric_labels[i]),
             x = metric_labels[i], y = "Density", color = legend_label) +
        theme_bw() + theme(legend.position = "right")

    } else if (has_replicates) {
      p <- ggplot(combined_df,
                  aes(x = !!sym(metric), fill = .data[[group_col]], color = .data[[group_col]])) +
        geom_density(alpha = 0.3) +
        scale_color_manual(values = colors) +
        scale_fill_manual(values = colors) +
        labs(title = paste0(title_prefix, ": ", metric_labels[i]),
             x = metric_labels[i], y = "Density", color = legend_label) +
        theme_bw() + theme(legend.position = "right")

    } else {
      p <- ggplot(combined_df,
                  aes(x = .data[[group_col]], y = !!sym(metric), fill = .data[[group_col]])) +
        geom_col() +
        scale_fill_manual(values = colors) +
        labs(title = paste0(title_prefix, ": ", metric_labels[i]),
             x = legend_label, y = metric_labels[i]) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none")
    }

    plot_name <- paste0(title_prefix, "_", metric)
    save_plot(p, plot_name)
  }
}

# Reference: plotting.R line 617
create_comparison_plots_internal <- function(NullHillValues, FilteredHillValues,
                                              plotname, outputfile,
                                              downSampleByPop, downSampleBySuperPop) {
  plotname_base <- if (!is.null(plotname)) {
    sub("\\.pdf$", "", plotname)
  } else {
    sub("\\.RData$", "", basename(outputfile))
  }

  if (downSampleByPop) {
    message("\nCreating comparison plots for population splits...")
    create_null_vs_filtered_comparison(
      NullHillValues$SplitByPopulation, FilteredHillValues$SplitByPopulation,
      "by_population", save_plot,
      numSampleSets = if (is.list(NullHillValues$SplitByPopulation) &&
                              !is.null(NullHillValues$SplitByPopulation$metadata))
        NullHillValues$SplitByPopulation$metadata$numSampleSets else NULL,
      numVariantSets = if (is.list(NullHillValues$SplitByPopulation) &&
                              !is.null(NullHillValues$SplitByPopulation$metadata))
        NullHillValues$SplitByPopulation$metadata$numVariantSets else NULL
    )
  }

  if (downSampleBySuperPop) {
    message("\nCreating comparison plots for superpopulation splits...")
    create_null_vs_filtered_comparison(
      NullHillValues$SplitBySuperPopulation, FilteredHillValues$SplitBySuperPopulation,
      "by_superpopulation", save_plot,
      numSampleSets = if (is.list(NullHillValues$SplitBySuperPopulation) &&
                              !is.null(NullHillValues$SplitBySuperPopulation$metadata))
        NullHillValues$SplitBySuperPopulation$metadata$numSampleSets else NULL,
      numVariantSets = if (is.list(NullHillValues$SplitBySuperPopulation) &&
                              !is.null(NullHillValues$SplitBySuperPopulation$metadata))
        NullHillValues$SplitBySuperPopulation$metadata$numVariantSets else NULL
    )
  }
}

# Reference: plotting.R line 681
create_null_vs_filtered_comparison <- function(null_data, filtered_data, level, save_plot,
                                              numSampleSets = NULL, numVariantSets = NULL) {

  message(paste0("  Creating null vs filtered comparison plots for ", level))
  message("    Statistical tests performed only on beta diversity; alpha and gamma are for visualization")

  if (!is.null(null_data$metadata)) {
    null_data <- null_data[names(null_data) != "metadata"]
  }

  populations <- intersect(names(null_data), names(filtered_data))
  message(paste0("    Common populations: ", paste(populations, collapse = ", ")))

  if (length(populations) == 0) {
    message("    No common populations found - skipping comparison plots")
    return(NULL)
  }

  metrics <- c("TD_beta", "TD_alpha", "TD_gamma")
  metric_labels <- c("Beta Diversity", "Alpha Diversity", "Gamma Diversity")
  metrics_to_test <- c("TD_beta")
  n_tests <- length(populations) * length(metrics_to_test)
  message(paste0("    Total tests for Bonferroni correction: ", n_tests, " (beta diversity only)"))

  plots_created <- 0
  pvalue_results <- list()

  for (metric_idx in seq_along(metrics)) {
    metric <- metrics[metric_idx]

    for (pop in populations) {
      null_pop_data <- null_data[[pop]]
      filtered_pop_data <- filtered_data[[pop]]

      if (is.null(null_pop_data) || is.null(filtered_pop_data)) {
        message(paste0("    Skipping ", pop, " - null or missing data"))
        next
      }

      null_df <- tryCatch(as.data.frame(null_pop_data), error = function(e) NULL)
      filtered_df <- tryCatch(as.data.frame(filtered_pop_data), error = function(e) NULL)

      if (is.null(null_df) || is.null(filtered_df)) next
      if (!metric %in% colnames(null_df) || !metric %in% colnames(filtered_df)) next

      is_point_estimate <- nrow(filtered_df) == 1
      test_this_metric <- metric %in% metrics_to_test
      p_value <- NA
      p_value_bonferroni <- NA
      filtered_mean <- mean(filtered_df[[metric]], na.rm = TRUE)

      if (test_this_metric) {
        null_values_raw <- null_df[[metric]]
        filtered_values <- filtered_df[[metric]]

        has_nested <- !is.null(numSampleSets) && !is.null(numVariantSets) &&
                      numSampleSets > 1 && numVariantSets > 1
        expected_total <- if (has_nested) numSampleSets * numVariantSets else 0

        if (has_nested && length(null_values_raw) == expected_total &&
            length(filtered_values) == numSampleSets) {

          null_matrix <- matrix(null_values_raw, nrow = numSampleSets,
                                ncol = numVariantSets, byrow = TRUE)
          p_values <- numeric(numSampleSets)
          for (si in seq_len(numSampleSets)) {
            p_values[si] <- calculate_two_sided_pvalue(filtered_values[si], null_matrix[si, ])
          }
          n_sig <- sum(p_values < 0.05)
          p_value <- mean(p_values)
          null_values <- null_values_raw

          pvalue_results[[length(pvalue_results) + 1]] <- data.frame(
            Level = level, Population = pop, Metric = metric,
            P_Value = p_value,
            Num_Significant = n_sig,
            Prop_Significant = n_sig / numSampleSets,
            stringsAsFactors = FALSE
          )

        } else if (is_point_estimate) {
          filtered_value <- filtered_values[1]
          null_values <- null_values_raw
          p_value <- calculate_two_sided_pvalue(filtered_value, null_values)
          message(sprintf("      Point estimate test: filtered=%.4f vs %d null replicates, p=%.4f",
                         filtered_value, length(null_values), p_value))
          pvalue_results[[length(pvalue_results) + 1]] <- data.frame(
            Level = level, Population = pop, Metric = metric,
            P_Value = p_value,
            Num_Significant = NA_integer_,
            Prop_Significant = NA_real_,
            stringsAsFactors = FALSE
          )
        } else {
          stop(sprintf("Unexpected data structure for %s %s %s", level, pop, metric))
        }

        p_value_bonferroni <- min(p_value * n_tests, 1.0)
      } else {
        null_values <- null_df[[metric]]
      }

      p <- create_histogram_plot(null_values, filtered_df, metric, pop)

      annotations <- create_plot_annotations(test_this_metric, p_value, p_value_bonferroni,
                                              pvalue_results, n_tests)

      p <- p +
        labs(title = paste0(pop, " (", level, "): ", metric_labels[metric_idx]),
             subtitle = annotations$p_label,
             x = metric_labels[metric_idx],
             caption = annotations$caption) +
        theme_bw() +
        theme(
          legend.position = "top",
          axis.title.y.left = element_text(color = "#0072B2"),
          axis.text.y.left = element_text(color = "#0072B2"),
          axis.title.y.right = element_text(color = "#D55E00"),
          axis.text.y.right = element_text(color = "#D55E00")
        )

      plot_name <- paste0(pop, "_", level, "_", metric, "_Comparison")
      save_plot(p, plot_name)
      plots_created <- plots_created + 1
    }
  }

  message(paste0("    Created ", plots_created, " comparison plots for ", level))

  if (length(pvalue_results) > 0) {
    return(dplyr::bind_rows(pvalue_results))
  } else {
    return(NULL)
  }
}

# Reference: plotting.R line 896
plot_pairwise_beta_heatmaps <- function(filtered_pw, null_pw_list, level, save_plot) {

  if (is.null(filtered_pw) || is.null(null_pw_list) || length(null_pw_list) == 0) {
    message(paste0("  Skipping pairwise heatmaps for ", level, " (no data)"))
    return(NULL)
  }
  message(paste0("  Creating pairwise beta diversity heatmaps for ", level, "..."))

  extract_pairs <- function(pw) {
    if (is.list(pw) && all(sapply(pw, is.list))) {
      site1 <- c(); site2 <- c(); beta <- c()
      for (pop in pw) {
        if (is.list(pop) && all(c("TD_beta", "site1", "site2") %in% names(pop))) {
          site1 <- c(site1, as.character(pop$site1))
          site2 <- c(site2, as.character(pop$site2))
          beta <- c(beta, as.numeric(pop$TD_beta))
        }
      }
      df <- data.frame(site1 = site1, site2 = site2, beta = beta, stringsAsFactors = FALSE)
      if (nrow(df) == 0) return(NULL)
      return(df)
    } else if (is.list(pw) && all(c("TD_beta", "site1", "site2") %in% names(pw))) {
      df <- data.frame(site1 = as.character(pw$site1),
                       site2 = as.character(pw$site2),
                       beta = as.numeric(pw$TD_beta),
                       stringsAsFactors = FALSE)
      if (nrow(df) == 0) return(NULL)
      return(df)
    } else {
      return(NULL)
    }
  }

  build_matrix <- function(df) {
    sites <- sort(unique(c(df$site1, df$site2)))
    n <- length(sites)
    mat <- matrix(NA_real_, nrow = n, ncol = n, dimnames = list(sites, sites))
    for (i in seq_len(nrow(df))) {
      a <- df$site1[i]; b <- df$site2[i]; v <- df$beta[i]
      prev <- mat[a, b]
      if (is.na(prev)) {
        mat[a, b] <- v; mat[b, a] <- v
      } else {
        mat[a, b] <- mean(c(prev, v), na.rm = TRUE); mat[b, a] <- mat[a, b]
      }
    }
    diag(mat) <- 1
    return(mat)
  }

  filtered_df <- extract_pairs(filtered_pw)
  null_dfs <- lapply(null_pw_list, extract_pairs)

  if (is.null(filtered_df) || nrow(filtered_df) == 0 ||
      any(sapply(null_dfs, function(df) is.null(df) || nrow(df) == 0))) {
    message(paste0("  Skipping pairwise heatmaps for ", level, ": no valid data"))
    return(NULL)
  }

  filtered_mat <- build_matrix(filtered_df)
  null_mats <- lapply(null_dfs, build_matrix)
  null_array <- simplify2array(null_mats)
  null_mat <- apply(null_array, 1:2, mean, na.rm = TRUE)
  diff_mat <- filtered_mat - null_mat

  mat_to_long <- function(mat, value_name) {
    df <- as.data.frame(as.table(mat))
    colnames(df) <- c("Group1", "Group2", value_name)
    return(df)
  }

  filtered_long <- mat_to_long(filtered_mat, "BetaDiversity")
  null_long <- mat_to_long(null_mat, "BetaDiversity")
  diff_long <- mat_to_long(diff_mat, "Difference")

  # 1. Filtered beta heatmap
  p1 <- ggplot(filtered_long, aes(x = Group1, y = Group2, fill = BetaDiversity)) +
    geom_tile() +
    scale_fill_gradientn(colors = pals::ocean.balance(100),
                         name = "Beta Diversity\n(Filtered)") +
    labs(title = paste0("Pairwise Beta Diversity: Filtered Data (", level, ")"),
         x = "", y = "") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          legend.position = "right")
  save_plot(p1, paste0(level, "_PairwiseBeta_Filtered"))

  # 2. Null beta heatmap
  p2 <- ggplot(null_long, aes(x = Group1, y = Group2, fill = BetaDiversity)) +
    geom_tile() +
    scale_fill_gradientn(colors = pals::ocean.balance(100),
                         name = "Beta Diversity\n(Null Mean)") +
    labs(title = paste0("Pairwise Beta Diversity: Null Distribution Mean (", level, ")"),
         x = "", y = "") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          legend.position = "right")
  save_plot(p2, paste0(level, "_PairwiseBeta_Null"))

  # 3. Difference heatmap
  half_range <- max(abs(diff_long$Difference), na.rm = TRUE)
  if (is.infinite(half_range) || half_range == 0) half_range <- 1
  p3 <- ggplot(diff_long, aes(x = Group1, y = Group2, fill = Difference)) +
    geom_tile() +
    scale_fill_gradientn(colors = pals::ocean.balance(100),
                         name = "Difference\n(Filtered - Null)",
                         limits = c(-half_range, half_range)) +
    labs(title = paste0("Pairwise Beta Diversity Difference (", level, ")"),
         subtitle = "Positive = Filtered > Null, Negative = Filtered < Null",
         x = "", y = "") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          legend.position = "right")
  save_plot(p3, paste0(level, "_PairwiseBeta_Difference"))

  message(paste0("    Created 3 pairwise heatmap plots for ", level))
  invisible(NULL)
}

# Reference: plotting.R line 1021
# (already vendored in the original plots.R — kept byte-faithful)
compare_null_distributions_by_population <- function(null_by_pop, metric = "TD_beta",
                                                      n_perm = 10000,
                                                      p_adjust_method = "bonferroni") {
  pop_names <- names(null_by_pop)
  all_vals <- lapply(null_by_pop, function(x) {
    if (is.data.frame(x) && metric %in% colnames(x)) return(x[[metric]])
    if (is.vector(x)) return(x)
    return(NA)
  })
  valid_idx <- sapply(all_vals, function(x) !all(is.na(x)))
  pop_names <- pop_names[valid_idx]
  all_vals <- all_vals[valid_idx]
  combined <- data.frame(
    value = unlist(all_vals),
    Population = rep(pop_names, times = sapply(all_vals, length))
  )
  kw <- kruskal.test(value ~ Population, data = combined)
  global_p <- kw$p.value
  pairwise_p <- matrix(NA, nrow = length(pop_names), ncol = length(pop_names),
                       dimnames = list(pop_names, pop_names))
  summary_table <- data.frame(Pop1 = character(), Pop2 = character(),
                              p_value = numeric(), stringsAsFactors = FALSE)
  if (global_p < 0.05) {
    for (i in seq_along(pop_names)) {
      for (j in seq_along(pop_names)) {
        if (i < j) {
          x <- all_vals[[i]]; y <- all_vals[[j]]
          obs_diff <- abs(mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE))
          pooled <- c(x, y)
          n_x <- length(x)
          perm_diffs <- replicate(n_perm, {
            perm <- sample(pooled)
            abs(mean(perm[1:n_x], na.rm = TRUE) - mean(perm[-(1:n_x)], na.rm = TRUE))
          })
          p_val <- (sum(perm_diffs >= obs_diff) + 1) / (n_perm + 1)
          pairwise_p[i, j] <- p_val
          summary_table <- rbind(summary_table,
                                 data.frame(Pop1 = pop_names[i], Pop2 = pop_names[j],
                                            p_value = p_val))
        }
      }
    }
    summary_table$p_adj <- p.adjust(summary_table$p_value, method = p_adjust_method)
    for (k in seq_len(nrow(summary_table))) {
      i <- which(pop_names == summary_table$Pop1[k])
      j <- which(pop_names == summary_table$Pop2[k])
      pairwise_p[i, j] <- summary_table$p_adj[k]
    }
  }
  return(list(global_test = kw, pairwise_p_matrix = pairwise_p, summary_table = summary_table))
}

# END vendored reference code
# =============================================================================

# -----------------------------------------------------------------------------
# Driver
# -----------------------------------------------------------------------------

plot_files <- character(0)
csv_files  <- character(0)
problems   <- character(0)

plot_counter <- 0
save_plot <- function(p, plot_name) {
  plot_counter <<- plot_counter + 1
  safe_name <- gsub("[^A-Za-z0-9_-]", "_", plot_name)
  fname <- file.path(outdir,
                     sprintf("%s_%03d_%s_mqc.png", out_prefix, plot_counter, safe_name))
  tryCatch({
    ggsave(fname, plot = p, width = 12, height = 8, dpi = 150, bg = "white")
    message("  saved: ", fname)
    plot_files <<- c(plot_files, fname)
  }, error = function(e) {
    msg <- paste0("could not write ", fname, ": ", conditionMessage(e))
    message("  WARNING: ", msg)
    problems <<- c(problems, msg)
  })
  invisible(NULL)
}

# Lightweight has_usable_split: same as has_valid_split (the reference function name)
has_usable_split <- has_valid_split

plot_density_family <- function(data_list, label) {
  tryCatch({
    create_density_plots_internal(data_list, label, save_plot)
  }, error = function(e) {
    msg <- paste0("density plots for ", label, " failed: ", conditionMessage(e))
    message("  WARNING: ", msg)
    problems <<- c(problems, msg)
  })
}

run_null_test <- function(split_list, level, label) {
  for (metric in c("TD_beta", "TD_alpha", "TD_gamma")) {
    res <- tryCatch(
      compare_null_distributions_by_population(split_list, metric = metric,
                                                n_perm = n_permutations),
      error = function(e) {
        msg <- paste0("null test (", level, ", ", metric, ") failed: ", conditionMessage(e))
        message("  WARNING: ", msg)
        problems <<- c(problems, msg)
        NULL
      })
    if (is.null(res)) next
    st <- res$summary_table
    if (is.data.frame(st) && nrow(st) > 0) {
      csv <- file.path(outdir, paste0(out_prefix, "_nulltest_", level, "_", metric, ".csv"))
      tryCatch({
        write.csv(st, csv, row.names = FALSE)
        message("  saved: ", csv)
        csv_files <<- c(csv_files, csv)
      }, error = function(e) {
        msg <- paste0("could not write ", csv, ": ", conditionMessage(e))
        message("  WARNING: ", msg)
        problems <<- c(problems, msg)
      })
    }
    message(sprintf("  %-16s %-9s global Kruskal-Wallis p = %s%s",
                    level, metric,
                    format(res$global_test$p.value, digits = 4),
                    if (nrow(st) > 0) sprintf("  (%d pairwise test(s))", nrow(st)) else "  (no pairwise tests)"))
    if (nrow(st) > 0) print(st, row.names = FALSE)
  }
}

# Helper to extract pairwise data from a split list (SplitByPopulation /
# SplitBySuperPopulation) — each element is a list with TD_beta, site1, site2.
extract_pairwise_from_split <- function(split_list) {
  # split_list is a named list; we need the Pairwise*ByVariant element.
  # For the split outputs: SplitByPopulation[[1]] is a list of population dfs,
  # and PairwiseSplitByPopulation[[1]] is a list of pairwise dfs.
  # The aggregate stores pairwise data in PairwiseSplitByPopulation and
  # PairwiseTotalData (overall).
  return(NULL)  # placeholder — the pairwise data comes from the Pairwise* fields
}

# -----------------------------------------------------------------------------
# Read and validate inputs
# -----------------------------------------------------------------------------

message("=== plots.R: Hill diversity plots ===")
message("  null input         : ", null_file)
if (has_filtered) message("  filtered input     : ", filtered_file)
message("  output prefix      : ", out_prefix)
message("  output directory   : ", outdir)
message("  by population      : ", want_pop)
message("  by superpopulation : ", want_superpop)
message("  has filtered arm   : ", has_filtered)
message("  permutations       : ", n_permutations)
message("  seed               : ", seed, " (permutation draws only)")
message("  q (reported)       : ", q_hill)
message("")

null_hill <- tryCatch(readRDS(null_file), error = function(e) {
  stop("Could not read the null rds: ", conditionMessage(e))
})
if (!is.list(null_hill)) {
  stop("Null input is not a hill.R result (expected a list, got ", class(null_hill)[1], ")")
}
if (is.null(null_hill$TotalData)) {
  stop("Null input has no TotalData; present names: ",
       paste(names(null_hill), collapse = ", "))
}

filtered_hill <- NULL
if (has_filtered) {
  filtered_hill <- tryCatch(readRDS(filtered_file), error = function(e) {
    stop("Could not read the filtered rds: ", conditionMessage(e))
  })
  if (!is.list(filtered_hill)) {
    stop("Filtered input is not a hill.R result (expected a list, got ",
         class(filtered_hill)[1], ")")
  }
  message("  Filtered arm loaded successfully.")
}

# -----------------------------------------------------------------------------
# 1. Null-only density plots and tests (unchanged from original plots.R)
# -----------------------------------------------------------------------------

set.seed(seed)

message("--- Null distributions: overall ---")
plot_density_family(list(Overall = null_hill$TotalData), "Null_Overall")

if (want_pop && has_usable_split(null_hill$SplitByPopulation)) {
  message("--- Null distributions: by population ---")
  plot_density_family(null_hill$SplitByPopulation, "Null_ByPopulation")
} else if (want_pop) {
  message("--- Null distributions: by population -- skipped (empty or missing)")
}

if (want_superpop && has_usable_split(null_hill$SplitBySuperPopulation)) {
  message("--- Null distributions: by superpopulation ---")
  plot_density_family(null_hill$SplitBySuperPopulation, "Null_BySuperPopulation")
} else if (want_superpop) {
  message("--- Null distributions: by superpopulation -- skipped (empty or missing)")
}

message("")
message("--- Null distribution tests ---")

if (want_pop && has_usable_split(null_hill$SplitByPopulation)) {
  run_null_test(null_hill$SplitByPopulation, "population", "population")
}
if (want_superpop && has_usable_split(null_hill$SplitBySuperPopulation)) {
  run_null_test(null_hill$SplitBySuperPopulation, "superpopulation", "superpopulation")
}
if ((want_pop && has_usable_split(null_hill$SplitByPopulation)) ||
    (want_superpop && has_usable_split(null_hill$SplitBySuperPopulation))) {
  message("(global test is Kruskal-Wallis; pairwise tests run only when it is significant at p < 0.05)")
}

# -----------------------------------------------------------------------------
# 2. Filtered vs null comparison plots (NEW)
# -----------------------------------------------------------------------------

if (has_filtered) {
  message("")
  message("--- Filtered vs null comparison plots ---")

  # --- Null vs filtered overlay plots (density + vertical line) ---
  if (!is.null(null_hill$TotalData) && !is.null(filtered_hill$TotalData)) {
    # TotalData is a flat data.frame (not split by population).
    # Overlay: null as density, filtered as vertical lines.
    message("  Total: overall overlay")

    metrics <- c("TD_beta", "TD_alpha", "TD_gamma")
    metric_labels <- c("Beta Diversity", "Alpha Diversity", "Gamma Diversity")

    for (mi in seq_along(metrics)) {
      m <- metrics[mi]
      if (!m %in% colnames(null_hill$TotalData) || !m %in% colnames(filtered_hill$TotalData)) next

      null_vals <- na.omit(null_hill$TotalData[[m]])
      filt_vals <- filtered_hill$TotalData[[m]]

      # Build overlay plot
      p <- ggplot() +
        geom_density(data = data.frame(x = null_vals),
                     aes(x = x), fill = "#0072B2", alpha = 0.4, color = "black") +
        geom_vline(data = data.frame(x = filt_vals),
                   aes(xintercept = x),
                   color = "#D55E00", linewidth = 1.2, linetype = "dashed") +
        labs(title = paste0("Total Data: ", metric_labels[mi]),
             x = metric_labels[mi], y = "Density") +
        theme_bw()

      save_plot(p, paste0("Overall_", m, "_NullVsFiltered"))
    }
  }

  # --- By-population overlay ---
  if (want_pop && has_usable_split(null_hill$SplitByPopulation) &&
      has_usable_split(filtered_hill$SplitByPopulation)) {
    message("  Population-level overlays")
    pop_names <- names(null_hill$SplitByPopulation)
    for (pop in pop_names) {
      if (!pop %in% names(filtered_hill$SplitByPopulation)) next
      null_pop <- null_hill$SplitByPopulation[[pop]]
      filt_pop <- filtered_hill$SplitByPopulation[[pop]]
      if (!is.data.frame(null_pop) || !is.data.frame(filt_pop)) next

      for (m in c("TD_beta", "TD_alpha", "TD_gamma")) {
        if (!m %in% colnames(null_pop) || !m %in% colnames(filt_pop)) next
        p <- ggplot() +
          geom_density(data = data.frame(x = na.omit(null_pop[[m]])),
                       aes(x = x), fill = "#0072B2", alpha = 0.4, color = "black") +
          geom_vline(data = data.frame(x = filt_pop[[m]]),
                     aes(xintercept = x),
                     color = "#D55E00", linewidth = 1.2, linetype = "dashed") +
          labs(title = paste0(pop, ": ", m), x = m, y = "Density") +
          theme_bw()
        save_plot(p, paste0("Pop_", pop, "_", m, "_NullVsFiltered"))
      }
    }
  }

  # --- By-superpopulation overlay ---
  if (want_superpop && has_usable_split(null_hill$SplitBySuperPopulation) &&
      has_usable_split(filtered_hill$SplitBySuperPopulation)) {
    message("  Superpopulation-level overlays")
    sp_names <- names(null_hill$SplitBySuperPopulation)
    for (sp in sp_names) {
      if (!sp %in% names(filtered_hill$SplitBySuperPopulation)) next
      null_sp <- null_hill$SplitBySuperPopulation[[sp]]
      filt_sp <- filtered_hill$SplitBySuperPopulation[[sp]]
      if (!is.data.frame(null_sp) || !is.data.frame(filt_sp)) next

      for (m in c("TD_beta", "TD_alpha", "TD_gamma")) {
        if (!m %in% colnames(null_sp) || !m %in% colnames(filt_sp)) next
        p <- ggplot() +
          geom_density(data = data.frame(x = na.omit(null_sp[[m]])),
                       aes(x = x), fill = "#0072B2", alpha = 0.4, color = "black") +
          geom_vline(data = data.frame(x = filt_sp[[m]]),
                     aes(xintercept = x),
                     color = "#D55E00", linewidth = 1.2, linetype = "dashed") +
          labs(title = paste0(sp, ": ", m), x = m, y = "Density") +
          theme_bw()
        save_plot(p, paste0("SuperPop_", sp, "_", m, "_NullVsFiltered"))
      }
    }
  }

  # --- Pairwise beta diversity heatmaps ---
  # Both arms carry PairwiseTotalData (600x9 df) and
  # PairwiseSplitByPopulation[[1]], PairwiseSplitBySuperPopulation[[1]].
  message("")
  message("--- Pairwise beta diversity heatmaps ---")

  draw_pairwise_heatmaps <- function(null_pw, filtered_pw, level) {
    if (is.null(null_pw) || is.null(filtered_pw)) return(NULL)
    # null_pw is a 600x9 df with columns: sampleSet, variantSet, replicate,
    #   metric, site1, site2, beta (or TD_beta).
    # filtered_pw is a 600x9 df (but effectively 1 replicate per pair = 6 rows).
    # plot_pairwise_beta_heatmaps expects:
    #   filtered_pw  = single-pair result: list(TD_beta=..., site1=..., site2=...)
    #   null_pw_list = LIST OF single-pair results, one per replicate (e.g. 100 elements),
    #                  each: list(TD_beta=scalar, site1=..., site2=...)
    # The original code bundled all 6 pairs and 100 replicates into two big lists,
    # which made filtered a 4x4 matrix and each null a 2x2 matrix — hence
    # "non-conformable arrays" (4x4 minus 2x2).  Fix: iterate over site pairs,
    # and for each pair build the correct single-pair structures.

    # Identify unique (site1, site2) pairs
    pair_df <- unique(data.frame(site1 = null_pw$site1, site2 = null_pw$site2,
                                 stringsAsFactors = FALSE))
    if (nrow(pair_df) == 0) {
      message(paste0("  Skipping heatmaps for ", level, ": no pairwise data"))
      return(NULL)
    }

    for (pi in seq_len(nrow(pair_df))) {
      s1 <- pair_df$site1[pi]; s2 <- pair_df$site2[pi]

      # Single-pair null replicates: 100 list elements, each with ONE scalar beta
      null_mask <- null_pw$site1 == s1 & null_pw$site2 == s2
      null_beta_vals <- null_pw$TD_beta[null_mask]
      null_pw_list <- lapply(seq_along(null_beta_vals), function(r) {
        list(TD_beta = null_beta_vals[r], site1 = s1, site2 = s2)
      })

      # Single-pair filtered beta: ONE scalar (or a single value vector)
      filt_mask <- filtered_pw$site1 == s1 & filtered_pw$site2 == s2
      filt_beta <- filtered_pw$TD_beta[filt_mask]
      filt_pw <- list(TD_beta = filt_beta, site1 = s1, site2 = s2)

      plot_pairwise_beta_heatmaps(filt_pw, null_pw_list, level, save_plot)
    }
  }

  # Overall (TotalData) pairwise
  if (!is.null(null_hill$PairwiseTotalData) && !is.null(filtered_hill$PairwiseTotalData)) {
    message("  Pairwise heatmaps: overall (TotalData)")
    draw_pairwise_heatmaps(null_hill$PairwiseTotalData, filtered_hill$PairwiseTotalData,
                           "TotalData")
  }

  # By population
  if (want_pop && has_usable_split(null_hill$PairwiseSplitByPopulation) &&
      has_usable_split(filtered_hill$PairwiseSplitByPopulation)) {
    message("  Pairwise heatmaps: by population")
    pop_names <- names(null_hill$PairwiseSplitByPopulation)
    for (pop in pop_names) {
      if (!pop %in% names(filtered_hill$PairwiseSplitByPopulation)) next
      np <- null_hill$PairwiseSplitByPopulation[[pop]]
      fp <- filtered_hill$PairwiseSplitByPopulation[[pop]]
      draw_pairwise_heatmaps(np, fp, paste0("SplitByPopulation_", pop))
    }
  }

  # By superpopulation
  if (want_superpop && has_usable_split(null_hill$PairwiseSplitBySuperPopulation) &&
      has_usable_split(filtered_hill$PairwiseSplitBySuperPopulation)) {
    message("  Pairwise heatmaps: by superpopulation")
    sp_names <- names(null_hill$PairwiseSplitBySuperPopulation)
    for (sp in sp_names) {
      if (!sp %in% names(filtered_hill$PairwiseSplitBySuperPopulation)) next
      np <- null_hill$PairwiseSplitBySuperPopulation[[sp]]
      fp <- filtered_hill$PairwiseSplitBySuperPopulation[[sp]]
      draw_pairwise_heatmaps(np, fp, paste0("SplitBySuperPopulation_", sp))
    }
  }
} else {
  message("")
  message("Filtered arm not provided — skipping comparison and heatmap plots.")
}

# -----------------------------------------------------------------------------
# Metadata
# -----------------------------------------------------------------------------

metadata <- list(
  script = "plots.R",
  null_file = null_file,
  filtered_file = if (has_filtered) filtered_file else NULL,
  qHillNumber = q_hill,
  nPermutations = n_permutations,
  seed = seed,
  evaluatePopulation = want_pop,
  evaluateSuperPopulation = want_superpop,
  hasFilteredArm = has_filtered,
  levelsPlotted = c("overall",
                    if (want_pop && has_usable_split(null_hill$SplitByPopulation)) "by_population",
                    if (want_superpop && has_usable_split(null_hill$SplitBySuperPopulation)) "by_superpopulation"),
  plotFiles = plot_files,
  csvFiles = csv_files,
  problems = problems,
  runFinished = format(Sys.time())
)
meta_file <- file.path(outdir, paste0(out_prefix, "_plots_metadata.rds"))
saveRDS(metadata, meta_file)

message("")
message("=== plots.R: done ===")
message("  plots written : ", length(plot_files))
message("  CSVs written  : ", length(csv_files))
if (length(problems) > 0) {
  message("  problems      : ", length(problems), " (see the warnings above)")
}
message("  metadata      : ", meta_file)
