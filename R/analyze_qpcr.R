#' Analyze Real-time PCR Data with Delta-Delta CT Method or Direct CT Comparison
#' 
#' This function analyzes qPCR data using either the delta-delta CT method (with a reference target)
#' or direct CT comparison method. It provides comprehensive statistical analysis and visualization
#' of the results.
#' 
#' @param data A dataframe containing PCR data
#' @param ct_col Name of column containing Ct values (default: "Cq")
#' @param sample_col Name of column containing sample IDs (default: "Sample")
#' @param target_col Name of column containing probe/target names (default: "Target")
#' @param target Target probe to analyze
#' @param reference_sample Name of the reference sample (e.g. "Control")
#' @param reference_target Name of the reference target/probe (optional, for delta-delta CT method)
#' @param method Analysis method ("delta_delta_ct" or "direct_ct")
#' @param undetermined_value Value to use for "Undetermined" results (default = 45)
#' @param pos_control Name of positive control sample (optional)
#' @param neg_control Name of negative control sample (optional)
#' @param omit_col Name of column indicating which rows to omit
#' 
#' @return A list containing:
#' \itemize{
#'   \item target: The analyzed target name
#'   \item reference_target: The reference target used (if applicable)
#'   \item reference_sample: The reference sample used
#'   \item ddct_data: Raw data with calculated values
#'   \item group_stats: Statistical summaries by group
#'   \item contrasts: Results of statistical comparisons
#'   \item model: Linear model results
#'   \item model_summary: Summary statistics of the model
#'   \item plots: List of generated plots
#'   \item table: Formatted results table
#'   \item simplified_analysis: Results from simplified t-test analysis
#' }
#' 
#' @examples
#' # Load example data
#' data(qpcr_example)
#' 
#' # Example 1: Delta-Delta CT method (with reference target)
#' results_ddct <- analyze_pcr(
#'   data = qpcr_example,
#'   target = "Target_Gene",
#'   reference_sample = "Control",
#'   reference_target = "Reference_Gene",
#'   method = "delta_delta_ct"
#' )
#' 
#' # View results
#' results_ddct$table  # View formatted results table
#' plot(results_ddct$plots$fold_change)  # View fold change plot
#' plot(results_ddct$plots$combined)  # View all plots combined
#' 
#' # Example 2: Direct CT method (without reference target)
#' results_direct <- analyze_pcr(
#'   data = qpcr_example,
#'   target = "Target_Gene",
#'   reference_sample = "Control",
#'   method = "direct_ct"
#' )
#' 
#' # Compare results between methods
#' results_ddct$table
#' results_direct$table
#' 
#' @references 
#' Yuan, J.S., Reed, A., Chen, F. et al. (2006) Statistical analysis of real-time PCR data. 
#' BMC Bioinformatics 7, 85. \doi{10.1186/1471-2105-7-85}
#' 
#' Livak, K.J., Schmittgen, T.D. (2001) Analysis of relative gene expression data using 
#' real-time quantitative PCR and the 2(-Delta Delta C(T)) Method. Methods 25(4), 402-8.
#' \doi{10.1006/meth.2001.1262}
#' 
#' @seealso 
#' \itemize{
#'   \item \url{https://doi.org/10.1186/1471-2105-7-85} for the statistical methods paper
#'   \item \url{https://doi.org/10.1006/meth.2001.1262} for the original delta-delta CT method paper
#' }
#' 
#' @importFrom tidyverse %>%
#' @importFrom dplyr filter mutate group_by ungroup summarise select left_join rename deframe
#' @importFrom broom tidy
#' @importFrom cli cli_progress_bar cli_progress_update
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_point position_jitter scale_fill_manual theme_minimal theme element_text labs coord_trans scale_y_continuous geom_hline geom_bar
#' @importFrom cowplot theme_minimal_grid
#' @importFrom gridExtra grid.arrange
#' @importFrom stats pt qt sd t.test lm
#' @importFrom knitr kable
#' 
#' @export
analyze_pcr <- function(data, 
                       ct_col = "Cq",
                       sample_col = "Sample",
                       target_col = "Target",
                       target,
                       reference_sample,
                       reference_target = NULL,
                       method = "delta_delta_ct",
                       undetermined_value = 45,
                       pos_control = NULL,
                       neg_control = "NTC",
                       omit_col = "Omit") {
  
  cli_progress_bar("Analyzing PCR data", total = 7)
  
  # Remove omitted samples and controls
  data <- data %>%
    filter(!.data[[omit_col]]) %>%
    filter(.data[[sample_col]] != neg_control)
  
  # Filter data based on method
  if (method == "delta_delta_ct" && !is.null(reference_target)) {
    data <- data %>%
      filter(.data[[target_col]] %in% c(target, reference_target))
  } else {
    data <- data %>%
      filter(.data[[target_col]] == target)
  }
  
  # Convert Cq to numeric and handle Undetermined values
  data <- data %>%
    mutate(across(all_of(ct_col), ~{
      case_when(
        .x == "Undetermined" ~ undetermined_value,
        TRUE ~ as.numeric(as.character(.x))
      )
    }))
  
  cli_progress_update()
  
  # Process data based on method
  if (method == "delta_delta_ct" && !is.null(reference_target)) {
    # Delta-Delta CT method (with reference target)
    processed_data <- data %>%
      group_by(.data[[sample_col]]) %>%
      mutate(
        ref_ct = mean(.data[[ct_col]][.data[[target_col]] == reference_target]),
        dCT = .data[[ct_col]] - ref_ct
      ) %>%
      ungroup() %>%
      filter(.data[[target_col]] == target)
  } else {
    # Direct CT method (without reference target)
    ref_ct <- data %>%
      filter(.data[[sample_col]] == reference_sample) %>%
      summarise(mean_ct = mean(.data[[ct_col]])) %>%
      pull(mean_ct)
    
    processed_data <- data %>%
      mutate(
        dCT = .data[[ct_col]] - ref_ct,    # Calculate ΔCT relative to reference sample
        fold_change = 2^(dCT)              # Changed sign for direct CT method
      )
  }
  
  cli_progress_update()
  
  # For delta-delta CT method, calculate fold changes as before
  if (method == "delta_delta_ct" && !is.null(reference_target)) {
    ref_dct <- processed_data %>%
      filter(.data[[sample_col]] == reference_sample) %>%
      summarise(ref_dct = mean(dCT))
    
    ddct_data <- processed_data %>%
      mutate(
        ddCT = dCT - ref_dct$ref_dct,
        fold_change = 2^(-ddCT)
      )
  } else {
    # For direct CT, we already have fold changes
    ddct_data <- processed_data
  }
  
  cli_progress_update()
  
  # Statistical analysis
  model <- lm(fold_change ~ get(sample_col), data = ddct_data)
  model_summary <- summary(model)
  
  # Calculate group statistics
  group_stats <- ddct_data %>%
    group_by(.data[[sample_col]]) %>%
    summarise(
      mean_fc = mean(fold_change),
      sd_fc = sd(fold_change),
      n = n(),
      se = sd_fc / sqrt(n),
      ci_lower = mean_fc - qt(0.975, n-1) * se,
      ci_upper = mean_fc + qt(0.975, n-1) * se,
      .groups = "drop"
    )
  
  cli_progress_update()
  
  # Perform pairwise comparisons
  contrasts <- ddct_data %>%
    filter(.data[[sample_col]] != reference_sample) %>%
    group_by(.data[[sample_col]]) %>%
    summarise(
      t_test = list(tidy(t.test(
        fold_change,
        ddct_data$fold_change[ddct_data[[sample_col]] == reference_sample]
      ))),
      .groups = "drop"
    ) %>%
    unnest(t_test)
  
  cli_progress_update()
  
  # Calculate simplified t-test only for delta-delta CT method
  simplified_ttest <- if(method == "delta_delta_ct" && !is.null(reference_target)) {
    tryCatch({
      # Calculate mean CT differences for each sample
      ct_differences <- data %>%
        group_by(.data[[sample_col]], .data[[target_col]]) %>%
        summarise(
          mean_ct = mean(.data[[ct_col]]),
          sd_ct = sd(.data[[ct_col]]),
          n = n(),
          .groups = "drop"
        ) %>%
        pivot_wider(
          names_from = .data[[target_col]],
          values_from = c(mean_ct, sd_ct, n)
        ) %>%
        mutate(
          delta_ct = get(paste0("mean_ct_", target)) - get(paste0("mean_ct_", reference_target)),
          # Calculate pooled standard deviation
          sd_delta = sqrt(
            (get(paste0("sd_ct_", target))^2 / get(paste0("n_", target))) +
            (get(paste0("sd_ct_", reference_target))^2 / get(paste0("n_", reference_target)))
          )
        )
      
      # Get reference sample's delta CT
      reference_delta <- ct_differences %>%
        filter(.data[[sample_col]] == reference_sample) %>%
        pull(delta_ct)
      
      # Perform t-tests as described in the paper
      t_test_results <- ct_differences %>%
        filter(.data[[sample_col]] != reference_sample) %>%
        mutate(
          delta_delta_ct = delta_ct - reference_delta,
          fold_change = 2^(-delta_delta_ct),
          # Calculate t-statistic using pooled standard deviation
          t_stat = (delta_ct - reference_delta) / sd_delta,
          df = get(paste0("n_", target)) + get(paste0("n_", reference_target)) - 2,
          p_value = 2 * pt(-abs(t_stat), df)
        )
      
      list(
        differences = ct_differences,
        t_tests = t_test_results %>% 
          select(.data[[sample_col]], t_stat, p_value, delta_delta_ct, fold_change)
      )
    }, error = function(e) {
      message("\nNote: Simplified t-test analysis failed with error: ", e$message)
      return(NULL)
    })
  } else {
    # For direct CT method, calculate simplified t-test differently
    tryCatch({
      # Calculate t-test directly on CT values
      ct_stats <- data %>%
        group_by(.data[[sample_col]]) %>%
        summarise(
          mean_ct = mean(.data[[ct_col]]),
          sd_ct = sd(.data[[ct_col]]),
          n = n(),
          .groups = "drop"
        )
      
      # Get reference sample stats
      ref_stats <- ct_stats %>%
        filter(.data[[sample_col]] == reference_sample)
      
      # Calculate t-tests
      t_test_results <- ct_stats %>%
        filter(.data[[sample_col]] != reference_sample) %>%
        mutate(
          t_stat = (mean_ct - ref_stats$mean_ct) / 
                   sqrt((sd_ct^2/n) + (ref_stats$sd_ct^2/ref_stats$n)),
          df = n + ref_stats$n - 2,
          p_value = 2 * pt(-abs(t_stat), df),
          fold_change = 2^(-(mean_ct - ref_stats$mean_ct))
        )
      
      list(
        differences = ct_stats,
        t_tests = t_test_results %>%
          select(.data[[sample_col]], t_stat, p_value, fold_change)
      )
    }, error = function(e) {
      message("\nNote: Direct CT t-test analysis failed with error: ", e$message)
      return(NULL)
    })
  }
  
  # Generate plots first
  cli_progress_update()
  
  # Calculate significance levels first
  plot_data <- group_stats %>%
    left_join(contrasts %>% select(!!sym(sample_col), p.value), by = sample_col) %>%
    mutate(sig_level = case_when(
      .data[[sample_col]] == reference_sample ~ "Reference",
      p.value < 0.001 ~ "p < 0.001",
      p.value < 0.01 ~ "p < 0.01",
      p.value < 0.05 ~ "p < 0.05",
      TRUE ~ "Not significant"
    ))
  
  # Create a lookup table for significance levels
  sig_levels <- plot_data %>%
    select(.data[[sample_col]], sig_level) %>%
    deframe()
  
  # 1. Raw Ct values boxplot
  raw_ct_plot <- ggplot(data %>%
                        mutate(sig_level = sig_levels[.data[[sample_col]]]), 
                        aes(x = .data[[sample_col]], y = .data[[ct_col]], fill = sig_level)) +
    geom_boxplot(alpha = 0.8) +
    geom_point(position = position_jitter(width = 0.1), alpha = 0.5, color = "black") +
    scale_fill_manual(
      name = "Significance",
      values = c(
        "Reference" = "grey80",
        "p < 0.001" = "#0072B2",
        "p < 0.01" = "#009E73",
        "p < 0.05" = "#56B4E9",
        "Not significant" = "#E69F00"
    )) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Raw Ct Values by Sample",
         x = "Sample",
         y = "Ct Value")
  
  # 2. Fold change plot with error bars
  fold_change_plot <- ggplot(ddct_data %>%
                            mutate(sig_level = sig_levels[.data[[sample_col]]]),
                            aes(x = .data[[sample_col]], y = fold_change, fill = sig_level)) +
    # Add boxplots with significance coloring
    geom_boxplot(alpha = 0.8) +
    # Add individual points
    geom_point(position = position_jitter(width = 0.1),
               alpha = 0.5, color = "black") +
    coord_trans(y = "log10") +
    scale_y_continuous(
      breaks = c(1, 2, 5, 10, 20, 50, 100),
      limits = c(1, max(ddct_data$fold_change) * 1.1)
    ) +
    scale_fill_manual(
      name = "Significance",
      values = c(
        "Reference" = "grey80",
        "p < 0.001" = "#0072B2",
        "p < 0.01" = "#009E73",
        "p < 0.05" = "#56B4E9",
        "Not significant" = "#E69F00"
    )) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste("Fold Change Relative to", reference_sample),
         subtitle = paste("Target:", target),
         x = "Sample",
         y = "Fold Change (log scale)") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.5)
  
  # 3. Delta CT values plot
  dct_plot <- ggplot(processed_data %>%
                     mutate(sig_level = sig_levels[.data[[sample_col]]]),
                     aes(x = .data[[sample_col]], y = dCT, fill = sig_level)) +
    geom_boxplot(alpha = 0.8) +
    geom_point(position = position_jitter(width = 0.1), alpha = 0.5, color = "black") +
    scale_fill_manual(
      name = "Significance",
      values = c(
        "Reference" = "grey80",
        "p < 0.001" = "#0072B2",
        "p < 0.01" = "#009E73",
        "p < 0.05" = "#56B4E9",
        "Not significant" = "#E69F00"
    )) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Delta CT Values by Sample",
         subtitle = paste("Target:", target, 
                         if(!is.null(reference_target)) paste("normalized to", reference_target) else ""),
         x = "Sample",
         y = "Delta CT")
  
  # 4. Significance plot
  significance_plot <- ggplot(contrasts, 
                            aes(x = .data[[sample_col]], y = -log10(p.value))) +
    geom_bar(stat = "identity", fill = "orange", alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Statistical Significance",
         subtitle = paste("Comparison to", reference_sample),
         x = "Sample",
         y = "-log10(p-value)")
  
  # Combine plots
  combined_plots <- gridExtra::grid.arrange(
    raw_ct_plot, fold_change_plot, 
    dct_plot, significance_plot,
    ncol = 2
  )
  
  # Print only the fold change plot
  print(fold_change_plot)
  
  # Create results list with all components
  results <- list(
    target = target,
    reference_target = reference_target,
    reference_sample = reference_sample,
    ddct_data = ddct_data,
    group_stats = group_stats,
    contrasts = contrasts,
    model = model,
    model_summary = list(
      r_squared = model_summary$r.squared,
      adj_r_squared = model_summary$adj.r.squared,
      f_statistic = model_summary$fstatistic
    ),
    plots = list(
      raw_ct = raw_ct_plot,
      fold_change = fold_change_plot,
      delta_ct = dct_plot,
      significance = significance_plot,
      combined = combined_plots
    )
  )
  
  # Add simplified analysis to results
  results$simplified_analysis <- simplified_ttest
  
  # Print summary tables
  cat("\nAnalysis for target:", target, "\n")
  cat("Method:", ifelse(method == "delta_delta_ct" && !is.null(reference_target), 
                       "Delta-Delta CT", "Direct CT"), "\n")
  if (!is.null(reference_target)) {
    cat("Reference target:", reference_target, "\n")
  }
  cat("Reference sample:", reference_sample, "\n\n")
  
  # Format and print the main results table
  results_table <- group_stats %>%
    left_join(
      contrasts %>% 
        select(!!sym(sample_col), statistic, p.value) %>%
        rename(pvalue_std = p.value),
      by = sample_col
    )
  
  # Add simplified t-test results if available
  if (!is.null(simplified_ttest) && !is.null(simplified_ttest$t_tests)) {
    results_table <- results_table %>%
      left_join(
        simplified_ttest$t_tests %>%
          select(!!sym(sample_col), t_stat, p_value) %>%
          rename(pvalue_simp = p_value),
        by = sample_col
      )
  }
  
  # Store the full table in results before any formatting for display
  results$full_table <- results_table
  
  # Format table for display
  display_table <- results_table %>%
    mutate(
      mean_fc = round(mean_fc, 2),
      sd_fc = round(sd_fc, 2),
      ci_lower = round(ci_lower, 2),
      ci_upper = round(ci_upper, 2),
      statistic = round(statistic, 2),
      t_stat = if("t_stat" %in% names(.)) round(t_stat, 2) else NA,
      pvalue = format(pvalue_std, scientific = TRUE, digits = 2),
      significance = case_when(
        is.na(pvalue_std) ~ "",
        pvalue_std < 0.001 ~ "***",
        pvalue_std < 0.01 ~ "**",
        pvalue_std < 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      simp_pvalue = if("pvalue_simp" %in% names(.)) 
                      format(pvalue_simp, scientific = TRUE, digits = 2) 
                    else NA,
      simp_significance = case_when(
        !("pvalue_simp" %in% names(.)) ~ "",
        is.na(pvalue_simp) ~ "",
        pvalue_simp < 0.001 ~ "***",
        pvalue_simp < 0.01 ~ "**",
        pvalue_simp < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    ) %>%
    select(!!sym(sample_col), mean_fc, sd_fc, ci_lower, ci_upper, 
           statistic, pvalue, significance, 
           t_stat, simp_pvalue, simp_significance)
  
  # Store the formatted table
  results$table <- display_table
  
  # Print the formatted table
  print(knitr::kable(
    display_table,
    col.names = c("Sample", "mean_fc", "sd_fc", "ci_lower", "ci_upper", 
                  "t.stat", "pvalue", "sig", 
                  "simp.t.stat", "simp.pvalue", "simp.sig"),
    align = c('l', rep('r', 4), 'r', 'r', 'c', 'r', 'r', 'c')
  ))
  
  return(results)
} 