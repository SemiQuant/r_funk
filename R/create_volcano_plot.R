#' @importFrom ggplot2 ggplot aes geom_point xlab ylab ggtitle theme_bw theme geom_hline geom_vline annotate scale_y_continuous layer_scales expansion element_text element_blank
#' @importFrom ggrepel geom_text_repel
#' @importFrom DESeq2 results resultsNames
#' @export
create_volcano_plot <- function(dds, result_name = NULL, title = NULL, condition_a = NULL, 
                                condition_b = NULL, p_value_threshold = 0.05, effect_size_threshold = 1,
                                label_significant = TRUE, n_labels = 20, point_alpha = 0.4,
                                point_size = 1.75, label_size = 2) {
  
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required for label placement. Please install it using install.packages('ggrepel')")
  }
  
  # Handle different types of result specification
  if (is.list(result_name)) {
    # Case for list contrast (e.g., M vs R from M vs U - R vs U)
    deseq_result <- results(dds, contrast = result_name)
    if (is.null(title)) {
      title <- "Contrast Comparison"
    }
    if (is.null(condition_a) || is.null(condition_b)) {
      # Try to extract conditions from the contrast names
      contrast_names <- unlist(result_name)
      parts1 <- strsplit(contrast_names[1], "_vs_|_VS_")[[1]]
      parts2 <- strsplit(contrast_names[2], "_vs_|_VS_")[[1]]
      condition_a <- tail(parts1, 1)  # Last part of first contrast
      condition_b <- tail(parts2, 1)  # Last part of second contrast
    }
  } else {
    # Original behavior for name-based results
    if (is.null(result_name)) {
      result_names <- resultsNames(dds)
      if (length(result_names) < 2) {
        stop("No second result name available. Please specify result_name explicitly.")
      }
      result_name <- result_names[2]
    }
    deseq_result <- results(dds, name = result_name)
    if (is.null(title)) {
      title <- result_name
    }
    if (is.null(condition_a) || is.null(condition_b)) {
      parts <- strsplit(result_name, "_vs_|_VS_")[[1]]
      if (length(parts) >= 2) {
        if (is.null(condition_b)) condition_b <- parts[1]
        if (is.null(condition_a)) condition_a <- parts[2]
      } else {
        condition_a <- "Condition A"
        condition_b <- "Condition B"
      }
    }
  }
  
  # Create plot data frame
  volc_plot_data <- data.frame(
    ID = rownames(deseq_result),
    P = deseq_result$pvalue,
    P.adj = deseq_result$padj,
    EffectSize = deseq_result$log2FoldChange,
    threshold = FALSE
  )
  
  # Remove NA values
  volc_plot_data <- volc_plot_data[complete.cases(volc_plot_data), ]
  
  # Set threshold factor
  volc_plot_data$threshold <- as.factor(
    volc_plot_data$P.adj <= p_value_threshold & 
      abs(volc_plot_data$EffectSize) >= effect_size_threshold
  )
  
  # Add labels
  volc_plot_data$lab <- ""
  if (label_significant) {
    sig_genes <- volc_plot_data[
      volc_plot_data$P.adj <= p_value_threshold & 
        abs(volc_plot_data$EffectSize) >= effect_size_threshold, 
    ]
    sig_genes <- sig_genes[order(abs(sig_genes$EffectSize), decreasing = TRUE), ]
    top_genes <- head(sig_genes$ID, n_labels)
    volc_plot_data$lab <- ifelse(volc_plot_data$ID %in% top_genes,
                                 as.character(volc_plot_data$ID), "")
  }
  
  # Handle zero p-values
  zeros <- any(volc_plot_data$P == 0)
  min_nonzero_p <- min(volc_plot_data$P[volc_plot_data$P > 0], na.rm = TRUE)
  volc_plot_data$P <- replace(volc_plot_data$P, volc_plot_data$P == 0, min_nonzero_p)
  
  # Calculate plot limits
  y_max <- max(-log10(min_nonzero_p), na.rm = TRUE)
  plot_ceiling <- -log10(min_nonzero_p)
  plot_ceiling <- y_max * 1.2
  annotation_y <- y_max * 1.1
  
  # Add plotY column
  volc_plot_data$plotY <- pmin(-log10(volc_plot_data$P), plot_ceiling)
  
  # Create plot
  p <- ggplot(data = volc_plot_data,
              aes(x = EffectSize, y = plotY, colour = threshold, text = ID)) +
    geom_point(alpha = point_alpha, size = point_size) +
    xlab("Effect Size (log2FC)") +
    ylab("-log10 p-value") +
    ggtitle(title) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    # Add annotations
    annotate("text",
             x = -max(abs(volc_plot_data$EffectSize))/2,
             y = annotation_y,
             label = paste("Upregulated in", condition_a),
             hjust = 0.5, size = 3) +
    annotate("text",
             x = max(abs(volc_plot_data$EffectSize))/2,
             y = annotation_y,
             label = paste("Upregulated in", condition_b),
             hjust = 0.5, size = 3) +
    # Add arrows
    annotate("segment",
             x = -effect_size_threshold,
             xend = -max(abs(volc_plot_data$EffectSize))/2,
             y = annotation_y * 0.98,
             yend = annotation_y * 0.98,
             arrow = arrow(length = unit(0.2, "cm"))) +
    annotate("segment",
             x = effect_size_threshold,
             xend = max(abs(volc_plot_data$EffectSize))/2,
             y = annotation_y * 0.98,
             yend = annotation_y * 0.98,
             arrow = arrow(length = unit(0.2, "cm"))) +
    scale_y_continuous(limits = c(0, plot_ceiling),
                       expand = expansion(mult = c(0.05, 0.05)))
  
  # Add zero p-value annotation if needed
  if (zeros) {
    p <- p + annotate("text",
                      x = 0, y = plot_ceiling,
                      label = "p ~ 0", size = 3)
  }
  
  # Add labels if requested
  if (label_significant) {
    p <- p + ggrepel::geom_text_repel(
      aes(label = lab),
      size = label_size,
      colour = "dark grey",
      max.overlaps = Inf,
      min.segment.length = 0,
      box.padding = 0.5,
      segment.color = "grey50",
      segment.size = 0.2,
      seed = 42
    )
  }
  
  # Add threshold lines
  p <- p +
    geom_hline(yintercept = -log10(p_value_threshold),
               colour = "grey", linetype = "longdash") +
    geom_vline(xintercept = c(-1 * effect_size_threshold, effect_size_threshold),
               colour = "grey", linetype = "longdash") +
    geom_vline(xintercept = 0, colour = "grey",
               linetype = "solid", alpha = 0.5)
  
  return(p)
}

#' Create volcano plots for all comparisons in a DESeq2 object
#'
#' @param dds DESeqDataSet object
#' @param ... Additional arguments passed to create_volcano_plot function
#'
#' @return A list of ggplot objects, one for each comparison
#' @examples
#' \dontrun{
#' # For a single comparison:
#' plot <- create_volcano_plot(dds, "treatment_vs_control")
#'
#' # For custom conditions:
#' plot <- create_volcano_plot(dds, 
#'                           result_name = "treatment_vs_control",
#'                           condition_a = "Control",
#'                           condition_b = "Treatment")
#'
#' # For all comparisons:
#' all_plots <- create_all_volcano_plots(dds)
#' }
#' @export
create_all_volcano_plots <- function(dds, ...) {
  result_names <- resultsNames(dds)
  plots <- lapply(result_names, function(name) {
    create_volcano_plot(dds, result_name = name, ...)
  })
  names(plots) <- result_names
  return(plots)
}

# Example usage:
# For a single comparison:
# plot <- create_volcano_plot(dds, "treatment_vs_control")

# For custom conditions:
# plot <- create_volcano_plot(dds, 
#                           result_name = "treatment_vs_control",
#                           condition_a = "Control",
#                           condition_b = "Treatment")

# For all comparisons:
# all_plots <- create_all_volcano_plots(dds)