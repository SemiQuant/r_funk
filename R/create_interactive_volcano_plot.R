#' Create an Interactive Volcano Plot
#' 
#' @param dds A DESeqDataSet object
#' @param result_name Name of the result to plot
#' @param title Custom title for the plot
#' @param condition_a Name of condition A
#' @param condition_b Name of condition B
#' @param p_value_threshold P-value threshold for significance (default: 0.05)
#' @param effect_size_threshold Log2 fold change threshold (default: 1)
#' @param n_labels Number of top genes to label (default: 20)
#' @param point_alpha Transparency of points (0-1) (default: 0.4)
#' 
#' @return A highcharter object
#' @export
#'
#' @importFrom highcharter highchart hc_add_series hc_xAxis hc_yAxis hc_title hc_tooltip hc_plotOptions hc_annotations
#' @importFrom DESeq2 results resultsNames
#' 
create_interactive_volcano_plot <- function(dds, 
                                          result_name = NULL, 
                                          title = NULL, 
                                          condition_a = NULL,
                                          condition_b = NULL,
                                          p_value_threshold = 0.05,
                                          effect_size_threshold = 1,
                                          n_labels = 20,
                                          point_alpha = 0.4) {
  
  # Get results
  deseq_result <- results(dds, name = result_name)
  
  # Set title if not provided
  if (is.null(title)) {
    title <- result_name
  }
  
  # Set conditions if not provided
  if (is.null(condition_a) || is.null(condition_b)) {
    parts <- strsplit(result_name, "_vs_|_VS_")[[1]]
    if (length(parts) >= 2) {
      condition_b <- parts[1]
      condition_a <- parts[2]
    } else {
      condition_a <- "Condition A"
      condition_b <- "Condition B"
    }
  }
  
  # Create plot data frame
  volc_data <- data.frame(
    gene = rownames(deseq_result),
    p_value = deseq_result$pvalue,
    p_adj = deseq_result$padj,
    log2FC = deseq_result$log2FoldChange,
    stringsAsFactors = FALSE
  )
  
  # Remove NA values and add categories
  volc_data <- volc_data[complete.cases(volc_data), ]
  volc_data$neg_log10_p <- -log10(volc_data$p_value)
  
  # Add categories
  volc_data$category <- "Not Significant"
  up_b_idx <- volc_data$p_adj <= p_value_threshold & volc_data$log2FC >= effect_size_threshold
  up_a_idx <- volc_data$p_adj <= p_value_threshold & volc_data$log2FC <= -effect_size_threshold
  
  volc_data$category[up_b_idx] <- paste("Upregulated in", condition_b)
  volc_data$category[up_a_idx] <- paste("Upregulated in", condition_a)
  
  # Calculate y-axis limits for annotation placement
  y_max <- max(-log10(volc_data$p_value), na.rm = TRUE)
  annotation_y <- y_max * 1.1
  x_max <- max(abs(volc_data$log2FC), na.rm = TRUE)
  
  # Create the base plot
  hc <- highchart() %>%
    hc_chart(type = "scatter", zoomType = "xy") %>%
    hc_title(text = title) %>%
    hc_xAxis(
      title = list(text = "Effect Size (log2FC)"),
      plotLines = list(
        list(value = 0, color = "#666666", width = 1),
        list(value = effect_size_threshold, color = "#666666", width = 1, dashStyle = "dash"),
        list(value = -effect_size_threshold, color = "#666666", width = 1, dashStyle = "dash")
      )
    ) %>%
    hc_yAxis(
      title = list(text = "-log10 p-value"),
      plotLines = list(
        list(value = -log10(p_value_threshold), color = "#666666", width = 1, dashStyle = "dash")
      ),
      max = annotation_y * 1.2
    ) %>%
    hc_tooltip(
      useHTML = TRUE,
      headerFormat = "",
      pointFormat = "<b>Gene:</b> {point.gene}<br/><b>log2FC:</b> {point.x:.2f}<br/><b>-log10 p-value:</b> {point.y:.2f}<br/><b>p-value:</b> {point.pvalue:.2e}"
    )
  
  # Add annotations
  hc <- hc %>%
    hc_annotations(list(
      labels = list(
        list(
          point = list(x = -x_max/2, y = annotation_y),
          text = paste("Upregulated in", condition_a),
          style = list(fontSize = "12px")
        ),
        list(
          point = list(x = x_max/2, y = annotation_y),
          text = paste("Upregulated in", condition_b),
          style = list(fontSize = "12px")
        )
      ),
      shapes = list(
        # Left arrow
        list(
          type = "path",
          points = list(
            list(x = -effect_size_threshold, y = annotation_y * 0.98),
            list(x = -x_max/2, y = annotation_y * 0.98)
          ),
          markerEnd = "arrow"
        ),
        # Right arrow
        list(
          type = "path",
          points = list(
            list(x = effect_size_threshold, y = annotation_y * 0.98),
            list(x = x_max/2, y = annotation_y * 0.98)
          ),
          markerEnd = "arrow"
        )
      )
    ))
  
  # Plot each category separately with ggplot2-like colors
  # Not significant points
  ns_data <- volc_data[volc_data$category == "Not Significant", ]
  hc <- hc %>% hc_add_series(
    data = list_parse(data.frame(
      x = ns_data$log2FC,
      y = ns_data$neg_log10_p,
      gene = ns_data$gene,
      pvalue = ns_data$p_value
    )),
    name = "Not Significant",
    color = paste0("rgba(150, 150, 150, ", point_alpha, ")"),  # grey
    marker = list(radius = 4, symbol = "circle", lineWidth = 0)
  )
  
  # Upregulated in condition B
  up_b_data <- volc_data[up_b_idx, ]
  hc <- hc %>% hc_add_series(
    data = list_parse(data.frame(
      x = up_b_data$log2FC,
      y = up_b_data$neg_log10_p,
      gene = up_b_data$gene,
      pvalue = up_b_data$p_value
    )),
    name = paste("Upregulated in", condition_b),
    color = paste0("rgba(227, 74, 51, ", point_alpha, ")"),  # ggplot2 red
    marker = list(radius = 4, symbol = "circle", lineWidth = 0)
  )
  
  # Upregulated in condition A
  up_a_data <- volc_data[up_a_idx, ]
  hc <- hc %>% hc_add_series(
    data = list_parse(data.frame(
      x = up_a_data$log2FC,
      y = up_a_data$neg_log10_p,
      gene = up_a_data$gene,
      pvalue = up_a_data$p_value
    )),
    name = paste("Upregulated in", condition_a),
    color = paste0("rgba(51, 124, 183, ", point_alpha, ")"),  # ggplot2 blue
    marker = list(radius = 4, symbol = "circle", lineWidth = 0)
  )
  
  # Add labels for top genes
  if (n_labels > 0) {
    sig_genes <- volc_data[volc_data$p_adj <= p_value_threshold & 
                          abs(volc_data$log2FC) >= effect_size_threshold, ]
    sig_genes <- sig_genes[order(abs(sig_genes$log2FC), decreasing = TRUE), ]
    top_genes <- head(sig_genes$gene, n_labels)
    
    label_data <- volc_data[volc_data$gene %in% top_genes, ]
    hc <- hc %>% hc_add_series(
      data = list_parse(data.frame(
        x = label_data$log2FC,
        y = label_data$neg_log10_p,
        name = label_data$gene,
        pvalue = label_data$p_value
      )),
      type = "scatter",
      name = "Top Genes",
      dataLabels = list(enabled = TRUE, format = "{point.name}"),
      marker = list(radius = 4),
      showInLegend = FALSE
    )
  }
  
  # Add hover effects
  hc <- hc %>% hc_plotOptions(
    scatter = list(
      marker = list(
        states = list(
          hover = list(enabled = TRUE, lineColor = "rgb(100,100,100)")
        )
      )
    )
  )
  
  return(hc)
}

create_interactive_volcano_plot(dds_rna, "Condition_R_vs_U", point_alpha = 0.3)
