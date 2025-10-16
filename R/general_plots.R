#' Create Volcano Plot
#'
#' Creates a volcano plot using ggplot2 with customizable columns, annotations, and threshold lines
#'
#' @param data A data frame containing differential expression results
#' @param id_col Column name for gene/feature IDs (default: "ID")
#' @param pvalue_col Column name for p-values (default: "pvalue")
#' @param padj_col Column name for adjusted p-values (default: "padj")
#' @param logfc_col Column name for log2 fold change (default: "log2FoldChange")
#' @param padj_threshold Significance threshold for adjusted p-value (default: 0.05)
#' @param effect_size_threshold Log2 fold change threshold for significance (default: 2)
#' @param num_annotate Number of top genes to annotate (default: 10)
#' @param title Plot title (default: NULL)
#' @param make_plot If TRUE, return ggplot object; if FALSE, return data frame (default: TRUE)
#'
#' @return A ggplot object or data frame
#' @export
#'
#' @examples
#' # Create toy data
#' toy_data <- create_toy_volcano_data()
#' 
#' # Basic volcano plot
#' volcano_plot(toy_data)
#' 
#' # Custom column names and thresholds
#' volcano_plot(toy_data, padj_threshold = 0.01, effect_size_threshold = 1.5)
#' 
#' # Annotate more genes
#' volcano_plot(toy_data, num_annotate = 20)
#' 
#' # Get data instead of plot
#' volcano_data <- volcano_plot(toy_data, make_plot = FALSE)
volcano_plot <- function(data, 
                         id_col = "ID",
                         pvalue_col = "pvalue",
                         padj_col = "padj",
                         logfc_col = "log2FoldChange",
                         padj_threshold = 0.05,
                         effect_size_threshold = 2,
                         num_annotate = 10,
                         title = NULL,
                         make_plot = TRUE) {
  
  # Check required columns
  required_cols <- c(id_col, pvalue_col, padj_col, logfc_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Create standardized data frame for plotting
  volc.plot <- data.frame(
    ID = data[[id_col]],
    P = data[[pvalue_col]],
    P.adj = data[[padj_col]],
    EffectSize = data[[logfc_col]],
    stringsAsFactors = FALSE
  )
  
  # Remove rows with NA p-values
  volc.plot <- volc.plot[!is.na(volc.plot$P), ]
  
  # Calculate threshold based on both p-value and effect size
  volc.plot$threshold <- as.factor(volc.plot$P.adj <= padj_threshold & 
                                     abs(volc.plot$EffectSize) >= effect_size_threshold)
  
  # Determine which genes to label
  # Get the effect size threshold for top N genes
  effect_sizeThresh_2 <- abs(volc.plot[volc.plot$P.adj <= padj_threshold, ]$EffectSize)
  if (length(effect_sizeThresh_2) > 0) {
    effect_sizeThresh_2 <- effect_sizeThresh_2[order(abs(effect_sizeThresh_2), decreasing = TRUE)]
    if (length(effect_sizeThresh_2) >= num_annotate) {
      effect_sizeThresh_2 <- effect_sizeThresh_2[num_annotate]
    } else {
      effect_sizeThresh_2 <- min(effect_sizeThresh_2)
    }
  } else {
    effect_sizeThresh_2 <- Inf  # No genes to annotate
  }
  
  # Create labels for top genes
  volc.plot$lab <- ifelse(volc.plot$P.adj <= padj_threshold & 
                            abs(volc.plot$EffectSize) >= effect_sizeThresh_2, 
                          as.character(volc.plot$ID), "")
  
  # Create the plot
  volc_DE_1 <- ggplot2::ggplot(data = volc.plot,
                               ggplot2::aes(x = EffectSize, y = -log10(P),
                                            colour = threshold, text = ID)) +
    ggplot2::geom_point(alpha = 0.4, size = 1.75) +
    ggplot2::xlab("Effect Size") + 
    ggplot2::ylab("-log10 p-value") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::geom_text(ggplot2::aes(label = volc.plot$lab),
                       hjust = 0, vjust = 0, colour = "dark grey", size = 2, nudge_x = 0.2)
  
  # Add threshold lines if there are significant genes
  if (any(volc.plot$P.adj <= padj_threshold, na.rm = TRUE)) {
    volc_DE_1 <- volc_DE_1 + 
      ggplot2::geom_hline(yintercept = -log10(min(volc.plot$P[volc.plot$P.adj <= padj_threshold], na.rm = TRUE)), 
                         colour = "grey", linetype = "longdash")
  }
  
  volc_DE_1 <- volc_DE_1 + 
    ggplot2::geom_vline(xintercept = c(-1 * effect_size_threshold, effect_size_threshold), 
                       colour = "grey", linetype = "longdash")
  
  # Add title if provided
  if (!is.null(title)) {
    volc_DE_1 <- volc_DE_1 + ggplot2::ggtitle(title)
  }
  
  if (make_plot)
    return(volc_DE_1)
  else
    return(volc.plot)
}


#' Create Interactive Enrichment Plot
#'
#' Creates an interactive scatter plot for pathway/enrichment analysis results
#'
#' @param data A data frame containing enrichment results
#' @param description_col Column name for pathway/term descriptions (default: "Description")
#' @param generatio_col Column name for gene ratio (default: "GeneRatio")
#' @param pvalue_col Column name for p-values or adjusted p-values (default: "p.adjust")
#' @param count_col Column name for gene counts (default: "Count")
#' @param title Plot title (default: "Enrichment Analysis")
#' @param top_n Show only top N pathways (default: NULL shows all)
#' @param sort_by Column to sort pathways by (default: same as pvalue_col)
#' @param dark_theme Use dark theme for plot (default: TRUE)
#'
#' @return A plotly object
#' @export
#'
#' @examples
#' # Create toy data
#' toy_enrichment <- create_toy_enrichment_data()
#' 
#' # Basic enrichment plot
#' enrichment_plot(toy_enrichment)
#' 
#' # Show top 10 pathways
#' enrichment_plot(toy_enrichment, top_n = 10)
#' 
#' # Custom columns
#' enrichment_plot(toy_enrichment, pvalue_col = "pvalue", count_col = "Count")
enrichment_plot <- function(data,
                           description_col = "Description",
                           generatio_col = "GeneRatio",
                           pvalue_col = "p.adjust",
                           count_col = "Count",
                           title = "Enrichment Analysis",
                           top_n = NULL,
                           sort_by = NULL,
                           dark_theme = TRUE) {
  
  # Check required columns
  required_cols <- c(description_col, generatio_col, pvalue_col, count_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Create standardized data frame for plotting
  plot_data <- data.frame(
    Description = data[[description_col]],
    GeneRatio = data[[generatio_col]],
    p.adjust = data[[pvalue_col]],
    Count = data[[count_col]]
  )
  
  # Convert GeneRatio to numeric if it's a character (e.g., "5/100")
  if (is.character(plot_data$GeneRatio)) {
    plot_data$GeneRatio <- sapply(strsplit(plot_data$GeneRatio, "/"), function(x) {
      as.numeric(x[1]) / as.numeric(x[2])
    })
  }
  
  # Ensure Count is numeric
  plot_data$Count <- as.numeric(plot_data$Count)
  
  # Ensure p.adjust is numeric
  plot_data$p.adjust <- as.numeric(plot_data$p.adjust)
  
  # Sort data
  if (is.null(sort_by)) {
    sort_by <- "p.adjust"
  }
  plot_data <- plot_data[order(plot_data[[sort_by]]), ]
  
  # Select top N if specified
  if (!is.null(top_n)) {
    plot_data <- head(plot_data, top_n)
  }
  
  # Keep description as character for now
  plot_data$Description <- as.character(plot_data$Description)
  
  # Set theme colors
  if (dark_theme) {
    bg_color <- "#3A3A39"
    text_color <- "#FEFBE8"
    grid_color <- "rgba(254, 251, 232, 0.2)"
  } else {
    bg_color <- "#FFFFFF"
    text_color <- "#000000"
    grid_color <- "rgba(0, 0, 0, 0.1)"
  }
  
  # Create plot
  p <- plot_data %>%
    plotly::plot_ly(
      y = ~Description,
      x = ~GeneRatio,
      type = 'scatter',
      mode = 'markers',
      color = ~p.adjust,
      size = ~Count,
      sizes = c(10, 100),
      text = paste0(
        "Pathway: ", plot_data$Description,
        '<br>Gene Ratio: ', round(plot_data$GeneRatio, 3),
        '<br>Adjusted p-value: ', formatC(plot_data$p.adjust, format = "e", digits = 2),
        '<br>Count: ', plot_data$Count
      ),
      hoverinfo = 'text'
    ) %>%
    plotly::layout(
      title = list(text = title, font = list(color = text_color)),
      xaxis = list(
        title = "Gene Ratio",
        tickcolor = text_color,
        gridcolor = grid_color,
        linecolor = text_color
      ),
      yaxis = list(
        title = "",
        categoryorder = "array",
        categoryarray = rev(plot_data$Description),
        tickcolor = text_color,
        gridcolor = grid_color,
        linecolor = text_color
      ),
      plot_bgcolor = bg_color,
      paper_bgcolor = bg_color,
      font = list(color = text_color)
    )
  
  return(p)
}


#' Create Toy Volcano Plot Data
#'
#' Generate example differential expression data for testing volcano plots
#'
#' @param n_genes Number of genes to simulate (default: 1000)
#' @param n_sig Number of significant genes (default: 100)
#' @param seed Random seed for reproducibility (default: 123)
#'
#' @return A data frame with columns: ID, pvalue, padj, log2FoldChange
#' @export
#'
#' @examples
#' toy_data <- create_toy_volcano_data()
#' head(toy_data)
create_toy_volcano_data <- function(n_genes = 1000, n_sig = 100, seed = 123) {
  set.seed(seed)
  
  # Generate gene IDs
  gene_ids <- paste0("Gene_", sprintf("%04d", 1:n_genes))
  
  # Generate p-values (some significant, some not)
  pvalues <- c(
    runif(n_sig, 0, 0.01),  # Significant genes
    runif(n_genes - n_sig, 0.05, 1)  # Non-significant genes
  )
  pvalues <- sample(pvalues)  # Shuffle
  
  # Generate adjusted p-values (slightly higher than p-values)
  padj <- pmin(pvalues * 1.5, 1)
  
  # Generate log2 fold changes (bigger changes for significant genes)
  log2FC <- rnorm(n_genes, mean = 0, sd = 1)
  sig_indices <- which(padj < 0.05)
  log2FC[sig_indices] <- log2FC[sig_indices] + sample(c(-2, 2), length(sig_indices), replace = TRUE)
  
  # Create data frame
  data <- data.frame(
    ID = gene_ids,
    pvalue = pvalues,
    padj = padj,
    log2FoldChange = log2FC,
    stringsAsFactors = FALSE
  )
  
  return(data)
}


#' Create Toy Enrichment Data
#'
#' Generate example pathway enrichment data for testing enrichment plots
#'
#' @param n_pathways Number of pathways to simulate (default: 20)
#' @param seed Random seed for reproducibility (default: 123)
#'
#' @return A data frame with columns: Description, GeneRatio, p.adjust, Count
#' @export
#'
#' @examples
#' toy_enrichment <- create_toy_enrichment_data()
#' head(toy_enrichment)
create_toy_enrichment_data <- function(n_pathways = 20, seed = 123) {
  set.seed(seed)
  
  # Example pathway names
  pathway_types <- c(
    "Immune response",
    "Cell cycle regulation",
    "Metabolic process",
    "Signal transduction",
    "DNA repair",
    "Protein synthesis",
    "Apoptosis",
    "Cell differentiation",
    "Transport",
    "Transcription regulation"
  )
  
  # Generate pathway descriptions
  descriptions <- paste(
    sample(pathway_types, n_pathways, replace = TRUE),
    "-",
    sample(1:100, n_pathways, replace = TRUE)
  )
  
  # Generate gene counts
  counts <- sample(5:100, n_pathways, replace = TRUE)
  
  # Generate p-values
  p.adjust <- sort(runif(n_pathways, 1e-10, 0.05))
  
  # Generate gene ratios as fractions (convert to character format like "5/100")
  gene_ratios <- paste0(
    counts,
    "/",
    sample(100:500, n_pathways, replace = TRUE)
  )
  
  # Create data frame
  data <- data.frame(
    Description = descriptions,
    GeneRatio = gene_ratios,
    p.adjust = p.adjust,
    Count = counts,
    stringsAsFactors = FALSE
  )
  
  return(data)
}

