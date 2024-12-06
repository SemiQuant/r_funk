#' Create Heatmaps from Gene Signatures
#'
#' This function creates heatmaps from gene signatures, allowing for visualization of
#' expression patterns across different conditions with significance annotations.
#'
#' @param vst_data A DESeq2 vst transformed object containing expression data
#' @param dds A DESeq2 object containing experimental design information
#' @param result_names Character vector of result names to use for significance annotation
#' @param signatures_file Optional path to an Excel file containing gene signatures.
#'        If NULL (default), uses the package's built-in gene signatures dataset
#' @param signature_id Optional specific signature ID to plot (default: NULL, plots all)
#' @param n_signatures Optional number of signatures to plot (default: NULL, plots all)
#' @param padj_threshold Adjusted p-value threshold for significance (default: 0.05)
#' @param scale Character indicating if the values should be centered and scaled in either
#'        the "row" direction or "column" direction, or "none" (default: "row")
#' @param condition_column Name of the condition column in colData(dds) (default: "Condition")
#'
#' @return A list of pheatmap objects
#' @export
#'
#' @importFrom pheatmap pheatmap
#' @importFrom readxl read_excel
#' @importFrom dplyr filter distinct pull select
#' @importFrom SummarizedExperiment assay colData
#' @importFrom DESeq2 results
#'
#' @examples
#' \dontrun{
#' create_signature_heatmaps(vst_data, dds, 
#'                          result_names = c("condition_treated_vs_control"),
#'                          "path/to/signatures.xlsx")
#' }
create_signature_heatmaps <- function(vst_data, dds, result_names, 
                                    signatures_file = NULL,
                                    signature_id = NULL, n_signatures = NULL,
                                    padj_threshold = 0.05, scale = "row",
                                    condition_column = "Condition") {
    # Get gene signatures
    if (is.null(signatures_file)) {
        gene_symbols_df <- gene_signatures
    } else {
        gene_symbols_df <- readxl::read_excel(signatures_file)
    }
    
    # Filter by signature_id if provided
    if (!is.null(signature_id)) {
        gene_symbols_df <- gene_symbols_df %>% 
            dplyr::filter(ID %in% signature_id)
    }
    
    # Get significant genes for each comparison
    sig_results <- list()
    for(result_name in result_names) {
        res <- DESeq2::results(dds, name = result_name)
        sig_results[[result_name]] <- rownames(res)[which(res$padj < padj_threshold)]
    }
    
    # Create row annotation for significant genes
    create_row_annotation <- function(genes, sig_list) {
        annotation_df <- data.frame(row.names = genes)
        for(result_name in names(sig_list)) {
            annotation_df[[result_name]] <- factor(genes %in% sig_list[[result_name]], 
                                                 levels = c(TRUE, FALSE))
        }
        return(annotation_df)
    }
    
    # Get expression data for signature genes
    hmap_data <- SummarizedExperiment::assay(vst_data) %>% 
        as.data.frame() %>%
        tibble::rownames_to_column("gene_id") %>%
        dplyr::filter(gene_id %in% gene_symbols_df$Symbol)
    
    # Get condition information from colData
    col_data <- SummarizedExperiment::colData(dds) %>% 
        as.data.frame()
    
    # Create annotation columns for the heatmap
    annotation_col <- data.frame(
        Condition = col_data[[condition_column]],
        row.names = rownames(col_data)
    )
    
    # Verify column names match
    if (!all(rownames(annotation_col) == colnames(hmap_data[-1]))) {
        stop("Sample names in vst_data do not match sample names in dds")
    }
    
    # Create color schemes
    ann_colors <- list(
        Condition = setNames(RColorBrewer::brewer.pal(
            n = min(length(unique(col_data[[condition_column]])), 8), 
            name = "Dark2"
        ), unique(col_data[[condition_column]]))
    )
    
    # Add significance colors
    for(result_name in result_names) {
        ann_colors[[result_name]] <- c("TRUE" = "red", "FALSE" = "white")
    }
    
    # Set row names for heatmap
    rownames(hmap_data) <- hmap_data$gene_id
    hmap_data <- hmap_data[,-1] # Remove gene_id column
    
    # Get unique combinations of ID
    unique_ids <- unique(gene_symbols_df$ID)
    
    # Limit number of signatures if specified
    if (!is.null(n_signatures)) {
        unique_ids <- unique_ids[1:min(n_signatures, length(unique_ids))]
    }
    
    # Create a list to store heatmaps
    heatmap_list <- list()
    
    # Create a heatmap for each signature
    for(current_id in unique_ids) {
        # Filter genes for current signature
        current_genes <- gene_symbols_df %>%
            dplyr::filter(ID == current_id) %>%
            dplyr::pull(Symbol)
        
        # Subset heatmap data
        current_hmap_data <- hmap_data[current_genes, ]
        current_hmap_data <- current_hmap_data[complete.cases(current_hmap_data), ]
        
        # Skip if no genes found
        if (nrow(current_hmap_data) == 0) next
        
        # Create row annotation for current genes
        row_annotation <- create_row_annotation(rownames(current_hmap_data), sig_results)
        
        # Create title
        plot_title <- sprintf("Gene Signature: %s", current_id)
        
        # Create heatmap
        heatmap_list[[current_id]] <- pheatmap::pheatmap(
            current_hmap_data,
            scale = scale,
            show_rownames = TRUE,
            annotation_col = annotation_col,
            annotation_row = row_annotation,
            annotation_colors = ann_colors,
            fontsize_row = 8,
            cluster_cols = TRUE,
            cluster_rows = TRUE,
            main = plot_title,
            clustering_distance_rows = "correlation",
            clustering_distance_cols = "correlation",
            silent = TRUE
        )
    }
    
    return(heatmap_list)
}

#' List Available Gene Signatures
#'
#' This function displays all available gene signatures in the package's built-in dataset
#' or from a provided Excel file.
#'
#' @param signatures_file Optional path to an Excel file containing gene signatures.
#'        If NULL (default), uses the package's built-in gene signatures dataset
#'
#' @return A data frame containing unique signature IDs and the number of genes in each signature
#' @export
#'
#' @examples
#' \dontrun{
#' list_signatures()
#' list_signatures("path/to/signatures.xlsx")
#' }
list_signatures <- function(signatures_file = NULL) {
    # Get gene signatures
    if (is.null(signatures_file)) {
        gene_symbols_df <- gene_signatures
    } else {
        gene_symbols_df <- readxl::read_excel(signatures_file)
    }
    
    # Return unique signatures with gene counts
    unique_sigs <- gene_symbols_df %>%
        dplyr::group_by(ID) %>%
        dplyr::summarise(
            n_genes = dplyr::n(),
            genes = paste(Symbol, collapse = ", ")
        )
    
    return(unique_sigs)
}
