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
        ann_colors[[result_name]] <- c("TRUE" = "#5bf5a5", "FALSE" = "white")
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
    
    # Set up the plotting device for multiple plots
    n_plots <- length(unique_ids)
    if(n_plots > 0) {
        # Calculate grid layout
        n_rows <- ceiling(sqrt(n_plots))
        n_cols <- ceiling(n_plots/n_rows)
        
        # Set up the plotting layout
        graphics::par(mfrow = c(n_rows, n_cols))
    }
    
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
        if (nrow(current_hmap_data) == 0) {
            warning(sprintf("No valid genes found for signature %s", current_id))
            next
        }
        
        # Check for infinite values and replace with NAs
        # Convert to matrix first to handle infinite values properly
        current_hmap_data <- as.matrix(current_hmap_data)
        current_hmap_data[!is.finite(current_hmap_data)] <- NA
        
        # Remove rows with all NAs
        current_hmap_data <- current_hmap_data[rowSums(!is.na(current_hmap_data)) > 0, , drop = FALSE]
        
        # Skip if no valid data remains
        if (nrow(current_hmap_data) == 0) {
            warning(sprintf("No valid data remains for signature %s after removing NA/Inf values", current_id))
            next
        }

        # Create row annotation for current genes
        row_annotation <- create_row_annotation(rownames(current_hmap_data), sig_results)
        
        # Create title
        plot_title <- sprintf("Gene Signature: %s", current_id)
        
        # Handle color palette for conditions
        n_conditions <- length(unique(col_data[[condition_column]]))
        color_palette <- if (n_conditions < 3) {
            # For 2 or fewer conditions, use a simple color pair
            c("#E41A1C", "#377EB8")[1:n_conditions]
        } else {
            RColorBrewer::brewer.pal(
                n = min(n_conditions, 8),
                name = "Dark2"
            )
        }
        
        # Update annotation colors
        ann_colors <- list(
            Condition = setNames(color_palette, 
                               unique(col_data[[condition_column]]))
        )
        
        for(result_name in result_names) {
            ann_colors[[result_name]] <- c("TRUE" = "#5bf5a5", "FALSE" = "white")
        }

        # Create heatmap with tryCatch
        tryCatch({
            # Prepare data for clustering
            clustering_data <- current_hmap_data
            
            # Scale the data if requested (before clustering)
            if (scale == "row") {
                clustering_data <- t(scale(t(clustering_data)))
            } else if (scale == "column") {
                clustering_data <- scale(clustering_data)
            }
            
            # Handle NAs in correlation calculation
            row_cor <- cor(t(clustering_data), use = "pairwise.complete.obs")
            col_cor <- cor(clustering_data, use = "pairwise.complete.obs")
            
            # Replace any remaining NAs with 0 correlation
            row_cor[is.na(row_cor)] <- 0
            col_cor[is.na(col_cor)] <- 0
            
            # Convert correlations to distances
            row_dist <- as.dist(1 - row_cor)
            col_dist <- as.dist(1 - col_cor)
            
            # Create hierarchical clustering objects
            row_hclust <- hclust(row_dist, method = "complete")
            col_hclust <- hclust(col_dist, method = "complete")
            
            # Create heatmap with pre-calculated clustering
            hmap <- pheatmap::pheatmap(
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
                clustering_method = "complete",
                clustering = list(
                    Gene = row_hclust,
                    Sample = col_hclust
                )
            )
            
            # Store the heatmap
            heatmap_list[[current_id]] <- hmap
        }, error = function(e) {
            warning(sprintf("Failed to create clustered heatmap for signature %s: %s", 
                          current_id, e$message))
            
            # Try without clustering if distance calculation fails
            tryCatch({
                hmap <- pheatmap::pheatmap(
                    current_hmap_data,
                    scale = scale,
                    show_rownames = TRUE,
                    annotation_col = annotation_col,
                    annotation_row = row_annotation,
                    annotation_colors = ann_colors,
                    fontsize_row = 8,
                    cluster_cols = FALSE,
                    cluster_rows = FALSE,
                    main = paste(plot_title, "(unclustered)"),
                )
                
                # Store the unclustered heatmap
                heatmap_list[[current_id]] <- hmap
            }, error = function(e2) {
                warning(sprintf("Failed to create unclustered heatmap for signature %s: %s", 
                              current_id, e2$message))
            })
        })
    }
    
    # Reset the plotting parameters
    graphics::par(mfrow = c(1, 1))
    
    return(invisible(heatmap_list))
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
