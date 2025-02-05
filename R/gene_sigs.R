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
#' @param lfc_threshold Log2 fold change threshold for filtering (default: NULL)
#' @param prefilter_padj Adjusted p-value threshold for pre-filtering (default: NULL)
#' @param prefilter_lfc Log2 fold change threshold for pre-filtering (default: NULL)
#' @param scale Character indicating if the values should be centered and scaled in either
#'        the "row" direction or "column" direction, or "none" (default: "row")
#' @param condition_column Name of the condition column in colData(dds) (default: "Condition")
#' @param perform_gsea Logical indicating whether to perform GSEA analysis (default: FALSE)
#' @param gsea_nperm Deprecated. Number of permutations is now automatically determined by fgseaMultilevel
#'
#' @return A list of pheatmap objects
#' @export
#'
#' @importFrom pheatmap pheatmap
#' @importFrom readxl read_excel
#' @importFrom dplyr filter distinct pull select
#' @importFrom SummarizedExperiment assay colData
#' @importFrom DESeq2 results
#' @importFrom fgsea fgsea
#' @importFrom stats p.adjust
#' @importFrom DT datatable
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
                                    condition_column = "Condition",
                                    lfc_threshold = NULL,
                                    prefilter_padj = NULL,
                                    prefilter_lfc = NULL,
                                    perform_gsea = FALSE,
                                    gsea_nperm = 1000) {
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
    
    # Get significant genes for each comparison and apply prefiltering if specified
    sig_results <- list()
    filtered_genes <- NULL
    
    for(result_name in result_names) {
        res <- DESeq2::results(dds, name = result_name)
        
        # Basic significance filtering
        sig_genes <- rownames(res)[which(res$padj < padj_threshold)]
        
        # Add LFC threshold if specified
        if (!is.null(lfc_threshold)) {
            sig_genes <- rownames(res)[which(res$padj < padj_threshold & 
                                           abs(res$log2FoldChange) > lfc_threshold)]
        }
        
        sig_results[[result_name]] <- sig_genes
        
        # Apply prefiltering if specified
        if (!is.null(prefilter_padj) || !is.null(prefilter_lfc)) {
            prefiltered <- rownames(res)
            
            if (!is.null(prefilter_padj)) {
                prefiltered <- rownames(res)[which(res$padj < prefilter_padj)]
            }
            
            if (!is.null(prefilter_lfc)) {
                prefiltered <- intersect(
                    prefiltered,
                    rownames(res)[which(abs(res$log2FoldChange) > prefilter_lfc)]
                )
            }
            
            if (is.null(filtered_genes)) {
                filtered_genes <- prefiltered
            } else {
                filtered_genes <- union(filtered_genes, prefiltered)
            }
        }
    }
    
    # Apply prefiltering to gene_symbols_df if filtered_genes is not NULL
    if (!is.null(filtered_genes)) {
        gene_symbols_df <- gene_symbols_df %>%
            dplyr::filter(Symbol %in% filtered_genes)
        
        if (nrow(gene_symbols_df) == 0) {
            stop("No genes remain after prefiltering. Consider relaxing the prefiltering thresholds.")
        }
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

        # Create row annotation for current genes only if no prefiltering is applied
        row_annotation <- if (is.null(prefilter_padj) && is.null(prefilter_lfc)) {
            create_row_annotation(rownames(current_hmap_data), sig_results)
        } else {
            NULL
        }
        
        # Create title
        plot_title <- sprintf("Gene Signature: %s", current_id)
        
        # Handle color palette for conditions
        n_conditions <- length(unique(col_data[[condition_column]]))
        color_palette <- if (n_conditions <= 2) {
            # For 2 or fewer conditions, use a simple color pair
            c("#E41A1C", "#377EB8")[1:n_conditions]
        } else if (n_conditions <= 8) {
            # For 3-8 conditions, use Dark2
            RColorBrewer::brewer.pal(n_conditions, "Dark2")
        } else {
            # For more than 8 conditions, interpolate colors
            colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(n_conditions)
        }
        
        # Update annotation colors
        ann_colors <- list(
            Condition = setNames(color_palette, 
                               unique(col_data[[condition_column]]))
        )
        
        # Add significance colors only if no prefiltering
        if (is.null(prefilter_padj) && is.null(prefilter_lfc)) {
            for(result_name in result_names) {
                ann_colors[[result_name]] <- c("TRUE" = "#5bf5a5", "FALSE" = "white")
            }
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
                annotation_row = row_annotation,  # Will be NULL if prefiltering is applied
                annotation_colors = ann_colors,
                fontsize_row = 8,
                cluster_cols = TRUE,
                cluster_rows = TRUE,
                main = plot_title
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
    
    # After creating sig_results, add GSEA analysis
    gsea_results <- list()
    gsea_formatted <- NULL
    if (perform_gsea) {
        for(result_name in result_names) {
            message(sprintf("Running GSEA for %s", result_name))
            res <- DESeq2::results(dds, name = result_name)
            
            # Create ranked gene list for GSEA using only log2FoldChange
            gene_ranks <- tryCatch({
                res_df <- as.data.frame(res)
                # Remove any rows with NA values in log2FoldChange
                res_df <- res_df[!is.na(res_df$log2FoldChange) & 
                                is.finite(res_df$log2FoldChange), ]
                
                # Use log2FoldChange as ranking metric
                ranks <- res_df$log2FoldChange
                names(ranks) <- rownames(res_df)
                
                # Sort ranks
                sort(ranks, decreasing = TRUE)
            }, error = function(e) {
                warning(sprintf("Failed to create gene ranks for %s: %s", 
                              result_name, e$message))
                return(NULL)
            })
            
            if (is.null(gene_ranks)) next
            
            # Debug info
            message(sprintf("Number of ranked genes: %d", length(gene_ranks)))
            
            # Create gene sets more efficiently
            gene_sets <- split(gene_symbols_df$Symbol, gene_symbols_df$ID)
            
            # Debug info
            message(sprintf("Number of gene sets: %d", length(gene_sets)))
            message(sprintf("Gene set sizes: %s", 
                          paste(sapply(gene_sets, length), collapse = ", ")))
            
            # Ensure gene symbols match
            common_genes <- intersect(names(gene_ranks), unlist(gene_sets))
            if (length(common_genes) == 0) {
                warning(sprintf("No overlap between ranked genes and signature genes for %s", 
                              result_name))
                next
            }
            
            message(sprintf("Number of common genes: %d", length(common_genes)))
            
            # Run GSEA with memory-efficient settings
            gsea_result <- tryCatch({
                fgsea::fgseaMultilevel(
                    pathways = gene_sets,
                    stats = gene_ranks,
                    minSize = 1,
                    maxSize = 500,
                    eps = 0.0,
                    nPermSimple = 1000
                )
            }, error = function(e) {
                warning(sprintf("GSEA failed for %s: %s", result_name, e$message))
                return(NULL)
            })
            
            if (!is.null(gsea_result) && nrow(gsea_result) > 0) {
                # Add comparison name to results
                gsea_result <- cbind(Comparison = result_name, gsea_result)
                
                # Store raw results
                gsea_results[[result_name]] <- gsea_result
                
                # Append to combined formatted results
                if (is.null(gsea_formatted)) {
                    gsea_formatted <- gsea_result
                } else {
                    gsea_formatted <- rbind(gsea_formatted, gsea_result)
                }
                
                message(sprintf("GSEA completed for %s with %d results", 
                              result_name, nrow(gsea_result)))
            } else {
                warning(sprintf("No significant GSEA results for %s", result_name))
            }
        }
        
        # Format the combined results if any exist
        if (!is.null(gsea_formatted) && nrow(gsea_formatted) > 0) {
            gsea_formatted <- gsea_formatted %>%
                dplyr::mutate(
                    pval = formatC(pval, format = "e", digits = 2),
                    padj = formatC(padj, format = "e", digits = 2),
                    log2err = round(log2err, 3),
                    ES = round(ES, 3),
                    NES = round(NES, 3),
                    leadingEdge = sapply(leadingEdge, paste, collapse = ", ")
                )
            
            # Create interactive datatable
            gsea_table <- DT::datatable(gsea_formatted,
                options = list(
                    pageLength = 10,
                    scrollX = TRUE,
                    dom = 'Bfrtip'
                ),
                filter = 'top',
                class = 'cell-border stripe'
            ) %>%
                DT::formatStyle(
                    columns = c('ES', 'NES'),
                    backgroundColor = DT::styleInterval(
                        c(0),
                        c('#FFC4C4', '#C4FFC4')
                    )
                )
        } else {
            gsea_table <- NULL
        }
    }

    # Return results
    result <- list(
        heatmaps = heatmap_list,
        gsea_results = if (perform_gsea) gsea_results else NULL,
        gsea_table = if (perform_gsea) gsea_table else NULL
    )

    return(invisible(result))
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
