#' Perform KEGG Pathway Enrichment Analysis
#'
#' This function performs KEGG pathway enrichment analysis on a set of
#' differentially expressed genes and creates visualization plots.
#'
#' @param results_df A data frame containing differential expression results
#' @param title Character string for plot titles and file names
#' @param p_cutoff Numeric value for adjusted p-value cutoff (default: 0.05)
#' @param fc_cutoff Numeric value for log2 fold change cutoff (default: 0)
#' @param genome Character string specifying the genome ("mmu" for mouse or "hsa" for human)
#' @param id_col Character string specifying the column name containing gene identifiers
#' @param cat_show Numeric value specifying number of categories to show in the plot (default: 20)
#'
#' @return A data frame containing KEGG enrichment results or a message if no
#'         significant pathways are found
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Filters significant genes based on p-value and fold change cutoffs
#'   \item Converts gene symbols to ENTREZ IDs
#'   \item Performs KEGG pathway enrichment analysis
#'   \item Creates visualization plots for significant pathways
#' }
#'
#' @import org.Mm.eg.db
#' @import clusterProfiler
#' @import dplyr
#' @import org.Hs.eg.db
#' @importFrom ggplot2 ggsave
#'
#' @examples
#' results_df <- data.frame(
#'   gene_id = c("Trp53", "Brca1"),
#'   padj = c(0.01, 0.03),
#'   log2FoldChange = c(1.5, -2.0)
#' )
#' kegg_results <- perform_kegg_analysis(results_df, "Treatment vs Control")
#'
#' @export

perform_kegg_analysis <- function(results_df, title, p_cutoff = 0.05, fc_cutoff = 0, 
                                genome = "mmu", cat_show = 20) {
    # Input validation
    if (!is.data.frame(results_df)) {
        stop("results_df must be a data frame")
    }
    if (!is.character(title) || length(title) != 1) {
        stop("title must be a single character string")
    }
    
    # Set up genome database
    if (genome == "mmu") {
        gen <- org.Mm.eg.db
    } else if (genome == "hsa") {
        gen <- org.Hs.eg.db
    } else {
        stop("genome must be either 'mmu' for mouse or 'hsa' for human")
    }

    # Identify ID column
    possible_cols <- c("gene_id", "peak_id", "site_id")
    id_col <- possible_cols[possible_cols %in% names(results_df)][1]
    
    if (is.na(id_col)) {
        warning("No valid ID column found in results dataframe")
        return(data.frame(message = "No valid ID column found"))
    }

    # Extract significant genes with error handling
    sig_genes <- try({
        results_df %>%
            filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
            filter(padj < p_cutoff & abs(log2FoldChange) > fc_cutoff) %>%
            pull(!!sym(id_col))
    })
    
    if (inherits(sig_genes, "try-error") || length(sig_genes) == 0) {
        message("No significant genes found matching the criteria")
        return(data.frame(message = "No significant genes found"))
    }
    
    if (length(sig_genes) < 10) {
        warning("Few significant genes found (", length(sig_genes), "). Results may not be meaningful.")
    }
  
    # Convert gene symbols to ENTREZ IDs with error handling
    entrez_ids <- try({
        bitr(sig_genes, 
             fromType = "SYMBOL",
             toType = "ENTREZID",
             OrgDb = gen)
    })
    
    if (inherits(entrez_ids, "try-error") || nrow(entrez_ids) == 0) {
        message("Failed to convert gene symbols to ENTREZ IDs")
        return(data.frame(message = "Gene ID conversion failed"))
    }
    
    # Perform KEGG enrichment analysis
    kk <- try({
        enrichKEGG(gene = entrez_ids$ENTREZID,
                  organism = genome,
                  pvalueCutoff = 0.05)
    })

    if (inherits(kk, "try-error") || nrow(kk@result) == 0) {
        message("No significant KEGG pathways found")
        return(data.frame(message = "No significant pathways"))
    }
    
    # Clean pathway descriptions
    kk@result$Description <- gsub(" - .*$", "", kk@result$Description)
    
    # Convert results to data frame
    results <- as.data.frame(kk)
    results$Comparison <- title
    
    # Create and save dotplot
    plot_file <- paste0("kegg_dotplot_", make.names(title), ".pdf")
    tryCatch({
        p <- dotplot(kk, showCategory = cat_show, title = title)
        # ggsave(filename = plot_file, 
        #        plot = p, 
        #        width = 10, 
        #        height = 8)
        # message("Plot saved to: ", plot_file)
        print(p)
    }, error = function(e) {
        warning("Failed to create or save plot: ", e$message)
    })
    return(results)
}
