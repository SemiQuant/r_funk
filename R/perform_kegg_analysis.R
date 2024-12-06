#' Perform KEGG Pathway Enrichment Analysis
#'
#' This function performs KEGG pathway enrichment analysis on a set of
#' differentially expressed genes and creates visualization plots.
#'
#' @param results_df A data frame containing differential expression results
#' @param title Character string for plot titles and file names
#' @param p_cutoff Numeric value for adjusted p-value cutoff (default: 0.05)
#' @param fc_cutoff Numeric value for log2 fold change cutoff (default: 0)
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

perform_kegg_analysis <- function(results_df, title, p_cutoff = 0.05, fc_cutoff = 0, genome = "mmu") {
  # Get significant genes - check for all possible ID columns
  possible_cols <- c("gene_id", "peak_id", "site_id")
  id_col <- possible_cols[possible_cols %in% names(results_df)][1]
  
  if (is.na(id_col)) {
    warning("No valid ID column found in results dataframe")
    return(data.frame(message = "No valid ID column found"))
  }
  
  sig_genes <- results_df %>%
    filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
    filter(padj < p_cutoff & abs(log2FoldChange) > fc_cutoff) %>%
    pull(!!sym(id_col))
  
  # Convert gene symbols to ENTREZ IDs
  if (genome == "mmu") {
    genome <- org.Mm.eg.db
  } else if (genome == "hsa") {
    genome <- org.Hs.eg.db
  }
  
  entrez_ids <- mapIds(genome,
                      keys = sig_genes,
                      keytype = "SYMBOL",
                      column = "ENTREZID",
                      multiVals = "first")
  
  # Remove NAs
  entrez_ids <- entrez_ids[!is.na(entrez_ids)]
  
  # Perform KEGG enrichment analysis
  kk <- enrichKEGG(gene = entrez_ids,
                   organism = genome,
                   pvalueCutoff = 0.05)
  
  # Create plots if there are results
  if (nrow(as.data.frame(kk)) > 0) {
    # Calculate term similarity matrix
    kk <- pairwise_termsim(kk)

    kk@result$Description <- gsub(" - .*$", "", kk@result$Description)
    
    # Dot plot
    p1 <- dotplot(kk, title = paste0("KEGG Pathways - ", title))
    print(p1)
    
    # Return results table
    return(as.data.frame(kk))
  } else {
    return(data.frame(message = "No significant KEGG pathways found"))
  }
}
