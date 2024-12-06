#' Perform GO Term Enrichment Analysis
#'
#' This function performs Gene Ontology (GO) enrichment analysis on a set of
#' differentially expressed genes and creates visualization plots.
#'
#' @param results_df A data frame containing differential expression results
#' @param title Character string for plot titles and file names
#' @param p_cutoff Numeric value for adjusted p-value cutoff (default: 0.05)
#' @param fc_cutoff Numeric value for log2 fold change cutoff (default: 0)
#' @param ont Character string specifying GO ontology: "BP" (Biological Process),
#'        "MF" (Molecular Function), or "CC" (Cellular Component) (default: "BP")
#' @param genome Character string specifying genome: "mmu" (mouse) or "hsa" (human)
#'        (default: "mmu")
#'
#' @return A data frame containing GO enrichment results or a message if no
#'         significant terms are found
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Filters significant genes based on p-value and fold change cutoffs
#'   \item Converts gene symbols to ENTREZ IDs
#'   \item Performs GO term enrichment analysis
#'   \item Creates visualization plots for significant terms
#' }
#'
#' @import org.Mm.eg.db
#' @import org.Hs.eg.db
#' @import clusterProfiler
#' @import dplyr
#' @importFrom AnnotationDbi mapIds
#'
#' @examples
#' results_df <- data.frame(
#'   gene_id = c("Trp53", "Brca1"),
#'   padj = c(0.01, 0.03),
#'   log2FoldChange = c(1.5, -2.0)
#' )
#' go_results <- perform_go_analysis(results_df, "Treatment vs Control", ont = "BP")
#'
#' @export

perform_go_analysis <- function(results_df, title, p_cutoff = 0.05, fc_cutoff = 0,
                              ont = "BP", genome = "mmu") {
  # Add error handling for package availability
  required_packages <- c("clusterProfiler", "dplyr")
  if (genome == "mmu") {
    required_packages <- c(required_packages, "org.Mm.eg.db")
  } else if (genome == "hsa") {
    required_packages <- c(required_packages, "org.Hs.eg.db")
  }
  
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not available. Please install it using BiocManager::install()"))
    }
  }
  
  # Add tryCatch for the main analysis
  tryCatch({
    # Get significant genes - check for all possible ID columns
    possible_cols <- c("gene_id", "peak_id", "site_id")
    id_col <- possible_cols[possible_cols %in% names(results_df)][1]
    
    if (is.na(id_col)) {
      warning("No valid ID column found in results dataframe")
      return(data.frame(message = "No valid ID column found"))
    }
    
    # Check for required columns
    required_cols <- c("padj", "log2FoldChange")
    missing_cols <- required_cols[!required_cols %in% names(results_df)]
    if (length(missing_cols) > 0) {
      stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
    }
    
    sig_genes <- results_df %>%
      dplyr::filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
      dplyr::filter(padj < p_cutoff & abs(log2FoldChange) > fc_cutoff) %>%
      dplyr::pull(!!sym(id_col))
    
    # Set up genome database
    if (genome == "mmu") {
      genome_db <- org.Mm.eg.db
      org <- "mmu"
    } else if (genome == "hsa") {
      genome_db <- org.Hs.eg.db
      org <- "hsa"
    } else {
      stop("Unsupported genome. Please use 'mmu' for mouse or 'hsa' for human.")
    }
    
    # Convert gene symbols to ENTREZ IDs
    tryCatch({
      entrez_ids <- mapIds(genome_db,
                          keys = sig_genes,
                          keytype = "SYMBOL",
                          column = "ENTREZID",
                          multiVals = "first")
    }, error = function(e) {
      stop(paste("Error converting gene symbols to ENTREZ IDs:", e$message))
    })
    
    # Remove NAs and provide more informative message
    na_count <- sum(is.na(entrez_ids))
    if (na_count > 0) {
      warning(sprintf("%d genes could not be mapped to ENTREZ IDs", na_count))
    }
    
    entrez_ids <- entrez_ids[!is.na(entrez_ids)]
    
    if (length(entrez_ids) == 0) {
      return(data.frame(message = "No valid ENTREZ IDs found"))
    }
    
    # Perform GO enrichment analysis
    ego <- enrichGO(gene = entrez_ids,
                   OrgDb = genome_db,
                   ont = ont,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   readable = TRUE)
    
    # Create plots if there are results
    if (nrow(as.data.frame(ego)) > 0) {
      # Calculate term similarity matrix
      ego <- pairwise_termsim(ego)
      
      # Dot plot
      p1 <- dotplot(ego, 
                   title = paste0("GO Terms (", ont, ") - ", title),
                   showCategory = 20)
      print(p1)
      
      # Return results table
      return(as.data.frame(ego))
    } else {
      return(data.frame(message = paste0("No significant GO terms found for ontology: ", ont)))
    }
  }, error = function(e) {
    message("Error in GO analysis: ", e$message)
    return(data.frame(message = paste("Error in GO analysis:", e$message)))
  })
} 