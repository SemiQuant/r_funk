#' Perform MSigDB Gene Set Enrichment Analysis
#' 
#' This function performs Gene Set Enrichment Analysis (GSEA) using MSigDB gene sets
#' and creates visualization plots. Uses the multilevel GSEA algorithm for more accurate
#' p-values and better computational efficiency.
#' 
#' @param results_df A data frame containing differential expression results with
#'        required columns: gene_symbol (or gene_id) and stat
#' @param title Character string for plot titles and file names
#' @param p_cutoff Numeric value for adjusted p-value cutoff (default: 0.05)
#' @param genome Character string specifying genome: "mmu" (mouse) or "hsa" (human)
#'        (default: "mmu")
#' @param category Character string specifying MSigDB category (default: "H" for Hallmark)
#' @param min_size Minimum size of gene sets to analyze (default: 15)
#' @param max_size Maximum size of gene sets to analyze (default: 500)
#' @param eps Double value, boundary for calculating the p-value (default: 0)
#' @param nPermSimple Number of permutations for the simple algorithm (default: 1000)
#' @param ranking_method Character string specifying the ranking method to use:
#'        "auto" (automatically choose based on available data),
#'        "stat" (use raw test statistic),
#'        "logfc" (use log fold change),
#'        "logfc_pval" (combine log fold change with p-value significance),
#'        "combined" (combine p-value and effect size),
#'        "expression" (incorporate expression level) (default: "auto")
#' @param handle_ties Logical indicating whether to automatically handle ties
#'        by using more sophisticated ranking metrics (default: TRUE)
#' 
#' @return A list containing:
#' \itemize{
#'   \item results: A data frame containing GSEA results with columns:
#'     \itemize{
#'       \item pathway: Name of the gene set
#'       \item pval: Nominal p-value
#'       \item padj: Adjusted p-value (BH method)
#'       \item ES: Enrichment score
#'       \item NES: Normalized enrichment score
#'       \item nMoreExtreme: Number of permutations with more extreme score
#'       \item size: Size of the gene set
#'       \item leadingEdge: Leading edge genes
#'     }
#'   \item plot: GSEA table visualization for significant pathways
#'   \item status: Character string indicating analysis status
#' }
#' 
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Validates input data and parameters
#'   \item Creates a ranked gene list based on the test statistic
#'   \item Retrieves MSigDB gene sets for the specified genome and category
#'   \item Performs multilevel GSEA analysis
#'   \item Creates visualization for significant pathways
#' }
#' 
#' The multilevel GSEA algorithm provides more accurate p-values than standard GSEA,
#' especially for pathways with strong enrichment scores. It automatically determines
#' the optimal number of permutations needed for each pathway.
#' 
#' @import msigdbr
#' @import fgsea
#' @import dplyr
#' @import tibble
#' @import ggplot2
#' 
#' @examples
#' results_df <- data.frame(
#'   gene_id = c("Trp53", "Brca1", "Myc"),
#'   padj = c(0.01, 0.03, 0.02),
#'   stat = c(1.5, -2.0, 2.5)
#' )
#' msig_results <- perform_mSig_analysis(results_df, "Treatment vs Control")
#' 
#' @export
perform_mSig_analysis <- function(results_df, title, p_cutoff = 0.05, 
                                genome = "mmu", category = "H",
                                min_size = 15, max_size = 500,
                                eps = 0, nPermSimple = 1000,
                                ranking_method = "auto",
                                handle_ties = TRUE) {
    # Input validation with informative error messages
    if (!is.data.frame(results_df)) {
        stop("results_df must be a data frame")
    }
    if (!is.character(title) || length(title) != 1) {
        stop("title must be a single character string")
    }
    if (!is.numeric(p_cutoff) || p_cutoff <= 0 || p_cutoff > 1) {
        stop("p_cutoff must be a number between 0 and 1")
    }
    if (!genome %in% c("mmu", "hsa")) {
        stop("genome must be either 'mmu' for mouse or 'hsa' for human")
    }
    if (!ranking_method %in% c("auto", "stat", "logfc", "logfc_pval", "combined", "expression")) {
        stop("ranking_method must be one of: 'auto', 'stat', 'logfc', 'logfc_pval', 'combined', 'expression'")
    }
    
    # Check for required packages
    required_packages <- c("msigdbr", "fgsea", "dplyr", "tibble", "ggplot2")
    missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
    if (length(missing_packages) > 0) {
        stop("Missing required packages: ", paste(missing_packages, collapse = ", "),
             ". Please install them using BiocManager::install()")
    }
    
    # Create empty result template
    empty_result <- list(
        results = data.frame(),
        plot = NULL,
        status = "initialized",
        warnings = character(),
        ranking_info = list()
    )
    
    # Wrap the main analysis in tryCatch for robust error handling
    tryCatch({
        # Check for required columns and handle gene_id/gene_symbol
        if ("gene_id" %in% names(results_df)) {
            results_df <- results_df %>%
                rename(gene_symbol = gene_id)
        } else if (!"gene_symbol" %in% names(results_df)) {
            stop("Missing required column: gene_id or gene_symbol")
        }
        
        # Check for required columns based on ranking method
        if (ranking_method == "logfc") {
            if (!"log2FoldChange" %in% names(results_df)) {
                stop("Missing required column: log2FoldChange (needed for logfc ranking method)")
            }
        } else if (ranking_method == "logfc_pval") {
            if (!all(c("log2FoldChange", "padj") %in% names(results_df))) {
                stop("Missing required columns: log2FoldChange and padj (needed for logfc_pval ranking method)")
            }
        } else if (!"stat" %in% names(results_df)) {
            stop("Missing required column: stat (needed for stat-based ranking methods)")
        }
        
        # Check for additional ranking metrics
        available_metrics <- intersect(
            c("log2FoldChange", "padj", "pvalue", "stat", "baseMean"),
            names(results_df)
        )
        
        # Function to create ranking statistic based on method
        create_ranking <- function(df, method = "auto") {
            if (method == "auto") {
                # Choose best method based on available metrics
                if (all(c("padj", "stat", "baseMean") %in% names(df))) {
                    method <- "expression"
                } else if (all(c("padj", "stat") %in% names(df))) {
                    method <- "combined"
                } else {
                    method <- "stat"
                }
            }
            
            ranking <- switch(method,
                "stat" = df$stat,
                "logfc" = {
                    if (!"log2FoldChange" %in% names(df)) {
                        stop("log2FoldChange column required for logfc ranking method")
                    }
                    # For logfc method, we'll use log2FoldChange directly
                    df$log2FoldChange
                },
                "logfc_pval" = {
                    if (!all(c("log2FoldChange", "padj") %in% names(df))) {
                        stop("Both log2FoldChange and padj columns required for logfc_pval ranking method")
                    }
                    # Combine log fold change with signed -log10(padj)
                    df$log2FoldChange * -log10(pmax(df$padj, 1e-10))
                },
                "combined" = {
                    # Combine p-value and effect size
                    -log10(pmax(df$padj, 1e-10)) * sign(df$stat)
                },
                "expression" = {
                    # Incorporate expression level
                    -log10(pmax(df$padj, 1e-10)) * sign(df$stat) * 
                        log2(pmax(df$baseMean, 1) + 1)
                },
                stop("Invalid ranking method")
            )
            
            # Add small random noise to break remaining ties if requested
            if (handle_ties) {
                set.seed(42)  # for reproducibility
                ranking <- ranking + rnorm(length(ranking), 0, sd(ranking, na.rm = TRUE) / 1e6)
            }
            
            return(ranking)
        }
        
        # Create initial ranking to check for ties
        initial_ranking <- if (ranking_method == "logfc") {
            results_df$log2FoldChange
        } else if (ranking_method == "logfc_pval") {
            results_df$log2FoldChange * -log10(pmax(results_df$padj, 1e-10))
        } else {
            results_df$stat
        }
        
        tie_counts <- table(initial_ranking)
        tied_values <- tie_counts[tie_counts > 1]
        tie_percentage <- (sum(tied_values * (tied_values - 1)) / length(initial_ranking)) * 100
        
        # Determine if we need to use a more sophisticated ranking
        if (handle_ties && tie_percentage > 25) {
            message("High proportion of ties detected. Using more sophisticated ranking method...")
            results_df$ranking <- create_ranking(results_df, method = ranking_method)
        } else {
            results_df$ranking <- initial_ranking
        }
        
        # Create ranked gene list
        ranked_genes <- results_df %>%
            filter(!is.na(ranking)) %>%
            arrange(desc(ranking)) %>%
            dplyr::select(gene_symbol, ranking) %>%
            deframe()
        
        if (length(ranked_genes) == 0) {
            return(modifyList(empty_result, list(
                status = "No valid genes found after filtering"
            )))
        }
        
        # Store ranking information
        ranking_info <- list(
            method_used = if(handle_ties && tie_percentage > 25) ranking_method else "stat",
            initial_ties = tie_percentage,
            final_ties = (sum(duplicated(ranked_genes)) / length(ranked_genes)) * 100
        )
        
        # Determine score type
        score_type <- if (all(ranked_genes > 0)) {
            message("All statistics are positive. Using scoreType = 'pos' for better results.")
            "pos"
        } else if (all(ranked_genes < 0)) {
            message("All statistics are negative. Using scoreType = 'neg' for better results.")
            "neg"
        } else {
            "std"
        }
        
        # Get MSigDB gene sets
        m_df <- msigdbr(species = if(genome == "mmu") "Mus musculus" else "Homo sapiens",
                       category = category)
        pathways <- split(m_df$gene_symbol, m_df$gs_name)
        
        if (length(pathways) == 0) {
            return(modifyList(empty_result, list(
                status = "No pathways found for the specified category",
                ranking_info = ranking_info
            )))
        }
        
        # Perform multilevel GSEA
        message("Running multilevel GSEA analysis...")
        
        fgsea_results <- fgsea::fgseaMultilevel(
            pathways = pathways,
            stats = ranked_genes,
            minSize = min_size,
            maxSize = max_size,
            eps = eps,
            nPermSimple = nPermSimple,
            scoreType = score_type
        )
        
        if (nrow(fgsea_results) == 0) {
            return(modifyList(empty_result, list(
                status = "No enriched pathways found",
                ranking_info = ranking_info
            )))
        }
        
        # Format results
        fgsea_results <- fgsea_results %>%
            as_tibble() %>%
            arrange(padj) %>%
            mutate(
                across(where(is.numeric), \(x) round(x, 3)),
                leadingEdge = sapply(leadingEdge, paste, collapse = ", ")
            )
        
        # Get significant pathways
        sig_results <- fgsea_results %>%
            filter(padj < p_cutoff)
        
        # Create GSEA table plot for significant pathways
        plot <- NULL
        if (nrow(sig_results) > 0) {
            message(sprintf("Creating GSEA table for %d significant pathways...", nrow(sig_results)))
            
            plot <- plotGseaTable(
                pathways = pathways[sig_results$pathway],
                stats = ranked_genes,
                fgseaRes = sig_results,
                gseaParam = 0.5
            )
        } else {
            message("No significant pathways found at the specified p-value cutoff")
        }
        
        # Return results with ranking information
        return(list(
            results = fgsea_results,
            plot = plot,
            status = "success",
            ranking_info = ranking_info
        ))
        
    }, error = function(e) {
        warning(sprintf("GSEA analysis failed: %s", e$message))
        return(modifyList(empty_result, list(
            status = paste("Error in GSEA analysis:", e$message)
        )))
    })
}
