#' Gene Signatures Dataset
#'
#' A dataset containing gene signatures information
#'
#' @format A data frame containing gene signatures
#' @source Internal gene signatures compilation
"gene_signatures" 

#' Example qPCR Data
#'
#' A dataset containing example qPCR data with control and treatment samples
#'
#' @format A data frame with 12 rows and 4 variables:
#' \describe{
#'   \item{Sample}{Sample name (Control or Treatment)}
#'   \item{Target}{Target gene name (Target_Gene or Reference_Gene)}
#'   \item{Cq}{Quantification cycle (Cq) values}
#'   \item{Omit}{Logical indicating whether to omit the data point}
#' }
"example_qpcr_data"

#' Example Gene Expression Data
#'
#' A dataset containing example differential expression results
#'
#' @format A data frame with 1000 rows and 6 variables:
#' \describe{
#'   \item{gene_symbol}{Gene identifier}
#'   \item{logFC}{Log2 fold change}
#'   \item{P.Value}{Raw p-value}
#'   \item{padj}{Adjusted p-value}
#'   \item{log2FoldChange}{Log2 fold change (alternative format)}
#'   \item{stat}{Test statistic}
#' }
"example_gene_expr"

#' Example Gene List
#'
#' A vector of significant genes from the example gene expression data
#'
#' @format A character vector containing gene symbols
"example_gene_list" 