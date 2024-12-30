#' Gene Signatures Dataset
#'
#' A dataset containing gene signatures information
#'
#' @format A data frame containing gene signatures
#' @source Internal gene signatures compilation
"gene_signatures" 

#' Example qPCR Data
#'
#' A dataset containing example qPCR results for demonstrating the analyze_pcr function.
#' The data includes Ct values for a target gene and a reference gene across multiple samples.
#'
#' @format A data frame with 24 rows and 5 columns:
#' \describe{
#'   \item{Sample}{Sample identifier}
#'   \item{Target}{Target gene or reference gene name}
#'   \item{Cq}{Ct value from qPCR}
#'   \item{Omit}{Logical indicating whether to omit the sample}
#'   \item{Replicate}{Technical replicate number}
#' }
#'
#' @examples
#' data(qpcr_example)
#' head(qpcr_example)
"qpcr_example" 