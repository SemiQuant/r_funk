#' Detect Outliers in RNA-seq Data
#'
#' This function identifies potential outlier samples in RNA-seq data using multiple methods:
#' MAD-based distances, Cook's distances, and inter-sample correlations.
#'
#' @param dds A DESeqDataSet object containing RNA-seq data
#' @param n_mad Numeric value specifying the number of MADs for the distance threshold (default: 3)
#' @param cook_threshold Numeric value for Cook's distance quantile threshold (default: 0.99)
#'
#' @return A list containing:
#' \itemize{
#'   \item mad_outliers: Sample names identified as outliers based on MAD distances
#'   \item cook_outliers: Sample names identified as outliers based on Cook's distances
#'   \item correlation_outliers: Sample names identified as outliers based on correlation
#'   \item all_outliers: Combined unique list of all identified outliers
#'   \item metrics: List of computed metrics used for outlier detection
#' }
#'
#' @import DESeq2
#' @import stats
#'
#' @examples
#' # Assuming you have a DESeqDataSet object 'dds'
#' outliers <- detect_outliers(dds)
#' print(outliers$all_outliers)
#'
#' @export
detect_outliers <- function(dds, n_mad = 3, cook_threshold = 0.99) {
  # VST transform
  vst_data <- vst(dds, blind = TRUE)
  
  # MAD-based detection
  sample_dists <- dist(t(assay(vst_data)))
  dist_matrix <- as.matrix(sample_dists)
  median_dists <- apply(dist_matrix, 1, median)
  mad_dists <- mad(median_dists)
  mad_outliers <- which(abs(median_dists - median(median_dists)) > n_mad * mad_dists)
  
  # Cook's distance
  cook_scores <- assays(dds)[["cooks"]]
  cook_outliers <- which(colSums(cook_scores > qchisq(cook_threshold, df = 2)) > nrow(dds) * 0.01)
  
  # Inter-sample correlation
  cor_matrix <- cor(assay(vst_data))
  mean_cors <- rowMeans(cor_matrix)
  iqr_outliers <- which(mean_cors < (quantile(mean_cors, 0.25) - 1.5 * IQR(mean_cors)))
  
  # Combined results
  list(
    mad_outliers = names(mad_outliers),
    cook_outliers = names(cook_outliers),
    correlation_outliers = colnames(dds)[iqr_outliers],
    all_outliers = unique(c(names(mad_outliers), 
                            names(cook_outliers), 
                            colnames(dds)[iqr_outliers])),
    metrics = list(
      mad_distances = median_dists,
      mad_threshold = median(median_dists) + n_mad * mad_dists,
      cook_scores = cook_scores,
      mean_correlations = mean_cors
    )
  )
}