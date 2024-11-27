#' Detect Outliers in RNA-seq Data with Optional Group-wise Comparisons
#'
#' This function identifies potential outlier samples in RNA-seq data using multiple methods.
#' It can perform analysis either across all samples or within experimental groups.
#'
#' @param dds A DESeqDataSet object containing RNA-seq data
#' @param group_col Optional character specifying the column name in colData(dds) for group information.
#'        This should match one of the column names in colData(dds). For example, if your
#'        experimental design includes a treatment column named "condition", you would set
#'        group_col = "condition". If NULL (default), analysis is performed across all samples.
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
#' # Without groups (analyze all samples together)
#' outliers <- detect_outliers(dds)
#'
#' # With groups (analyze within treatment groups)
#' outliers <- detect_outliers(dds, group_col = "condition")
#'
#' @export
detect_outliers <- function(dds, group_col = NULL, n_mad = 3, cook_threshold = 0.99) {
  # VST transform
  vst_data <- vst(dds, blind = TRUE)
  
  # Initialize vectors to store outliers
  mad_outliers <- character(0)
  correlation_outliers <- character(0)
  
  if (is.null(group_col)) {
    # Original behavior - analyze all samples together
    
    # MAD-based detection
    sample_dists <- dist(t(assay(vst_data)))
    dist_matrix <- as.matrix(sample_dists)
    median_dists <- apply(dist_matrix, 1, median)
    mad_dists <- mad(median_dists)
    mad_outlier_idx <- which(abs(median_dists - median(median_dists)) > n_mad * mad_dists)
    mad_outliers <- names(mad_outlier_idx)
    
    # Inter-sample correlation
    cor_matrix <- cor(assay(vst_data))
    mean_cors <- rowMeans(cor_matrix)
    iqr_outliers <- which(mean_cors < (quantile(mean_cors, 0.25) - 1.5 * IQR(mean_cors)))
    correlation_outliers <- colnames(dds)[iqr_outliers]
    
  } else {
    # Group-wise analysis
    
    # Verify group_col exists in colData
    if (!group_col %in% names(colData(dds))) {
      stop(sprintf("Column '%s' not found in colData(dds). Available columns: %s",
                   group_col, paste(names(colData(dds)), collapse = ", ")))
    }
    
    # Get group information
    groups <- colData(dds)[[group_col]]
    
    # Calculate MAD-based distances within each group
    for (group in unique(groups)) {
      # Get samples for this group
      group_samples <- which(groups == group)
      
      if (length(group_samples) > 2) {  # Need at least 3 samples for meaningful MAD
        # Calculate distances within this group
        group_data <- assay(vst_data)[, group_samples]
        sample_dists <- dist(t(group_data))
        dist_matrix <- as.matrix(sample_dists)
        
        # Calculate median distances and MAD within group
        median_dists <- apply(dist_matrix, 1, median)
        mad_dists <- mad(median_dists)
        
        # Identify outliers within this group
        group_outliers <- which(abs(median_dists - median(median_dists)) > n_mad * mad_dists)
        
        # Add to overall outlier list
        if (length(group_outliers) > 0) {
          mad_outliers <- c(mad_outliers, 
                            colnames(dds)[group_samples[group_outliers]])
        }
        
        # Correlation analysis within group
        cor_matrix <- cor(group_data)
        mean_cors <- rowMeans(cor_matrix)
        iqr_outliers <- which(mean_cors < (quantile(mean_cors, 0.25) - 1.5 * IQR(mean_cors)))
        
        if (length(iqr_outliers) > 0) {
          correlation_outliers <- c(correlation_outliers,
                                    colnames(dds)[group_samples[iqr_outliers]])
        }
      } else {
        warning(sprintf("Group '%s' has fewer than 3 samples (%d). Skipping outlier detection for this group.",
                        group, length(group_samples)))
      }
    }
  }
  
  # Cook's distance calculation (same for both grouped and ungrouped)
  cook_scores <- assays(dds)[["cooks"]]
  cook_outliers <- which(colSums(cook_scores > qchisq(cook_threshold, df = 2)) > 
                           nrow(dds) * 0.01)
  
  # Combined results
  list(
    mad_outliers = mad_outliers,
    cook_outliers = names(cook_outliers),
    correlation_outliers = correlation_outliers,
    all_outliers = unique(c(mad_outliers, 
                            names(cook_outliers), 
                            correlation_outliers)),
    metrics = list(
      cook_scores = cook_scores
    )
  )
}