% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gene_sigs.R
\name{create_signature_heatmaps}
\alias{create_signature_heatmaps}
\title{Create Heatmaps from Gene Signatures}
\usage{
create_signature_heatmaps(
  vst_data,
  dds,
  result_names,
  signatures_file = NULL,
  signature_id = NULL,
  n_signatures = NULL,
  padj_threshold = 0.05,
  scale = "row",
  condition_column = "Condition",
  lfc_threshold = NULL,
  prefilter_padj = NULL,
  prefilter_lfc = NULL,
  perform_gsea = FALSE,
  gsea_nperm = 1000
)
}
\arguments{
\item{vst_data}{A DESeq2 vst transformed object containing expression data}

\item{dds}{A DESeq2 object containing experimental design information}

\item{result_names}{Character vector of result names to use for significance annotation}

\item{signatures_file}{Optional path to an Excel file containing gene signatures.
If NULL (default), uses the package's built-in gene signatures dataset}

\item{signature_id}{Optional specific signature ID to plot (default: NULL, plots all)}

\item{n_signatures}{Optional number of signatures to plot (default: NULL, plots all)}

\item{padj_threshold}{Adjusted p-value threshold for significance (default: 0.05)}

\item{scale}{Character indicating if the values should be centered and scaled in either
the "row" direction or "column" direction, or "none" (default: "row")}

\item{condition_column}{Name of the condition column in colData(dds) (default: "Condition")}

\item{lfc_threshold}{Log2 fold change threshold for filtering (default: NULL)}

\item{prefilter_padj}{Adjusted p-value threshold for pre-filtering (default: NULL)}

\item{prefilter_lfc}{Log2 fold change threshold for pre-filtering (default: NULL)}

\item{perform_gsea}{Logical indicating whether to perform GSEA analysis (default: FALSE)}

\item{gsea_nperm}{Deprecated. Number of permutations is now automatically determined by fgseaMultilevel}
}
\value{
A list of pheatmap objects
}
\description{
This function creates heatmaps from gene signatures, allowing for visualization of
expression patterns across different conditions with significance annotations.
}
\examples{
\dontrun{
create_signature_heatmaps(vst_data, dds, 
                         result_names = c("condition_treated_vs_control"),
                         "path/to/signatures.xlsx")
}
}
