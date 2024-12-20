% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/perform_kegg_analysis.R
\name{perform_kegg_analysis}
\alias{perform_kegg_analysis}
\title{Perform KEGG Pathway Enrichment Analysis}
\usage{
perform_kegg_analysis(
  results_df,
  title,
  p_cutoff = 0.05,
  fc_cutoff = 0,
  genome = "mmu",
  cat_show = 20
)
}
\arguments{
\item{results_df}{A data frame containing differential expression results}

\item{title}{Character string for plot titles and file names}

\item{p_cutoff}{Numeric value for adjusted p-value cutoff (default: 0.05)}

\item{fc_cutoff}{Numeric value for log2 fold change cutoff (default: 0)}

\item{genome}{Character string specifying the genome ("mmu" for mouse or "hsa" for human)}

\item{cat_show}{Numeric value specifying number of categories to show in the plot (default: 20)}

\item{id_col}{Character string specifying the column name containing gene identifiers}
}
\value{
A data frame containing KEGG enrichment results or a message if no
significant pathways are found
}
\description{
This function performs KEGG pathway enrichment analysis on a set of
differentially expressed genes and creates visualization plots.
}
\details{
The function performs the following steps:
\itemize{
\item Filters significant genes based on p-value and fold change cutoffs
\item Converts gene symbols to ENTREZ IDs
\item Performs KEGG pathway enrichment analysis
\item Creates visualization plots for significant pathways
}
}
\examples{
results_df <- data.frame(
  gene_id = c("Trp53", "Brca1"),
  padj = c(0.01, 0.03),
  log2FoldChange = c(1.5, -2.0)
)
kegg_results <- perform_kegg_analysis(results_df, "Treatment vs Control")

}
