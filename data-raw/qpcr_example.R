# Create example datasets for the package

# Example qPCR data
example_qpcr_data <- data.frame(
  Sample = rep(c("Control", "Treatment"), each = 6),
  Target = rep(c("Target_Gene", "Reference_Gene"), each = 3, times = 2),
  Cq = c(
    rnorm(3, 25, 0.3), # Control - Target
    rnorm(3, 20, 0.3), # Control - Reference
    rnorm(3, 24, 0.3), # Treatment - Target
    rnorm(3, 20, 0.3)  # Treatment - Reference
  ),
  Omit = FALSE
)

# Example gene expression data
set.seed(123)
example_gene_expr <- data.frame(
  gene_symbol = paste0("Gene", 1:1000),
  logFC = rnorm(1000),
  P.Value = runif(1000),
  padj = p.adjust(runif(1000), method = "BH"),
  log2FoldChange = rnorm(1000),
  stat = rnorm(1000)
)
rownames(example_gene_expr) <- example_gene_expr$gene_symbol

# Example gene list for enrichment analysis
example_gene_list <- example_gene_expr$gene_symbol[example_gene_expr$padj < 0.05 & abs(example_gene_expr$logFC) > 1]

# Save the data
usethis::use_data(example_qpcr_data, overwrite = TRUE)
usethis::use_data(example_gene_expr, overwrite = TRUE)
usethis::use_data(example_gene_list, overwrite = TRUE) 