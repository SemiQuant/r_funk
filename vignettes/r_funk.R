## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE,
  fig.width = 7,
  fig.height = 5
)


## -----------------------------------------------------------------------------
require(Rfunk)


## -----------------------------------------------------------------------------
# Create minimal example data
test_data <- data.frame(
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

# Run analysis
results <- analyze_pcr(
  data = test_data,
  ct_col = "Cq",
  sample_col = "Sample",
  target_col = "Target",
  target = "Target_Gene",
  reference_sample = "Control",
  reference_target = "Reference_Gene"
)


## -----------------------------------------------------------------------------
# Create example data with outliers
set.seed(123)
data <- data.frame(
  values = c(rnorm(97, mean = 0, sd = 1), c(5, 6, 7))  # 97 normal points and 3 outliers
)

# Detect outliers
outliers <- detect_outliers(data$values)
print(paste("Number of outliers detected:", sum(outliers)))


## -----------------------------------------------------------------------------
# Example gene expression data
data("example_gene_expr")
data("example_gene_list")

# Perform GO analysis
go_results <- perform_go_analysis(
  gene_list = example_gene_list,
  organism = "hsapiens",
  ont = "BP"
)

# Perform KEGG analysis
kegg_results <- perform_kegg_analysis(
  gene_list = example_gene_list,
  organism = "hsa"
)

# Perform MSigDB analysis
msig_results <- perform_mSig_analysis(
  gene_list = example_gene_list,
  collection = "H"
)


## -----------------------------------------------------------------------------
# Create a basic volcano plot
volcano_plot <- create_volcano_plot(
  data = example_gene_expr,
  x_col = "logFC",
  y_col = "P.Value",
  labels = rownames(example_gene_expr)
)
print(volcano_plot)

# Create an interactive volcano plot
interactive_volcano <- create_interactive_volcano_plot(
  data = example_gene_expr,
  x_col = "logFC",
  y_col = "P.Value",
  labels = rownames(example_gene_expr)
)
interactive_volcano

# Create an interactive datatable
dt <- create_interactive_datatable(
  data = example_gene_expr,
  rownames = TRUE,
  filter = "top",
  download_buttons = TRUE
)
dt


## -----------------------------------------------------------------------------
sessionInfo()

