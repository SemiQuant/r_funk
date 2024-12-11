# Rfunk

A collection of useful R functions for bioinformatics analysis, focusing on enhanced visualization tools for differential expression analysis.

## Installation

You can install the development version of r_funk from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("SemiQuant/r_funk")
```

## Usage

### Static Volcano Plots

```r
library(r_funk)

# Create a single volcano plot for a specific comparison
plot <- create_volcano_plot(dds, "treatment_vs_control")

# Create volcano plots for all comparisons in your DESeq2 object
all_plots <- create_all_volcano_plots(dds)
```

### Interactive Volcano Plots

```r
# Create an interactive volcano plot
interactive_plot <- create_interactive_volcano_plot(dds, "treatment_vs_control")

# The interactive plot allows you to:
# - Hover over points to see gene names and statistics
# - Zoom in/out
# - Pan across the plot
# - Select and highlight specific genes
```

### Enrichment Analysis

```r
# Perform GO enrichment analysis
go_results <- run_go_analysis(deg_list, 
                             organism = "human",
                             ont = "BP")  # BP, MF, or CC

# Run KEGG pathway analysis
kegg_results <- run_kegg_analysis(deg_list, 
                                 organism = "hsa")

# Visualize enrichment results
plot_enrichment(go_results, type = "GO")
plot_enrichment(kegg_results, type = "KEGG")
```

### Gene Signature Analysis

```r
# List available gene signatures
signatures <- list_signatures()  # Use built-in signatures
signatures <- list_signatures("path/to/signatures.xlsx")  # Or use custom signatures

# Create heatmaps for gene signatures with optional GSEA
results <- create_signature_heatmaps(
    vst_data,
    dds,
    result_names = c("treatment_vs_control"),
    signature_id = "IFNb – Kotov",  # Optional: specify signature
    padj_threshold = 0.05,
    scale = "row",
    perform_gsea = TRUE  # Enable GSEA analysis
)

# Access results
heatmaps <- results$heatmaps
gsea_results <- results$gsea_results  # Contains GSEA statistics if perform_gsea = TRUE

# GSEA results include:
# - Normalized Enrichment Score (NES)
# - P-values and adjusted p-values
# - Leading edge genes
# - Size of gene sets
```
