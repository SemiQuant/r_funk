# r_funk

## Installation

You can install the development version of Rfunk from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("SemiQuant/Rfunk")
```

## Usage

```r
require(r_funk)

# Create a volcano plot
plot <- create_volcano_plot(dds, "treatment_vs_control")

# Create volcano plots for all comparisons
all_plots <- create_all_volcano_plots(dds)
```
