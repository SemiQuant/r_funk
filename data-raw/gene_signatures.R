# Script to prepare gene signatures dataset
gene_signatures <- readxl::read_excel("data-raw/gene_signatures.xlsx")

# Clean up any potential issues (NA values, etc)
gene_signatures <- as.data.frame(gene_signatures)

# Save as internal data
usethis::use_data(gene_signatures, internal = FALSE, overwrite = TRUE)

