# Script to prepare gene signatures dataset
gene_signatures <- readxl::read_excel("data-raw/gene_signatures.xlsx")

# Save as internal data
usethis::use_data(gene_signatures, internal = FALSE, overwrite = TRUE) 