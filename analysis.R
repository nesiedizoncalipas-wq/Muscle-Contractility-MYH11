# Run installation only once on your computer:
install.packages(c("rentrez", "ggplot2"))

library(rentrez)
library(ggplot2)

# CHANGE THESE FOUR VALUES to your assignment:
gene_symbol <- "MYH11"
transcript_accession <- "NM_002474.3"
clinvar_id <- "1467492"
variant_protein <- "" # Use "" for splice variants

for (x in c("data", "results", "figures")) {
  dir.create(x, showWarnings = FALSE)
}