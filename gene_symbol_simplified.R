setwd("D:/MSc/GNBF6010/Manuscript_writing")
# Load required R packages
library(AnnotationDbi)
library(org.Hs.eg.db)
library(tidyverse)
library(dplyr)
# Load result data frame
output_dir <- "results/deseq2/"
# deg_result_file <- file.path(output_dir, "DEG_results_Interaction_Celltype_Treatment.csv")
# deg_result_file <- file.path(output_dir, "DEG_results_Mu_vs_Wt_Control.csv")
# deg_result_file <- file.path(output_dir, "")
# deg_result_file <- file.path(output_dir, "DEG_results_Sor_vs_Ctl_in_MU.csv")
# deg_result_file <- file.path(output_dir, "DEG_results_Sor_vs_Ctl_in_WT.csv")
deg_result_file <- file.path(output_dir, "DEG_results_Sor_vs_Ctl_MainEffect.csv")
if (file.exists(deg_result_file)) {
  df_interaction <- read.csv(deg_result_file)
  print(paste("Loaded:", deg_result_file))
} else {
  if (!exists("df_interaction")) {
    stop("Error: df_interaction not found. Please load or generate it first.")
  }
  print("Using existing df_interaction from R environment.")
}

print("Original data frame (preview):")
head(df_interaction)

# Prepare Ensembl IDs for conversion
ensembl_ids_with_version <- df_interaction$GeneID
ensembl_ids_no_version <- gsub("\\..*$", "", ensembl_ids_with_version)
print("Ensembl ID without version (preview):")
head(ensembl_ids_no_version)

# Convert Ensembl IDs to Gene Symbols
print("Starting ID conversion (Ensembl to Symbol)...")
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = ensembl_ids_no_version,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")
print("Converted Gene Symbol (preview):")
head(gene_symbols)

# Add Gene Symbol to data frame
df_interaction$GeneSymbol <- gene_symbols
df_interaction <- df_interaction %>% dplyr::select(GeneID, GeneSymbol, everything())
print("Data frame with Gene Symbol (preview):")
head(df_interaction)

# Save data frame with Gene Symbol
output_file_with_symbol <- file.path(output_dir, "DEG_results_Sor_vs_Ctl_MainEffect_with_Symbols.csv")
print(paste("Saving results with Gene Symbol to:", output_file_with_symbol))
write.csv(df_interaction,
          file = output_file_with_symbol,
          row.names = FALSE)
print("File saved.")