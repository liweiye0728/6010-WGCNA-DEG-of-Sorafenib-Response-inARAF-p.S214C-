setwd("D:/MSc/GNBF6010/Manuscript_writing")
library(tidyverse)

input_dir <- "results/deseq2"
output_dir_genelists <- "results/deg_gene_lists"
if (!dir.exists(output_dir_genelists)) {
  dir.create(output_dir_genelists, recursive = TRUE)
  print(paste("Created directory:", output_dir_genelists))
}

padj_threshold <- 0.05
log2fc_threshold <- 1

process_and_save_sig_genes <- function(file_path, comparison_name, padj_cutoff, lfc_cutoff, output_dir) {
  print(paste("Processing:", comparison_name, "from file:", basename(file_path)))
  
  if (!file.exists(file_path)) {
    warning(paste("File not found, skipping:", file_path))
    return(NULL)
  }
  
  deg_results_df <- read.csv(file_path)
  sig_genes_df <- deg_results_df %>%
    filter(!is.na(padj) & padj < padj_cutoff)
  
  if (nrow(sig_genes_df) == 0) {
    print(paste("No significant genes found in", comparison_name, "(padj <", padj_cutoff, ")"))
    return(NULL)
  }
  
  up_genes_df <- sig_genes_df %>%
    filter(log2FoldChange > lfc_cutoff) %>%
    arrange(padj)
  
  down_genes_df <- sig_genes_df %>%
    filter(log2FoldChange < -lfc_cutoff) %>%
    arrange(padj)
  
  print(paste("Found significant upregulated genes:", nrow(up_genes_df)))
  print(paste("Found significant downregulated genes:", nrow(down_genes_df)))
  
  prepare_output_df <- function(df) {
    df %>%
      mutate(Identifier = ifelse(is.na(GeneSymbol) | GeneSymbol == "", GeneID, GeneSymbol)) %>%
      select(Identifier, GeneID, GeneSymbol, log2FoldChange, padj)
  }
  
  if (nrow(up_genes_df) > 0) {
    output_up_df <- prepare_output_df(up_genes_df)
    output_up_file <- file.path(output_dir, paste0("Significant_Up_", comparison_name, ".csv"))
    write.csv(output_up_df, file = output_up_file, row.names = FALSE)
    print(paste("Saved upregulated gene list to:", output_up_file))
  } else {
    print("No significant upregulated genes to save.")
  }
  
  if (nrow(down_genes_df) > 0) {
    output_down_df <- prepare_output_df(down_genes_df)
    output_down_file <- file.path(output_dir, paste0("Significant_Down_", comparison_name, ".csv"))
    write.csv(output_down_df, file = output_down_file, row.names = FALSE)
    print(paste("Saved downregulated gene list to:", output_down_file))
  } else {
    print("No significant downregulated genes to save.")
  }
}

comparisons <- list(
  "Mu_vs_Wt_Control" = "DEG_results_Mu_vs_Wt_Control_with_Symbols.csv",
  "Sor_vs_Ctl_in_WT" = "DEG_results_Sor_vs_Ctl_in_WT_with_Symbols.csv",
  "Sor_vs_Ctl_in_MU" = "DEG_results_Sor_vs_Ctl_in_MU_with_Symbols.csv",
  "Interaction_Celltype_Treatment" = "DEG_results_Interaction_Celltype_Treatment_with_Symbols.csv"
)

for (comp_name in names(comparisons)) {
  file_name <- comparisons[[comp_name]]
  full_file_path <- file.path(input_dir, file_name)
  process_and_save_sig_genes(full_file_path, comp_name, padj_threshold, log2fc_threshold, output_dir_genelists)
  cat("---\n")
}

print("Extraction and saving of significant gene lists for all comparisons completed.")
