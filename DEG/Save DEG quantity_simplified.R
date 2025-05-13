setwd("D:/MSc/GNBF6010/Manuscript_writing")
library(tidyverse)

padj_threshold <- 0.05
log2fc_threshold <- 1

calculate_sig_counts <- function(result_object, padj_cutoff, lfc_cutoff) {
  res_df <- as.data.frame(result_object)
  res_sig <- subset(res_df, !is.na(padj) & padj < padj_cutoff)
  res_sig_up <- subset(res_sig, log2FoldChange > lfc_cutoff)
  res_sig_down <- subset(res_sig, log2FoldChange < -lfc_cutoff)
  return(list(
    Total_Significant = nrow(res_sig),
    Upregulated = nrow(res_sig_up),
    Downregulated = nrow(res_sig_down)
  ))
}

print("Calculating significant gene counts...")

counts_mu_vs_wt <- calculate_sig_counts(res_mu_vs_wt_ordered, padj_threshold, log2fc_threshold)
print("Mutant vs Wildtype (Control):")
print(counts_mu_vs_wt)

counts_sor_vs_ctl_wt <- calculate_sig_counts(res_sor_vs_ctl_wt_ordered, padj_threshold, log2fc_threshold)
print("Drug vs Control (WT):")
print(counts_sor_vs_ctl_wt)

counts_sor_vs_ctl_mu <- calculate_sig_counts(res_sor_vs_ctl_mu_ordered, padj_threshold, log2fc_threshold)
print("Drug vs Control (MU):")
print(counts_sor_vs_ctl_mu)

counts_interaction <- calculate_sig_counts(res_interaction_ordered, padj_threshold, log2fc_threshold)
print("Interaction:")
print(counts_interaction)

print("Creating summary dataframe...")
summary_df <- data.frame(
  Comparison = c("Mu_vs_Wt_Control",
                 "Sor_vs_Ctl_in_WT",
                 "Sor_vs_Ctl_in_MU",
                 "Interaction"),
  Total_Significant = c(counts_mu_vs_wt$Total_Significant,
                        counts_sor_vs_ctl_wt$Total_Significant,
                        counts_sor_vs_ctl_mu$Total_Significant,
                        counts_interaction$Total_Significant),
  Upregulated = c(counts_mu_vs_wt$Upregulated,
                  counts_sor_vs_ctl_wt$Upregulated,
                  counts_sor_vs_ctl_mu$Upregulated,
                  counts_interaction$Upregulated),
  Downregulated = c(counts_mu_vs_wt$Downregulated,
                    counts_sor_vs_ctl_wt$Downregulated,
                    counts_sor_vs_ctl_mu$Downregulated,
                    counts_interaction$Downregulated)
)

print("DEG count summary:")
print(summary_df)

output_dir <- "results/deg_summary"
output_file <- file.path(output_dir, "DEG_summary_counts.csv")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  print(paste("Created directory:", output_dir))
}

print(paste("Saving summary results to:", output_file))
write.csv(summary_df, file = output_file, row.names = FALSE)

print("DEG count summary successfully saved as CSV file.")
