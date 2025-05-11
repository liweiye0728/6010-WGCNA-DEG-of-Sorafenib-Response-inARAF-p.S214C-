setwd("D:/MSc/GNBF6010/Manuscript_writing")
# Define output directory
output_dir <- "results/deseq2"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  print(paste("Created directory:", output_dir))
}

# Save results as CSV

# Comparison 1: Mutant vs Wild-type (Control)
print("Saving: Mutant vs Wild-type (Control) results...")
df_mu_vs_wt <- as.data.frame(res_mu_vs_wt_ordered)
df_mu_vs_wt <- tibble::rownames_to_column(df_mu_vs_wt, "GeneID")
write.csv(df_mu_vs_wt,
          file = file.path(output_dir, "DEG_results_Mu_vs_Wt_Control.csv"),
          row.names = FALSE)
print("Saved: DEG_results_Mu_vs_Wt_Control.csv")

# Comparison 2: Sorafenib vs Control (Wild-type)
print("Saving: Sorafenib vs Control (Wild-type) results...")
df_sor_vs_ctl_wt <- as.data.frame(res_sor_vs_ctl_wt_ordered)
df_sor_vs_ctl_wt <- tibble::rownames_to_column(df_sor_vs_ctl_wt, "GeneID")
write.csv(df_sor_vs_ctl_wt,
          file = file.path(output_dir, "DEG_results_Sor_vs_Ctl_in_WT.csv"),
          row.names = FALSE)
print("Saved: DEG_results_Sor_vs_Ctl_in_WT.csv")

# Comparison 3: Sorafenib vs Control (Mutant)
print("Saving: Sorafenib vs Control (Mutant) results...")
df_sor_vs_ctl_mu <- as.data.frame(res_sor_vs_ctl_mu_ordered)
df_sor_vs_ctl_mu <- tibble::rownames_to_column(df_sor_vs_ctl_mu, "GeneID")
write.csv(df_sor_vs_ctl_mu,
          file = file.path(output_dir, "DEG_results_Sor_vs_Ctl_in_MU.csv"),
          row.names = FALSE)
print("Saved: DEG_results_Sor_vs_Ctl_in_MU.csv")

# Comparison 4: Interaction effect
print("Saving: Interaction results...")
df_interaction <- as.data.frame(res_interaction_ordered)
df_interaction <- tibble::rownames_to_column(df_interaction, "GeneID")
write.csv(df_interaction,
          file = file.path(output_dir, "DEG_results_Interaction_Celltype_Treatment.csv"),
          row.names = FALSE)
print("Saved: DEG_results_Interaction_Celltype_Treatment.csv")

# Optional: Sorafenib vs Control (Main effect)
print("Saving: Sorafenib vs Control (Main effect) results...")
df_sor_vs_ctl_main <- as.data.frame(res_sor_vs_ctl_wt_ordered)
df_sor_vs_ctl_main <- tibble::rownames_to_column(df_sor_vs_ctl_main, "GeneID")
write.csv(df_sor_vs_ctl_main,
          file = file.path(output_dir, "DEG_results_Sor_vs_Ctl_MainEffect.csv"),
          row.names = FALSE)
print("Saved: DEG_results_Sor_vs_Ctl_MainEffect.csv")

print("All DEG results successfully saved as CSV files.")