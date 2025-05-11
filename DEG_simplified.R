setwd("D:/MSc/GNBF6010/Manuscript_writing")

# --- 1. Load Required R Packages ---
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)


sample_metadata <- read.csv("results/sample_metadata_raw.csv",
                            header = TRUE,
                            check.names = FALSE,
                            row.names = 1)
count_matrix <- read.csv("results/count_matrix.csv",
                         header = TRUE,      
                         check.names = FALSE,  
                         row.names = 1)        
count_matrix <- as.matrix(count_matrix)
sample_metadata$cell_type <- factor(sample_metadata$cell_type)
sample_metadata$treatment <- factor(sample_metadata$treatment)

sample_metadata$cell_type <- factor(sample_metadata$cell_type, levels = c("wt", "mu"))
sample_metadata$treatment <- factor(sample_metadata$treatment, levels = c("Ctl", "Sor"))


# --- 2. Create DESeqDataSet Object ---
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_metadata,
                              design = ~ cell_type + treatment + cell_type:treatment)

# --- 3. Filter Lowly Expressed Genes ---
keep <- rowSums(counts(dds) >= 10) >= 4
dds <- dds[keep,]
print(paste("Number of genes after filtering:", nrow(dds)))

# --- 4. Run DESeq2 Analysis ---
print("Running DESeq()...")
dds <- DESeq(dds)
print("DESeq() completed.")

# View model coefficients
print("Model coefficients (resultsNames):")
resultsNames(dds)

# --- 5. Extract Differential Expression Results ---
padj_threshold <- 0.05
log2fc_threshold <- 1  # Updated threshold for fold change cutoff

# 5a. Comparison 1: ARAF mutant vs wild-type under control condition
print("Extracting Comparison 1: Mutant vs Wild-Type under Control condition")
res_mu_vs_wt <- results(dds, name="cell_type_mu_vs_wt", alpha=padj_threshold)
res_mu_vs_wt_ordered <- res_mu_vs_wt[order(res_mu_vs_wt$padj),]
n_sig_mu_vs_wt <- sum(res_mu_vs_wt_ordered$padj < padj_threshold & abs(res_mu_vs_wt_ordered$log2FoldChange) >= log2fc_threshold, na.rm=TRUE)
print(paste("Significant DEGs (padj <", padj_threshold, ", |log2FC| >=", log2fc_threshold, "):", n_sig_mu_vs_wt))

# 5b. Comparison 2: Sorafenib vs Control in Wild-Type cells
print("Extracting Comparison 2: Sorafenib vs Control in Wild-Type cells")
res_sor_vs_ctl_wt <- results(dds, name="treatment_Sor_vs_Ctl", alpha=padj_threshold)
res_sor_vs_ctl_wt_ordered <- res_sor_vs_ctl_wt[order(res_sor_vs_ctl_wt$padj),]
n_sig_sor_vs_ctl_wt <- sum(res_sor_vs_ctl_wt_ordered$padj < padj_threshold & abs(res_sor_vs_ctl_wt_ordered$log2FoldChange) >= log2fc_threshold, na.rm=TRUE)
print(paste("Significant DEGs (padj <", padj_threshold, ", |log2FC| >=", log2fc_threshold, "):", n_sig_sor_vs_ctl_wt))

# 5c. Comparison 3: Sorafenib vs Control in Mutant cells
print("Extracting Comparison 3: Sorafenib vs Control in Mutant cells")
res_sor_vs_ctl_mu <- results(dds, contrast=list(c("treatment_Sor_vs_Ctl", "cell_typemu.treatmentSor")), alpha=padj_threshold)
res_sor_vs_ctl_mu_ordered <- res_sor_vs_ctl_mu[order(res_sor_vs_ctl_mu$padj),]
n_sig_sor_vs_ctl_mu <- sum(res_sor_vs_ctl_mu_ordered$padj < padj_threshold & abs(res_sor_vs_ctl_mu_ordered$log2FoldChange) >= log2fc_threshold, na.rm=TRUE)
print(paste("Significant DEGs (padj <", padj_threshold, ", |log2FC| >=", log2fc_threshold, "):", n_sig_sor_vs_ctl_mu))

# 5d. Comparison 4: Interaction effect
print("Extracting Comparison 4: Interaction effect (Mutant vs Wild-Type under treatment)")
res_interaction <- results(dds, name="cell_typemu.treatmentSor", alpha=padj_threshold)
res_interaction_ordered <- res_interaction[order(res_interaction$padj),]
n_sig_interaction <- sum(res_interaction_ordered$padj < padj_threshold & abs(res_interaction_ordered$log2FoldChange) >= log2fc_threshold, na.rm=TRUE)
print(paste("Significant interaction genes (padj <", padj_threshold, ", |log2FC| >=", log2fc_threshold, "):", n_sig_interaction))

# --- 6. Visualization ---

# 6a. Volcano Plot for Interaction Effect
print("Plotting Volcano Plot for Interaction results...")
res_interaction_df <- as.data.frame(res_interaction_ordered) %>%
  rownames_to_column(var = "gene_id") %>%
  mutate(
    significance = case_when(
      padj < padj_threshold & log2FoldChange >= log2fc_threshold  ~ "Upregulated",
      padj < padj_threshold & log2FoldChange <= -log2fc_threshold ~ "Downregulated",
      TRUE                                                       ~ "Not Significant"
    ),
    significance = factor(significance, levels = c("Upregulated", "Downregulated", "Not Significant"))
  )

ggplot(res_interaction_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significance), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  labs(title = "Volcano Plot: Interaction (Cell Type * Treatment)",
       x = "Log2 Fold Change (Interaction Effect)",
       y = "-Log10 Adjusted P-value") +
  theme_bw() +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "grey50")

# 6b. Heatmap of Significant Interaction Genes
sig_interaction_genes <- subset(res_interaction_ordered, padj < 0.01 & abs(log2FoldChange) >= log2fc_threshold)
sig_gene_ids <- rownames(sig_interaction_genes)

if(length(sig_gene_ids) > 1) {
  print(paste("Plotting heatmap for", length(sig_gene_ids), "significant interaction genes..."))
  vsd <- vst(dds, blind = FALSE)
  vst_counts_sig <- assay(vsd)[sig_gene_ids, ]
  df_col <- as.data.frame(colData(dds)[,c("cell_type","treatment")])
  
  pheatmap(vst_counts_sig,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_rownames = (length(sig_gene_ids) < 50),
           show_colnames = TRUE,
           annotation_col = df_col,
           scale = "row",
           main = "Heatmap of Top Interaction Genes (VST values, row-scaled)",
           fontsize_row = 8,
           border_color = NA)
} else {
  print("No sufficient significant interaction genes for heatmap plotting.")
}

print("DEG analysis completed.")
