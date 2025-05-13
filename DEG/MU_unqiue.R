library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(AnnotationDbi)
library(ggplot2)

input_dir <- "results/deseq2"
output_dir_genelists <- "results/deg_gene_lists"
enrichment_dir <- "results/enrichment_MU_Unique_Down"
plot_dir <- "results/plots"
if (!dir.exists(output_dir_genelists)) dir.create(output_dir_genelists, recursive = TRUE)
if (!dir.exists(enrichment_dir)) dir.create(enrichment_dir, recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

padj_threshold <- 0.05
log2fc_threshold <- 1 

if (!exists("sig_genes_wt") || !exists("sig_genes_mu")) {
  print("sig_genes_wt or sig_genes_mu not found, attempting to reload from files...")
  file_wt <- file.path(input_dir, "DEG_results_Sor_vs_Ctl_in_WT_with_Symbols.csv")
  file_mu <- file.path(input_dir, "DEG_results_Sor_vs_Ctl_in_MU_with_Symbols.csv")
  if (!file.exists(file_wt) || !file.exists(file_mu)) {
    stop("Error: DEG result files not found, cannot proceed.")
  }
  deg_wt <- read.csv(file_wt)
  deg_mu <- read.csv(file_mu)
  
  filter_sig_genes <- function(df, padj_cutoff, lfc_cutoff) {
    df %>%
      filter(!is.na(padj) & padj < padj_cutoff & abs(log2FoldChange) > lfc_cutoff) %>%
      mutate(Direction = ifelse(log2FoldChange > lfc_cutoff, "Upregulated", "Downregulated")) %>%
      mutate(Identifier = ifelse(is.na(GeneSymbol) | GeneSymbol == "", GeneID, GeneSymbol)) %>%
      dplyr::select(Identifier, GeneID, GeneSymbol, log2FoldChange, padj, Direction)
  }
  sig_genes_wt <- filter_sig_genes(deg_wt, padj_threshold, log2fc_threshold)
  sig_genes_mu <- filter_sig_genes(deg_mu, padj_threshold, log2fc_threshold)
  print("Significant genes reloaded and filtered.")
}

genes_sig_wt_ids <- sig_genes_wt$Identifier
genes_sig_mu_ids <- sig_genes_mu$Identifier

unique_genes_mu_ids <- setdiff(genes_sig_mu_ids, genes_sig_wt_ids)
print(paste("Found", length(unique_genes_mu_ids), "MU-specific responsive genes."))

if (length(unique_genes_mu_ids) == 0) {
  stop("Error: No MU-specific responsive genes found, cannot proceed with analysis.")
}

unique_genes_mu_df <- sig_genes_mu %>%
  filter(Identifier %in% unique_genes_mu_ids)

unique_mu_down <- unique_genes_mu_df %>%
  filter(Direction == "Downregulated") %>%
  arrange(padj)

print(paste("Number of MU-specific downregulated genes:", nrow(unique_mu_down)))

if (nrow(unique_mu_down) > 0) {
  output_down_file <- file.path(output_dir_genelists, "Unique_MU_Sor_Responsive_Down.csv")
  write.csv(unique_mu_down %>% dplyr::select(Identifier, GeneID, GeneSymbol, log2FoldChange, padj),
            file = output_down_file, row.names = FALSE)
  print(paste("MU-specific downregulated gene list saved to:", output_down_file))
} else {
  stop("Error: No MU-specific downregulated genes found, cannot proceed with enrichment analysis.")
}

ensembl_ids_unique_mu_down <- unique_mu_down$GeneID
ensembl_no_version_unique_mu_down <- gsub("\\..*$", "", ensembl_ids_unique_mu_down)

print("Converting MU-specific downregulated gene Ensembl IDs to Entrez IDs...")
entrez_ids_unique_mu_down <- mapIds(org.Hs.eg.db,
                                    keys = ensembl_no_version_unique_mu_down,
                                    column = "ENTREZID",
                                    keytype = "ENSEMBL",
                                    multiVals = "first")

entrez_ids_unique_mu_down_clean <- unique(entrez_ids_unique_mu_down[!is.na(entrez_ids_unique_mu_down)])
print(paste("Obtained", length(entrez_ids_unique_mu_down_clean), "unique Entrez IDs for enrichment analysis."))

if (length(entrez_ids_unique_mu_down_clean) < 10) {
  warning("Warning: Too few MU-specific downregulated genes, enrichment analysis may be insignificant or unstable.")
}

if (length(entrez_ids_unique_mu_down_clean) >= 10) {
  print("Starting GO BP enrichment analysis for MU-specific downregulated genes...")
  ego_bp_mu_unique_down <- enrichGO(gene = entrez_ids_unique_mu_down_clean,
                                    OrgDb = org.Hs.eg.db,
                                    keyType = 'ENTREZID',
                                    ont = "BP",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff = 0.05,
                                    qvalueCutoff = 0.20,
                                    readable = TRUE)
  
  if (!is.null(ego_bp_mu_unique_down) && nrow(ego_bp_mu_unique_down@result) > 0) {
    print("GO BP enrichment results (partial):")
    print(head(ego_bp_mu_unique_down@result[, c("ID", "Description", "p.adjust", "qvalue", "geneID", "Count")], 20))
    write.csv(as.data.frame(ego_bp_mu_unique_down@result),
              file = file.path(enrichment_dir, "GO_BP_enrichment_Unique_MU_Down_DEG.csv"),
              row.names = FALSE)
    print(paste("GO BP results saved to:", enrichment_dir))
    
    go_dotplot_mu_down <- dotplot(ego_bp_mu_unique_down, showCategory = 20, title = "GO BP Enrichment (Unique MU Downregulated DEG)")
    print(go_dotplot_mu_down)
    ggsave(file.path(plot_dir, "GO_BP_dotplot_Unique_MU_Down_DEG.png"), plot = go_dotplot_mu_down, width = 10, height = 8)
    
    go_barplot_mu_down <- barplot(ego_bp_mu_unique_down, showCategory = 15, title = "GO BP Enrichment (Unique MU Downregulated DEG)")
    print(go_barplot_mu_down)
    ggsave(file.path(plot_dir, "GO_BP_barplot_Unique_MU_Down_DEG.png"), plot = go_barplot_mu_down, width = 10, height = 8)
  } else {
    print("No significant GO BP enrichment found for MU-specific downregulated genes.")
  }
  
  print("Starting KEGG pathway enrichment analysis for MU-specific downregulated genes...")
  ekegg_mu_unique_down <- enrichKEGG(gene = entrez_ids_unique_mu_down_clean,
                                     organism = 'hsa',
                                     keyType = 'ncbi-geneid',
                                     pAdjustMethod = "BH",
                                     pvalueCutoff = 0.05,
                                     qvalueCutoff = 0.20)
  
  if (!is.null(ekegg_mu_unique_down) && nrow(ekegg_mu_unique_down@result) > 0) {
    ekegg_mu_unique_down_readable <- setReadable(ekegg_mu_unique_down, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    print("KEGG enrichment results (partial):")
    print(head(ekegg_mu_unique_down_readable@result[, c("ID", "Description", "p.adjust", "qvalue", "geneID", "Count")], 20))
    write.csv(as.data.frame(ekegg_mu_unique_down_readable@result),
              file = file.path(enrichment_dir, "KEGG_enrichment_Unique_MU_Down_DEG.csv"),
              row.names = FALSE)
    print(paste("KEGG results saved to:", enrichment_dir))
    
    kegg_dotplot_mu_down <- dotplot(ekegg_mu_unique_down_readable, showCategory = 20, title = "KEGG Pathway Enrichment (Unique MU Downregulated DEG)")
    print(kegg_dotplot_mu_down)
    ggsave(file.path(plot_dir, "KEGG_dotplot_Unique_MU_Down_DEG.png"), plot = kegg_dotplot_mu_down, width = 10, height = 8)
    
    kegg_barplot_mu_down <- barplot(ekegg_mu_unique_down_readable, showCategory = 15, title = "KEGG Pathway Enrichment (Unique MU Downregulated DEG)")
    print(kegg_barplot_mu_down)
    ggsave(file.path(plot_dir, "KEGG_barplot_Unique_MU_Down_DEG.png"), plot = kegg_barplot_mu_down, width = 10, height = 8)
  } else {
    print("No significant KEGG pathways enriched for MU-specific downregulated genes.")
  }
  
  print("Checking for specific pathways/processes in MU-specific downregulated gene enrichment results...")
  
  if (!is.null(ego_bp_mu_unique_down) && nrow(ego_bp_mu_unique_down@result) > 0) {
    go_results_mu_unique_down_df <- as.data.frame(ego_bp_mu_unique_down@result)
    dna_rep_go <- go_results_mu_unique_down_df[grepl("DNA replication", go_results_mu_unique_down_df$Description, ignore.case = TRUE), ]
    cell_cycle_go <- go_results_mu_unique_down_df[grepl("cell cycle|mitotic", go_results_mu_unique_down_df$Description, ignore.case = TRUE), ]
    mapk_pi3k_go <- go_results_mu_unique_down_df[grepl("MAPK cascade|Ras protein signal|ERK1 and ERK2 cascade|PI3K|phosphatidylinositol 3-kinase signaling", go_results_mu_unique_down_df$Description, ignore.case = TRUE), ]
    
    print("GO BP terms related to DNA replication:")
    print(dna_rep_go[, c("ID", "Description", "p.adjust", "qvalue", "Count")])
    print("GO BP terms related to Cell cycle:")
    print(cell_cycle_go[, c("ID", "Description", "p.adjust", "qvalue", "Count")])
    print("GO BP terms related to MAPK/ERK/PI3K:")
    print(mapk_pi3k_go[, c("ID", "Description", "p.adjust", "qvalue", "Count")])
  }
  
  if (!is.null(ekegg_mu_unique_down) && nrow(ekegg_mu_unique_down@result) > 0) {
    kegg_results_mu_unique_down_df <- as.data.frame(ekegg_mu_unique_down_readable@result)
    dna_rep_kegg <- kegg_results_mu_unique_down_df[grepl("DNA replication", kegg_results_mu_unique_down_df$Description, ignore.case = TRUE), ]
    cell_cycle_kegg <- kegg_results_mu_unique_down_df[grepl("Cell cycle", kegg_results_mu_unique_down_df$Description, ignore.case = TRUE), ]
    mapk_kegg <- kegg_results_mu_unique_down_df[grepl("MAPK signaling pathway", kegg_results_mu_unique_down_df$Description, ignore.case = TRUE), ]
    pi3k_kegg <- kegg_results_mu_unique_down_df[grepl("PI3K-Akt signaling pathway", kegg_results_mu_unique_down_df$Description, ignore.case = TRUE), ]
    
    print("KEGG pathways related to DNA replication:")
    print(dna_rep_kegg[, c("ID", "Description", "p.adjust", "qvalue", "Count")])
    print("KEGG pathways related to Cell cycle:")
    print(cell_cycle_kegg[, c("ID", "Description", "p.adjust", "qvalue", "Count")])
    print("KEGG pathways related to MAPK signaling:")
    print(mapk_kegg[, c("ID", "Description", "p.adjust", "qvalue", "Count")])
    print("KEGG pathways related to PI3K-Akt signaling:")
    print(pi3k_kegg[, c("ID", "Description", "p.adjust", "qvalue", "Count")])
  }
} else {
  print("Too few MU-specific downregulated genes, skipping enrichment analysis.")
}

print("MU-specific downregulated gene enrichment analysis completed.")
