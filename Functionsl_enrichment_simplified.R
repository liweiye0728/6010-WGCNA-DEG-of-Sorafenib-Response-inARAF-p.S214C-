setwd("D:/MSc/GNBF6010/Manuscript_writing")
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(AnnotationDbi)
library(tidyverse)
library(WGCNA)

target_module <- "blue"
print(paste("Selected target module:", target_module))

all_gene_ids_in_datExpr <- colnames(datExpr)
module_genes <- all_gene_ids_in_datExpr[moduleColors == target_module]
print(paste("Module", target_module, "contains", length(module_genes), "genes"))

print("Checking gene ID format and converting if necessary...")
input_gene_ids <- module_genes
is_ensembl <- grepl("^ENSG", input_gene_ids[1])

if (is_ensembl) {
  print("Input IDs appear to be Ensembl IDs, converting to Entrez IDs...")
  ensembl_no_version <- gsub("\\..*$", "", input_gene_ids)
  entrez_ids <- mapIds(org.Hs.eg.db,
                       keys = ensembl_no_version,
                       column = "ENTREZID",
                       keytype = "ENSEMBL",
                       multiVals = "first")
  entrez_ids <- entrez_ids[!is.na(entrez_ids)]
  entrez_ids <- unique(entrez_ids)
  print(paste("Successfully converted to", length(entrez_ids), "unique Entrez IDs"))
  target_gene_list <- entrez_ids
} else {
  print("Input IDs are not Ensembl format, assuming they are Entrez IDs or Gene Symbols.")
  target_gene_list <- unique(input_gene_ids[!is.na(input_gene_ids)])
}

if (length(target_gene_list) == 0) {
  stop("Error: Target gene list is empty, cannot perform enrichment analysis.")
}

print(paste("Starting GO enrichment analysis for module", target_module, "..."))
ego <- enrichGO(gene = target_gene_list,
                OrgDb = org.Hs.eg.db,
                keyType = 'ENTREZID',
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.20)

if (is.null(ego) || nrow(ego@result) == 0) {
  print(paste("Module", target_module, "has no significant GO BP enrichment terms."))
} else {
  print(paste("Module", target_module, "found", nrow(ego@result), "significant GO BP enrichment terms."))
  print("GO BP enrichment results (partial):")
  print(head(ego@result[, c("ID", "Description", "p.adjust", "qvalue", "geneID", "Count")]))
  
  write.csv(as.data.frame(ego@result),
            file = paste0("results/enrichment/GO_BP_enrichment_", target_module, ".csv"),
            row.names = FALSE)
  
  print("Plotting GO enrichment results...")
  if (nrow(ego@result) > 0) {
    print(dotplot(ego, showCategory = 20, title = paste("GO BP Enrichment -", target_module, "Module")))
    ggsave(paste0("results/plots/GO_BP_dotplot_", target_module, ".png"), width = 10, height = 8)
    
    print(barplot(ego, showCategory = 15, title = paste("GO BP Enrichment -", target_module, "Module")))
    ggsave(paste0("results/plots/GO_BP_barplot_", target_module, ".png"), width = 10, height = 8)
    
    ego_readable <- setReadable(ego, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    if (nrow(ego_readable@result) > 0) {
      print(cnetplot(ego_readable, categorySize="pvalue", foldChange=NULL, showCategory = 5))
      ggsave(paste0("results/plots/GO_BP_cnetplot_", target_module, ".png"), width = 12, height = 10)
    }
  }
}

print(paste("Starting KEGG pathway enrichment analysis for module", target_module, "..."))
ekegg <- enrichKEGG(gene = target_gene_list,
                    organism = 'hsa',
                    keyType = 'ncbi-geneid',
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.20)

if (is.null(ekegg) || nrow(ekegg@result) == 0) {
  print(paste("Module", target_module, "has no significant KEGG pathways."))
} else {
  print(paste("Module", target_module, "found", nrow(ekegg@result), "significant KEGG pathways."))
  ekegg_readable <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  print("KEGG enrichment results (partial):")
  print(head(ekegg_readable@result[, c("ID", "Description", "p.adjust", "qvalue", "geneID", "Count")]))
  
  write.csv(as.data.frame(ekegg_readable@result),
            file = paste0("results/enrichment/KEGG_enrichment_", target_module, ".csv"),
            row.names = FALSE)
  
  print("Plotting KEGG enrichment results...")
  if (nrow(ekegg@result) > 0) {
    print(dotplot(ekegg_readable, showCategory = 20, title = paste("KEGG Pathway Enrichment -", target_module, "Module")))
    ggsave(paste0("results/plots/KEGG_dotplot_", target_module, ".png"), width = 10, height = 8)
    
    print(barplot(ekegg_readable, showCategory = 15, title = paste("KEGG Pathway Enrichment -", target_module, "Module")))
    ggsave(paste0("results/plots/KEGG_barplot_", target_module, ".png"), width = 10, height = 8)
  }
  
  library(pathview)
  kegg_pathway_id <- ekegg@result$ID[1]
  if (!is.na(kegg_pathway_id)) {
    print(paste("Generating KEGG pathway map:", kegg_pathway_id))
    pv.out <- pathview(gene.data = target_gene_list,
                       pathway.id = kegg_pathway_id,
                       species = "hsa",
                       limit = list(gene=max(abs(geneList)), cpd=1))
    print(paste("Pathway map saved as", kegg_pathway_id, ".pathview.png/.xml"))
  }
}

print("Functional enrichment analysis completed.")
print("Next steps: result integration and interpretation.")