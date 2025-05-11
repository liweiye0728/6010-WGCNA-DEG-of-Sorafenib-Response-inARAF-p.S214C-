setwd("D:/MSc/GNBF6010/Manuscript_writing")
# BiocManager::install(c("WGCNA", "preprocessCore", "impute"))
library(DESeq2)
library(WGCNA)
library(tidyverse)
# install.packages("dynamicTreeCut")

print("Performing VST transformation (blind = FALSE)...")
vsd <- vst(dds, blind = FALSE)
print("VST transformation completed.")

datExpr0 <- t(assay(vsd))
print(paste("Initial expression matrix dimensions (samples x genes):", dim(datExpr0)[1], "x", dim(datExpr0)[2]))

print("Checking gene and sample validity (goodSamplesGenes)...")
gsg <- goodSamplesGenes(datExpr0, verbose = 3)

if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) {
    printFlush(paste("Removing low-quality/zero-variance genes:", paste(colnames(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  }
  if (sum(!gsg$goodSamples) > 0) {
    printFlush(paste("Removing low-quality samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  }
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
  print(paste("Expression matrix dimensions after removing problematic genes/samples:", dim(datExpr0)[1], "x", dim(datExpr0)[2]))
} else {
  print("All genes and samples passed goodSamplesGenes check.")
}

print("Performing gene filtering (based on variance)...")
gene_variances <- apply(datExpr0, 2, var)
n_genes_total <- ncol(datExpr0)
n_genes_to_keep <- floor(n_genes_total * 0.75)

print(paste("Total genes:", n_genes_total))
print(paste("Planned to keep top variance genes:", n_genes_to_keep))

keep_gene_indices <- order(gene_variances, decreasing = TRUE)[1:n_genes_to_keep]
datExpr <- datExpr0[, keep_gene_indices]

print(paste("Final expression matrix dimensions after variance filtering (samples x genes):", dim(datExpr)[1], "x", dim(datExpr)[2]))

print("Performing sample clustering to detect outliers...")
sampleTree <- hclust(dist(datExpr), method = "average")

par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

print("WGCNA data preparation completed. Final expression matrix 'datExpr' generated.")
print("Next step is selecting soft thresholding power.")