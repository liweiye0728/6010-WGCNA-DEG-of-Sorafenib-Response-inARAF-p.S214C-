setwd("D:/MSc/GNBF6010/Manuscript_writing")
# Load required R packages
library(WGCNA)
library(tidyverse)

# Enable multi-threading
allowWGCNAThreads()

# Check data dimensions
print(paste("Expression matrix dimensions (samples x genes):", nrow(datExpr), "x", ncol(datExpr)))
print("Sample metadata:")
print(sample_metadata)

# Select soft thresholding power
print("Selecting soft threshold (pickSoftThreshold)...")
powers <- c(c(1:10), seq(from = 12, to = 30, by = 2))
sft <- pickSoftThreshold(datExpr,
                         powerVector = powers,
                         verbose = 5,
                         networkType = "signed hybrid")

# Plot results for power selection
par(mfrow = c(1, 2))
cex1 = 0.9
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
     main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.85, col = "red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
par(mfrow = c(1, 1))

# Set soft power
softPower <- sft$powerEstimate
if (is.na(softPower)) {
  softPower <- 9
  print(paste("Automatic estimation failed, manually selected power =", softPower))
} else {
  print(paste("Using automatically estimated power =", softPower))
}

# Construct network and identify modules
print("Constructing network and identifying modules (blockwiseModules)...")
net_type <- "signed hybrid"
tom_type <- "signed"
merge_cut_height <- 0.25
min_module_size <- 30
net <- blockwiseModules(datExpr,
                        power = softPower,
                        networkType = net_type,
                        TOMType = tom_type,
                        minModuleSize = min_module_size,
                        reassignThreshold = 0,
                        mergeCutHeight = merge_cut_height,
                        numericLabels = TRUE,
                        pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "results/WGCNA_TOM",
                        verbose = 3,
                        maxBlockSize = ncol(datExpr) + 100)
print("blockwiseModules completed.")
print("Module counts and sizes:")
print(table(net$colors))

# Convert numeric labels to colors
moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)
print("Module color labels example:")
head(moduleColors)

# Visualize gene dendrogram and module colors
print("Plotting gene dendrogram and module colors...")
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Relate modules to traits
print("Relating modules to traits...")
nSamples <- nrow(datExpr)
traitData <- data.frame(
  is_mutant = ifelse(sample_metadata[rownames(datExpr), "cell_type"] == "mu", 1, 0),
  is_treated = ifelse(sample_metadata[rownames(datExpr), "treatment"] == "Sor", 1, 0)
)
rownames(traitData) <- rownames(datExpr)
print("Prepared trait data (traitData):")
head(traitData)

# Calculate module eigengenes (MEs)
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

# Calculate module-trait correlations
moduleTraitCor <- cor(MEs, traitData, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# Visualize module-trait correlation heatmap
print("Plotting module-trait correlation heatmap...")
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traitData),
               yLabels = names(MEs),
               ySymbols = gsub("ME", "", names(MEs)),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               main = "Module-trait relationships")

# Identify hub genes in a target module
target_module <- "blue"
target_trait <- "is_mutant"
print(paste("Analyzing module:", target_module, "with trait:", target_trait, "..."))

# Calculate module membership (kME)
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
modNames <- substring(names(MEs), 3)
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")

# Calculate gene significance (GS) for the target trait
geneTraitSignificance <- as.data.frame(cor(datExpr, traitData[, target_trait], use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", target_trait, sep = "")
names(GSPvalue) <- paste("p.GS.", target_trait, sep = "")

# Extract genes in the target module
module_genes_indices <- which(moduleColors == target_module)
module_gene_ids <- colnames(datExpr)[module_genes_indices]

# Get MM and GS values for module genes
module_MM_values <- geneModuleMembership[module_genes_indices, paste0("MM", target_module)]
module_GS_values <- geneTraitSignificance[module_genes_indices, paste0("GS.", target_trait)]

# Identify hub genes
hub_gene_order <- order(abs(module_MM_values), decreasing = TRUE)
hub_genes <- module_gene_ids[hub_gene_order]
print(paste("Hub genes in module", target_module, "(sorted by kME):"))
head(hub_genes, 20)

# WGCNA analysis completed
print("WGCNA analysis main steps completed.")
print("Next steps: functional enrichment analysis and result integration.")
