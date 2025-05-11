setwd("D:/MSc/GNBF6010/Manuscript_writing")
# Objective: Exploratory Data Analysis
# Before conducting complex statistical analysis, preliminarily examine the data quality, the relationships between samples, batch effects, etc

sample_metadata <- read.csv("results/sample_metadata_raw.csv",
                            header = TRUE,
                            check.names = FALSE,
                            row.names = 1)
count_matrix <- read.csv("results/count_matrix.csv",
                            header = TRUE,      
                            check.names = FALSE,  
                            row.names = 1)        
count_matrix <- as.matrix(count_matrix)
# # Validate consistency
# if(all(colnames(count_matrix) == rownames(sample_metadata))) {
#   print("Validation successful: Count matrix columns match metadata rows.")
# } else {
#   stop("Error: Count matrix columns do not match metadata rows!")
# }
# BiocManager::install("DESeq2")
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_metadata,
                              design = ~ 1)

vsd <- vst(dds, blind = TRUE)
print("VST transformed data (first few rows):")
head(assay(vsd))

pcaData <- plotPCA(vsd, intgroup = c("cell_type", "treatment"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

print("Plotting PCA...")
ggplot(pcaData, aes(x = PC1, y = PC2, color = cell_type, shape = treatment)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA Plot of Samples (VST transformed data)") +
  theme_bw() +
  scale_color_brewer(palette="Set1") +
  theme(plot.title = element_text(hjust = 0.5))

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
# rownames(sampleDistMatrix) <- paste(vsd$cell_type, vsd$treatment, vsd$sample_id, sep="-")
rownames(sampleDistMatrix) <- paste(vsd$sample_id)
colnames(sampleDistMatrix) <- NULL

colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

print("Plotting sample distance heatmap...")
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         main = "Heatmap of Sample Distances (VST transformed data)",
         fontsize_row = 8,
         fontsize_col = 8)

print("EDA analysis completed. Please check the generated PCA plot and heatmap.")