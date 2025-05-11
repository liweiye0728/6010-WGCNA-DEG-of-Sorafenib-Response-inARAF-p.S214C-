setwd("D:/MSc/GNBF6010/Manuscript_writing")
# Load required package
library(tidyverse)

# Define file path and sample names
featurecounts_file <- "results/counts_new/featureCounts_gene_counts.txt"
sample_names <- c("ARAFmu_Ctl_1", "ARAFmu_Ctl_2", "ARAFmu_sor_1", "ARAFmu_sor_2",
                  "ARAFwt_Ctl_1", "ARAFwt_Ctl_2", "ARAFwt_sor_1", "ARAFwt_sor_2")

# Read featureCounts output
fc_results <- read.table(featurecounts_file, header = TRUE, sep = "\t", comment.char = "#")
print("Raw featureCounts head:")
head(fc_results)
print(paste("Raw data dimensions:", dim(fc_results)[1], "rows", dim(fc_results)[2], "columns"))

# Extract and clean count matrix
count_data <- fc_results[, 7:ncol(fc_results)]
cleaned_colnames <- gsub("^results\\.alignment_new\\.", "", colnames(count_data))
cleaned_colnames <- gsub("_sorted.bam$", "", basename(colnames(count_data)))
colnames(count_data) <- cleaned_colnames
rownames(count_data) <- fc_results$Geneid

# Check for duplicate gene IDs
if(any(duplicated(rownames(count_data)))) {
  warning("Duplicate gene IDs found in count matrix!")
}

# Verify and order columns
print("Cleaned column names:")
print(colnames(count_data))
if(!all(colnames(count_data) %in% sample_names) || !all(sample_names %in% colnames(count_data))) {
  warning("Column names do not fully match sample_names!")
  print("Defined sample names:")
  print(sample_names)
} else {
  count_data <- count_data[, sample_names]
  print("Columns ordered by sample_names.")
}

# Convert to matrix
count_matrix <- as.matrix(count_data)
print("Final count matrix head:")
head(count_matrix)
print(paste("Final count matrix dimensions:", dim(count_matrix)[1], "rows (genes)", dim(count_matrix)[2], "columns (samples)"))

# Create sample metadata
sample_metadata <- data.frame(
  sample_id = sample_names,
  cell_type = factor(gsub("_(Ctl|sor)_\\d$", "", gsub("^ARAF", "", sample_names))),
  treatment = factor(gsub("^.*_(Ctl|sor)_\\d$", "\\1", sample_names)),
  row.names = sample_names
)

# Check factor levels
print("Cell type factor levels:")
levels(sample_metadata$cell_type)
print("Treatment factor levels:")
levels(sample_metadata$treatment)
print("Sample metadata:")
print(sample_metadata)

# Validate consistency
if(all(colnames(count_matrix) == rownames(sample_metadata))) {
  print("Validation successful: Count matrix columns match metadata rows.")
} else {
  stop("Error: Count matrix columns do not match metadata rows!")
}