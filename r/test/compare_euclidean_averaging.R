

# Load all TensorOmics helper functions
source("r/tensoromics_functions.R")
library(readr)
library(ggplot2)

# This script compares two normalization strategies for gene expression data:
# 1. Mean-after-log: Log2(x+1) transformation is applied before averaging replicates.
# 2. Mean-before-log: Replicates are averaged before Log2(x+1) transformation.
# For each method, the Euclidean distance from each gene to its orthogroup centroid is calculated.
# The plot shows the density distribution of these distances for both methods, allowing you to visually compare how normalization affects gene-to-centroid distances.


# === Example of full normalization pipeline ===

# Load raw expression data
input_file <- "material/kallisto_sex_data.tsv"  # Update the input file path if necessary
output_file1 <- "results/normalization_after_log.tsv"  
output_file2 <- "results/normalization_before_log.tsv"  

df <- read.table(input_file, header = TRUE, sep = "\t")
# Prepare matrix for processing (removing the gene ID column)
gene_ids <- df[[1]]              # Save gene identifiers
col_names <- colnames(df)[-1]    
df_matrix <- as.matrix(df[,-1])   # Convert the expression values into a matrix

# === Diagnose data quality before normalization ===
cat("Analyzing data quality before normalization...\n")
diagnostics <- tox_diagnose_data_quality(df_matrix)

# Clean data if there are problems
if (diagnostics$problems$na_count > 0 || diagnostics$problems$inf_count > 0 || diagnostics$problems$nan_count > 0) {
  cat("Data contains problematic values. Cleaning data...\n")
  
  # Clean the data using our cleaning function
  df_matrix_clean <- tox_clean_data_for_normalization(
    df_matrix,
    remove_all_zero_genes = TRUE,
    na_strategy = "impute_mean",  # Impute NA with gene means instead of removing genes
    min_expression_threshold = 0.0,  # Don't set a threshold for TPM data
    convert_small_to_zero = FALSE    # Preserve small TPM values
  )
  
  # Update gene_ids to match cleaned matrix
  gene_ids <- gene_ids[1:nrow(df_matrix_clean)]
  
  cat("Data cleaning completed. Proceeding with normalization...\n")
} else {
  df_matrix_clean <- df_matrix
  cat("Data is clean. Proceeding with normalization...\n")
}


# === Apply normalization steps: mean-after-log ===
normalized_matrix_std <- tox_normalize_by_std_dev(df_matrix_clean)    # Normalize by standard deviation
normalized_matrix_qtl <- tox_quantile_normalization(normalized_matrix_std)  # Apply quantile normalization
normalized_matrix_log <- tox_log2_transformation(normalized_matrix_qtl)     # Log2(x+1) transformation
averaged_df <- tox_calculate_tissue_averages(normalized_matrix_log)         # Average replicates by tissue

# Create final normalized dataframe with averages (after log)
normalized_df_after <- data.frame(gene_id = gene_ids, averaged_df)
write.table(normalized_df_after, output_file1, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Normalization pipeline completed successfully!\n")
cat("Output saved to:", output_file1, "\n")
cat("Final dimensions:", nrow(normalized_df_after), "genes x", ncol(normalized_df_after)-1, "conditions\n")

# === Apply normalization steps: mean-before-log ===
averaged_df <- tox_calculate_tissue_averages(normalized_matrix_qtl)         # Average replicates by tissue
# Ensure averaged_df is a numeric matrix
if (is.list(averaged_df)) {
  averaged_df <- do.call(cbind, averaged_df)
}
averaged_df <- as.matrix(averaged_df)
normalized_matrix_log <- tox_log2_transformation(averaged_df)     # Log2(x+1) transformation

normalized_df_before <- data.frame(gene_id = gene_ids, normalized_matrix_log)
write.table(normalized_df_before, output_file2, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Normalization pipeline completed successfully!\n")
cat("Output saved to:", output_file2, "\n")
cat("Final dimensions:", nrow(normalized_df_before), "genes x", ncol(normalized_df_before)-1, "conditions\n")


# === Calculate Euclidean distance to centroids for both normalization methods ===
centroids_file <- "material/centroids_orthologs_only.tsv"
orthogroups_file <- "material/filtered_families.tsv"

# Read centroids and prepare matrix
centroids_df <- read_tsv(centroids_file, show_col_types = FALSE)
centroids_matrix <- as.matrix(centroids_df[,-1])
centroids_matrix <- t(centroids_matrix) # genes x dimensions

# Read gene-to-family mapping
orthogroups_df <- read_tsv(orthogroups_file, show_col_types = FALSE)
gene_ids_centroids <- centroids_df[[1]]
orthogroup_to_index <- setNames(1:nrow(centroids_df), gene_ids_centroids)

# Create mapping gene_id -> family_index for genes in the data
gene_to_family <- rep(0, length(gene_ids))
names(gene_to_family) <- gene_ids
species_columns <- colnames(orthogroups_df)[colnames(orthogroups_df) != "Orthogroup"]
for (i in 1:nrow(orthogroups_df)) {
  orthogroup_id <- orthogroups_df$Orthogroup[i]
  if (!orthogroup_id %in% names(orthogroup_to_index)) next
  family_index <- orthogroup_to_index[orthogroup_id]
  for (species in species_columns) {
    target_genes <- orthogroups_df[[species]][i]
    if (is.na(target_genes) || target_genes == "") next
    genes_list <- trimws(strsplit(target_genes, ",")[[1]])
    for (gene_id in genes_list) {
      if (gene_id %in% gene_ids && gene_to_family[gene_id] == 0) {
        gene_to_family[gene_id] <- family_index
      }
    }
  }
}
gene_to_family <- as.integer(gene_to_family)

# Prepare matrices for Fortran (genes x conditions)
mat_after <- as.matrix(normalized_df_after[,-1])
mat_before <- as.matrix(normalized_df_before[,-1])

# Transpose for Fortran column-major format
mat_after_t <- t(mat_after)
mat_before_t <- t(mat_before)
centroids_t <- centroids_matrix

d <- nrow(mat_after_t)

# Calculate distances to centroids
distances_after <- tox_distance_to_centroid(as.vector(mat_after_t), as.vector(centroids_t), gene_to_family, d)
distances_before <- tox_distance_to_centroid(as.vector(mat_before_t), as.vector(centroids_t), gene_to_family, d)

# Filter genes with assigned families
has_family <- gene_to_family > 0
distances_after <- distances_after[has_family]
distances_before <- distances_before[has_family]
gene_ids_fam <- gene_ids[has_family]

# Plot comparison of distances to centroids
# df_dist <- data.frame(
#   gene_id = rep(gene_ids_fam, 2),
#   distance = c(distances_after, distances_before),
#   method = rep(c("Mean-after-log", "Mean-before-log"), each = length(gene_ids_fam))
# )

# p2 <- ggplot(df_dist, aes(x = distance, fill = method)) +
#   geom_density(alpha = 0.5) +
#   labs(title = "Euclidean Distance to Centroid: Mean-after-log vs Mean-before-log",
#        x = "Distance to Centroid", y = "Density") +
#   theme_minimal()

# # Save the plot in results/
# plot_file <- "results/compare_centroid_distances.png"
# ggsave(plot_file, plot = p2, width = 7, height = 5)
# cat("Plot saved to:", plot_file, "\n")
# print(p2)


# Summary statistics
cat("Distance to centroids (mean-after-log): mean=", mean(distances_after), " sd=", sd(distances_after), "\n")
cat("Distance to centroids (mean-before-log): mean=", mean(distances_before), " sd=", sd(distances_before), "\n")

# === Gen-by-gen difference and boxplot for centroid distances ===
# Only for genes with assigned family
distance_diff <- distances_after - distances_before
# df_dist_diff <- data.frame(
#   gene_id = gene_ids_fam,
#   diff = distance_diff
# )

# p_box_dist <- ggplot(df_dist_diff, aes(y = diff)) +
#   geom_boxplot(fill = "skyblue", alpha = 0.7) +
#   labs(title = "Gen-by-gen difference: Distance to Centroid (After - Before)",
#        y = "Difference in Distance (After - Before)") +
#   theme_minimal()

# # Save and print boxplot
# plot_file_box_dist <- "results/boxplot_centroid_distance_diff.png"
# ggsave(plot_file_box_dist, plot = p_box_dist, width = 7, height = 6)
# cat("Boxplot of centroid distance differences saved to:", plot_file_box_dist, "\n")
# print(p_box_dist)
# === Calculate tissue versatility for both normalization strategies ===
# Each gene's expression vector is its values across 13 tissues (columns)


# For mean-after-log
versatility_after <- tox_calculate_tissue_versatility(
  t(mat_after),
  vector_selection = 1:ncol(t(mat_after)),
  axis_selection = 1:nrow(t(mat_after))
)$tissue_versatilities
# For mean-before-log
versatility_before <- tox_calculate_tissue_versatility(
  t(mat_before),
  vector_selection = 1:ncol(t(mat_before)),
  axis_selection = 1:nrow(t(mat_before))
)$tissue_versatilities

# # Plot comparison of tissue versatility
# df_vers <- data.frame(
#   gene_id = rep(gene_ids, 2),
#   versatility = c(versatility_after, versatility_before),
#   method = rep(c("Mean-after-log", "Mean-before-log"), each = length(gene_ids))
# )

# p3 <- ggplot(df_vers, aes(x = versatility, fill = method)) +
#   geom_density(alpha = 0.5) +
#   labs(title = "Tissue Versatility: Mean-after-log vs Mean-before-log",
#        x = "Tissue Versatility", y = "Density") +
#   theme_minimal()

# # Save the plot in results/
# plot_file_vers <- "results/compare_tissue_versatility.png"
# ggsave(plot_file_vers, plot = p3, width = 7, height = 5)
# cat("Tissue versatility plot saved to:", plot_file_vers, "\n")
# print(p3)

# Summary statistics for tissue versatility
cat("Tissue versatility (mean-after-log): mean=", mean(versatility_after), " sd=", sd(versatility_after), "\n")
cat("Tissue versatility (mean-before-log): mean=", mean(versatility_before), " sd=", sd(versatility_before), "\n")

# === Gen-by-gen difference and boxplot for tissue versatility ===
versatility_diff <- versatility_after - versatility_before
# df_vers_diff <- data.frame(
#   gene_id = gene_ids,
#   diff = versatility_diff
# )

# p_box_vers <- ggplot(df_vers_diff, aes(y = diff)) +
#   geom_boxplot(fill = "orange", alpha = 0.7) +
#   labs(title = "Gen-by-gen difference: Tissue Versatility (After - Before)",
#        y = "Difference in Versatility (After - Before)") +
#   theme_minimal()

# # Save and print boxplot
# plot_file_box_vers <- "results/boxplot_tissue_versatility_diff.png"
# ggsave(plot_file_box_vers, plot = p_box_vers, width = 7, height = 6)
# cat("Boxplot of tissue versatility differences saved to:", plot_file_box_vers, "\n")
# print(p_box_vers)


# Test t para distancia al centroide
t_test_dist <- t.test(distance_diff)
cat("T-test para distancia al centroide:\n")
print(t_test_dist)

# Test Wilcoxon para distancia al centroide
wilcox_dist <- wilcox.test(distance_diff)
cat("Wilcoxon test para distancia al centroide:\n")
print(wilcox_dist)

# Test t para versatilidad
t_test_vers <- t.test(versatility_diff)
cat("T-test para versatilidad:\n")
print(t_test_vers)

# Test Wilcoxon para versatilidad
wilcox_vers <- wilcox.test(versatility_diff)
cat("Wilcoxon test para versatilidad:\n")
print(wilcox_vers)