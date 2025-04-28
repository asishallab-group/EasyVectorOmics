# Load all TensorOmics helper functions
source("r/tensoromics_functions.R")

# === Example of full normalization pipeline ===

# Load raw expression data
input_file <- "material/tpm.tsv"  # Update the input file path if necessary
output_file <- "results/normalization.tsv"  # Output file path for normalized data

df <- read.table(input_file, header = TRUE, sep = "\t")

# Prepare matrix for processing (removing the gene ID column)
gene_ids <- df[[1]]              # Save gene identifiers
col_names <- colnames(df)[-1]    
df_matrix <- as.matrix(df[,-1])   # Convert the expression values into a matrix

# === Apply normalization steps sequentially ===
normalized_matrix_std <- normalize_by_std_dev(df_matrix)    # Normalize by standard deviation
normalized_matrix_qtl <- quantile_normalization(normalized_matrix_std)  # Apply quantile normalization
normalized_matrix_log <- log2_transformation(normalized_matrix_qtl)     # Log2(x+1) transformation
averaged_df <- calculate_tissue_averages(normalized_matrix_log)         # Average replicates by tissue

# === Calculate fold changes between specified groups ===
fc_df <- calculate_fc_by_patterns(
  df = averaged_df,
  control_pattern = "dietM",           # Define control pattern
  condition_patterns = c("dietP")       # Define conditions pattern(s)
)

# === Save the resulting normalized and transformed data ===
normalized_df <- data.frame(gene_id = gene_ids, fc_df)
write.table(normalized_df, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
