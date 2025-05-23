library(dplyr)
library(tidyr)
library(preprocessCore)

# Function to read data
read_data <- function(input_file) {
  cat("Loading file:", input_file, "\n")
  df <- read.table(input_file, header = TRUE, sep = "\t")
  
  # Always take the first column as GeneID
  gene_id_column <- df[[1]]
  rownames(df) <- gene_id_column
  
  # Remove the first column from the dataframe
  df <- df[, -1]
  
  return(df)
}

normalize_by_std_dev <- function(df) {
  # Calculate sqrt(mean(x^2)) for each gene
  std_dev <- apply(df, 1, function(x) {
    sqrt(mean(x^2, na.rm = TRUE))
  })
  
  # If the standard deviation is 0, set it to 1 to avoid division by zero
  std_dev[std_dev == 0] <- 1
  
  # Normalize TPM values by dividing each gene by its standard deviation
  normalized_df <- sweep(df, 1, std_dev, FUN = "/")
  
  return(normalized_df)
}


# Function to perform quantile normalization
quantile_normalization <- function(df) {
  # Save original column names before converting to a matrix
  original_colnames <- colnames(df)
  original_rownames <- rownames(df)
  
  # Perform quantile normalization
  normalized_data <- normalize.quantiles(as.matrix(df))
#   print(normalized_data)
  
  # Convert the normalized data back to a dataframe
  normalized_quantile_df <- as.data.frame(normalized_data)
  
#   # Restore the original column names
  colnames(normalized_quantile_df) <- original_colnames
#   rownames(normalized_quantile_df) <- original_rownames
  normalized_quantile_df$gene_id <- original_rownames
  normalized_quantile_df <- normalized_quantile_df[, c("gene_id", original_colnames)]


  return(normalized_quantile_df)
}

# Function to apply log2 transformation
log2_transformation <- function(df) {
  # Apply log2(x + 1) transformation, ignoring NaN values
  log_transformed_df <- df
  log_transformed_df[, -1] <- log2(df[, -1] + 1)
  
  return(log_transformed_df)
}

# Function to calculate log fold change
calculate_log2fc_from_log <- function(df, gene_col = "gene_id") {
  # Ensure numeric columns
  numeric_cols <- setdiff(colnames(df), gene_col)
  df[numeric_cols] <- lapply(df[numeric_cols], as.numeric)
  
  # Detect tissues from column names
  tissues <- unique(gsub("_diet[MP]", "", numeric_cols))
  
  # Calculate logFC
  log_fc_cols <- c()  # to track created columns
  for (t in tissues) {
    col_M <- paste0(t, "_dietM")
    col_P <- paste0(t, "_dietP")
    
    if (col_M %in% colnames(df) && col_P %in% colnames(df)) {
      log_fc_col <- paste0(t, "_logFC")
      df[[log_fc_col]] <- df[[col_P]] - df[[col_M]]
      log_fc_cols <- c(log_fc_cols, log_fc_col)
    }
  }
  
  # Return only gene_id + logFC columns
  return(df[, c(gene_col, log_fc_cols)])
}


# Function to calculate tissue averages, grouping replicates
calculate_tissue_averages <- function(df) {
  tissue_columns <- colnames(df)
  tissue_groups <- list()

  # Group columns by tissue type (e.g., 'heart', 'testis') by removing the replicate part
  for (col in tissue_columns) {
    # Extract tissue name from 'heart_rep1' -> 'heart'
    tissue_name <- paste(strsplit(col, "_")[[1]][1:(length(strsplit(col, "_")[[1]]) - 1)], collapse = "_")

    if (!(tissue_name %in% names(tissue_groups))) {
      tissue_groups[[tissue_name]] <- c()
    }
    tissue_groups[[tissue_name]] <- c(tissue_groups[[tissue_name]], col)
  }

  # Create a new dataframe to store tissue averages
  averaged_df <- data.frame(row.names = rownames(df))

  # For each tissue group, calculate the average of the replicates
  for (tissue_name in names(tissue_groups)) {
    columns <- tissue_groups[[tissue_name]]
    
    # Ensure the columns are treated as a dataframe for rowMeans
    tissue_data <- df[, columns, drop = FALSE]  # Prevent reduction to vector if only one column
    
    # Calculate the average across the replicates
    averaged_df[[tissue_name]] <- rowMeans(tissue_data, na.rm = TRUE)
  }

  return(averaged_df)
}

# Main function to normalize TPM and save the result
normalize_tpm <- function(input_file, output_file) {
  df <- read_data(input_file)
  if (is.null(df)) {
    return(NULL)
  }
  
  print(head(df))

  # Normalize by standard deviation
  normalized_df <- normalize_by_std_dev(df)
  cat("Normalized Data (by Standard Deviation):\n")
  print(head(normalized_df))  # Print the normalized data
  
  # Apply quantile normalization
  quantile_normalized_df <- quantile_normalization(normalized_df)
  cat("Quantile Normalized Data:\n")
  print(head(quantile_normalized_df))  # Print the quantile normalized data
  
#   Apply log2 transformation
  cat("Log transformation:\n")

  log_transformed_df <- log2_transformation(quantile_normalized_df)
  print(head(log_transformed_df))

    # Calculate log fold change
  cat("Log fold change:\n")
  log_fc_df <- calculate_log2fc_from_log(log_transformed_df)
  print(head(log_fc_df))

  # Calculate tissue averages
  # averaged_df <- calculate_tissue_averages(log_transformed_df)
  # cat("Averaged Tissue Data:\n")
  # print(averaged_df)  # Print the averaged tissue data
  
  # Write the result to a new file
  write.table(log_fc_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Averaged tissue values saved:", output_file, "\n")
}

# Example usage
input_file <- "material/tpm.tsv"  # Adjust the input file path
output_file <- "results/normalization.tsv"  # Output file path
normalize_tpm(input_file, output_file)