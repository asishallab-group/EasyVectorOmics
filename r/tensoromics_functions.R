# === Load the shared library ===
dyn.load("build/tox_normalization.so")



#' Normalize gene expression values by standard deviation
#'
#' This function wraps the Fortran subroutine `normalize_by_std_dev`
#' to normalize each gene's expression across tissues by the standard deviation.
#'
#' @param input_matrix A numeric matrix with genes as rows and tissues as columns.
#' @return A normalized numeric matrix with the same dimensions and names as the input.
#' @details
#' - The input matrix is flattened into a column-major vector.
#' - Calls a Fortran routine that normalizes each gene across tissues.
#' - Restores the original row and column names after normalization.
#'
#' @examples
#' normalized_matrix <- normalize_by_std_dev(input_matrix)
normalize_by_std_dev <- function(input_matrix) {
  n_genes <- nrow(input_matrix)  # Number of genes (rows)
  n_tissues <- ncol(input_matrix)  # Number of tissues (columns)

  # print(head(input_matrix))

  # Prepare the input vector (flatten matrix column-major) and allocate output space
  input_vector <- as.numeric(as.vector(input_matrix))
  output_vector <- numeric(n_genes * n_tissues)

  # Call the Fortran subroutine
  result <- .Fortran("normalize_by_std_dev_r_py",
               as.integer(n_genes),
               as.integer(n_tissues),
               input_vector,
               output_vector)

  matrix(result[[4]], nrow = n_genes, ncol = n_tissues,
         dimnames = dimnames(input_matrix))

}


#' Quantile normalization of gene expression values
#'
#' This function wraps the Fortran subroutine `quantile_normalization`
#' to apply quantile normalization across all tissue columns.
#'
#' @param input_matrix A numeric matrix with genes as rows and tissues as columns.
#' @return A quantile-normalized numeric matrix with the same dimensions and names as the input.
#' @details
#' - The input matrix is flattened into a column-major vector.
#' - The Fortran subroutine sorts and aligns the distributions across tissues.
#' - After normalization, the original row and column names are restored.
#'
#' @examples
#' normalized_matrix <- quantile_normalization(input_matrix)
quantile_normalization <- function(input_matrix) {
  n_genes <- nrow(input_matrix)
  n_tissues <- ncol(input_matrix)

  # Flatten input matrix (column-major) and preallocate vectors
  input_vector <- as.numeric(as.vector(input_matrix))
  output_vector <- numeric(n_genes * n_tissues)
  temp_col <- numeric(n_genes)
  rank_means <- numeric(n_genes)
  perm <- integer(n_genes)

  # Estimar tamaño máximo para la pila (según pseudocódigo: log2(n) + 10)
  max_stack <- as.integer(ceiling(log2(n_genes)) + 10)
  stack_left <- integer(max_stack)
  stack_right <- integer(max_stack)

  # Fortran interop: asegurar tipos
  storage.mode(input_vector) <- "double"
  storage.mode(output_vector) <- "double"
  storage.mode(temp_col) <- "double"
  storage.mode(rank_means) <- "double"
  storage.mode(perm) <- "integer"
  storage.mode(stack_left) <- "integer"
  storage.mode(stack_right) <- "integer"

  # Initialize first stack entry manually
  # stack_left[1] <- 1L
  # stack_right[1] <- n_genes

  result <- .Fortran("quantile_normalization_r_py",
    as.integer(n_genes),
    as.integer(n_tissues),
    input_vector,
    output_vector,
    temp_col,
    rank_means,
    perm,
    stack_left,
    stack_right,
    as.integer(max_stack)
  )

  matrix(result[[4]], nrow = n_genes, ncol = n_tissues,
         dimnames = dimnames(input_matrix))
}




# #' Apply log2(x + 1) transformation to gene expression values
# #'
# #' This function wraps the Fortran subroutine `log2_transformation`
# #' to apply a log2(x + 1) transformation to each element in the input matrix.
# #'
# #' @param input_matrix A numeric matrix with genes as rows and tissues as columns.
# #' @return A numeric matrix with log2-transformed expression values, 
# #' preserving the same dimensions and names as the input.
# #' @details
# #' - The input matrix is flattened into a column-major vector.
# #' - The Fortran subroutine applies a log2(x+1) transformation to each value.
# #' - After transformation, the original row and column names are restored.
# #'
# #' @examples
# #' log_matrix <- log2_transformation(input_matrix)
# log2_transformation <- function(input_matrix) {
#   n_genes <- nrow(input_matrix)  # Number of genes (rows)
#   n_tissues <- ncol(input_matrix)  # Number of tissues (columns)

#   # # Save column and row names for later restoration
#   # col_names <- colnames(input_matrix)
#   # row_names <- rownames(input_matrix)

#   # Prepare the input vector (flatten matrix column-major) and allocate output space
#   input_vector <- as.numeric(as.vector(input_matrix))
#   output_vector <- numeric(n_genes * n_tissues)

#   # Call the Fortran subroutine
#   result <- .Fortran("log2_transformation",
#                as.integer(n_genes),
#                as.integer(n_tissues),
#                input_vector,
#                output_vector)

#   # Reconstruct the transformed matrix
#   matrix(result[[4]], nrow = n_genes, ncol = n_tissues,
#   dimnames = dimnames(input_matrix))

#   # # Restore row and column names
#   # colnames(normalized_matrix) <- col_names
#   # rownames(normalized_matrix) <- row_names

#   # return(normalized_matrix)
# }

# #' Parse tissue group name from column name
# #'
# #' Helper function to extract the tissue group from a column name.
# #' If the column name ends with an underscore followed by a number (e.g., a replicate),
# #' it removes the replicate identifier and returns the base group name.
# #'
# #' @param colname A string with the column name to parse.
# #' @return A string representing the parsed tissue group name.
# #' @examples
# #' parse_tissue_group("muscle_dietM_1") # returns "muscle_dietM"
# #' parse_tissue_group("brain_dietP")    # returns "brain_dietP"
# parse_tissue_group <- function(colname) {
#   parts <- strsplit(colname, "_")[[1]]
  
#   if (length(parts) >= 2 && grepl("^[0-9]+$", parts[length(parts)])) {
#     # If last part is a number (replicate), remove it
#     return(paste(parts[1:(length(parts)-1)], collapse = "_"))
#   } else {
#     # Otherwise, return full name
#     return(colname)
#   }
# }


# #' Calculate average expression across replicates for each tissue group
# #'
# #' This function wraps the Fortran subroutine `calc_tiss_avg`
# #' to compute the mean expression value for replicates grouped by tissue.
# #'
# #' @param df A data frame or matrix with genes as rows and tissue replicates as columns.
# #' @return A data frame with genes as rows and averaged tissues as columns.
# #' @details
# #' - Replicate columns are grouped based on their parsed tissue group names.
# #' - The input matrix is sorted according to groups before being processed.
# #' - Calls a Fortran subroutine that calculates averages within each group.
# #' - The result restores the original gene IDs as row names.
# #'
# #' @examples
# #' averaged_df <- calculate_tissue_averages(df)
# calculate_tissue_averages <- function(df) {
#   n_genes <- nrow(df)      # Number of genes (rows)
#   n_columns <- ncol(df)    # Number of columns (tissues)

#   # --- Parse all column names to find their corresponding tissue group ---
#   tissue_groups <- sapply(colnames(df), parse_tissue_group)

#   # --- Identify unique tissue groups ---
#   unique_groups <- unique(tissue_groups)
#   n_groups <- length(unique_groups)

#   # --- Initialize mapping for groups ---
#   group_starts <- integer(n_groups)
#   group_counts <- integer(n_groups)

#   # --- Sort the dataframe by tissue group name ---
#   df_sorted <- df[, order(tissue_groups)]
#   sorted_tissue_groups <- tissue_groups[order(tissue_groups)]

#   current_group <- sorted_tissue_groups[1]
#   group_starts[1] <- 1
#   group_counts[1] <- 1
#   group_idx <- 1

#   # --- Build group_starts and group_counts arrays ---
#   for (i in 2:length(sorted_tissue_groups)) {
#     if (sorted_tissue_groups[i] == current_group) {
#       group_counts[group_idx] <- group_counts[group_idx] + 1
#     } else
#       {
#       group_idx <- group_idx + 1
#       group_starts[group_idx] <- i
#       group_counts[group_idx] <- 1
#       current_group <- sorted_tissue_groups[i]
#     }
#   }

#   # --- Prepare input vector and allocate output space ---
#   input_vector <- as.numeric(as.vector(as.matrix(df_sorted)))
#   output_vector <- numeric(n_genes * n_groups)

#   # --- Call the Fortran subroutine to calculate averages ---
#   result <- .Fortran("calc_tiss_avg",
#                as.integer(n_genes),
#                as.integer(n_groups),
#                as.integer(group_starts),
#                as.integer(group_counts),
#                as.numeric(input_vector),
#                as.numeric(output_vector))

#   # --- Reconstruct output matrix ---
#   output_matrix <- matrix(result[[6]], nrow = n_genes, ncol = n_groups)

#   # --- Restore column and row names ---
#   colnames(output_matrix) <- unique_groups
#   rownames(output_matrix) <- rownames(df)

#   return(as.data.frame(output_matrix))
# }

# #' Prepare control and condition column indices based on naming patterns
# #'
# #' This helper function searches for columns in the input dataframe that match
# #' specified control and condition patterns. It builds the mapping necessary 
# #' to calculate fold changes.
# #'
# #' @param df A data frame with expression data, genes as rows and tissues/conditions as columns.
# #' @param control_pattern A string pattern used to identify control columns.
# #' @param condition_patterns A character vector with patterns used to identify condition columns.
# #' @return A list containing:
# #' \describe{
# #'   \item{control_cols}{Indices of control columns}
# #'   \item{condition_cols}{Indices of condition columns}
# #'   \item{condition_labels}{Column names for the resulting fold changes}
# #' }
# #'
# #' @examples
# #' indices_info <- prepare_indices_by_patterns(df, control_pattern = "dietM", condition_patterns = c("dietP"))
# prepare_indices_by_patterns <- function(df, control_pattern, condition_patterns) {
#   colnames_df <- colnames(df)

#   control_cols <- integer(0)         # Will store control column indices
#   condition_cols <- integer(0)       # Will store condition column indices
#   condition_labels <- character(0)   # Labels for resulting logFC columns

#   # --- Find all control columns matching control_pattern ---
#   control_candidates <- grep(control_pattern, colnames_df, value = TRUE)

#   if (length(control_candidates) == 0) {
#     stop("No control columns found matching pattern: ", control_pattern)
#   }

#   # --- For each control, find associated condition columns ---
#   for (control in control_candidates) {
#     tissue_prefix <- sub(paste0("_", control_pattern), "", control)

#     for (cond_pattern in condition_patterns) {
#       pattern_to_match <- paste0("^", tissue_prefix, "_", cond_pattern)
#       matching_conditions <- grep(pattern_to_match, colnames_df, value = TRUE)

#       for (condition in matching_conditions) {
#         control_cols <- c(control_cols, which(colnames_df == control))
#         condition_cols <- c(condition_cols, which(colnames_df == condition))
#         condition_labels <- c(condition_labels, paste0(condition, "_logFC"))
#       }
#     }
#   }

#   if (length(control_cols) == 0) {
#     stop("No valid control-condition pairs found!")
#   }

#   return(list(control_cols = control_cols, condition_cols = condition_cols, condition_labels = condition_labels))
# }

# #' Calculate log2 fold changes based on control and condition patterns
# #'
# #' This function wraps the Fortran subroutine `calc_fchange`
# #' to calculate the fold changes between control and condition columns.
# #'
# #' @param df A data frame with genes as rows and tissues/conditions as columns.
# #' @param control_pattern A string pattern to detect control columns.
# #' @param condition_patterns A character vector with patterns to detect condition columns.
# #' @return A data frame with genes as rows and log2 fold change values as columns.
# #'
# #' @examples
# #' fc_df <- calculate_fc_by_patterns(df, control_pattern = "dietM", condition_patterns = c("dietP"))
# calculate_fc_by_patterns <- function(df, control_pattern, condition_patterns) {
#   n_genes <- nrow(df)       # Number of genes
#   n_columns <- ncol(df)     # Number of columns (conditions)

#   # --- Identify control and condition columns ---
#   indices_info <- prepare_indices_by_patterns(df, control_pattern, condition_patterns)
#   control_cols <- indices_info$control_cols
#   condition_cols <- indices_info$condition_cols
#   condition_labels <- indices_info$condition_labels
  
#   n_pairs <- length(control_cols)

#   # --- Prepare input and output vectors ---
#   input_vector <- as.numeric(as.vector(as.matrix(df)))
#   output_vector <- numeric(n_genes * n_pairs)

#   print(control_cols)
#   print(condition_cols)
#   print(head(input_vector))
#   # --- Call Fortran subroutine to calculate fold changes ---
#   result <- .Fortran("calc_fchange",
#                as.integer(n_genes),
#                as.integer(n_pairs),
#                as.integer(control_cols),
#                as.integer(condition_cols),
#                input_vector,
#                output_vector)

#   # --- Reconstruct the fold change matrix ---
#   output_matrix <- matrix(result[[6]], nrow = n_genes, ncol = n_pairs)
#   colnames(output_matrix) <- condition_labels
#   rownames(output_matrix) <- rownames(df)

#   return(as.data.frame(output_matrix))
# }

