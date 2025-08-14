# === Load the shared library ===
dyn.load("build/libtensor-omics.so")



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

  # Validate input data before calling Fortran
  if (any(is.na(input_vector))) {
    stop("Error: Input matrix contains NA values. Found ", sum(is.na(input_vector)), " NA values.")
  }
  if (any(is.infinite(input_vector))) {
    stop("Error: Input matrix contains infinite values. Found ", sum(is.infinite(input_vector)), " infinite values.")
  }
  if (any(is.nan(input_vector))) {
    stop("Error: Input matrix contains NaN values. Found ", sum(is.nan(input_vector)), " NaN values.")
  }

  # Call the Fortran subroutine
  result <- .Fortran("normalize_by_std_dev_r",
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

  result <- .Fortran("quantile_normalization_r",
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




#' Apply log2(x + 1) transformation to gene expression values
#'
#' This function wraps the Fortran subroutine `log2_transformation`
#' to apply a log2(x + 1) transformation to each element in the input matrix.
#'
#' @param input_matrix A numeric matrix with genes as rows and tissues as columns.
#' @return A numeric matrix with log2-transformed expression values, 
#' preserving the same dimensions and names as the input.
#' @details
#' - The input matrix is flattened into a column-major vector.
#' - The Fortran subroutine applies a log2(x+1) transformation to each value.
#' - After transformation, the original row and column names are restored.
#'
#' @examples
#' log_matrix <- log2_transformation(input_matrix)
log2_transformation <- function(input_matrix) {
  n_genes <- nrow(input_matrix)  # Number of genes (rows)
  n_tissues <- ncol(input_matrix)  # Number of tissues (columns)

  # # Save column and row names for later restoration
  # col_names <- colnames(input_matrix)
  # row_names <- rownames(input_matrix)

  # Prepare the input vector (flatten matrix column-major) and allocate output space
  input_vector <- as.numeric(as.vector(input_matrix))
  output_vector <- numeric(n_genes * n_tissues)

  # Call the Fortran subroutine
  result <- .Fortran("log2_transformation_r",
               as.integer(n_genes),
               as.integer(n_tissues),
               input_vector,
               output_vector)

  # Reconstruct the transformed matrix
  matrix(result[[4]], nrow = n_genes, ncol = n_tissues,
  dimnames = dimnames(input_matrix))

  # # Restore row and column names
  # colnames(normalized_matrix) <- col_names
  # rownames(normalized_matrix) <- row_names

  # return(normalized_matrix)
}

#' Parse tissue group name from column name
#'
#' Helper function to extract the tissue group from a column name.
#' Handles various replicate naming patterns:
#' - "muscle_dietM_1" -> "muscle_dietM" 
#' - "Adipose_rep1" -> "Adipose"
#' - "Brain_rep2" -> "Brain"
#' - "tissue_condition_rep3" -> "tissue_condition"
#'
#' @param colname A string with the column name to parse.
#' @return A string representing the parsed tissue group name.
#' @examples
#' parse_tissue_group("muscle_dietM_1") # returns "muscle_dietM"
#' parse_tissue_group("Adipose_rep1")   # returns "Adipose"
#' parse_tissue_group("brain_dietP")    # returns "brain_dietP"
parse_tissue_group <- function(colname) {
  parts <- strsplit(colname, "_")[[1]]
  
  # Pattern 1: ends with just a number (e.g., "muscle_dietM_1")
  if (length(parts) >= 2 && grepl("^[0-9]+$", parts[length(parts)])) {
    return(paste(parts[1:(length(parts)-1)], collapse = "_"))
  }
  
  # Pattern 2: ends with "rep" followed by number (e.g., "Adipose_rep1")
  if (length(parts) >= 2 && grepl("^rep[0-9]+$", parts[length(parts)])) {
    return(paste(parts[1:(length(parts)-1)], collapse = "_"))
  }
  
  # Pattern 3: ends with "replicate" followed by number (e.g., "Tissue_replicate1")
  if (length(parts) >= 2 && grepl("^replicate[0-9]+$", parts[length(parts)])) {
    return(paste(parts[1:(length(parts)-1)], collapse = "_"))
  }
  
  # If no pattern matches, return full name
  return(colname)
}


#' Calculate average expression across replicates for each tissue group
#'
#' This function wraps the Fortran subroutine `calc_tiss_avg`
#' to compute the mean expression value for replicates grouped by tissue.
#'
#' @param df A data frame or matrix with genes as rows and tissue replicates as columns.
#' @return A data frame with genes as rows and averaged tissues as columns.
#' @details
#' - Replicate columns are grouped based on their parsed tissue group names.
#' - The input matrix is sorted according to groups before being processed.
#' - Calls a Fortran subroutine that calculates averages within each group.
#' - The result restores the original gene IDs as row names.
#'
#' @examples
#' averaged_df <- calculate_tissue_averages(df)
calculate_tissue_averages <- function(df) {
  n_genes <- nrow(df)      # Number of genes (rows)
  n_columns <- ncol(df)    # Number of columns (tissues)

  # --- Parse all column names to find their corresponding tissue group ---
  tissue_groups <- sapply(colnames(df), parse_tissue_group)

  # --- Identify unique tissue groups ---
  unique_groups <- unique(tissue_groups)
  n_groups <- length(unique_groups)

  # --- Initialize mapping for groups ---
  group_starts <- integer(n_groups)
  group_counts <- integer(n_groups)

  # --- Sort the dataframe by tissue group name ---
  df_sorted <- df[, order(tissue_groups)]
  sorted_tissue_groups <- tissue_groups[order(tissue_groups)]

  current_group <- sorted_tissue_groups[1]
  group_starts[1] <- 1
  group_counts[1] <- 1
  group_idx <- 1

  # --- Build group_starts and group_counts arrays ---
  for (i in 2:length(sorted_tissue_groups)) {
    if (sorted_tissue_groups[i] == current_group) {
      group_counts[group_idx] <- group_counts[group_idx] + 1
    } else
      {
      group_idx <- group_idx + 1
      group_starts[group_idx] <- i
      group_counts[group_idx] <- 1
      current_group <- sorted_tissue_groups[i]
    }
  }

  # --- Prepare input vector and allocate output space ---
  input_vector <- as.numeric(as.vector(as.matrix(df_sorted)))
  output_vector <- numeric(n_genes * n_groups)

  # --- Call the Fortran subroutine to calculate averages ---
  result <- .Fortran("calc_tiss_avg_r",
               as.integer(n_genes),
               as.integer(n_groups),
               as.integer(group_starts),
               as.integer(group_counts),
               as.numeric(input_vector),
               as.numeric(output_vector))

  # --- Reconstruct output matrix ---
  output_matrix <- matrix(result[[6]], nrow = n_genes, ncol = n_groups)

  # --- Restore column and row names ---
  colnames(output_matrix) <- unique_groups
  rownames(output_matrix) <- rownames(df)

  return(as.data.frame(output_matrix))
}

#' Prepare control and condition column indices based on naming patterns
#'
#' This helper function searches for columns in the input dataframe that match
#' specified control and condition patterns. It builds the mapping necessary 
#' to calculate fold changes.
#'
#' @param df A data frame with expression data, genes as rows and tissues/conditions as columns.
#' @param control_pattern A string pattern used to identify control columns.
#' @param condition_patterns A character vector with patterns used to identify condition columns.
#' @return A list containing:
#' \describe{
#'   \item{control_cols}{Indices of control columns}
#'   \item{condition_cols}{Indices of condition columns}
#'   \item{condition_labels}{Column names for the resulting fold changes}
#' }
#'
#' @examples
#' indices_info <- prepare_indices_by_patterns(df, control_pattern = "dietM", condition_patterns = c("dietP"))
prepare_indices_by_patterns <- function(df, control_pattern, condition_patterns) {
  colnames_df <- colnames(df)

  control_cols <- integer(0)         # Will store control column indices
  condition_cols <- integer(0)       # Will store condition column indices
  condition_labels <- character(0)   # Labels for resulting logFC columns

  # --- Find all control columns matching control_pattern ---
  control_candidates <- grep(control_pattern, colnames_df, value = TRUE)

  if (length(control_candidates) == 0) {
    stop("No control columns found matching pattern: ", control_pattern)
  }

  # --- For each control, find associated condition columns ---
  for (control in control_candidates) {
    tissue_prefix <- sub(paste0("_", control_pattern), "", control)

    for (cond_pattern in condition_patterns) {
      pattern_to_match <- paste0("^", tissue_prefix, "_", cond_pattern)
      matching_conditions <- grep(pattern_to_match, colnames_df, value = TRUE)

      for (condition in matching_conditions) {
        control_cols <- c(control_cols, which(colnames_df == control))
        condition_cols <- c(condition_cols, which(colnames_df == condition))
        condition_labels <- c(condition_labels, paste0(condition, "_logFC"))
      }
    }
  }

  if (length(control_cols) == 0) {
    stop("No valid control-condition pairs found!")
  }

  return(list(control_cols = control_cols, condition_cols = condition_cols, condition_labels = condition_labels))
}

#' Calculate log2 fold changes based on control and condition patterns
#'
#' This function wraps the Fortran subroutine `calc_fchange`
#' to calculate the fold changes between control and condition columns.
#'
#' @param df A data frame with genes as rows and tissues/conditions as columns.
#' @param control_pattern A string pattern to detect control columns.
#' @param condition_patterns A character vector with patterns to detect condition columns.
#' @return A data frame with genes as rows and log2 fold change values as columns.
#'
#' @examples
#' fc_df <- calculate_fc_by_patterns(df, control_pattern = "dietM", condition_patterns = c("dietP"))
calculate_fc_by_patterns <- function(df, control_pattern, condition_patterns) {
  n_genes <- nrow(df)       # Number of genes
  n_columns <- ncol(df)     # Number of columns (conditions)

  # --- Identify control and condition columns ---
  indices_info <- prepare_indices_by_patterns(df, control_pattern, condition_patterns)
  control_cols <- indices_info$control_cols
  condition_cols <- indices_info$condition_cols
  condition_labels <- indices_info$condition_labels
  
  n_pairs <- length(control_cols)

  # --- Prepare input and output vectors ---
  input_vector <- as.numeric(as.vector(as.matrix(df)))
  output_vector <- numeric(n_genes * n_pairs)

  print(control_cols)
  print(condition_cols)
  print(head(input_vector))
  # --- Call Fortran subroutine to calculate fold changes ---
  result <- .Fortran("calc_fchange_r",
               as.integer(n_genes),
               as.integer(n_columns),   # Pass n_cols as required by Fortran
               as.integer(n_pairs),
               as.integer(control_cols),
               as.integer(condition_cols),
               input_vector,
               output_vector)

  # --- Reconstruct the fold change matrix ---
  output_matrix <- matrix(result[[7]], nrow = n_genes, ncol = n_pairs)
  colnames(output_matrix) <- condition_labels
  rownames(output_matrix) <- rownames(df)

  return(as.data.frame(output_matrix))
}

#' Diagnose data quality issues in gene expression matrix
#'
#' This function examines the input matrix for common data quality issues
#' that could cause problems in downstream analysis.
#'
#' @param input_matrix A numeric matrix with genes as rows and tissues as columns.
#' @param show_details Logical indicating whether to show detailed information.
#' @return A list with diagnostic information about the data quality.
#' @examples
#' diagnostics <- diagnose_data_quality(input_matrix)
diagnose_data_quality <- function(input_matrix, show_details = TRUE) {
  n_genes <- nrow(input_matrix)
  n_tissues <- ncol(input_matrix)
  total_values <- n_genes * n_tissues
  
  # Check for different types of problematic values
  na_count <- sum(is.na(input_matrix))
  inf_count <- sum(is.infinite(input_matrix))
  nan_count <- sum(is.nan(input_matrix))
  zero_count <- sum(input_matrix == 0, na.rm = TRUE)
  negative_count <- sum(input_matrix < 0, na.rm = TRUE)
  
  # Find problematic genes (rows with issues)
  genes_with_na <- which(apply(input_matrix, 1, function(x) any(is.na(x))))
  genes_with_inf <- which(apply(input_matrix, 1, function(x) any(is.infinite(x))))
  genes_with_nan <- which(apply(input_matrix, 1, function(x) any(is.nan(x))))
  genes_all_zero <- which(apply(input_matrix, 1, function(x) all(x == 0, na.rm = TRUE)))
  
  # Summary statistics
  if (na_count == 0 && inf_count == 0 && nan_count == 0) {
    min_val <- min(input_matrix, na.rm = TRUE)
    max_val <- max(input_matrix, na.rm = TRUE)
    mean_val <- mean(input_matrix, na.rm = TRUE)
  } else {
    min_val <- NA
    max_val <- NA
    mean_val <- NA
  }
  
  diagnostics <- list(
    dimensions = c(genes = n_genes, tissues = n_tissues, total_values = total_values),
    problems = list(
      na_count = na_count,
      inf_count = inf_count,
      nan_count = nan_count,
      zero_count = zero_count,
      negative_count = negative_count
    ),
    problematic_genes = list(
      genes_with_na = genes_with_na,
      genes_with_inf = genes_with_inf,
      genes_with_nan = genes_with_nan,
      genes_all_zero = genes_all_zero
    ),
    statistics = list(
      min_val = min_val,
      max_val = max_val,
      mean_val = mean_val
    )
  )
  
  if (show_details) {
    cat("=== DATA QUALITY DIAGNOSTICS ===\n")
    cat("Matrix dimensions:", n_genes, "genes x", n_tissues, "tissues (", total_values, "total values)\n\n")
    
    cat("Problem summary:\n")
    cat("  - NA values:", na_count, "(", round(100*na_count/total_values, 2), "%)\n")
    cat("  - Infinite values:", inf_count, "(", round(100*inf_count/total_values, 2), "%)\n")
    cat("  - NaN values:", nan_count, "(", round(100*nan_count/total_values, 2), "%)\n")
    cat("  - Zero values:", zero_count, "(", round(100*zero_count/total_values, 2), "%)\n")
    cat("  - Negative values:", negative_count, "(", round(100*negative_count/total_values, 2), "%)\n\n")
    
    cat("Problematic genes:\n")
    cat("  - Genes with NA:", length(genes_with_na), "\n")
    cat("  - Genes with Inf:", length(genes_with_inf), "\n") 
    cat("  - Genes with NaN:", length(genes_with_nan), "\n")
    cat("  - Genes all zero:", length(genes_all_zero), "\n\n")
    
    if (na_count == 0 && inf_count == 0 && nan_count == 0) {
      cat("Data range:\n")
      cat("  - Min value:", min_val, "\n")
      cat("  - Max value:", max_val, "\n")
      cat("  - Mean value:", mean_val, "\n\n")
    }
    
    # Show examples of problematic genes
    if (length(genes_with_na) > 0) {
      cat("First few genes with NA values:\n")
      print(head(genes_with_na, 5))
      cat("\n")
    }
    
    if (length(genes_with_inf) > 0) {
      cat("First few genes with infinite values:\n")
      print(head(genes_with_inf, 5))
      cat("\n")
    }
  }
  
  return(invisible(diagnostics))
}

#' Clean data by removing or imputing problematic values
#'
#' This function handles NA, NaN, Inf values and genes that are all zeros
#' to prepare data for Fortran normalization routines.
#'
#' @param df_matrix A numeric matrix with genes as rows and tissues as columns
#' @param remove_all_zero_genes Logical, whether to remove genes that are all zeros
#' @param na_strategy Strategy for handling NA values: "remove_genes", "remove_samples", "impute_zero", "impute_mean"
#' @param min_expression_threshold Minimum expression value to consider (values below this become 0)
#' @return A cleaned matrix ready for normalization
clean_data_for_normalization <- function(df_matrix, 
                                        remove_all_zero_genes = TRUE,
                                        na_strategy = "remove_genes",
                                        min_expression_threshold = 0.0,  # Changed default to 0.0
                                        convert_small_to_zero = FALSE) {   # New parameter to control this behavior
  
  cat("=== CLEANING DATA FOR NORMALIZATION ===\n")
  original_dims <- dim(df_matrix)
  cat("Original dimensions:", original_dims[1], "genes x", original_dims[2], "tissues\n")
  
  # Step 1: Handle very small values (only if explicitly requested)
  if (convert_small_to_zero && min_expression_threshold > 0.0) {
    small_values <- df_matrix > 0 & df_matrix < min_expression_threshold
    if (sum(small_values, na.rm = TRUE) > 0) {
      cat("Converting", sum(small_values, na.rm = TRUE), "values <", min_expression_threshold, "to zero\n")
      df_matrix[small_values] <- 0
    }
  } else {
    cat("Preserving all small values (convert_small_to_zero = FALSE)\n")
  }
  
  # Step 2: Handle infinite values
  inf_values <- is.infinite(df_matrix)
  if (sum(inf_values, na.rm = TRUE) > 0) {
    cat("WARNING: Converting", sum(inf_values, na.rm = TRUE), "infinite values to NA\n")
    df_matrix[inf_values] <- NA
  }
  
  # Step 3: Handle NaN values
  nan_values <- is.nan(df_matrix)
  if (sum(nan_values, na.rm = TRUE) > 0) {
    cat("WARNING: Converting", sum(nan_values, na.rm = TRUE), "NaN values to NA\n")
    df_matrix[nan_values] <- NA
  }
  
  # Step 4: Handle NA values according to strategy
  na_count <- sum(is.na(df_matrix))
  if (na_count > 0) {
    cat("Handling", na_count, "NA values using strategy:", na_strategy, "\n")
    
    if (na_strategy == "remove_genes") {
      # Remove genes with any NA values
      genes_with_na <- apply(df_matrix, 1, function(x) any(is.na(x)))
      df_matrix <- df_matrix[!genes_with_na, , drop = FALSE]
      cat("Removed", sum(genes_with_na), "genes with NA values\n")
      
    } else if (na_strategy == "remove_samples") {
      # Remove samples/tissues with any NA values
      samples_with_na <- apply(df_matrix, 2, function(x) any(is.na(x)))
      df_matrix <- df_matrix[, !samples_with_na, drop = FALSE]
      cat("Removed", sum(samples_with_na), "samples with NA values\n")
      
    } else if (na_strategy == "impute_zero") {
      # Replace NA with 0
      df_matrix[is.na(df_matrix)] <- 0
      cat("Imputed", na_count, "NA values with zero\n")
      
    } else if (na_strategy == "impute_mean") {
      # Replace NA with gene mean (row-wise)
      for (i in 1:nrow(df_matrix)) {
        na_positions <- is.na(df_matrix[i, ])
        if (any(na_positions)) {
          gene_mean <- mean(df_matrix[i, ], na.rm = TRUE)
          if (is.finite(gene_mean)) {
            df_matrix[i, na_positions] <- gene_mean
          } else {
            df_matrix[i, na_positions] <- 0  # If all values are NA, use 0
          }
        }
      }
      cat("Imputed", na_count, "NA values with gene means\n")
      
    } else if (na_strategy == "smart_impute") {
      # More sophisticated strategy: remove genes with >50% NA, impute the rest
      na_threshold <- 0.5  # Remove genes with more than 50% NA values
      
      genes_with_many_na <- apply(df_matrix, 1, function(x) {
        sum(is.na(x)) / length(x) > na_threshold
      })
      
      if (sum(genes_with_many_na) > 0) {
        df_matrix <- df_matrix[!genes_with_many_na, , drop = FALSE]
        cat("Removed", sum(genes_with_many_na), "genes with >", na_threshold*100, "% NA values\n")
      }
      
      # Impute remaining NA values with gene means
      remaining_na <- sum(is.na(df_matrix))
      if (remaining_na > 0) {
        for (i in 1:nrow(df_matrix)) {
          na_positions <- is.na(df_matrix[i, ])
          if (any(na_positions)) {
            gene_mean <- mean(df_matrix[i, ], na.rm = TRUE)
            if (is.finite(gene_mean)) {
              df_matrix[i, na_positions] <- gene_mean
            } else {
              df_matrix[i, na_positions] <- 0
            }
          }
        }
        cat("Imputed", remaining_na, "remaining NA values with gene means\n")
      }
    }
  }
  
  # Step 5: Handle all-zero genes
  if (remove_all_zero_genes) {
    all_zero_genes <- apply(df_matrix, 1, function(x) all(x == 0, na.rm = TRUE))
    if (sum(all_zero_genes) > 0) {
      df_matrix <- df_matrix[!all_zero_genes, , drop = FALSE]
      cat("Removed", sum(all_zero_genes), "genes with all zero values\n")
    }
  }
  
  # Final validation
  final_dims <- dim(df_matrix)
  cat("Final dimensions:", final_dims[1], "genes x", final_dims[2], "tissues\n")
  cat("Genes removed:", original_dims[1] - final_dims[1], "\n")
  cat("Samples removed:", original_dims[2] - final_dims[2], "\n")
  
  # Check for remaining problematic values
  remaining_na <- sum(is.na(df_matrix))
  remaining_inf <- sum(is.infinite(df_matrix))
  remaining_nan <- sum(is.nan(df_matrix))
  
  if (remaining_na > 0 || remaining_inf > 0 || remaining_nan > 0) {
    cat("WARNING: Still have problematic values:\n")
    cat("  NA:", remaining_na, "\n")
    cat("  Inf:", remaining_inf, "\n") 
    cat("  NaN:", remaining_nan, "\n")
  } else {
    cat("✓ Data is clean and ready for Fortran normalization\n")
  }
  
  return(df_matrix)
}

#' Calculate signed clock hand angle between two normalized vectors
#' 
#' @param v1 First normalized vector in RAP space
#' @param v2 Second normalized vector in RAP space
#' @param selected_axes_for_signed Indices of 3 axes for directionality (ignored if dim <= 3)
#' @return Signed angle in radians [-π, π]
tox_clock_hand_angle_between_vectors <- function(v1, v2, selected_axes_for_signed = c(1, 2, 3)) {
  n_dims <- length(v1)
  if (length(v2) != n_dims) {
    stop("Vectors must have the same dimension")
  }
  
  # Ensure selected_axes_for_signed has exactly 3 elements
  if (length(selected_axes_for_signed) != 3) {
    stop("selected_axes_for_signed must have exactly 3 elements")
  }
  
  # Call Fortran wrapper
  result <- .Fortran("clock_hand_angle_between_vectors_r",
                    v1 = as.double(v1),
                    v2 = as.double(v2),
                    n_dims = as.integer(n_dims),
                    signed_angle = as.double(0),
                    selected_axes_for_signed = as.integer(selected_axes_for_signed))
  
  return(result$signed_angle)
}

#' Calculate signed clock hand angles for multiple vector pairs
#' 
#' @param origins Matrix of origin vectors (n_dims x n_vecs)
#' @param targets Matrix of target vectors (n_dims x n_vecs)
#' @param vecs_selection_mask Logical vector indicating which pairs to compute
#' @param selected_axes_for_signed Indices of 3 axes for directionality
#' @return Vector of signed angles in radians [-π, π]
tox_clock_hand_angles_for_shift_vectors <- function(origins, targets, 
                                               vecs_selection_mask = NULL,
                                               selected_axes_for_signed = c(1, 2, 3)) {
  if (!is.matrix(origins) || !is.matrix(targets)) {
    stop("origins and targets must be matrices")
  }
  
  if (!identical(dim(origins), dim(targets))) {
    stop("origins and targets must have the same dimensions")
  }
  
  n_dims <- nrow(origins)
  n_vecs <- ncol(origins)
  
  # Default selection mask (all TRUE)
  if (is.null(vecs_selection_mask)) {
    vecs_selection_mask <- rep(TRUE, n_vecs)
  }
  
  if (length(vecs_selection_mask) != n_vecs) {
    stop("vecs_selection_mask length must equal number of vector pairs")
  }
  
  n_selected_vecs <- sum(vecs_selection_mask)
  
  # Ensure selected_axes_for_signed has exactly 3 elements
  if (length(selected_axes_for_signed) != 3) {
    stop("selected_axes_for_signed must have exactly 3 elements")
  }
  
  # Call Fortran wrapper
  result <- .Fortran("clock_hand_angles_for_shift_vectors_r",
                    origins = as.double(origins),
                    targets = as.double(targets),
                    n_dims = as.integer(n_dims),
                    n_vecs = as.integer(n_vecs),
                    vecs_selection_mask = as.logical(vecs_selection_mask),
                    n_selected_vecs = as.integer(n_selected_vecs),
                    selected_axes_for_signed = as.integer(selected_axes_for_signed),
                    signed_angles = as.double(rep(0, n_selected_vecs)))
  
  return(result$signed_angles)
}

#' Calculate relative axis contributions from a shift vector
#'
#' This function wraps the Fortran subroutine `relative_axes_changes_from_shift_vector_r`
#' to compute the relative axis contributions for a given shift vector in RAP space.
#'
#' @param shift_vector Numeric vector representing the shift in RAP space.
#' @return Numeric vector of relative axis contributions (sums to 1).
relative_axes_changes_from_shift_vector <- function(shift_vector) {
  n_dims <- length(shift_vector)
  contrib <- numeric(n_dims)
  result <- .Fortran("relative_axes_changes_from_shift_vector_r",
                     as.double(shift_vector),
                     as.integer(n_dims),
                     contrib)
  return(result[[3]])
}

#' Calculate relative axis contributions from an expression vector
#'
#' This function wraps the Fortran subroutine `relative_axes_expression_from_expression_vector_r`
#' to compute the relative axis contributions for a given expression vector in RAP space.
#'
#' @param expression_vector Numeric vector representing the expression in RAP space.
#' @return Numeric vector of relative axis contributions (sums to 1).
relative_axes_expression_from_expression_vector <- function(expression_vector) {
  n_dims <- length(expression_vector)
  contrib <- numeric(n_dims)
  result <- .Fortran("relative_axes_expression_from_expression_vector_r",
                     as.double(expression_vector),
                     as.integer(n_dims),
                     contrib)
  return(result[[3]])
}
