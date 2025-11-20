library(Rcpp)

# Get absolute path to build directory containing the compiled Fortran library

lib_path <- normalizePath("build")

# Set up compilation flags for linking with Fortran library
Sys.setenv(PKG_LIBS = paste0("-Wl,-rpath,", lib_path, " -L", lib_path, " -ltensor-omics -lgfortran"))

# Compile and load all TensorOmics Rcpp wrapper functions (includes error_handling.cpp)
sourceCpp("rcpp/tensoromics_functions.cpp", env = .GlobalEnv)

cat("✓ TensorOmics Rcpp functions loaded successfully\n")

source("rcpp/error_handling.R")

# ===================================================================
# EUCLIDEAN DISTANCE FUNCTIONS
# ===================================================================

#' Calculate Euclidean distance between two vectors
#' 
#' Computes the Euclidean distance between two vectors of the same dimension.
#' This function automatically checks for errors and throws informative exceptions.
#' 
#' @param vec1 First vector (numeric)
#' @param vec2 Second vector (numeric, same length as vec1)
#' 
#' @return Numeric value representing the Euclidean distance between the vectors
#' 
tox_euclidean_distance <- function(vec1, vec2) {
  # Input validation
  validate_numeric_vector(vec1, "vec1")
  validate_numeric_vector(vec2, "vec2")
  validate_same_length(vec1, vec2, "vec1", "vec2")
  validate_nonempty(vec1, "vec1")

  # Call Rcpp wrapper 
  return(tox_euclidean_distance_rcpp(as.numeric(vec1), as.numeric(vec2)))
}


#' Calculate distances from genes to their family centroids
#' 
#' Computes the Euclidean distance from each gene to its corresponding family centroid.
#' This function automatically checks for errors and throws informative exceptions.
#' 
#' @param genes Matrix of gene expression data (genes as columns, dimensions as rows)
#' @param centroids Matrix of family centroids (families as columns, dimensions as rows) 
#' @param gene_to_fam Integer vector mapping each gene to its family index (1-based)
#' @param d Integer number of dimensions
#' 
#' @return Numeric vector of distances from each gene to its family centroid
#' 
tox_distance_to_centroid <- function(genes, centroids, gene_to_fam, d) {
  #Convert to appropriate types
  genes <- as.numeric(genes)
  centroids <- as.numeric(centroids)
  gene_to_fam <- as.integer(gene_to_fam)
  d <- as.integer(d)

  # Input validation
  validate_positive_integer_scalar(d, "d")
  validate_divisible_length(genes, d, "genes")
  validate_divisible_length(centroids, d, "centroids")

  n_genes <- as.integer(length(genes) / d)
  n_families <- as.integer(length(centroids) / d)

  validate_gene_to_family(gene_to_fam, n_genes, n_families, "gene_to_fam")

  # Call Rcpp wrapper
  return(tox_distance_to_centroid_rcpp(genes, centroids, gene_to_fam, d))
}


#' Calculate Tissue Versatility
#' 
#' Computes normalized tissue versatility for selected expression vectors.
#' The metric is based on the angle between each gene expression vector and the space diagonal.
#' Versatility is normalized to [0, 1], where 0 means uniform expression and 1 means expression in only one axis.
#' This function automatically checks for errors and throws informative exceptions.
#' 
#' @param expression_vectors Matrix where each column is a gene expression vector (n_axes x n_vectors)
#' @param vector_selection Logical vector indicating which vectors to process (length n_vectors)
#' @param axis_selection Logical vector indicating which axes to include in calculation (length n_axes)
#' 
#' @return List containing:
#'   \item{tissue_versatilities}{Normalized tissue versatility values [0,1] for selected vectors}
#'   \item{tissue_angles_deg}{Angles in degrees [0,90] for selected vectors}
#'   \item{n_selected_vectors}{Number of vectors processed}
#'   \item{n_selected_axes}{Number of axes used in calculation}
#' 
tox_calculate_tissue_versatility <- function(expression_vectors, vector_selection, axis_selection) {
  #Input validation
  validate_numeric_matrix(expression_vectors, "expression_vectors")

  # Ensure selectors have expected lengths
  validate_logical_or_index_vector(vector_selection, expected_length = ncol(expression_vectors), name = "vector_selection")
  validate_logical_or_index_vector(axis_selection, expected_length = nrow(expression_vectors), name = "axis_selection")

  #Convert to appropriate types for Rcpp
  if (is.numeric(vector_selection)) {
    vector_selection <- as.integer(as.logical(vector_selection))
  } else {
    vector_selection <- as.integer(vector_selection)
  }
  
  if (is.numeric(axis_selection)) {
    axis_selection <- as.integer(as.logical(axis_selection))
  } else {
    axis_selection <- as.integer(axis_selection)
  }

  # Call Rcpp wrapper
  result <- tox_calculate_tissue_versatility_rcpp(expression_vectors, vector_selection, axis_selection)
  
  # Check for errors
  if (result$ierr != 0) {
    check_err_code(result$ierr)
  }
  
  # Return structured result 
  return(list(
    tissue_versatilities = result$tissue_versatilities,
    tissue_angles_deg = result$tissue_angles_deg,
    n_selected_vectors = result$n_selected_vectors,
    n_selected_axes = result$n_selected_axes
  ))

}
# ===================================================================
# NORMALIZATION FUNCTIONS
# ===================================================================
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
#' normalized_matrix <- tox_normalize_by_std_dev(input_matrix)
tox_normalize_by_std_dev <- function(input_matrix) {
  # Validate input matrix values (NA / Inf / NaN)
  validate_numeric_matrix_values(input_matrix, "input_matrix")
  result <- tox_normalize_by_std_dev_rcpp(input_matrix)
  if (result$ierr != 0) {
    check_err_code(result$ierr)
  }

  return(matrix(result$output_vector, nrow = nrow(input_matrix), ncol = ncol(input_matrix), dimnames = dimnames(input_matrix)))
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
#' normalized_matrix <- tox_quantile_normalization(input_matrix)
tox_quantile_normalization <- function(input_matrix) {
  validate_matrix(input_matrix, "input_matrix")
  result <- tox_quantile_normalization_rcpp(input_matrix)

  if (result$ierr != 0) {
    check_err_code(result$ierr)
  }
 
  return(matrix(result$output_vector, nrow = n_genes, ncol = n_tissues, dimnames = dimnames(input_matrix)))
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
#' log_matrix <- tox_log2_transformation(input_matrix)
tox_log2_transformation <- function(input_matrix) {
  validate_matrix(input_matrix, "input_matrix")
  result <- tox_log2_transformation_rcpp(input_matrix)
  
  if (result$ierr != 0) {
    check_err_code(result$ierr)
  }
  
  return(matrix(result$output_vector, nrow = n_genes, ncol = n_tissues, dimnames = dimnames(input_matrix)))
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
#' averaged_df <- tox_calculate_tissue_averages(df)
tox_calculate_tissue_averages <- function(df) {
  validate_matrix(as.matrix(df), "df")
  result <- tox_calc_tiss_avg_rcpp(as.matrix(df))
  
  if (result$ierr != 0) {
    check_err_code(result$ierr)
  }
  n_genes <- nrow(df)
  n_cols <- length(result$output_vector) / n_genes
  output_matrix <- matrix(result$output_vector, nrow = n_genes, ncol = n_cols)
  return(as.data.frame(output_matrix))
}


#' Calculate log2 fold changes based on control and condition patterns
#' @param df A data frame with genes as rows and tissues/conditions as columns.
#' @param control_pattern A string pattern to detect control columns.
#' @param condition_patterns A character vector with patterns to detect condition columns.
tox_calculate_fc_by_patterns <- function(df, control_pattern, condition_patterns) {
  validate_matrix(as.matrix(df), "df")
  validate_string_scalar(control_pattern, "control_pattern")
  validate_character_vector(condition_patterns, "condition_patterns")

  result <- tox_calc_fchange_rcpp(as.matrix(df), control_pattern, condition_patterns)
  
  if (result$ierr != 0) {
    check_err_code(result$ierr)
  }

  n_genes <- nrow(df)
  n_conditions <- length(condition_patterns)
  output_matrix <- matrix(result$output_vector, nrow = n_genes, ncol = n_conditions)
  return(as.data.frame(output_matrix))
}


#' Complete normalization pipeline for gene expression data (up to log2(x+1))
#' @param input_matrix Numeric matrix (genes x tissues)
#' @param group_s Integer vector: start column index for each replicate group (1-based)
#' @param group_c Integer vector: number of columns per replicate group
tox_normalization_pipeline <- function(input_matrix, group_s, group_c) {
  validate_matrix(input_matrix, "input_matrix")
  group_s <- as.integer(group_s)
  group_c <- as.integer(group_c)
  validate_group_vectors(group_s, group_c, ncol(input_matrix))

  result <- tox_normalization_pipeline_rcpp(input_matrix, group_s, group_c)
  
  if (result$ierr != 0) {
    check_err_code(result$ierr)
  }

  return(matrix(result$buf_log, nrow = nrow(input_matrix), ncol = length(group_s)))
}
