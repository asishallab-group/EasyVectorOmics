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
# OUTLIER DETECTION FUNCTIONS 
# ===================================================================
#' Complete outlier detection workflow
#'
#' This function performs the complete outlier detection workflow:
#' 1. Computes family scaling factors using LOESS
#' 2. Calculates RDI values
#' 3. Identifies outliers based on percentile threshold
#'
#' @param distances Numeric vector of gene distances
#' @param gene_to_fam Integer vector mapping genes to family indices
#' @param n_families Integer number of families
#' @param percentile Percentile threshold for outlier detection (default: 95.0)
#' @return List with components:
#'   - is_outlier: Logical vector indicating outliers
#'   - loess_x: Family median distances
#'   - loess_y: Family standard deviations
#'   - loess_n: Number of genes used per family
tox_detect_outliers <- function(distances, gene_to_fam, n_families, percentile = 95.0) {
  # Input validation
  validate_numeric_vector(distances, "distances")
  n_genes <- as.integer(length(distances))

  # Call Rcpp wrapper
  result <- tox_detect_outliers_rcpp(distances, gene_to_fam, n_families, percentile)

  # Check for error (ierr returned by the Rcpp wrapper)
  check_err_code(result$ierr)

  # Return structured result
  return(list(
    is_outlier = result$is_outlier,
    loess_x = result$loess_x,
    loess_y = result$loess_y,
    loess_n = result$loess_n
  ))

}

#' Compute family scaling factors using LOESS smoothing
#'
#' This function calculates scaling factors for gene families based on distance distributions.
#' It uses LOESS smoothing to estimate the relationship between family median distances
#' and their standard deviations.
#'
#' @param distances Numeric vector of gene distances
#' @param gene_to_fam Integer vector mapping genes to family indices
#' @param n_families Integer number of families
#' @return List with components:
#'   - dscale: Scaling factors for each family
#'   - loess_x: Family median distances 
#'   - loess_y: Family standard deviations
#'   - indices_used: Number of genes used per family
tox_compute_family_scaling <- function(distances, gene_to_fam, n_families) {
  # Input validation
  validate_numeric_vector(distances, "distances")
  n_genes <- as.integer(length(distances))
  validate_length_equals_n(gene_to_fam, n_genes, "gene_to_fam")
  validate_index_bounds(gene_to_fam, low = 1, high = n_families, name = "gene_to_fam")

  # Call the Rcpp forwarder.
  result <- tox_compute_family_scaling_rcpp(distances, gene_to_fam, n_families)

  # Check for error
  if (result$ierr != 0) {
    check_err_code(result$ierr)
  }

  # Return structured result
  return(list(
    dscale = as.numeric(result$dscale),
    loess_x = as.numeric(result$loess_x),
    loess_y = as.numeric(result$loess_y),
    indices_used = as.integer(result$indices_used)
  ))
}

#' Compute family scaling factors using LOESS smoothing (Expert Version)
#'
#' Expert version of compute_family_scaling with user-provided work arrays.
#' This version requires pre-allocated work arrays for maximum performance and control.
#' Use this when you need fine-grained control over memory allocation or are calling
#' this function many times in a tight loop.
#'
#' @param distances Numeric vector of gene distances
#' @param gene_to_fam Integer vector mapping genes to family indices
#' @param n_families Integer number of families
#' @param perm_tmp Pre-allocated permutation array for sorting (n_genes)
#' @param stack_left_tmp Pre-allocated stack array for sorting (n_genes)
#' @param stack_right_tmp Pre-allocated stack array for sorting (n_genes)
#' @param family_distances Pre-allocated work array for family distances (n_genes)
#' @return List with components:
#'   - dscale: Scaling factors for each family
#'   - loess_x: Family median distances 
#'   - loess_y: Family standard deviations
#'   - indices_used: Number of genes used per family
#'   - perm_tmp: Final state of permutation array
#'   - stack_left_tmp: Final state of left stack array
#'   - stack_right_tmp: Final state of right stack array
#'   - family_distances: Final state of family distances array
tox_compute_family_scaling_expert <- function(distances, gene_to_fam, n_families,
                                              perm_tmp, stack_left_tmp, stack_right_tmp,
                                              family_distances) {
# Input validation
  validate_numeric_vector(distances, "distances")
  n_genes <- as.integer(length(distances))
  validate_length_equals_n(gene_to_fam, n_genes, "gene_to_fam")
  validate_index_bounds(gene_to_fam, low = 1, high = n_families, name = "gene_to_fam")
  validate_length_equals_n(perm_tmp, n_genes, "perm_tmp")
  validate_length_equals_n(stack_left_tmp, n_genes, "stack_left_tmp")
  validate_length_equals_n(stack_right_tmp, n_genes, "stack_right_tmp")
  validate_length_equals_n(family_distances, n_genes, "family_distances")

  # Call the Rcpp forwarder.
  result <- tox_compute_family_scaling_expert_rcpp(n_families, distances, gene_to_fam, perm_tmp, stack_left_tmp, stack_right_tmp,
    family_distances
  )

   # Check for error
  if (result$ierr != 0) {
    check_err_code(result$ierr)
  }
  # Return structured result
  return(list(
    dscale = result$dscale,
    loess_x = result$loess_x,
    loess_y = result$loess_y,
    indices_used = result$indices_used,
    perm_tmp     = result$perm_tmp,
    stack_left_tmp   = result$stack_left_tmp,
    stack_right_tmp  = result$stack_right_tmp,
    family_distances = result$family_distances

  ))
}

#' Compute Relative Distance Index (RDI) for genes
#'
#' This function calculates the Relative Distance Index for each gene,
#' which is the absolute distance divided by the family scaling factor.
#'
#' @param distances Numeric vector of gene distances
#' @param gene_to_fam Integer vector mapping genes to family indices
#' @param dscale Numeric vector of scaling factors for each family
#' @return List with components:
#'   - rdi: Relative Distance Index for each gene
#'   - sorted_rdi: RDI values sorted in ascending order
tox_compute_rdi <- function(distances, gene_to_fam, dscale) {
# Input validation
  validate_numeric_vector(distances, "distances")
  n_genes <- as.integer(length(distances))
  validate_length_equals_n(gene_to_fam, n_genes, "gene_to_fam")
  n_families <- as.integer(length(dscale))
  validate_index_bounds(gene_to_fam, low = 1, high = n_families, name = "gene_to_fam")

  # Call Rcpp forwarder
  result <- tox_compute_rdi_rcpp(distances, gene_to_fam, dscale)

  # Return 
  return(list(
    rdi = result$rdi,
    sorted_rdi = result$sorted_rdi
  ))
}

#' Identify outliers based on RDI percentiles
#'
#' This function identifies outliers by comparing each gene's RDI value
#' against a percentile threshold of the sorted RDI distribution.
#'
#' @param rdi Numeric vector of RDI values
#' @param percentile Percentile threshold (default: 95.0)
#' @return List with components:
#'   - is_outlier: Logical vector indicating outliers
#'   - threshold: The RDI threshold value used
tox_identify_outliers <- function(rdi, percentile = 95.0) {
  #Calculate dimensions
  n_genes <- as.integer(length(rdi))
  # Call the Rcpp forwarder
  result <- tox_identify_outliers_rcpp(rdi, percentile)
  # Return structured result
  return(list(
    is_outlier = result$is_outlier,
    threshold = result$threshold
  ))

}
