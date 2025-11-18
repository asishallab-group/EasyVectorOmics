library(Rcpp)

# Get absolute path to build directory containing the compiled Fortran library
lib_path <- normalizePath("build")

# Set up compilation flags for linking with Fortran library
Sys.setenv(PKG_LIBS = paste0("-Wl,-rpath,", lib_path, " -L", lib_path, " -ltensor-omics -lgfortran"))

# Compile and load all TensorOmics Rcpp wrapper functions (includes error_handling.cpp)
sourceCpp("rcpp/tensoromics_functions.cpp", env = .GlobalEnv)

cat("✓ TensorOmics Rcpp functions loaded successfully\n")
## Use the rcpp-local validators (we keep all edits inside rcpp/ per workspace policy)
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
  # Perform R-layer validation using rcpp-local validators, then forward to Rcpp.
  validate_numeric_vector(vec1, "vec1")
  validate_numeric_vector(vec2, "vec2")
  validate_equal_length(vec1, vec2, "vec1", "vec2")
  validate_nonempty_vector(vec1, "vec1")

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
  # R-layer validation in rcpp/ (kept here because r/ must not be changed)
  validate_numeric_vector(genes, "genes")
  validate_numeric_vector(centroids, "centroids")
  validate_positive_integer_scalar(d, "d")

   # Convert to appropriate types
  genes <- as.numeric(genes)
  centroids <- as.numeric(centroids)
  gene_to_fam <- as.integer(gene_to_fam)
  d <- as.integer(d)

#  # Validate flattened lengths are compatible with d
  validate_divisible_length(genes, d, "genes")
  validate_divisible_length(centroids, d, "centroids")
  # Calculate dimensions
  n_genes <- as.integer(length(genes) / d)
  n_families <- as.integer(length(centroids) / d)
  validate_gene_to_family(gene_to_fam, n_genes, n_families, "gene_to_fam")
  validate_length_equals_n(gene_to_fam, n_genes, "gene_to_fam")
  
  
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
  # R-layer validation (kept in rcpp/ to avoid touching r/)
  validate_numeric_matrix(expression_vectors, "expression_vectors")
  n_axes <- nrow(as.matrix(expression_vectors))
  n_vectors <- ncol(as.matrix(expression_vectors))

  # Ensure selection vectors have correct lengths and types
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

  validate_logical_vector(as.logical(vector_selection), "vector_selection", expected_length = n_vectors)
  validate_logical_vector(as.logical(axis_selection), "axis_selection", expected_length = n_axes)

  result <- tox_calculate_tissue_versatility_rcpp(as.matrix(expression_vectors), vector_selection, axis_selection)
  if (result$ierr != 0) check_err_code(result$ierr)

  return(list(
    tissue_versatilities = result$tissue_versatilities,
    tissue_angles_deg = result$tissue_angles_deg,
    n_selected_vectors = result$n_selected_vectors,
    n_selected_axes = result$n_selected_axes
  ))
}
# ===================================================================
# SHIFT VECTOR FIELD FUNCTIONS
# ===================================================================
#' Calculate Shift Vector Field 
#' Computes the shift vector field for each gene expression vector based on its family centroid.
#' The shift vector is defined as the difference between the gene expression vector and its corresponding family centroid,
#' starting at the expression vector and pointing to its family centroid.
#' This function automatically checks for errors and throws informative exceptions.
#'
#' @param expression_vectors: Matrix where each column is a gene expression vector (n_axes x n_vectors)
#' @param family_centroids: Matrix where each column is a family centroid vector (n_axes x n_families)
#' @param gene_to_centroid: Array mapping each gene to its corresponding family centroid ID in family_centroids (length n_vectors)
#' 
#' @return List containing:
#'   \item{shift_vectors}{The computed shift vectors for each gene expression vector}
#'

tox_compute_shift_vector_field <- function(expression_vectors, family_centroids, gene_to_centroid) {
  # R-layer validation (kept in rcpp/)
  validate_numeric_matrix(expression_vectors, "expression_vectors")
  validate_numeric_matrix(family_centroids, "family_centroids")
  # Ensure matching axes (rows)
  validate_matching_rows(as.matrix(expression_vectors), as.matrix(family_centroids), "expression_vectors", "family_centroids")

  # gene_to_centroid should be integer vector with length equal to number of vectors
  gene_to_centroid <- as.integer(gene_to_centroid)
  n_vectors <- ncol(as.matrix(expression_vectors))
  validate_length_equals_n(gene_to_centroid, n_vectors, "gene_to_centroid")
  if (any(is.na(gene_to_centroid))) stop("`gene_to_centroid` must not contain NA values.")
  if (any(gene_to_centroid < 0L)) stop("`gene_to_centroid` must not contain negative indices.")

  result <- tox_compute_shift_vector_field_rcpp(as.matrix(expression_vectors), as.matrix(family_centroids), as.integer(gene_to_centroid))
  check_err_code(result$ierr)
  return(list(shift_vectors = result$shift_vectors))
}

# ===================================================================
# GENE CENTROIDS FUNCTIONS
# ===================================================================
#' Calculate Gene Centroids

#' Computes the centroids for each gene family based on the expression vectors of its member genes.
#' This function automatically checks for errors and throws informative exceptions.
#'
#' @param expression_vectors: Matrix where each column is a gene expression vector (n_axes x n_vectors)
#' @param gene_to_family: Array mapping each gene to its corresponding family ID (length n_vectors)
#' @param n_families: Total number of gene families
#' @param ortholog_set: Logical array indicating if a gene is part of a specific subset (e.g., orthologs)
#' @param mode: Character string indicating the mode of operation ('all' or 'ortho')
#'
#' @return List containing:
#'   \item{centroid_matrix}{The computed centroids for each gene family}
#'

 
tox_group_centroid <- function(expression_vectors, gene_to_family, n_families, ortholog_set, mode = 'all') {
  # R-layer validation (kept in rcpp/)
  validate_group_centroid_inputs(as.matrix(expression_vectors), gene_to_family, n_families, ortholog_set, mode)

  result <- tox_group_centroid_rcpp(as.matrix(expression_vectors), as.integer(gene_to_family), as.integer(n_families), as.integer(ortholog_set), as.integer(ifelse(mode == 'all', 1, 0)))
  check_err_code(result$ierr)
  return(result)
}

#' Compute the element-wise mean for a given set of gene expression vectors
#'
#' This function wraps the Fortran subroutine `mean_vector_r`
#' to compute the centroid (mean vector) for a selected set of genes.
#'
#' @param expression_vectors Numeric matrix (n_axes x n_genes) of gene expression vectors
#' @param gene_indices Integer vector of column indices of selected genes (1-based)
#'
#' @return Numeric vector of length n_axes representing the computed centroid
#'

tox_mean_vector <- function(expression_vectors, gene_indices) {
  # R-layer validation (kept in rcpp/)
  validate_mean_vector_inputs(as.matrix(expression_vectors), as.integer(gene_indices))

  result <- tox_mean_vector_rcpp(as.matrix(expression_vectors), as.integer(gene_indices))
  check_err_code(result$ierr)
  return(result)
}
