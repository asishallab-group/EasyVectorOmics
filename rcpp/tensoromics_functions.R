library(Rcpp)

# Get absolute path to build directory containing the compiled Fortran library
lib_path <- normalizePath("build")

# Set up compilation flags for linking with Fortran library
Sys.setenv(PKG_LIBS = paste0("-Wl,-rpath,", lib_path, " -L", lib_path, " -ltensor-omics -lgfortran"))

# Compile and load all TensorOmics Rcpp wrapper functions (includes error_handling.cpp)
sourceCpp("rcpp/tensoromics_functions.cpp", env = .GlobalEnv)

cat("✓ TensorOmics Rcpp functions loaded successfully\n")

source("r/error_handling.R")


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
  if (!is.numeric(vec1) || !is.numeric(vec2)) {
    stop("Both vectors must be numeric")
  }
  if (length(vec1) != length(vec2)) {
    stop("Vectors must have the same length")
  }
  if (length(vec1) == 0) {
    stop("Vectors cannot be empty")
  }
  
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
  # Input validation
  if (!is.numeric(genes) || !is.numeric(centroids)) {
    stop("genes and centroids must be numeric")
  }
  if (!is.numeric(gene_to_fam) && !is.integer(gene_to_fam)) {
    stop("gene_to_fam must be numeric or integer")
  }
  if (!is.numeric(d) && !is.integer(d)) {
    stop("d must be numeric or integer")
  }
  
  # Convert to appropriate types
  genes <- as.numeric(genes)
  centroids <- as.numeric(centroids)
  gene_to_fam <- as.integer(gene_to_fam)
  d <- as.integer(d)
  
  # Calculate dimensions
  n_genes <- as.integer(length(genes) / d)
  n_families <- as.integer(length(centroids) / d)
  
  # Validate dimensions
  if (length(genes) %% d != 0) {
    stop("Length of genes must be divisible by d")
  }
  if (length(centroids) %% d != 0) {
    stop("Length of centroids must be divisible by d")
  }
  if (length(gene_to_fam) != n_genes) {
    stop("Length of gene_to_fam must equal number of genes")
  }
  if (any(gene_to_fam < 0)) {
    stop("gene_to_fam indices must be between 0 and n_families (0 = no family assignment)")
  }
  
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
  # Input validation
  if (!is.matrix(expression_vectors)) {
    stop("expression_vectors must be a matrix")
  }
  if (!is.logical(vector_selection) && !is.numeric(vector_selection)) {
    stop("vector_selection must be logical or numeric")
  }
  if (!is.logical(axis_selection) && !is.numeric(axis_selection)) {
    stop("axis_selection must be logical or numeric")
  }
  
  # Convert to appropriate types for Rcpp
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
  
  # Validate dimensions
  if (length(vector_selection) != ncol(expression_vectors)) {
    stop("vector_selection length must match number of columns in expression_vectors")
  }
  if (length(axis_selection) != nrow(expression_vectors)) {
    stop("axis_selection length must match number of rows in expression_vectors")
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
  # Input validation
  if (!is.matrix(expression_vectors)) {
    stop("expression_vectors must be a matrix")
  }

  if (!is.matrix(family_centroids)) {
    stop("family_centroids must be a matrix")
  }

  # Validate length of gene_to_centroid
  if (n_vectors != length(gene_to_centroid)) {
    stop("number of expression_vectors must be equal to length of gene_to_centroid")
  }

  # Validate dimensions
  if (n_axes_genes != n_axes_centroids) {
    stop("family_centroids must have the same number of axes as expression_vectors")
  }
  

  
  # Call Rcpp wrapper
  result <- tox_compute_shift_vector_field_rcpp(expression_vectors, family_centroids, gene_to_centroid)
  
  # Check for errors and throw informative messages
  check_err_code(result$ierr)
  
  return(list(
    shift_vectors = result$shift_vectors
  ))
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
  
  # 1) Validate inputs
  if (!is.matrix(expression_vectors) || !is.numeric(expression_vectors)) {
    stop("`expression_vectors` must be a numeric matrix.")
  }
 
  if (!is.integer(gene_to_family) || length(gene_to_family) != n_genes) {
    stop("`gene_to_family` must be an integer vector of length n_genes.")
  }
  if (!is.logical(ortholog_set) || length(ortholog_set) != n_genes) {
    stop("`ortholog_set` must be a logical vector of length n_genes.")
  }
  if (!mode %in% c('all', 'ortho')) {
    stop("`mode` must be either 'all' or 'ortho'.")
  }

 
  # 3) Call Rcpp wrapper
  result <- tox_group_centroid_rcpp(expression_vectors, gene_to_family, n_families, ortholog_set, mode = 'all')
  
  # Check for errors and throw informative messages
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
  # Validate inputs
  if (!is.matrix(expression_vectors) || !is.numeric(expression_vectors)) {
    stop("`expression_vectors` must be a numeric matrix.")
  }
 
  if (!is.integer(gene_indices) || any(gene_indices < 1) || any(gene_indices > n_genes)) {
    stop("`gene_indices` must be integer indices between 1 and n_genes.")
  }

 
  # Call Rcpp wrapper
  result <- tox_mean_vector_rcpp(expression_vectors, gene_indices)
  
  # Check for errors and throw informative messages
  check_err_code(result$ierr)
  

  return(result)
}
