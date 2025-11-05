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
# OUTLIER DETECTION FUNCTIONS 
# ===================================================================
#' Detect Gene Outliers
#' 
#' Identifies gene outliers based on relative distance indices (RDI) and LOESS-smoothed
#' family scaling factors.
#' The function computes family-wise median and standard deviation of gene-to-centroid distances,
#' applies LOESS smoothing, and marks genes whose relative distance exceeds a specified percentile threshold as outliers.
#' This function automatically checks for errors and throws informative exceptions.
#' 
#' @param distances Numeric vector of gene-to-centroid distances (length = n_genes)
#' @param gene_to_fam Integer vector mapping each gene to a family (0 = no family, >0 = family index)
#' @param n_families Integer scalar representing the total number of families
#' @param percentile Numeric percentile threshold for outlier detection (default = 95.0)
#' 
#' @return List containing:
#'   \item{is_outlier}{Logical vector indicating whether each gene is an outlier}
#'   \item{loess_x}{Numeric vector of median distances used for LOESS smoothing}
#'   \item{loess_y}{Numeric vector of standard deviations used for LOESS smoothing}
#'   \item{loess_n}{Integer vector with number of genes per family used in LOESS computation}
#' 
tox_detect_outliers <- function(distances, gene_to_fam, n_families, percentile = 95.0) {
  # Input validation
  if (!is.numeric(distances)) {
    stop("distances must be numeric")
  }
  if (!is.numeric(gene_to_fam) && !is.integer(gene_to_fam)) {
    stop("gene_to_fam must be numeric or integer")
  }
  if (!is.numeric(n_families) && !is.integer(n_families)) {
    stop("n_families must be numeric or integer")
  }
  if (!is.numeric(percentile)) {
    stop("percentile must be numeric")
  }
  
  # Convert to appropriate types for Rcpp
  distances   <- as.numeric(distances)
  gene_to_fam <- as.integer(gene_to_fam)
  n_families  <- as.integer(n_families)
  percentile  <- as.numeric(percentile)
  
  # Calculate dimensions
  n_genes <- as.integer(length(distances))

  # Validate dimensions
  if (length(gene_to_fam) != n_genes) {
    stop("length of gene_to_fam must match number of distances")
  }
  if (length(distances) == 0) {
    stop("distances cannot be empty")
  }
  if (any(!is.finite(distances))) {
    stop("distances contains non-finite values")
  }
  if (n_families < 1L) {
    stop("n_families must be >= 1")
  }
  if (percentile < 0 || percentile > 100) {
    stop("percentile must be between 0 and 100")
  }
  
  # Call Rcpp wrapper
  result <- tox_detect_outliers_rcpp(distances, gene_to_fam, n_families, percentile)
  
  # Check for backend error
  if (result$ierr != 0) {
    check_err_code(result$ierr)
  }
  
  # Return structured result
  return(list(
    is_outlier = result$is_outlier,
    loess_x    = result$loess_x,
    loess_y    = result$loess_y,
    loess_n    = result$loess_n
  ))
}

#' Compute per-family scaling factors and LOESS reference points
#'
#' Computes a scaling factor for each gene family by calculating the median
#' and standard deviation of gene-to-centroid distances within each family,
#' then applying LOESS smoothing to these summary values. 
#' This function automatically checks for errors and throws informative exceptions.
#'
#' @param distances Numeric vector of distances for each gene to its family centroid.
#' @param gene_to_fam Integer vector mapping each gene to a family index (0 = no family, >0 = family index).
#'        Must have the same length as `distances`.
#' @param n_families Integer number of families (>= 1).
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{dscale}{Numeric vector of scaling factors per family}
#'     \item{loess_x}{Numeric vector of x coordinates used as reference points for LOESS (median distances)}
#'     \item{loess_y}{Numeric vector of y coordinates used as reference points for LOESS (standard deviations)}
#'     \item{indices_used}{Integer vector of family indices that were included in the calculation}
#'
#'   }
#'
tox_compute_family_scaling <- function(distances, gene_to_fam, n_families) {

  # Input validation
  if (!is.numeric(distances)) {
    stop("distances must be numeric")
  }
  if (!is.numeric(gene_to_fam) && !is.integer(gene_to_fam)) {
    stop("gene_to_fam must be numeric or integer")
  }
  if (!is.numeric(n_families) && !is.integer(n_families)) {
    stop("n_families must be numeric or integer")
  }
  
  # Convert to appropriate types
  distances <- as.numeric(distances)
  gene_to_fam <- as.integer(gene_to_fam)
  n_families <- as.integer(n_families)

  # Dimensions 
  n_genes <- as.integer(length(distances))

  # Validate dimensions 
  if (n_genes == 0L) {
    stop("distances cannot be empty")  
  }
  if (any(!is.finite(distances))) {
    stop("distances contains non-finite values")  
  }
  if (length(gene_to_fam) != n_genes) {
    stop("length of gene_to_fam must match number of distances")  
  }
  if (n_families < 1L) {
    stop("n_families must be >= 1")  
  }
  if (any(gene_to_fam < 0L) || any(gene_to_fam > n_families)) {
    stop("Invalid input: gene_to_fam indices must be between 0 and n_families (0 = no family assignment)")  
  }
  
  # Call the Rcpp forwarder.
  result <- tox_compute_family_scaling_rcpp( distances, gene_to_fam, n_families)

  # Check for backend error
  if (result$ierr != 0) {
    check_err_code(result$ierr)
  }

  # Return structured result
  return(list(
    dscale = as.numeric(result$dscale),
    loess_x = as.numeric(result$loess_x),
    loess_y = as.numeric(result$loess_y),
    indices_used = as.integer(result$indices_used)
  ))}

#' Expert: Compute per-family scaling with caller-provided work arrays
#'
#' Expert version of `tox_compute_family_scaling`. This variant accepts
#' caller-allocated workspace vectors for sorting/permutation and family
#' distance storage. Use this when you need fine-grained control over
#' memory allocation or to avoid repeated allocations in tight loops.
#'
#' @param distances Numeric vector of distances for each gene to its family centroid.
#' @param gene_to_fam Integer vector mapping each gene to a family index (0 = no family, >0 = family index).
#' @param n_families Integer number of families (>= 1).
#' @param perm_tmp Integer permutation vector for sorting (length == length(distances)).
#' @param stack_left_tmp Integer stack vector for sorting (length == length(distances)).
#' @param stack_right_tmp Integer stack vector for sorting (length == length(distances)).
#' @param family_distances Numeric work vector for sorting (length == length(distances)).
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{dscale}{Numeric vector of scaling factors per family}
#'     \item{loess_x}{Numeric vector of LOESS x reference points (median distances)}
#'     \item{loess_y}{Numeric vector of LOESS y reference points (stddev)}
#'     \item{indices_used}{Integer vector of family indices included in the calculation}
#'   }
tox_compute_family_scaling_expert <- function(distances, gene_to_fam, n_families,
                                              perm_tmp, stack_left_tmp, stack_right_tmp,
                                              family_distances) {

  # Input validation
  if (!is.numeric(distances)) {
    stop("distances must be numeric")
  }
  if (!is.numeric(gene_to_fam) && !is.integer(gene_to_fam)) {
    stop("gene_to_fam must be numeric or integer")
  }
  if (!is.numeric(n_families) && !is.integer(n_families)) {
    stop("n_families must be numeric or integer")
  }
  if (!is.numeric(family_distances) && !is.integer(family_distances)) {
    stop("family_distances must be numeric")
  }
  if (!is.numeric(perm_tmp) && !is.integer(perm_tmp)) {
    stop("perm_tmp must be numeric or integer")
  }
  if (!is.numeric(stack_left_tmp) && !is.integer(stack_left_tmp)) {
    stop("stack_left_tmp must be numeric or integer")
  }
  if (!is.numeric(stack_right_tmp) && !is.integer(stack_right_tmp)) {
    stop("stack_right_tmp must be numeric or integer")
  }

  # Convert to appropriate types
  distances        <- as.numeric(distances)
  gene_to_fam      <- as.integer(gene_to_fam)
  n_families       <- as.integer(n_families)
  perm_tmp         <- as.integer(perm_tmp)
  stack_left_tmp   <- as.integer(stack_left_tmp)
  stack_right_tmp  <- as.integer(stack_right_tmp)
  family_distances <- as.numeric(family_distances)

  # Dimention
  n_genes <- as.integer(length(distances))

  # Validate Dimentions
  if (n_genes == 0L) {
    stop("distances cannot be empty") 
  }
  if (any(!is.finite(distances))) {
    stop("distances contains non-finite values")  
  }
  if (length(gene_to_fam) != n_genes) {
    stop("length of gene_to_fam must match number of distances")  
  }
  if (n_families < 1L) {
    stop("n_families must be >= 1")  
  }
  # Work arrays must match distances length
  if (length(perm_tmp) != n_genes ||
      length(stack_left_tmp) != n_genes ||
      length(stack_right_tmp) != n_genes ||
      length(family_distances) != n_genes) {
    stop("workspace arrays (perm_tmp, stack_left_tmp, stack_right_tmp, family_distances) must each have length equal to number of genes")  
  }
  if (any(gene_to_fam < 0L) || any(gene_to_fam > n_families)) {
    stop("Invalid input: gene_to_fam indices must be between 0 and n_families (0 = no family assignment)")  
  }

 
  # Call Rcpp expert forwarder
  result <- tox_compute_family_scaling_expert_rcpp( n_families, distances, gene_to_fam, perm_tmp, stack_left_tmp, stack_right_tmp,
    family_distances
  )

   # Check for backend error
  if (result$ierr != 0) {
    check_err_code(result$ierr)
  }

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

#Compute Relative Distance Index (RDI)
#'
#' Computes the Relative Distance Index (RDI) for each gene by dividing its
#' Euclidean distance to the family centroid by the scaling factor of its family.
#' The function validates inputs in R, then delegates computation to the
#' Rcpp/Fortran backend. This function automatically checks for errors and
#' throws informative exceptions.
#'
#' @param distances Numeric vector of distances for each gene to its centroid
#'   (length = n_genes).
#' @param gene_to_fam Integer vector mapping each gene to a family
#'   (0 = no family, >0 = family index); must have the same length as `distances`.
#' @param dscale Numeric vector of scaling factors per family
#'   (length = n_families; indices in `gene_to_fam` must be in [0, n_families]).
#'
#' @return A list containing:
#' \describe{
#'   \item{rdi}{Numeric vector of RDI values for each gene.}
#'   \item{sorted_rdi}{Numeric vector of RDI values sorted in ascending order.}
#' }
tox_compute_rdi <- function(distances, gene_to_fam, dscale) {
  # Input validation
  if (!is.numeric(distances)) {
    stop("distances must be numeric")
  }
  if (!is.numeric(gene_to_fam) && !is.integer(gene_to_fam)) {
    stop("gene_to_fam must be numeric or integer")
  }
  if (!is.numeric(dscale)) {
    stop("dscale must be numeric")
  }

  # Type conversion
  distances <- as.numeric(distances)
  gene_to_fam <- as.integer(gene_to_fam)
  dscale <- as.numeric(dscale)

  # Calculate Dimentions 
  n_genes <- as.integer(length(distances))
  n_families <- as.integer(length(dscale))
  
  # Validate Dimentions
  if (n_genes == 0L) {
    stop("distances cannot be empty")  
  }
  if (length(gene_to_fam) != n_genes) {
    stop("length of gene_to_fam must match number of distances") 
  }
  if (any(!is.finite(distances)) || any(!is.finite(dscale))) {
    stop("distances and dscale must contain finite numeric values")  
  }
  if (n_families < 1L) {
    stop("dscale must have length >= 1 (n_families >= 1)") 
  }
  # gene_to_fam indices must be within [0, n_families]
  if (any(gene_to_fam < 0L) || any(gene_to_fam > n_families)) {
    stop("Invalid input: gene_to_fam indices must be between 0 and n_families (0 = no family assignment)")
  }
 
  # Call Rcpp forwarder (no ierr check requested)
  result <- tox_compute_rdi_rcpp(distances, gene_to_fam, dscale)


  # Return RDI values (forwarder is expected to fill rdi and sorted_rdi)
  return(list(
    rdi = result$rdi,
    sorted_rdi = result$sorted_rdi
  ))
}
#' Identify outliers from RDI values
#'
#' Identifies genes whose Relative Distance Index (RDI) falls in the top
#' percentile. The function delegates computation to the Rcpp forwarder and
#' returns a logical vector and the numeric threshold used.
#'
#' @param rdi Numeric vector of RDI values per gene
#' @param percentile Numeric scalar percentile in [0,100] (default 95.0)
#'
#' @return A named list with:
#'   \describe{
#'     \item{is_outlier}{logical vector indicating which genes are outliers}
#'     \item{threshold}{numeric scalar threshold value used}
#'   }
tox_identify_outliers <- function(rdi, percentile = 95.0) {

  # Input validation
  if (!is.numeric(rdi)) {
    stop("rdi must be numeric")
  }
  if (!is.numeric(percentile)) {
    stop("percentile must be numeric")
  }

  # Type conversion
  rdi <- as.numeric(rdi)
  percentile <- as.numeric(percentile)

  # Calculate dimentions
  n_genes <- as.integer(length(rdi))

   # Validate dimentions
  if (n_genes == 0L) {
    stop("rdi cannot be empty") 
  }
  if (any(!is.finite(rdi))) {
    stop("rdi contains non-finite values") 
  }
  if (percentile < 0.0 || percentile > 100.0) {
    stop("percentile must be between 0 and 100") 
  }
  
  # Call forwarder
  result <- tox_identify_outliers_rcpp(rdi,percentile)

  return(list(
    is_outlier = result$is_outlier,
    threshold = result$threshold
  ))

}
