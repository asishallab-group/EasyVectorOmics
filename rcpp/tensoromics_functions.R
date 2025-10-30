
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
# OUTLIER DETECTION FUNCTIONS (Validation + direct Rcpp return)
# ===================================================================

#' Compute family scaling factors for outlier detection
#'
#' Computes per-family scaling (and LOESS helper arrays) used by the RDI pipeline.
#'
#' @param distances Numeric vector of per-gene distances (length = number of genes).
#' @param gene_to_fam Integer vector mapping each gene to a 1-based family ID (length = number of genes).
#' @return A list of computed scaling and LOESS arrays.

tox_compute_family_scaling <- function(distances, gene_to_fam) {
 
  # Input validation
  if (!is.numeric(distances)) {
    check_err_code(201)
  }
  if (!is.numeric(gene_to_fam) && !is.integer(gene_to_fam)) {
    check_err_code(201)
  }

  # Type conversion
  distances <- as.numeric(distances)
  gene_to_fam <- as.integer(gene_to_fam)
  
  n_genes <- as.integer(length(distances))
  # infer number of families from the maximum family index (0 is allowed meaning "no family")
  n_families <- as.integer(max(gene_to_fam, na.rm = TRUE))

  if (n_genes == 0) {
    check_err_code(202)
  }
  if (length(gene_to_fam) != n_genes) {
    check_err_code(206) 
  }
  if (n_families < 1) {
    check_err_code(201)
  }
  # allow 0 = no family; reject negatives or indices > n_families
  if (any(gene_to_fam < 0) || any(gene_to_fam > n_families)) {
    check_err_code(201)
  }

  # NaN / Inf checks per manual
  if (any(is.na(distances)) || any(!is.finite(distances))) {
    check_err_code(204)
  }

  # Allocate family-level outputs
  dscale <- numeric(n_families)
  loess_x <- numeric(n_families)
  loess_y <- numeric(n_families)
  indices_used <- integer(n_families)

  # Call the Rcpp implementation
    result <- tox_compute_family_scaling_rcpp(
      n_genes, n_families,
      distances, gene_to_fam,
      dscale, loess_x, loess_y, indices_used
    )

  # Check Fortran error code
  if (result$ierr != 0) {
    check_err_code(result$ierr)
  }
  
  return(result)
}

#' Compute family scaling (expert)
#'
#' Expert entrypoint that accepts pre-allocated per-gene work arrays. 
#' workspace buffers to avoid repeated allocations.
#'
#' @param distances Numeric vector of per-gene distances (length = number of genes).
#' Compute family scaling factors (expert variant)
#'
#' Expert entrypoint that accepts pre-allocated per-gene workspace buffers.
#' This mirrors the internal Fortran/C expert entrypoint and is useful when
#' callers want to reuse temporary buffers to avoid repeated allocations.
#'
#' @param distances Numeric vector of per-gene distances (length = number of genes).
#' @param gene_to_fam Integer vector mapping each gene to a 1-based family ID (0 = no family).
#' @param perm_tmp Integer vector; permutation workspace buffer (length = number of genes).
#' @param stack_left_tmp Integer vector; left-stack workspace buffer (length = number of genes).
#' @param stack_right_tmp Integer vector; right-stack workspace buffer (length = number of genes).
#' @param family_distances Numeric vector; per-gene working buffer for family-distance computations (length = number of genes).
#' @param n_families Optional integer; if NULL it will be inferred from `gene_to_fam`.
#'
#' @return A named list with elements: dscale (numeric, length = n_families),
#' loess_x (numeric), loess_y (numeric), indices_used (integer indices of used families),
#' and ierr (integer canonical error code). The function will call `check_err_code()`
#' on non-zero `ierr` before returning in higher-level wrappers.
tox_compute_family_scaling_expert <- function(distances, gene_to_fam, perm_tmp, stack_left_tmp, stack_right_tmp, family_distances, n_families = NULL) {

  # Input validation
  if (!is.numeric(distances)) {
    check_err_code(201)
  }
  if (!is.numeric(gene_to_fam) && !is.integer(gene_to_fam)) {
    check_err_code(201)
  }
  if (!is.integer(perm_tmp) && !is.numeric(perm_tmp)) {
    check_err_code(201)
  }
  if (!is.integer(stack_left_tmp) && !is.numeric(stack_left_tmp)) {
    check_err_code(201)
  }
  if (!is.integer(stack_right_tmp) && !is.numeric(stack_right_tmp)) {
    check_err_code(201)
  }
  if (!is.numeric(family_distances)) {
    check_err_code(201)
  }

# Type conversion
  distances <- as.numeric(distances)
  gene_to_fam <- as.integer(gene_to_fam)
  perm_tmp <- as.integer(perm_tmp)
  stack_left_tmp <- as.integer(stack_left_tmp)
  stack_right_tmp <- as.integer(stack_right_tmp)
  family_distances <- as.numeric(family_distances)

  n_genes <- as.integer(length(distances))
  if (n_genes == 0) check_err_code(202)

  if (is.null(n_families)) {
    n_families <- as.integer(max(gene_to_fam, na.rm = TRUE))
    if (is.na(n_families) || n_families < 1L) check_err_code(201)
  } else {
    n_families <- as.integer(n_families)
  }

# Dimension checks for workspace buffers
  if (length(gene_to_fam) != n_genes) {
    check_err_code(206)
  }
  if (length(perm_tmp) != n_genes) {
    check_err_code(203)
  }
  if (length(stack_left_tmp) != n_genes) {
    check_err_code(203)
  }
  if (length(stack_right_tmp) != n_genes) {
    check_err_code(203)
  }
  if (length(family_distances) != n_genes) {
    check_err_code(203)
  }

  # allow 0 = no family; reject negatives or indices > n_families
  if (any(gene_to_fam < 0) || any(gene_to_fam > n_families)) {
    check_err_code(201)
  }

  # NaN / Inf checks
  if (any(is.na(distances)) || any(!is.finite(distances))) check_err_code(204)
  if (any(is.na(family_distances)) || any(!is.finite(family_distances))) check_err_code(204)

  # Allocate family-level outputs
  dscale <- numeric(n_families)
  loess_x <- numeric(n_families)
  loess_y <- numeric(n_families)
  indices_used <- integer(n_families)

  # Call the Rcpp expert forwarder
  result <- tox_compute_family_scaling_expert_rcpp(
    n_genes, n_families,
    distances, gene_to_fam,
    dscale, loess_x, loess_y, indices_used,
    perm_tmp, stack_left_tmp, stack_right_tmp,
    family_distances
  )

  if (!is.null(result$ierr) && result$ierr != 0) {
    check_err_code(result$ierr)
  }

  if (!is.null(result$indices_used)) {
    idx <- as.integer(result$indices_used)
    idx <- idx[idx != 0L]
    result$indices_used <- idx
  }

  return(result)
}

#' Calculate Euclidean distance between two vectors
#'
#' Computes the Euclidean (L2) distance between two numeric vectors of the same length.
#'
#' @param vec1 Numeric vector (first vector)
#' @param vec2 Numeric vector (second vector, same length as vec1)
#'
#' @return Numeric scalar: the Euclidean distance between the two vectors
tox_euclidean_distance <- function(vec1, vec2) {

  # Input validation 
  if (!is.numeric(vec1) || !is.numeric(vec2)) {
    check_err_code(201)
  }
  if (length(vec1) != length(vec2)) {
    check_err_code(203)
  }
  if (length(vec1) == 0) {
    check_err_code(202)
  }

  # NaN / Inf checks
  if (any(is.na(vec1)) || any(!is.finite(vec1)) || any(is.na(vec2)) || any(!is.finite(vec2))) {
    check_err_code(204)
  }


  # Call Rcpp implementation 
   result <- tox_euclidean_distance_rcpp(as.numeric(vec1), as.numeric(vec2))
  return(result$distance)
}

#' Calculate tissue versatility (Rcpp wrapper)
#' Computes normalized tissue (axis) versatility and angles for selected expression vectors and axes.
#'
#' @param expression_vectors Numeric matrix of shape (n_axes × n_vectors) where each column is a gene expression vector
#' @param vector_selection Logical or numeric vector (length = n_vectors) indicating which vectors to process
#' @param axis_selection Logical or numeric vector (length = n_axes) indicating which axes (rows) to include
#'
#' @return A list with elements:
#'   - tissue_versatilities: numeric vector (length = n_selected_vectors)
#'   - tissue_angles_deg: numeric vector (length = n_selected_vectors)
#'   - n_selected_vectors: integer
#'   - n_selected_axes: integer
tox_calculate_tissue_versatility <- function(expression_vectors, vector_selection, axis_selection) {
  # Input validation
  if (!is.matrix(expression_vectors)) {
    check_err_code(201)
  }
  if (!is.logical(vector_selection) && !is.numeric(vector_selection)) {
    check_err_code(201)
  }
  if (!is.logical(axis_selection) && !is.numeric(axis_selection)) {
    check_err_code(201)
  }

  if (is.numeric(vector_selection)) {
    vector_selection <- as.logical(vector_selection)
  }
  if (is.numeric(axis_selection)) {
    axis_selection <- as.logical(axis_selection)
  }

  n_axes <- nrow(expression_vectors)
  n_vectors <- ncol(expression_vectors)

  if (n_axes == 0 || n_vectors == 0) {
    check_err_code(202)
  }

  if (length(vector_selection) != n_vectors) {
    check_err_code(206)
  }
  if (length(axis_selection) != n_axes) {
    check_err_code(206)
  }

  # Prepare integer selection vectors for the Rcpp call
  vec_sel_int <- as.integer(vector_selection)
  axis_sel_int <- as.integer(axis_selection)

  # NaN / Inf checks
  if (any(is.na(expression_vectors)) || any(!is.finite(expression_vectors))) {
    check_err_code(204)
  }

  # Call Rcpp implementation
  result <- tox_calculate_tissue_versatility_rcpp(as.matrix(expression_vectors), vec_sel_int, axis_sel_int)

  # Check and rethrow canonical errors
  if (!is.null(result$ierr) && result$ierr != 0) {
    check_err_code(result$ierr)
  }

  return(list(
    tissue_versatilities = result$tissue_versatilities,
    tissue_angles_deg = result$tissue_angles_deg,
    n_selected_vectors = result$n_selected_vectors,
    n_selected_axes = result$n_selected_axes
  ))
}

#' Compute Relative Distance Index (RDI) for each gene
#'
#' @param distances Numeric vector of per-gene distances (length = number of genes).
#' @param gene_to_fam Integer vector mapping each gene to a 1-based family ID (length = number of genes).
#' @param dscale Numeric vector of per-family scaling factors (length = number of families).
#'
#' @return A list with RDI outputs (e.g., rdi, sorted_rdi, permutation), as returned by the Rcpp implementation.
tox_compute_rdi <- function(distances, gene_to_fam, dscale) {
  # Input validation
  if (!is.numeric(distances)) {
    check_err_code(201)
  }
  if (!is.numeric(gene_to_fam) && !is.integer(gene_to_fam)) {
    check_err_code(201)
  }
  if (!is.numeric(dscale)) {
    check_err_code(201)
  }
# Type conversion
  distances <- as.numeric(distances)
  gene_to_fam <- as.integer(gene_to_fam)
  dscale <- as.numeric(dscale)

  n_genes <- as.integer(length(distances))
  n_families <- as.integer(length(dscale))

  if (n_genes == 0) {
    check_err_code(202)
  }
  if (length(gene_to_fam) != n_genes) {
    check_err_code(206)
  }
  if (n_families < 1) {
    check_err_code(201)
  }
  # allow 0 = no family; reject negatives or indices > n_families
  if (any(gene_to_fam < 0) || any(gene_to_fam > n_families)) {
    check_err_code(201)
  }
  # NaN / Inf checks
  if (any(is.na(distances)) || any(!is.finite(distances)) || any(is.na(dscale)) || any(!is.finite(dscale))) {
    check_err_code(204)
  }
# Allocate outputs
  rdi <- numeric(n_genes)
  sorted_rdi <- numeric(n_genes)
  perm <- integer(n_genes)
  stack_left <- integer(n_genes)
  stack_right <- integer(n_genes)

  # Call Rcpp wrapper (no ierr check - Fortran routine doesn't return it yet)
  result <- tox_compute_rdi_rcpp(
    n_genes, n_families,
    distances, gene_to_fam, dscale,
    rdi, sorted_rdi,
    perm, stack_left, stack_right
  )

  # Check Fortran error code
  if (!is.null(result$ierr) && result$ierr != 0) {
    check_err_code(result$ierr)
  }

  return(result)
}

#' Identify gene outliers based on RDI percentile threshold
#'
#' @param rdi Numeric vector of per-gene RDI values.
#' @param sorted_rdi Numeric vector; sorted copy of rdi used for thresholding (same length as rdi).
#' @param percentile Numeric scalar in [0, 100] specifying the percentile cutoff.
#'
#' @return Integer vector of 0/1 flags (or equivalent structure) indicating outliers, as returned by the Rcpp implementation.
tox_identify_outliers <- function(rdi, threshold=NULL, percentile=95.0) {
  # Input validation
  if (!is.numeric(rdi)) {
     check_err_code(201)
  }
  if (!is.numeric(sorted_rdi)) {
    check_err_code(201)
  }
  if (!is.numeric(percentile)) {
    check_err_code(201)
  }
# Type conversion
  rdi <- as.numeric(rdi)
  sorted_rdi <- as.numeric(sorted_rdi)
  percentile <- as.numeric(percentile)

  # NaN / Inf checks
  if (any(is.na(rdi)) || any(!is.finite(rdi)) || any(is.na(sorted_rdi)) || any(!is.finite(sorted_rdi))) {
    check_err_code(204)
  }

  n_genes <- as.integer(length(rdi))
  if (n_genes == 0) {
     check_err_code(202)
  }
  if (length(sorted_rdi) != n_genes) {
     check_err_code(203)
  }
  if (length(percentile) != 1L || is.na(percentile)) {
    check_err_code(201)
  }
  if (percentile < 0 || percentile > 100) {
    check_err_code(201)
  }

  is_outlier_int <- integer(n_genes)
# Call Rcpp wrapper
  result <- tox_identify_outliers_rcpp(
    n_genes,
    rdi, sorted_rdi,
    is_outlier_int,
    percentile
  )
  
  # Check Fortran error code
  if (result$ierr != 0) {
    check_err_code(result$ierr)
  }
  
  return(result)
}

#' Detect gene outliers using complete RDI pipeline
#'
#' Runs scaling -> RDI -> thresholding in one call.
#'
#' @param distances Numeric vector of per-gene distances (length = number of genes).
#' @param gene_to_fam Integer vector mapping each gene to a 1-based family ID (length = number of genes).
#' @param n_families Integer count of distinct families (>= 1).
#' @param percentile Numeric scalar in [0, 100] for the outlier cutoff.
#'
#' @return A list with pipeline outputs (e.g., RDI values, sorted RDI, and outlier flags), as returned by the Rcpp implementation.
tox_detect_outliers <- function(distances, gene_to_fam, n_families, percentile = 95.0) {
  # Input validation
  if (!is.numeric(distances)) {
     check_err_code(201)
  }
  if (!is.numeric(gene_to_fam) && !is.integer(gene_to_fam)) {
    check_err_code(201)
  }
  if (!is.numeric(n_families) && !is.integer(n_families)) {
    check_err_code(201)
  }
  if (!is.numeric(percentile)) {
    check_err_code(201)
  }
# Type conversion
  distances <- as.numeric(distances)
  gene_to_fam <- as.integer(gene_to_fam)
  n_families <- as.integer(n_families)
  percentile <- as.numeric(percentile)

  n_genes <- as.integer(length(distances))
  if (n_genes == 0) {
    check_err_code(202)
  }
  if (length(gene_to_fam) != n_genes) {
    check_err_code(206)
  }
  if (n_families < 1) {
    check_err_code(201)
  }
  if (any(gene_to_fam < 1) || any(gene_to_fam > n_families)) {
    check_err_code(201)
  }
  if (length(percentile) != 1L || is.na(percentile)) {
    check_err_code(201)
  }
  if (percentile < 0 || percentile > 100) {
    check_err_code(201)
  }
# Allocate outputs and work arrays
  work_array <- numeric(n_genes)
  perm <- integer(n_genes)
  stack_left <- integer(n_genes)
  stack_right <- integer(n_genes)
  is_outlier_int <- integer(n_genes)
  loess_x <- numeric(n_families)
  loess_y <- numeric(n_families)
  loess_n <- integer(n_families)
# Call Rcpp wrapper
  result <- tox_detect_outliers_rcpp(
    n_genes, n_families,
    distances, gene_to_fam,
    work_array,
    perm, stack_left, stack_right,
    is_outlier_int,
    loess_x, loess_y, loess_n,
    percentile
  )
  
  # Check Fortran error code
  if (result$ierr != 0) {
    check_err_code(result$ierr)
  }
  
  return(result)
}
