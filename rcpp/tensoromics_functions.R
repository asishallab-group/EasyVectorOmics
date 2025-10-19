
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
#' @param n_families Integer count of distinct families (>= 1).
#' @param expert Logical; if TRUE, uses the expert variant with extra work arrays.
#'
#' @return A list of computed scaling and LOESS arrays (structure as returned by the Rcpp implementation).
tox_compute_family_scaling <- function(distances, gene_to_fam, n_families, expert = FALSE) {
  # Input validation
  if (!is.numeric(distances)) {
    check_err_code(207)
  }
  if (!is.numeric(gene_to_fam) && !is.integer(gene_to_fam)) {
    check_err_code(207)
  }
  if (!is.numeric(n_families) && !is.integer(n_families)) {
    check_err_code(207)
  }
  if (!is.logical(expert)) {
     check_err_code(201) 
  }

  distances <- as.numeric(distances)
  gene_to_fam <- as.integer(gene_to_fam)
  n_families <- as.integer(n_families)
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
    check_err_code(209)
  }
  
  # Allocate outputs
  dscale <- numeric(n_families)
  loess_x <- numeric(n_families)
  loess_y <- numeric(n_families)
  indices_used <- integer(n_families)

  if (expert) {
    perm_tmp <- integer(n_genes)
    stack_left_tmp <- integer(n_genes)
    stack_right_tmp <- integer(n_genes)
    family_distances <- numeric(n_genes)

    result <- tox_compute_family_scaling_expert_rcpp(
      n_genes, n_families,
      distances, gene_to_fam,
      dscale, loess_x, loess_y, indices_used,
      perm_tmp, stack_left_tmp, stack_right_tmp,
      family_distances
    )
  } else {
    result <- tox_compute_family_scaling_rcpp(
      n_genes, n_families,
      distances, gene_to_fam,
      dscale, loess_x, loess_y, indices_used
    )
  }
  
  # Check Fortran error code
  if (result$ierr != 0) {
    check_err_code(result$ierr)
  }
  
  return(result)
}

#' Compute Relative Distance Index (RDI) for each gene
#'
#' @param distances Numeric vector of per-gene distances (length = number of genes).
#' @param gene_to_fam Integer vector mapping each gene to a 1-based family ID (length = number of genes).
#' @param dscale Numeric vector of per-family scaling factors (length = number of families).
#'
#' @return A list with RDI outputs (e.g., rdi, sorted_rdi, permutation), as returned by the Rcpp implementation.
tox_compute_rdi <- function(distances, gene_to_fam, dscale) {
  if (!is.numeric(distances)) {
    check_err_code(207)
  }
  if (!is.numeric(gene_to_fam) && !is.integer(gene_to_fam)) {
    check_err_code(207) 
  }
  if (!is.numeric(dscale)) {
    check_err_code(207) 
  }

  distances <- as.numeric(distances)
  gene_to_fam <- as.integer(gene_to_fam)
  dscale <- as.numeric(dscale)

  n_genes <- as.integer(length(distances))
  n_families <- as.integer(length(dscale))

  if (n_genes == 0) {
     check_err_code(212) 
  }
  if (length(gene_to_fam) != n_genes) {
    check_err_code(206)
  }
  if (n_families < 1) {
    check_err_code(201)
  }
  if (any(gene_to_fam < 1) || any(gene_to_fam > n_families)) {
    check_err_code(209)
  }

  rdi <- numeric(n_genes)
  sorted_rdi <- numeric(n_genes)
  perm <- integer(n_genes)
  stack_left <- integer(n_genes)
  stack_right <- integer(n_genes)

  # Call Rcpp wrapper (no ierr check - Fortran routine doesn't return it yet)
  return(tox_compute_rdi_rcpp(
    n_genes, n_families,
    distances, gene_to_fam, dscale,
    rdi, sorted_rdi,
    perm, stack_left, stack_right
  ))
}

#' Identify gene outliers based on RDI percentile threshold
#'
#' @param rdi Numeric vector of per-gene RDI values.
#' @param sorted_rdi Numeric vector; sorted copy of rdi used for thresholding (same length as rdi).
#' @param percentile Numeric scalar in [0, 100] specifying the percentile cutoff.
#'
#' @return Integer vector of 0/1 flags (or equivalent structure) indicating outliers, as returned by the Rcpp implementation.
tox_identify_outliers <- function(rdi, sorted_rdi, percentile = 95.0) {
  if (!is.numeric(rdi)) {
     check_err_code(207) 
  }
  if (!is.numeric(sorted_rdi)) {
    check_err_code(207)
  }
  if (!is.numeric(percentile)) {
    check_err_code(207)
  }

  rdi <- as.numeric(rdi)
  sorted_rdi <- as.numeric(sorted_rdi)
  percentile <- as.numeric(percentile)

  n_genes <- as.integer(length(rdi))
  if (n_genes == 0) {
     check_err_code(212) 
  }
  if (length(sorted_rdi) != n_genes) {
     check_err_code(214) 
  }
  if (length(percentile) != 1L || is.na(percentile)) {
    check_err_code(213) 
  }
  if (percentile < 0 || percentile > 100) {
    check_err_code(210)
  }

  is_outlier_int <- integer(n_genes)

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
  if (!is.numeric(distances)) {
     check_err_code(207) 
  }
  if (!is.numeric(gene_to_fam) && !is.integer(gene_to_fam)) {
    check_err_code(207)
  }
  if (!is.numeric(n_families) && !is.integer(n_families)) {
    check_err_code(207)
  }
  if (!is.numeric(percentile)) {
    check_err_code(207)
  }

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
    check_err_code(209)
  }
  if (length(percentile) != 1L || is.na(percentile)) {
    check_err_code(213)
  }
  if (percentile < 0 || percentile > 100) {
    check_err_code(210)
  }

  work_array <- numeric(n_genes)
  perm <- integer(n_genes)
  stack_left <- integer(n_genes)
  stack_right <- integer(n_genes)
  is_outlier_int <- integer(n_genes)
  loess_x <- numeric(n_families)
  loess_y <- numeric(n_families)
  loess_n <- integer(n_families)

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
