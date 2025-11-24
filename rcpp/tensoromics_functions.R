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
  
  check_err_code(result$ierr)

  
  # Return structured result 
  return(list(
    tissue_versatilities = result$tissue_versatilities,
    tissue_angles_deg = result$tissue_angles_deg,
    n_selected_vectors = result$n_selected_vectors,
    n_selected_axes = result$n_selected_axes
  ))

}


# ============================================================
#  1) Array metadata wrapper
# ============================================================

tox_get_array_metadata <- function(filename, max_dims = 5, with_clen = FALSE) {
  validate_array_deserialize_inputs(filename, max_dims)

  filename_ascii <- utf8ToInt(filename)
  result <- get_array_metadata_rcpp(
    filename_ascii = as.integer(filename_ascii),
    fn_len         = length(filename_ascii),
    dims_out_capacity = max_dims,
    with_clen      = with_clen
  )

  check_err_code(result$ierr)

  dims <- result$dims[seq_len(result$ndim)]
  if (isTRUE(with_clen) && !is.null(result$clen)) {
    return(list(
      dims = dims,
      ndim = result$ndim,
      clen = result$clen
    ))
  }

  list(
    dims = dims,
    ndim = result$ndim
  )
}


# ============================================================
#  2) Deserialization (int / real / char)
# ============================================================

tox_deserialize_int_array <- function(filename, max_dims = 5) {
  validate_array_deserialize_inputs(filename, max_dims)

  meta <- tox_get_array_metadata(filename, max_dims)
  total_size <- prod(meta$dims)

  filename_ascii <- utf8ToInt(filename)

  result <- tox_deserialize_int_rcpp(
    filename_ascii = as.integer(filename_ascii),
    fn_len         = length(filename_ascii),
    total          = as.integer(total_size)
  )

  check_err_code(result$ierr)

  array(result$int_arr[seq_len(total_size)], dim = meta$dims)
}


tox_deserialize_real_array <- function(filename, max_dims = 5) {
  validate_array_deserialize_inputs(filename, max_dims)

  meta <- tox_get_array_metadata(filename, max_dims)
  total_size <- prod(meta$dims)

  filename_ascii <- utf8ToInt(filename)

  result <- tox_deserialize_real_rcpp(
    filename_ascii = as.integer(filename_ascii),
    fn_len         = length(filename_ascii),
    total          = as.integer(total_size)
  )

  check_err_code(result$ierr)

  array(result$real_arr[seq_len(total_size)], dim = meta$dims)
}


tox_deserialize_char_array <- function(filename, max_dims = 5) {
  validate_array_deserialize_inputs(filename, max_dims)

  meta <- tox_get_array_metadata(filename, max_dims, with_clen = TRUE)
  actual_dims <- meta$dims
  clen        <- meta$clen
  total_array_size <- prod(actual_dims)

  filename_ascii <- utf8ToInt(filename)

  result <- tox_deserialize_char_flat_rcpp(
    clen           = as.integer(clen),
    total          = as.integer(total_array_size),
    filename_ascii = as.integer(filename_ascii),
    fn_len         = length(filename_ascii)
  )

  check_err_code(result$ierr)

  mat <- matrix(result$ascii_arr, nrow = clen)
  chars <- apply(mat, 2L, function(col) {
    rawToChar(as.raw(col[col > 0L]))
  })

  array(chars, dim = actual_dims[seq_len(meta$ndim)])
}


# ============================================================
#  3) Serialization (int / real / char)
# ============================================================

tox_serialize_int_array <- function(arr, filename) {
  validate_filename(filename)

  flat <- as.integer(arr)
  dims <- if (is.null(dim(arr))) {
    as.integer(length(arr))
  } else {
    as.integer(dim(arr))
  }
  ndim <- length(dims)

  filename_ascii <- utf8ToInt(filename)

  ierr <- tox_serialize_int_nd_rcpp(
    arr            = flat,
    dims           = dims,
    ndim           = as.integer(ndim),
    filename_ascii = as.integer(filename_ascii),
    fn_len         = length(filename_ascii)
  )

  check_err_code(ierr)
  invisible(NULL)
}


tox_serialize_real_array <- function(arr, filename) {
  validate_filename(filename)

  flat <- as.double(arr)
  dims <- if (is.null(dim(arr))) {
    as.integer(length(arr))
  } else {
    as.integer(dim(arr))
  }
  ndim <- length(dims)

  filename_ascii <- utf8ToInt(filename)

  ierr <- tox_serialize_real_nd_rcpp(
    arr            = flat,
    dims           = dims,
    ndim           = as.integer(ndim),
    filename_ascii = as.integer(filename_ascii),
    fn_len         = length(filename_ascii)
  )

  check_err_code(ierr)
  invisible(NULL)
}


tox_serialize_char_array <- function(arr, filename) {
  validate_filename(filename)
  validate_character_vector(arr)

  arr  <- as.array(arr)
  dims <- dim(arr)
  if (is.null(dims)) dims <- length(arr)

  clen <- max(nchar(arr, type = "chars"))

  mat <- matrix(0L, nrow = clen, ncol = length(arr))
  for (i in seq_along(arr)) {
    chars <- utf8ToInt(substr(arr[i], 1L, clen))
    mat[seq_along(chars), i] <- chars
  }

  ascii_arr <- as.integer(mat)
  ndim      <- length(dims)
  filename_ascii <- utf8ToInt(filename)

  ierr <- tox_serialize_char_flat_rcpp(
    ascii_arr      = ascii_arr,
    dims           = as.integer(dims),
    ndim           = as.integer(ndim),
    clen           = as.integer(clen),
    filename_ascii = as.integer(filename_ascii),
    fn_len         = length(filename_ascii)
  )

  check_err_code(ierr)
  invisible(NULL)
}


# ============================================================
#  4) KD-tree index (multidimensional) + spherical KD
# ============================================================

build_kd_index <- function(X, dim_order = NULL) {
  validate_numeric_matrix(X, "X")

  d <- nrow(X)
  n <- ncol(X)

  if (is.null(dim_order)) {
    dim_order <- seq_len(d)
  }

  validate_integer_vector(as.integer(dim_order), "dim_order", expected_length = d)

  result <- tox_build_kd_index_rcpp(
    points          = X,
    dimension_order = as.integer(dim_order)
  )

  check_err_code(result$ierr)

  return(result$kd_indices)
}


build_spherical_kd <- function(vectors) {
  validate_numeric_matrix(vectors, "vectors")

  result <- build_spherical_kd_rcpp(vectors)

  check_err_code(result$ierr)

  return(result)
}


# ============================================================
#  5) BST index + range query
# ============================================================

build_bst_index <- function(x) {
  validate_numeric_vector(x, "x")
  validate_nonempty_vector(x, "x")

  result <- build_bst_index_rcpp(x)

  check_err_code(result$ierr)

  return(result$ix)
}


bst_range_query <- function(x, ix, lo, hi) {
  validate_numeric_vector(x, "x")
  validate_integer_vector(as.integer(ix), "ix")
  validate_equal_length(x, ix, "x", "ix")

  result <- bst_range_query_rcpp(
    values = x,
    ix     = as.integer(ix),
    lo     = lo,
    hi     = hi
  )

  check_err_code(result$ierr)

  indices <- result$output_indices[seq_len(result$num_matches)]

  list(indices = indices, count = result$num_matches)
}


# ============================================================
#  6) which() helper via C backend
# ============================================================

tox_which <- function(mask, m_max = length(mask)) {
  validate_integer_vector(as.integer(mask), "mask")

  result <- tox_which_rcpp(
    mask = as.integer(mask),
    m_max = as.integer(m_max)
  )

  check_err_code(result$ierr)

  result$idx_out[seq_len(result$m_out)]
}


# ============================================================
#  7) LOESS smoothing (2D)
# ============================================================

tox_loess_smooth_2d <- function(x_ref, y_ref, x_query,
                                indices_used = NULL,
                                kernel_sigma,
                                kernel_cutoff) {

  validate_loess_smooth_2d_inputs(
    x_ref         = x_ref,
    y_ref         = y_ref,
    x_query       = x_query,
    indices_used  = indices_used,
    kernel_sigma  = kernel_sigma,
    kernel_cutoff = kernel_cutoff
  )

  if (is.null(dim(y_ref))) {
    y_ref <- matrix(y_ref, nrow = 1L)
  }
  # We only need the first row as numeric vector for Rcpp
  y_ref_vec <- as.numeric(y_ref[1L, , drop = TRUE])

  if (is.null(indices_used)) {
    indices_used <- seq_along(x_ref)
  }

  result <- tox_loess_smooth_2d_rcpp(
    x_ref        = x_ref,
    y_ref        = y_ref_vec,
    x_query      = x_query,
    indices_used = as.integer(indices_used),
    kernel_sigma = kernel_sigma,
    kernel_cutoff= kernel_cutoff
  )

  check_err_code(result$ierr)

  list(
    y_out          = result$y_out,
    smoothed_values = as.vector(result$y_out)
  )
}

