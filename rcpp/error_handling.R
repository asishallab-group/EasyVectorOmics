check_err_code <- function(ierr) {
  if (ierr == 0) return(invisible(NULL))
  msg <- switch(as.character(ierr),
    # I/O errors
    '101' = "Could not open file.",
    '102' = "Could not read magic number.",
    '103' = "Could not read type code.",
    '104' = "Could not read number of dimensions.",
    '105' = "Could not read array dimensions",
    '106' = "Could not read character length.",
    '107' = "Could not read array data.",
    '112' = "Could not write magic number",
    '113' = "Could not write type code",
    '114' = "Could not write number of dimensions",
    '115' = "Could not write dimensions",
    '116' = "Could not write character length",
    '117' = "Could not write array data",
    # ADD MORE HERE
    
    # FORMAT ERRORS
    '200' = "Invalid format detected.",
    '201' = "Invalid input provided.",
    '202' = "Empty input arrays provided.",
    '203' = "Dimension mismatch detected.",
    '204' = "NaN or Inf found in input data.",
    '205' = "Unsupported data type encountered.",
    '206' = "Array size mismatch detected",

    # MEMORY ERRORS
    '301' = "Memory allocation failed.",
    '302' = "Null pointer reference encountered.",

    # FORTRAN RUNTIME ERRORS
    '5002' = "Fortran runtime error: unit not open / not connected.",

    # Internal errors
    '9001' = "Internal error: unexpected state.",
    '9999' = "Unknown error.",
    paste("Unmapped error code:", ierr)
  )
  stop(msg)
}

## Centralized R-side validators (kept inside rcpp/ when r/ must not be changed)
## Messages include specific substrings the tests may assert on.

validate_numeric_matrix <- function(x, name = "matrix") {
  if (!is.matrix(x) || !is.numeric(x)) {
    stop(sprintf("`%s` must be a numeric matrix (must be numeric).", name))
  }
  invisible(NULL)
}

validate_integer_vector <- function(x, name = "vector", expected_length = NULL, bounds = NULL) {
  if (!is.integer(x)) stop(sprintf("`%s` must be an integer vector.", name))
  if (!is.null(expected_length) && length(x) != expected_length) stop(sprintf("`%s` must be of length %d.", name, expected_length))
  if (!is.null(bounds)) {
    if (!is.null(bounds$min) && any(x < bounds$min)) stop(sprintf("`%s` contains values < %d.", name, bounds$min))
    if (!is.null(bounds$max) && any(x > bounds$max)) stop(sprintf("`%s` contains values > %d.", name, bounds$max))
  }
  invisible(NULL)
}

validate_logical_vector <- function(x, name = "vector", expected_length = NULL) {
  if (!is.logical(x)) stop(sprintf("`%s` must be a logical vector.", name))
  if (!is.null(expected_length) && length(x) != expected_length) stop(sprintf("`%s` must be of length %d.", name, expected_length))
  invisible(NULL)
}

validate_mode <- function(mode, allowed = c('all', 'ortho')) {
  if (!mode %in% allowed) stop(sprintf("`mode` must be one of: %s.", paste(allowed, collapse = ", ")))
  invisible(NULL)
}

# Higher-level validators for specific wrappers
validate_group_centroid_inputs <- function(expression_vectors, gene_to_family, n_families, ortholog_set, mode = 'all') {
  validate_numeric_matrix(expression_vectors, "expression_vectors")
  n_axes <- nrow(expression_vectors)
  n_genes <- ncol(expression_vectors)
  validate_integer_vector(as.integer(gene_to_family), "gene_to_family", expected_length = n_genes)
  validate_logical_vector(as.logical(ortholog_set), "ortholog_set", expected_length = n_genes)
  if (!is.numeric(n_families) || length(n_families) != 1 || as.integer(n_families) < 1) stop("`n_families` must be a positive integer scalar.")
  validate_mode(mode)
  invisible(NULL)
}

validate_mean_vector_inputs <- function(expression_vectors, gene_indices) {
  validate_numeric_matrix(expression_vectors, "expression_vectors")
  n_genes <- ncol(expression_vectors)
  if (!is.integer(gene_indices)) stop("`gene_indices` must be an integer vector of 1-based column indices.")
  if (any(gene_indices < 1) || any(gene_indices > n_genes)) stop("`gene_indices` must contain indices between 1 and ncol(expression_vectors).")
  invisible(NULL)
}

# Additional generic validators
validate_numeric_vector <- function(x, name = "x") {
  if (!is.numeric(x)) stop(sprintf("`%s` must be numeric (must be numeric).", name))
  invisible(NULL)
}

validate_positive_integer_scalar <- function(x, name = "value") {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || as.integer(x) < 1) stop(sprintf("`%s` must be a positive integer scalar.", name))
  invisible(NULL)
}

validate_equal_length <- function(a, b, name_a = "a", name_b = "b") {
  if (length(a) != length(b)) stop(sprintf("`%s` length must match length of `%s` (same length).", name_a, name_b))
  invisible(NULL)
}

# Ensure vector is non-empty
validate_nonempty_vector <- function(x, name = "x") {
  if (length(x) == 0) stop(sprintf("`%s` must not be empty (cannot be empty).", name))
  invisible(NULL)
}

# Ensure length(x) is divisible by d (used when vectors represent flattened matrices)
validate_divisible_length <- function(x, d, name = "x") {
  if (!is.numeric(d) || length(d) != 1 || is.na(d) || as.integer(d) < 1) stop("`d` must be a positive integer scalar.")
  if (length(x) %% as.integer(d) != 0) stop(sprintf("Length of `%s` must be divisible by d (divisible by %d).", name, as.integer(d)))
  invisible(NULL)
}

# Validate integer index vector (1-based or 0-based depending on contract)
validate_index_vector <- function(idx, min_val, max_val, name = "indices") {
  if (!is.integer(idx)) stop(sprintf("`%s` must be an integer vector.", name))
  if (any(idx < min_val) || any(idx > max_val)) stop(sprintf("`%s` must contain indices between %d and %d (must be between %d and %d).", name, min_val, max_val, min_val, max_val))
  invisible(NULL)
}

# Validate that two matrices have the same number of rows
validate_matching_rows <- function(A, B, name_A = "A", name_B = "B") {
  if (!is.matrix(A) || !is.matrix(B)) stop(sprintf("`%s` and `%s` must both be matrices.", name_A, name_B))
  if (nrow(A) != nrow(B)) stop(sprintf("Number of rows in `%s` must match number of rows in `%s`.", name_A, name_B))
  invisible(NULL)
}

# Validate length equals given n (used to match gene_to_fam length to number of genes)
validate_length_equals_n <- function(x, n, name = "x") {
  if (length(x) != as.integer(n)) stop(sprintf("Length of %s must equal number of genes.", name))
  invisible(NULL)
}
