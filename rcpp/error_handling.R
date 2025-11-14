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
    # FORMAT ERRORS
    '200' = "Invalid format detected.",
    '201' = "Invalid input provided.",
    '202' = "Empty input arrays provided.",
    '203' = "Dimension mismatch detected.",
    '204' = "NaN or Inf found in input data.",
    '205' = "Unsupported data type encountered.",
    '206' = "Array size mismatch detected.",
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


## Centralized, reusable validators (message-based)
## These functions only validate and produce descriptive, test-friendly
## error messages. They should not perform heavy coercion; callers may
## still coerce types when appropriate.

validate_numeric_vector <- function(x, name = deparse(substitute(x))) {
  if (!is.numeric(x)) stop(sprintf("%s must be numeric", name))
  invisible(TRUE)
}

validate_nonempty <- function(x, name = deparse(substitute(x))) {
  if (length(x) == 0) stop(sprintf("%s cannot be empty", name))
  invisible(TRUE)
}

validate_same_length <- function(x, y, name_x = deparse(substitute(x)), name_y = deparse(substitute(y))) {
  if (length(x) != length(y)) stop(sprintf("%s and %s must have the same length", name_x, name_y))
  invisible(TRUE)
}

validate_divisible_length <- function(x, d, name = deparse(substitute(x))) {
  if (d <= 0) stop(sprintf("d must be a positive integer for %s", name))
  if (length(x) %% d != 0) stop(sprintf("Length of %s must be divisible by d", name))
  invisible(TRUE)
}

validate_numeric_matrix <- function(m, name = deparse(substitute(m))) {
  if (!is.matrix(m)) stop(sprintf("%s must be a matrix", name))
  if (!is.numeric(m)) stop(sprintf("%s must be numeric", name))
  invisible(TRUE)
}

validate_logical_or_index_vector <- function(v, expected_length = NULL, name = deparse(substitute(v))) {
  if (!(is.logical(v) || is.numeric(v) || is.integer(v))) stop(sprintf("%s must be logical or numeric", name))
  if (!is.null(expected_length) && length(v) != expected_length) stop(sprintf("%s length must match expected length %d", name, expected_length))
  invisible(TRUE)
}

validate_gene_to_family <- function(gene_to_fam, n_genes, n_families, name = deparse(substitute(gene_to_fam))) {
  if (!is.numeric(gene_to_fam) && !is.integer(gene_to_fam)) stop(sprintf("%s must be numeric or integer", name))
  if (length(gene_to_fam) != n_genes) stop(sprintf("Length of %s must equal number of genes (%d)", name, n_genes))
  if (any(is.na(gene_to_fam))) stop(sprintf("%s contains NA values", name))
  if (any(gene_to_fam < 0)) stop(sprintf("%s indices must be between 0 and %d (0 = no family assignment)", name, n_families))
  invisible(TRUE)
}

validate_index_bounds <- function(idx, low = 1, high = Inf, name = deparse(substitute(idx))) {
  if (!is.numeric(idx) && !is.integer(idx)) stop(sprintf("%s must be numeric or integer", name))
  if (any(idx < low, na.rm = TRUE) || any(idx > high, na.rm = TRUE)) stop(sprintf("Invalid input: %s indices must be between %s and %s", name, low, ifelse(is.infinite(high), "Inf", as.character(high))))
  invisible(TRUE)
}

## Small helpers for messages used in tests
validate_length_equals_n <- function(x, n, name = deparse(substitute(x))) {
  if (length(x) != n) stop(sprintf("%s length must equal %d", name, n))
  invisible(TRUE)
}

validate_positive_integer_scalar <- function(x, name = deparse(substitute(x))) {
  if (!(is.numeric(x) || is.integer(x))) stop(sprintf("%s must be numeric or integer", name))
  if (length(x) != 1 || is.na(x) || as.integer(x) <= 0) stop(sprintf("%s must be a positive integer scalar", name))
  invisible(TRUE)
}

