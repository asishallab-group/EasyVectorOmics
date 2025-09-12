dyn.load("build/libtensor-omics.so")

#' Check error code and throw informative error if needed
#' 
#' @param ierr Error code from Fortran routine
tox_errors <- function(ierr) {
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

strings_to_ascii_matrix <- function(arr, clen) {
  n <- length(arr)
  # Create a matrix of integers with dimensions clen x n
  mat <- matrix(0L, nrow = clen, ncol = n)
  for (i in seq_along(arr)) {
    chars <- utf8ToInt(substr(arr[i], 1, clen))
    mat[seq_along(chars), i] <- chars
  }
  mat  # Return the matrix directly (not as.integer)
}

ascii_matrix_to_strings <- function(ascii_mat, clen) {
  # ascii_mat: integer vector or matrix, possibly 1D if only one gene
  if (is.null(dim(ascii_mat))) {
    # 1D vector: convert to matrix with one column
    ascii_mat <- matrix(ascii_mat, nrow = clen)
  }
  n <- ncol(ascii_mat)
  strings <- character(n)
  for (i in seq_len(n)) {
    chars <- ascii_mat[, i]
    non_zero <- which(chars != 0)
    if (length(non_zero) > 0) {
      strings[i] <- intToUtf8(chars[1:max(non_zero)])
    } else {
      strings[i] <- ""
    }
  }
  strings
}

string_to_ascii <- function(s, len) {
  chars <- utf8ToInt(substr(s, 1, len))
  if (length(chars) < len) {
    chars <- c(chars, rep(0L, len - length(chars)))
  }
  chars
}

#' Read expression values from multiple files into a numeric matrix
#' @param file_list Character vector of filenames
#' @param gene_ids Character vector of gene IDs to match
#' @param n_header_rows Number of header rows to skip in each file
#' @param gene_col Column index (1-based) of gene IDs in each file
#' @param value_cols Vector of column indices (1-based) of expression values in each file
#' @param start_row Row index (1-based) to start reading data in each file
#' @param delimiter Delimiter used in the files (default: tab)
#' @param n_samples Total number of samples (rows) to read across all files
#' @return A list with elements:
#'   - expression_vectors: A numeric matrix of dimensions (n_samples x n_genes)
#'   - ierr: Integer error code (0 if successful)
read_expression_vectors <- function(file_list, gene_ids, 
                                    n_header_rows, gene_col, value_cols, start_row, delimiter = "\t", n_samples) {
    nfiles <- length(file_list)
    ngenes <- length(gene_ids)
    gene_len <- max(nchar(gene_ids))

    # Convert file names to ASCII integers (pad with spaces, not zeros)
    file_ascii_list <- lapply(file_list, function(f) {
    raw <- charToRaw(enc2native(f))          # native single-byte encoding
    chars <- as.integer(raw)
    c(chars, rep(32L, 256 - length(chars)))  # pad with spaces (32)
    })
    file_ascii <- do.call(cbind, file_ascii_list)

  
    # Convert gene IDs to ASCII matrix
    gene_ascii <- strings_to_ascii_matrix(gene_ids, gene_len)
    
    # Convert delimiter to ASCII
    delimiter_ascii <- utf8ToInt(delimiter)
    
    # Prepare output
    expression_vectors <- double(n_samples * ngenes)
    
    out <- .Fortran("read_expression_vectors_R",
        file_list_ascii = as.integer(file_ascii),
        file_list_len = as.integer(256),  # Fixed length of 256
        n_files = as.integer(nfiles),
        gene_ids_ascii = as.integer(gene_ascii),
        gene_ids_len = as.integer(gene_len),
        n_genes = as.integer(ngenes),
        expression_vectors = as.double(expression_vectors),
        n_samples = as.integer(n_samples),
        n_genes2 = as.integer(ngenes),
        n_header_rows = as.integer(n_header_rows),
        gene_col = as.integer(gene_col),
        value_cols = as.integer(value_cols),
        n_value_cols = as.integer(length(value_cols)),
        start_row = as.integer(start_row),
        ierr = 0L,
        delimiter_ascii = as.integer(delimiter_ascii),
        dlen = as.integer(length(delimiter_ascii))
    )

    tox_errors(out$ierr)
    
    list(
        expression_vectors = matrix(out$expression_vectors, nrow = n_samples, ncol = ngenes),
        ierr = out$ierr
    )
}

#' Read gene IDs from a single file
#' @param filename Character string of the filename
#' @param ngenes Number of gene IDs to read
#' @param gene_len Maximum length of each gene ID
#' @param n_header_rows Number of header rows to skip
#' @param gene_col Column index (1-based) of gene IDs in the file
#' @return A list with elements:
#'   - gene_ids: Character vector of gene IDs read from the file
#'   - ierr: Integer error code (0 if successful)
read_gene_ids_from_file <- function(filename, ngenes, gene_len, n_header_rows = 1, gene_col = 1) {
  # Convert filename to ASCII integers
  filename_ascii <- utf8ToInt(filename)
  filename_ascii <- c(filename_ascii, rep(0L, 256 - length(filename_ascii)))  # Pad to 256
  
  gene_ascii <- matrix(0L, nrow = gene_len, ncol = ngenes)
  ierr <- integer(1)
  
  out <- .Fortran("read_gene_ids_from_file_R",
    filename_ascii = as.integer(filename_ascii),
    fn_len = as.integer(256),  # Fixed length of 256
    gene_ids_ascii = as.integer(gene_ascii),
    gene_ids_len = as.integer(gene_len),
    n_genes = as.integer(ngenes),
    n_header_rows = as.integer(n_header_rows),
    gene_col = as.integer(gene_col),
    ierr = 0
  )
  tox_errors(out$ierr)
  
  list(
    gene_ids = ascii_matrix_to_strings(out$gene_ids_ascii, gene_len),
    ierr = out$ierr
  )
}

#' Read gene family assignments from a file
#' @param filename Character string of the filename
#' @param gene_ids Character vector of gene IDs to match
#' @param n_families Number of unique gene families expected
#' @param family_len Maximum length of each family ID
#' @return A list with elements:
#'   - family_ids: Character vector of family IDs read from the file
#'   - gene_to_fam: Integer vector mapping each gene to its family index
#'   - ierr: Integer error code (0 if successful)
#' Note: If genes are not found in the gene IDs list a message will be printed but NO error code is thrown
read_family_file <- function(filename, gene_ids, n_families, family_len) {
  ngenes <- length(gene_ids)
  gene_len <- max(nchar(gene_ids))
  
  # Convert filename to ASCII integers
  filename_ascii <- utf8ToInt(filename)
  filename_ascii <- c(filename_ascii, rep(0L, 256 - length(filename_ascii)))  # Pad to 256
  
  gene_ascii <- strings_to_ascii_matrix(gene_ids, gene_len)
  family_ascii <- matrix(0L, nrow = family_len, ncol = n_families)
  gene_to_fam <- integer(ngenes)
  ierr <- integer(1)
  
  out <- .Fortran("read_family_file_R",
    filename_ascii = as.integer(filename_ascii),
    fn_len = as.integer(256),  # Fixed length of 256
    gene_ids_ascii = as.integer(gene_ascii),
    gene_ids_len = as.integer(gene_len),
    n_genes = as.integer(ngenes),
    family_ids_ascii = as.integer(family_ascii),
    family_ids_len = as.integer(family_len),
    n_families = as.integer(n_families),
    gene_to_fam = as.integer(gene_to_fam),
    ierr = 0
  )
  tox_errors(out$ierr)
  list(
    family_ids = ascii_matrix_to_strings(matrix(out$family_ids_ascii, nrow = family_len, ncol = n_families)),
    gene_to_fam = out$gene_to_fam,
    ierr = out$ierr
  )
}

#' Filter out genes without family assignments
#' @param gene_ids Character vector of gene IDs
#' @param expression_vectors Numeric matrix of expression values (n_samples x n_genes)
#' @param gene_to_fam Integer vector mapping each gene to its family index (0 if unassigned)
#' @return A list with elements:    
#'   - gene_ids: Filtered character vector of gene IDs
#'   - expression_vectors: Filtered numeric matrix of expression values
#'   - gene_to_fam: Filtered integer vector mapping each gene to its family index
#'   - n_genes_kept: Number of genes kept after filtering
#'   - ierr: Integer error code (0 if successful)
filter_unassigned_genes <- function(gene_ids, expression_vectors, gene_to_fam) {
  ngenes <- length(gene_ids)
  gene_len <- max(nchar(gene_ids))
  n_samples <- nrow(expression_vectors)
  ierr <- integer(1)

  gene_ascii <- strings_to_ascii_matrix(gene_ids, gene_len)
  
  out <- .Fortran("filter_unassigned_genes_R",
    gene_ids_ascii = as.integer(gene_ascii),
    gene_ids_len = as.integer(gene_len),
    n_genes = as.integer(ngenes),
    expression_vectors_flat = as.double(as.vector(expression_vectors)),
    n_samples = as.integer(n_samples),
    gene_to_fam = as.integer(gene_to_fam),
    mask = logical(ngenes),
    n_genes_kept = integer(1),
    ierr = 0
  )
  tox_errors(out$ierr)
  
  # Apply the mask on the R side
  mask <- out$mask
  n_genes_kept <- out$n_genes_kept
  
  list(
    gene_ids = gene_ids[mask],
    expression_vectors = expression_vectors[, mask, drop = FALSE],
    gene_to_fam = gene_to_fam[mask],
    n_genes_kept = n_genes_kept,
    ierr = out$ierr
  )
}

# tox_data_validation.R
# R wrappers for Fortran validation routines
# Uses ASCII conversion helpers: strings_to_ascii_matrix, ascii_matrix_to_strings

validate_data_structure <- function(n_genes, n_families, n_samples, d,
                                    gene_ids, gene_family_ids,
                                    gene_to_fam, expression_vectors,
                                    family_centroids, shift_vectors) {
  gene_len <- max(nchar(gene_ids))
  fam_len  <- max(nchar(gene_family_ids))
  gene_ascii <- strings_to_ascii_matrix(gene_ids, gene_len)
  fam_ascii  <- strings_to_ascii_matrix(gene_family_ids, fam_len)
  ierr <- integer(1)

  out <- .Fortran("validate_data_structure_R",
                  as.integer(n_genes),
                  as.integer(n_families),
                  as.integer(n_samples),
                  as.integer(d),
                  as.integer(gene_ascii),
                  as.integer(gene_len),
                  as.integer(fam_ascii),
                  as.integer(fam_len),
                  as.integer(gene_to_fam),
                  as.double(expression_vectors),
                  as.double(family_centroids),
                  as.double(shift_vectors),
                  ierr = ierr)
  tox_errors(out$ierr)
  list(ierr = out$ierr)
}

validate_gene_to_family_mapping <- function(gene_to_fam, n_genes, n_families) {
  ierr <- integer(1)
  out <- .Fortran("validate_gene_to_family_mapping_R",
                  as.integer(gene_to_fam),
                  as.integer(n_genes),
                  as.integer(n_families),
                  ierr = ierr)
  tox_errors(out$ierr)
  list(ierr = out$ierr)
}

validate_expression_data <- function(expression_vectors, n_genes, n_samples, check_non_negative = TRUE) {
  ierr <- integer(1)
  out <- .Fortran("validate_expression_data_R",
                  as.double(expression_vectors),
                  as.integer(n_genes),
                  as.integer(n_samples),
                  as.logical(check_non_negative),
                  ierr = ierr)
  tox_errors(out$ierr)
  list(ierr = out$ierr)
}

validate_family_centroids <- function(family_centroids, n_families, d) {
  ierr <- integer(1)
  out <- .Fortran("validate_family_centroids_R",
                  as.double(family_centroids),
                  as.integer(n_families),
                  as.integer(d),
                  ierr = ierr)
  tox_errors(out$ierr)
  list(ierr = out$ierr)
}

# Update the validate_shift_vectors function to match Fortran signature
validate_shift_vectors <- function(shift_vectors, expression_vectors, family_centroids, gene_to_fam, d) {
  ierr <- integer(1)
  
  # Get dimensions
  n_genes <- ncol(expression_vectors)
  n_samples <- nrow(expression_vectors)
  n_families <- ncol(family_centroids)
  
  out <- .Fortran("validate_shift_vectors_R",
                  as.double(shift_vectors),
                  as.double(expression_vectors),
                  as.double(family_centroids),
                  as.integer(gene_to_fam),
                  as.integer(d),
                  as.integer(n_genes),
                  as.integer(n_samples),
                  as.integer(n_families),
                  ierr = ierr)
  tox_errors(out$ierr)
  list(ierr = out$ierr)
}

validate_gene_ids_uniqueness <- function(gene_ids, n_genes) {
  gene_len <- max(nchar(gene_ids))
  gene_ascii <- strings_to_ascii_matrix(gene_ids, gene_len)
  ierr <- integer(1)
  out <- .Fortran("validate_gene_ids_uniqueness_R",
                  as.integer(gene_ascii),
                  as.integer(gene_len),
                  as.integer(n_genes),
                  ierr = ierr)
  tox_errors(out$ierr)
  list(ierr = out$ierr)
}

validate_family_ids_uniqueness <- function(family_ids, n_families) {
  fam_len <- max(nchar(family_ids))
  fam_ascii <- strings_to_ascii_matrix(family_ids, fam_len)
  ierr <- integer(1)
  out <- .Fortran("validate_family_ids_uniqueness_R",
                  as.integer(fam_ascii),
                  as.integer(fam_len),
                  as.integer(n_families),
                  ierr = ierr)
  tox_errors(out$ierr)
  list(ierr = out$ierr)
}

validate_empty_strings <- function(strings, n) {
  str_len <- max(nchar(strings))
  str_ascii <- strings_to_ascii_matrix(strings, str_len)
  ierr <- integer(1)
  out <- .Fortran("validate_empty_strings_R",
                  as.integer(str_ascii),
                  as.integer(str_len),
                  as.integer(n),
                  ierr = ierr)
  tox_errors(out$ierr)
  list(ierr = out$ierr)
}

validate_all_data <- function(n_genes, n_families, n_samples, d,
                              gene_ids, gene_family_ids,
                              gene_to_fam, expression_vectors,
                              family_centroids, shift_vectors) {
  # Determine max lengths for string encoding
  gene_len <- max(nchar(gene_ids))
  fam_len  <- max(nchar(gene_family_ids))

  # Convert to ASCII matrices
  gene_ascii <- strings_to_ascii_matrix(gene_ids, gene_len)
  fam_ascii  <- strings_to_ascii_matrix(gene_family_ids, fam_len)

  ierr <- integer(1)

  out <- .Fortran("validate_all_data_R",
                  as.integer(n_genes),
                  as.integer(n_families),
                  as.integer(n_samples),
                  as.integer(d),
                  as.integer(gene_ascii),
                  as.integer(gene_len),
                  as.integer(fam_ascii),
                  as.integer(fam_len),
                  as.integer(gene_to_fam),
                  as.double(expression_vectors),
                  as.double(family_centroids),
                  as.double(shift_vectors),
                  ierr = ierr)

  tox_errors(out$ierr)
  list(ierr = out$ierr)
}

# ONLY FOR TEST PURPOSES
tox_compute_shift_vector_field <- function(expression_vectors, family_centroids, gene_to_centroid) {
  # Input validation
  if (!is.matrix(expression_vectors)) {
    stop("expression_vectors must be a matrix")
  }

  if (!is.matrix(family_centroids)) {
    stop("family_centroids must be a matrix")
  }
  
  # Dimensions and counts
  n_axes_genes <- nrow(expression_vectors)
  n_vectors <- ncol(expression_vectors)
  n_axes_centroids <- nrow(family_centroids)
  n_families <- ncol(family_centroids)
  
  # Validate length of gene_to_centroid
  if (n_vectors != length(gene_to_centroid)) {
    stop("number of expression_vectors must be equal to length of gene_to_centroid")
  }

  # Validate dimensions
  if (n_axes_genes != n_axes_centroids) {
    stop("family_centroids must have the same number of axes as expression_vectors")
  }
  
  # Prepare output arrays
  shift_vectors <- matrix(0.0, nrow = 2 *n_axes_genes, ncol = n_vectors)
  ierr <- as.integer(0)
  
  # Call Fortran wrapper
  result <- .Fortran("compute_shift_vector_field_r",
                     n_axes_genes = as.integer(n_axes_genes),
                     n_vectors = as.integer(n_vectors),
                     n_families = as.integer(n_families),
                     expression_vectors = as.double(expression_vectors),
                     family_centroids = as.double(family_centroids),
                     gene_to_centroid = as.integer(gene_to_centroid),
                     shift_vectors = as.double(shift_vectors),
                     ierr = ierr)
  
  # Check for errors and throw informative messages
  tox_errors(result$ierr)
  
  # Return structured result (no ierr since we checked for errors)
  return(list(
    shift_vectors = result$shift_vectors
  ))
}

tox_group_centroid <- function(expression_vectors, gene_to_family, n_families, ortholog_set, mode = 'all') {
  
  # 1) Validate inputs
  if (!is.matrix(expression_vectors) || !is.numeric(expression_vectors)) {
    stop("`expression_vectors` must be a numeric matrix.")
  }
  n_axes <- nrow(expression_vectors)
  n_genes <- ncol(expression_vectors)

  if (!is.integer(gene_to_family) || length(gene_to_family) != n_genes) {
    stop("`gene_to_family` must be an integer vector of length n_genes.")
  }
  if (!is.logical(ortholog_set) || length(ortholog_set) != n_genes) {
    stop("`ortholog_set` must be a logical vector of length n_genes.")
  }
  if (!mode %in% c('all', 'ortho')) {
    stop("`mode` must be either 'all' or 'ortho'.")
  }

  # 2) Prepare inputs/outputs for Fortran
  use_all_mode <- (mode == 'all')
  centroid_matrix_out <- matrix(0.0, nrow = n_axes, ncol = n_families)
  selected_indices_ws <- integer(n_genes) # Workspace buffer
  ierr <- as.integer(0)
  # 3) Call Fortran
  result <- .Fortran("group_centroid_r",
                     expression_vectors = as.double(expression_vectors),
                     n_axes = as.integer(n_axes),
                     n_genes = as.integer(n_genes),
                     gene_to_family = as.integer(gene_to_family),
                     num_families = as.integer(n_families),
                     centroid_matrix = centroid_matrix_out,
                     use_all_mode = as.logical(use_all_mode),
                     ortholog_set = as.logical(ortholog_set),
                     selected_indices = selected_indices_ws,
                     selected_indices_len = as.integer(n_genes),
                     ierr = ierr)
  
  # Check for errors and throw informative messages
  tox_errors(result$ierr)

  # 4) Return the populated output matrix (no ierr since we checked for errors)
  return(result$centroid_matrix)
}
