dyn.load("build/libtensor-omics.so")
source("r/tensoromics_functions.R")
source("r/error_handling.R")

strings_to_raw_matrix <- function(arr, clen) {
  n <- length(arr)
  # Create a matrix of raw bytes with dimensions clen x n
  mat <- matrix(raw(1), nrow = clen, ncol = n)
  for (i in seq_along(arr)) {
    # Convert string to raw bytes
    raw_bytes <- charToRaw(arr[i])
    length_bytes <- length(raw_bytes)
    if (length_bytes > 0) {
      # Copy raw bytes into matrix
      mat[1:min(length_bytes, clen), i] <- raw_bytes[1:min(length_bytes, clen)]
    }
    # Note: No null termination needed - Fortran handles this
  }
  mat  # Return the raw matrix
}

raw_matrix_to_strings <- function(raw_mat, clen) {
  # raw_mat: raw matrix with dimensions clen x n
  if (is.null(dim(raw_mat))) {
    # 1D vector: convert to matrix with one column
    raw_mat <- matrix(raw_mat, nrow = clen)
  }
  n <- ncol(raw_mat)
  strings <- character(n)
  for (i in seq_len(n)) {
    raw_vec <- raw_mat[, i]
    # Find first null byte (if any)
    null_pos <- which(raw_vec == as.raw(0))
    if (length(null_pos) > 0) {
      end_pos <- min(null_pos[1] - 1, length(raw_vec))
    } else {
      end_pos <- length(raw_vec)
    }
    if (end_pos > 0) {
      strings[i] <- rawToChar(raw_vec[1:end_pos])
    } else {
      strings[i] <- ""
    }
  }
  strings
}

string_to_raw <- function(s, len) {
  raw_bytes <- charToRaw(s)
  if (length(raw_bytes) < len) {
    # Pad with zeros if needed (Fortran will handle null termination)
    raw_bytes <- c(raw_bytes, raw(1))
    length(raw_bytes) <- len
    raw_bytes[is.na(raw_bytes)] <- as.raw(0)
  }
  raw_bytes[1:len]  # Ensure exact length
}

# Keep the old integer conversion functions for functions that still need them
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
                                    n_header_rows, gene_col, value_cols, delimiter = "\t", n_samples) {
    nfiles <- length(file_list)
    ngenes <- length(gene_ids)
    gene_len <- max(nchar(gene_ids)) + 1  # +1 for null terminator

    # Convert file names to raw bytes
    file_raw_list <- lapply(file_list, function(f) {
      raw_bytes <- charToRaw(enc2native(f))
      # Pad with zeros to fixed length of 256
      if (length(raw_bytes) < 256) {
        raw_bytes <- c(raw_bytes, raw(256 - length(raw_bytes)))
      }
      raw_bytes[1:256]
    })
    file_raw <- do.call(cbind, file_raw_list)

    # Convert gene IDs to raw matrix
    gene_raw <- strings_to_raw_matrix(gene_ids, gene_len)
    
    # Convert delimiter to raw
    delimiter_raw <- charToRaw(delimiter)
    
    # Prepare output
    expression_vectors <- double(n_samples * ngenes)
    
    out <- .Fortran("read_expression_vectors_R",
        file_list_raw = file_raw,              # Pass raw bytes directly
        file_list_len = as.integer(256),       # Fixed length of 256
        n_files = as.integer(nfiles),
        gene_ids_raw = gene_raw,               # Pass raw bytes directly
        gene_ids_len = as.integer(gene_len),
        n_genes = as.integer(ngenes),
        expression_vectors = as.double(expression_vectors),
        n_samples = as.integer(n_samples),
        n_header_rows = as.integer(n_header_rows),
        gene_col = as.integer(gene_col),
        value_cols = as.integer(value_cols),
        n_value_cols = as.integer(length(value_cols)),
        ierr = 0L,
        delimiter_raw = delimiter_raw,         # Pass raw bytes directly
        dlen = as.integer(length(delimiter_raw))
    )

    check_err_code(out$ierr)
    
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
  # Convert filename to raw bytes
  filename_raw <- charToRaw(filename)
  
  gene_raw <- matrix(raw(1), nrow = gene_len + 1, ncol = ngenes)  # +1 for null terminator
  ierr <- integer(1)
  
  out <- .Fortran("read_gene_ids_from_file_R",
    filename_raw = filename_raw,              # Pass raw bytes directly
    fn_len = as.integer(length(filename_raw)),
    gene_ids_raw = gene_raw,                  # Pass raw bytes directly
    gene_ids_len = as.integer(gene_len + 1),  # +1 for null terminator
    n_genes = as.integer(ngenes),
    n_header_rows = as.integer(n_header_rows),
    gene_col = as.integer(gene_col),
    ierr = 0
  )
  check_err_code(out$ierr)
  
  list(
    gene_ids = raw_matrix_to_strings(out$gene_ids_raw, gene_len + 1),
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
  gene_len <- max(nchar(gene_ids)) + 1  # +1 for null terminator
  
  # Convert filename to raw bytes
  filename_raw <- charToRaw(filename)
  
  gene_raw <- strings_to_raw_matrix(gene_ids, gene_len)
  family_raw <- matrix(raw(1), nrow = family_len + 1, ncol = n_families)  # +1 for null terminator
  gene_to_fam <- integer(ngenes)
  ierr <- integer(1)
  
  out <- .Fortran("read_family_file_R",
    filename_raw = filename_raw,              # Pass raw bytes directly
    fn_len = as.integer(length(filename_raw)),
    gene_ids_raw = gene_raw,                  # Pass raw bytes directly
    gene_ids_len = as.integer(gene_len),
    n_genes = as.integer(ngenes),
    family_ids_raw = family_raw,              # Pass raw bytes directly
    family_ids_len = as.integer(family_len + 1),  # +1 for null terminator
    n_families = as.integer(n_families),
    gene_to_fam = as.integer(gene_to_fam),
    ierr = 0
  )
  check_err_code(out$ierr)
  list(
    family_ids = raw_matrix_to_strings(matrix(out$family_ids_raw, nrow = family_len + 1, ncol = n_families)),
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
  gene_len <- max(nchar(gene_ids)) + 1  # +1 for null terminator
  n_samples <- nrow(expression_vectors)
  ierr <- integer(1)

  gene_raw <- strings_to_raw_matrix(gene_ids, gene_len)
  
  out <- .Fortran("filter_unassigned_genes_R",
    gene_ids_raw = gene_raw,                  # Pass raw bytes directly
    gene_ids_len = as.integer(gene_len),
    n_genes = as.integer(ngenes),
    expression_vectors_flat = as.double(as.vector(expression_vectors)),
    n_samples = as.integer(n_samples),
    gene_to_fam = as.integer(gene_to_fam),
    mask = logical(ngenes),
    n_genes_kept = integer(1),
    ierr = 0
  )
  check_err_code(out$ierr)
  
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

# R wrappers for Fortran validation routines
# Uses raw conversion helpers: strings_to_raw_matrix, raw_matrix_to_strings

validate_data_structure <- function(n_genes, n_families, n_samples, d,
                                    gene_ids, gene_family_ids,
                                    gene_to_fam, expression_vectors,
                                    family_centroids, shift_vectors) {
  gene_len <- max(nchar(gene_ids)) + 1  # +1 for null terminator
  fam_len  <- max(nchar(gene_family_ids)) + 1  # +1 for null terminator
  gene_raw <- strings_to_raw_matrix(gene_ids, gene_len)
  fam_raw  <- strings_to_raw_matrix(gene_family_ids, fam_len)
  ierr <- integer(1)

  out <- .Fortran("validate_data_structure_R",
                  as.integer(n_genes),
                  as.integer(n_families),
                  as.integer(n_samples),
                  gene_raw,                    # Pass raw bytes directly
                  as.integer(gene_len),
                  fam_raw,                     # Pass raw bytes directly
                  as.integer(fam_len),
                  as.integer(gene_to_fam),
                  as.double(expression_vectors),
                  as.double(family_centroids),
                  as.double(shift_vectors),
                  ierr = ierr)
  check_err_code(out$ierr)
  list(ierr = out$ierr)
}

validate_gene_to_family_mapping <- function(gene_to_fam, n_genes, n_families) {
  ierr <- integer(1)
  out <- .Fortran("validate_gene_to_family_mapping_R",
                  as.integer(gene_to_fam),
                  as.integer(n_genes),
                  as.integer(n_families),
                  ierr = ierr)
  check_err_code(out$ierr)
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
  check_err_code(out$ierr)
  list(ierr = out$ierr)
}

validate_family_centroids <- function(family_centroids, n_families, n_samples) {
  ierr <- integer(1)
  out <- .Fortran("validate_family_centroids_R",
                  as.double(family_centroids),
                  as.integer(n_families),
                  as.integer(n_samples),
                  ierr = ierr)
  check_err_code(out$ierr)
  list(ierr = out$ierr)
}

validate_shift_vectors <- function(shift_vectors, expression_vectors, family_centroids, gene_to_fam, n_samples) {
  ierr <- integer(1)
  
  # Get dimensions
  n_genes <- ncol(expression_vectors)
  n_families <- ncol(family_centroids)
  
  out <- .Fortran("validate_shift_vectors_R",
                  as.double(shift_vectors),
                  as.double(expression_vectors),
                  as.double(family_centroids),
                  as.integer(gene_to_fam),
                  as.integer(n_genes),
                  as.integer(n_samples),
                  as.integer(n_families),
                  ierr = ierr)
  check_err_code(out$ierr)
  list(ierr = out$ierr)
}

validate_gene_ids_uniqueness <- function(gene_ids, n_genes) {
  gene_len <- max(nchar(gene_ids)) + 1  # +1 for null terminator
  gene_raw <- strings_to_raw_matrix(gene_ids, gene_len)
  ierr <- integer(1)
  out <- .Fortran("validate_gene_ids_uniqueness_R",
                  gene_raw,                    # Pass raw bytes directly
                  as.integer(gene_len),
                  as.integer(n_genes),
                  ierr = ierr)
  check_err_code(out$ierr)
  list(ierr = out$ierr)
}

validate_family_ids_uniqueness <- function(family_ids, n_families) {
  fam_len <- max(nchar(family_ids)) + 1  # +1 for null terminator
  fam_raw <- strings_to_raw_matrix(family_ids, fam_len)
  ierr <- integer(1)
  out <- .Fortran("validate_family_ids_uniqueness_R",
                  fam_raw,                     # Pass raw bytes directly
                  as.integer(fam_len),
                  as.integer(n_families),
                  ierr = ierr)
  check_err_code(out$ierr)
  list(ierr = out$ierr)
}

validate_all_data <- function(n_genes, n_families, n_samples,
                              gene_ids, gene_family_ids,
                              gene_to_fam, expression_vectors,
                              family_centroids, shift_vectors) {
  # Determine max lengths for string encoding (+1 for null terminator)
  gene_len <- max(nchar(gene_ids)) + 1
  fam_len  <- max(nchar(gene_family_ids)) + 1

  # Convert to raw matrices
  gene_raw <- strings_to_raw_matrix(gene_ids, gene_len)
  fam_raw  <- strings_to_raw_matrix(gene_family_ids, fam_len)

  ierr <- integer(1)

  out <- .Fortran("validate_all_data_R",
                  as.integer(n_genes),
                  as.integer(n_families),
                  as.integer(n_samples),
                  gene_raw,                    # Pass raw bytes directly
                  as.integer(gene_len),
                  fam_raw,                     # Pass raw bytes directly
                  as.integer(fam_len),
                  as.integer(gene_to_fam),
                  as.double(expression_vectors),
                  as.double(family_centroids),
                  as.double(shift_vectors),
                  ierr = ierr)

  check_err_code(out$ierr)
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
  check_err_code(result$ierr)
  
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
  check_err_code(result$ierr)

  # 4) Return the populated output matrix (no ierr since we checked for errors)
  return(result$centroid_matrix)
}

tox_save_data_archive <- function(zip_filename,
                                 gene_ids = NULL, gene_ids_name = NULL,
                                 expression_vectors = NULL, expression_vectors_name = NULL,
                                 gene_to_fam = NULL, gene_to_fam_name = NULL,
                                 family_ids = NULL, family_ids_name = NULL,
                                 family_centroids = NULL, family_centroids_name = NULL,
                                 shift_vectors = NULL, shift_vectors_name = NULL) {
  
  # Flags initialization
  gene_ids_array_valid <- gene_ids_name_valid <- FALSE
  expression_vectors_array_valid <- expression_vectors_name_valid <- FALSE
  gene_to_fam_array_valid <- gene_to_fam_name_valid <- FALSE
  family_ids_array_valid <- family_ids_name_valid <- FALSE
  family_centroids_array_valid <- family_centroids_name_valid <- FALSE
  shift_vectors_array_valid <- shift_vectors_name_valid <- FALSE
  
  if (!is.character(zip_filename)) {
    stop("Type mismatch: zip_filename must be a string")
  }
  
  # Validations
  if (!is.null(gene_ids)) {
    if (length(dim(gene_ids)) == 1 || is.vector(gene_ids)) {
      gene_ids_array_valid <- TRUE
    } else {
      stop(paste("Gene IDs dimensions mismatch: Expected 1 but got", length(dim(gene_ids))))
    }
  }
  if (!is.null(gene_ids_name)) {
    if (is.character(gene_ids_name)) {
      gene_ids_name_valid <- TRUE
    } else {
      stop("Gene IDs name must be a string")
    }
  }
  if (xor(gene_ids_array_valid, gene_ids_name_valid)) {
    message("Gene IDs: Either provide array and filename or none. Skipping.")
  }
  
  if (!is.null(expression_vectors)) {
    if (length(dim(expression_vectors)) == 2) {
      expression_vectors_array_valid <- TRUE
    } else {
      stop(paste("Expression vectors dim mismatch: Expected 2 but got", length(dim(expression_vectors))))
    }
  }
  if (!is.null(expression_vectors_name)) {
    if (is.character(expression_vectors_name)) {
      expression_vectors_name_valid <- TRUE
    } else {
      stop("Expression vector name must be a string")
    }
  }
  if (xor(expression_vectors_array_valid, expression_vectors_name_valid)) {
    message("Expression vectors: Either provide array and filename or none. Skipping.")
  }
  
  if (!is.null(gene_to_fam)) {
    if (length(dim(gene_to_fam)) == 1 || is.vector(gene_to_fam)) {
      gene_to_fam_array_valid <- TRUE
    } else {
      stop(paste("Gene to family dim mismatch: Expected 1 but got", length(dim(gene_to_fam))))
    }
  }
  if (!is.null(gene_to_fam_name)) {
    if (is.character(gene_to_fam_name)) {
      gene_to_fam_name_valid <- TRUE
    } else {
      stop("Gene to family name must be a string")
    }
  }
  if (xor(gene_to_fam_array_valid, gene_to_fam_name_valid)) {
    message("Gene to family: Either provide array and filename or none. Skipping.")
  }
  
  if (!is.null(family_ids)) {
    if (length(dim(family_ids)) == 1 || is.vector(family_ids)) {
      family_ids_array_valid <- TRUE
    } else {
      stop(paste("Family IDs dim mismatch: Expected 1 but got", length(dim(family_ids))))
    }
  }
  if (!is.null(family_ids_name)) {
    if (is.character(family_ids_name)) {
      family_ids_name_valid <- TRUE
    } else {
      stop("Family IDs name must be a string")
    }
  }
  if (xor(family_ids_array_valid, family_ids_name_valid)) {
    message("Family IDs: Either provide array and filename or none. Skipping.")
  }
  
  if (!is.null(family_centroids)) {
    if (length(dim(family_centroids)) == 2) {
      family_centroids_array_valid <- TRUE
    } else {
      stop(paste("Family centroids dim mismatch: Expected 2 but got", length(dim(family_centroids))))
    }
  }
  if (!is.null(family_centroids_name)) {
    if (is.character(family_centroids_name)) {
      family_centroids_name_valid <- TRUE
    } else {
      stop("Family centroids name must be a string")
    }
  }
  if (xor(family_centroids_array_valid, family_centroids_name_valid)) {
    message("Family centroids: Either provide array and filename or none. Skipping.")
  }
  
  if (!is.null(shift_vectors)) {
    if (length(dim(shift_vectors)) == 2) {
      shift_vectors_array_valid <- TRUE
    } else {
      stop(paste("Shift vectors dim mismatch: Expected 2 but got", length(dim(shift_vectors))))
    }
  }
  if (!is.null(shift_vectors_name)) {
    if (is.character(shift_vectors_name)) {
      shift_vectors_name_valid <- TRUE
    } else {
      stop("Shift vectors name must be a string")
    }
  }
  if (xor(shift_vectors_array_valid, shift_vectors_name_valid)) {
    message("Shift vectors: Either provide array and filename or none. Skipping.")
  }
  
  # Manifest as list
  manifest_lines <- character()
  
  # Create temporary directory for files
  temp_dir <- tempfile()
  dir.create(temp_dir)
  
  tryCatch({
    # Write files to temporary directory
    if (gene_ids_array_valid && gene_ids_name_valid) {
      tox_serialize_char_array(gene_ids, file.path(temp_dir, gene_ids_name))
      message(paste("Wrote gene IDs to", gene_ids_name))
      manifest_lines <- c(manifest_lines, paste0("gene_ids=", gene_ids_name))
    }
    
    if (expression_vectors_array_valid && expression_vectors_name_valid) {
      tox_serialize_real_array(expression_vectors, file.path(temp_dir, expression_vectors_name))
      message(paste("Wrote expression vectors to", expression_vectors_name))
      manifest_lines <- c(manifest_lines, paste0("expression=", expression_vectors_name))
    }
    
    if (gene_to_fam_array_valid && gene_to_fam_name_valid) {
      tox_serialize_int_array(gene_to_fam, file.path(temp_dir, gene_to_fam_name))
      message(paste("Wrote gene to family to", gene_to_fam_name))
      manifest_lines <- c(manifest_lines, paste0("gene_to_family=", gene_to_fam_name))
    }
    
    if (family_ids_array_valid && family_ids_name_valid) {
      tox_serialize_char_array(family_ids, file.path(temp_dir, family_ids_name))
      message(paste("Wrote family IDs to", family_ids_name))
      manifest_lines <- c(manifest_lines, paste0("family_ids=", family_ids_name))
    }
    
    if (family_centroids_array_valid && family_centroids_name_valid) {
      tox_serialize_real_array(family_centroids, file.path(temp_dir, family_centroids_name))
      message(paste("Wrote family centroids to", family_centroids_name))
      manifest_lines <- c(manifest_lines, paste0("family_centroids=", family_centroids_name))
    }
    
    if (shift_vectors_array_valid && shift_vectors_name_valid) {
      tox_serialize_real_array(shift_vectors, file.path(temp_dir, shift_vectors_name))
      message(paste("Wrote shift vectors to", shift_vectors_name))
      manifest_lines <- c(manifest_lines, paste0("shift_vectors=", shift_vectors_name))
    }
    
    # Write manifest file
    writeLines(manifest_lines, file.path(temp_dir, "manifest.txt"))
    
    # Create zip archive with same compression level as Python (level 6)
    files_to_zip <- list.files(temp_dir, full.names = TRUE)
    zip::zip(zip_filename, files = files_to_zip, mode = "cherry-pick", 
             include_directories = FALSE, compression_level = 6)
    
  }, finally = {
    # Clean up temporary directory
    unlink(temp_dir, recursive = TRUE)
  })
}

tox_read_data_archive <- function(zip_filename,
                                 gene_ids = NULL,
                                 expression_vectors = NULL,
                                 gene_to_fam = NULL,
                                 family_ids = NULL,
                                 family_centroids = NULL,
                                 shift_vectors = NULL) {
  
  if (!is.character(zip_filename) || nchar(zip_filename) == 0) {
    stop("Zip name needs to be a non-empty string")
  }
  
  result <- list(
    gene_ids = NULL,
    expression_vectors = NULL,
    gene_to_fam = NULL,
    family_ids = NULL,
    family_centroids = NULL,
    shift_vectors = NULL
  )
  
  # Create temporary directory for extraction
  temp_dir <- tempfile()
  dir.create(temp_dir)
  
  tryCatch({
    # Extract all files
    extracted_files <- zip::unzip(zip_filename, exdir = temp_dir)
    
    # Read manifest file - use the actual file path
    manifest_path <- file.path(temp_dir, "manifest.txt")
    if (!file.exists(manifest_path)) {
      stop("Manifest file not found in archive")
    }
    
    manifest_lines <- readLines(manifest_path)
    
    # Parse manifest
    file_mapping <- list()
    for (line in manifest_lines) {
      parts <- strsplit(line, "=")[[1]]
      if (length(parts) == 2) {
        file_mapping[[parts[1]]] <- parts[2]
      }
    }
    
    # Extract and read files based on parameters
    if (!is.null(gene_ids) && !is.null(file_mapping[["gene_ids"]])) {
      filename <- file_mapping[["gene_ids"]]
      file_path <- file.path(temp_dir, filename)
      if (file.exists(file_path)) {
        result$gene_ids <- tox_deserialize_char_array(file_path)
        message(paste("Gene ids extracted from", filename))
      }
    }
    
    if (!is.null(expression_vectors) && !is.null(file_mapping[["expression"]])) {
      filename <- file_mapping[["expression"]]
      file_path <- file.path(temp_dir, filename)
      if (file.exists(file_path)) {
        result$expression_vectors <- tox_deserialize_real_array(file_path)
        message(paste("Expression Vectors extracted from", filename))
      }
    }
    
    if (!is.null(gene_to_fam) && !is.null(file_mapping[["gene_to_family"]])) {
      filename <- file_mapping[["gene_to_family"]]
      file_path <- file.path(temp_dir, filename)
      if (file.exists(file_path)) {
        result$gene_to_fam <- tox_deserialize_int_array(file_path)
        message(paste("Gene to family mapping extracted from", filename))
      }
    }
    
    if (!is.null(family_ids) && !is.null(file_mapping[["family_ids"]])) {
      filename <- file_mapping[["family_ids"]]
      file_path <- file.path(temp_dir, filename)
      if (file.exists(file_path)) {
        result$family_ids <- tox_deserialize_char_array(file_path)
        message(paste("Extracted family IDs from", filename))
      }
    }
    
    if (!is.null(family_centroids) && !is.null(file_mapping[["family_centroids"]])) {
      filename <- file_mapping[["family_centroids"]]
      file_path <- file.path(temp_dir, filename)
      if (file.exists(file_path)) {
        result$family_centroids <- tox_deserialize_real_array(file_path)
        message(paste("Extracted family centroids from", filename))
      }
    }
    
    if (!is.null(shift_vectors) && !is.null(file_mapping[["shift_vectors"]])) {
      filename <- file_mapping[["shift_vectors"]]
      file_path <- file.path(temp_dir, filename)
      if (file.exists(file_path)) {
        result$shift_vectors <- tox_deserialize_real_array(file_path)
        message(paste("Extracted shift vectors from", filename))
      }
    }
    
  }, finally = {
    # Clean up temporary directory
    unlink(temp_dir, recursive = TRUE)
  })
  
  return(result)
}