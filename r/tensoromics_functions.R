dyn.load("./build/libtensor-omics.so")

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

tox_get_array_metadata <- function(filename, max_dims = 5, with_clen = FALSE) {
  ascii <- utf8ToInt(filename)
  dims <- integer(max_dims)
  ndims <- integer(1)
  ierr <- integer(1)
  
  if (with_clen) {
    clen <- integer(1)
    res <- .Fortran("get_array_metadata_r",
                    as.integer(ascii),          # filename_ascii
                    as.integer(length(ascii)),  # fn_len
                    dims,                       # dims_out
                    ndims,                      # ndims
                    ierr,                       # ierr
                    clen)                       # clen (present)
    
    check_err_code(res[[5]])  # ierr
    return(list(
      dims = res[[3]][1:res[[4]]],  # dims_out[1:ndims]
      ndim = res[[4]],              # ndims
      clen = res[[6]]               # clen
    ))
    
  } else {
    res <- .Fortran("get_array_metadata_r",
                    as.integer(ascii),          # filename_ascii
                    as.integer(length(ascii)),  # fn_len
                    dims,                       # dims_out
                    ndims,                      # ndims
                    ierr)                       # ierr (kein clen)
    
    check_err_code(res[[5]])  # ierr
    return(list(
      dims = res[[3]][1:res[[4]]],  # dims_out[1:ndims]
      ndim = res[[4]]               # ndims
    ))
  }
}

# deserializes an integer array from a file, reads array dimensions first and then creates a proper array
# That is then being filled by fortran
tox_deserialize_int_array <- function(filename, max_dims = 5) {
    ascii <- utf8ToInt(filename)

    meta <- tox_get_array_metadata(filename, max_dims)
    total_size <- prod(meta$dims)

    flat <- integer(total_size)
    ndim <- integer(1)
    ierr <- integer(1)

    res <- .Fortran("deserialize_int_r",
                flat_arr = flat,
                arr_size = as.integer(total_size),
                filename_ascii = as.integer(ascii),
                fn_len = as.integer(length(ascii)),
                ierr = ierr)
    check_err_code(res$ierr)

    array(res$flat_arr[1:prod(meta$dims)], dim = meta$dims)
}

# Deserializes a real array from a file, reads array dimensions first and then creates a proper array
# That is then being filled by fortran
tox_deserialize_real_array <- function(filename, max_dims = 5) {
    ascii <- utf8ToInt(filename)

    meta <- tox_get_array_metadata(filename, max_dims)
    total_size <- prod(meta$dims)

    flat <- double(total_size)
    dims <- as.integer(meta$dims)
    ndim <- integer(1)
    ierr <- integer(1)

    res <- .Fortran("deserialize_real_flat_r",
                flat_arr = flat,
                arr_size = as.integer(total_size),
                filename_ascii = as.integer(ascii),
                fn_len = as.integer(length(ascii)),
                ierr = ierr)
    check_err_code(res$ierr)
    array(res$flat_arr[1:prod(meta$dims)], dim = meta$dims)
}

# Deserializes a character array from a file, reads array dimensions and character length first
# Then creates a proper array that is then being filled by fortran
# Note that the array needs to be translated back to characters
tox_deserialize_char_array <- function(filename, max_dims = 5) {
  ascii <- utf8ToInt(filename)
  dims <- integer(max_dims)
  ndim <- integer(1)
  clen <- integer(1)
  ierr <- integer(1)
  # Load metadata dimensions + clen
  meta <- tox_get_array_metadata(filename, max_dims, with_clen = TRUE)

  actual_dims <- meta$dims
  clen <- meta$clen
  total_array_size <- prod(actual_dims)
  cat("actual_dims:", actual_dims, "clen:", clen, "\n")

  ascii_arr <- integer(clen * total_array_size)

  res <- .Fortran("deserialize_char_flat_r",
    ascii_arr = ascii_arr,
    arr_size = as.integer(clen * total_array_size),
    filename_ascii = ascii,
    fn_len = as.integer(length(ascii)),
    ierr = ierr
  )
  check_err_code(res$ierr)
  # translate ASCII back to char
  mat <- matrix(res$ascii_arr, nrow = clen)
  chars <- apply(mat, 2, function(col) rawToChar(as.raw(col[col > 0])))

  array(chars, dim = meta$dims[1:meta$ndim])
}


# BASE R arrays are column-major just like fortran, so no serialization is needed for the array structure.
# Array can simply be passed with with as.integer()
tox_serialize_int_array <- function(arr, filename) {
  flat <- as.integer(arr)
  dims <- if (is.null(dim(arr))) {
    as.integer(length(arr))  # 1D-Vector
  } else {
    as.integer(dim(arr))
  }
  ndim <- as.integer(length(dims))
  ascii <- utf8ToInt(filename)
  ierr <- integer(1)

  res <- .Fortran("serialize_int_flat_r",
           arr = flat,
           array_size = length(flat),
           dims = dims,
           ndim = ndim,
           filename_ascii = as.integer(ascii),
           fn_len = as.integer(length(ascii)),
           ierr = ierr)
  check_err_code(res$ierr)
}

# BASE R arrays are column-major just like fortran, so no serialization is needed for the array structure.
# Array can simply be passed with with as.double() to pass it in a flat format.
tox_serialize_real_array <- function(arr, filename) {
  flat <- as.double(arr)

  dims <- if (is.null(dim(arr))) {
    as.integer(length(arr))  # 1D-Vector
  } else {
    as.integer(dim(arr))
  }

  ndim <- as.integer(length(dims))
  ascii <- utf8ToInt(filename)
  ierr <- integer(1)

  res <- .Fortran("serialize_real_flat_r",
           arr = flat,
           array_size = length(flat),
           dims = dims,
           ndim = ndim,
           filename_ascii = as.integer(ascii),
           fn_len = as.integer(length(ascii)),
           ierr = ierr)
  check_err_code(res$ierr)
}

# Serializes a character array to a file, encoding it as an integer matrix
# Each character is converted to its ASCII integer representation
# The matrix is then serialized with Fortran
tox_serialize_char_array <- function(arr, filename) {
  stopifnot(is.character(arr))
  arr <- as.array(arr)
  dims <- dim(arr)
  if (is.null(dims)) dims <- length(arr)
  clen <- max(nchar(arr, type = "chars"))
  ierr <- integer(1)

  # encode to integer matrix
  # Chars can not be passed via .Fortran directly
  mat <- matrix(0L, nrow = clen, ncol = length(arr))
  for (i in seq_along(arr)) {
    chars <- utf8ToInt(substr(arr[i], 1, clen))
    mat[seq_along(chars), i] <- chars
  }

  res <- .Fortran("serialize_char_flat_r",
    ascii_arr = as.integer(mat),
    array_size = length(mat),
    dims = as.integer(dims),
    ndim = as.integer(length(dims)),
    clen = as.integer(clen),
    filename_ascii = utf8ToInt(filename),
    fn_len = nchar(filename),
    ierr = ierr
  )
  check_err_code(res$ierr)
}