dyn.load("./build/libtensor-omics.so")

# returns the dimensions of an array stored in a file, works for integer and real
get_array_dims <- function(filename, max_dims = 5) {
  ascii <- utf8ToInt(filename)
  dims <- integer(max_dims)
  ndims <- integer(1)

  res <- .Fortran("get_array_dims_r",
                  filename_ascii = as.integer(ascii),
                  fn_len = as.integer(length(ascii)),
                  dims_out = dims,
                  ndims = ndims)

  res$dims_out[1:res$ndims]
}

#returns the dimensions, type code and character length of an array stored in a file; works for character arrays
get_array_metadata_chars <- function(filename, max_dims = 5) {
  ascii <- utf8ToInt(filename)
  dims <- integer(max_dims)
  ndims <- integer(1)
  typecode <- integer(1)
  clen <- integer(1)

  res <- .Fortran("get_array_metadata_chars_r",
                  filename_ascii = as.integer(ascii),
                  fn_len = as.integer(length(ascii)),
                  dims_out = dims,
                  ndims = ndims,
                  type_code_out = typecode,
                  clen_out = clen)

  list(dims = res$dims_out[1:res$ndims],
       type = res$type_code_out,
       clen = res$clen_out)
}

# deserializes an integer array from a file, reads array dimensions first and then creates a proper array
# That is then being filled by fortran
deserialize_int_array <- function(filename, max_dims = 5) {
    ascii <- utf8ToInt(filename)

    actual_dims <- get_array_dims(filename, max_dims)
    total_size <- prod(actual_dims)

    flat <- integer(total_size)
    dims <- integer(max_dims)
    ndim <- integer(1)

    res <- .Fortran("deserialize_int_r",
                flat_arr = flat,
                dims_out = dims,
                ndim_out = ndim,
                filename_ascii = as.integer(ascii),
                fn_len = as.integer(length(ascii)))

    actual_dims <- res$dims_out[1:res$ndim_out]
    array(res$flat_arr[1:prod(actual_dims)], dim = actual_dims)
}

# Deserializes a real array from a file, reads array dimensions first and then creates a proper array
# That is then being filled by fortran
deserialize_real_array <- function(filename, max_dims = 5) {
    ascii <- utf8ToInt(filename)

    actual_dims <- get_array_dims(filename, max_dims)
    total_size <- prod(actual_dims)

    flat <- double(total_size)
    dims <- integer(max_dims)
    ndim <- integer(1)

    res <- .Fortran("deserialize_real_flat_r",
                flat_arr = flat,
                dims_out = dims,
                ndim_out = ndim,
                filename_ascii = as.integer(ascii),
                fn_len = as.integer(length(ascii)))
    actual_dims <- res$dims_out[1:res$ndim_out]
    array(res$flat_arr[1:prod(actual_dims)], dim = actual_dims)
}

# Deserializes a character array from a file, reads array dimensions and character length first
# Then creates a proper array that is then being filled by fortran
# Note that the array needs to be translated back to characters
deserialize_char_array <- function(filename, max_dims = 5) {
  ascii <- utf8ToInt(filename)
  dims <- integer(max_dims)
  ndim <- integer(1)
  clen <- integer(1)

  # Load metadata dimensions + clen
  meta <- get_array_metadata_chars(filename, max_dims)

  actual_dims <- meta$dims
  clen <- meta$clen
  total <- prod(actual_dims)
  cat("actual_dims:", actual_dims, "clen:", clen, "\n")

  ascii_arr <- integer(clen * total)

  res <- .Fortran("deserialize_char_flat_r",
    ascii_arr = ascii_arr,
    arr_size = integer(clen * total),  # size is ignored in the Fortran code
    dims_out = actual_dims,
    ndim_out = ndim,
    clen_out = clen,
    filename_ascii = ascii,
    fn_len = length(ascii),
    ndim_actual = as.integer(length(actual_dims))
  )

  # translate ASCII back to char
  mat <- matrix(res$ascii_arr, nrow = clen)
  chars <- apply(mat, 2, function(col) rawToChar(as.raw(col[col > 0])))

  array(chars, dim = res$dims_out[1:res$ndim_out])
}


# BASE R arrays are column-major just like fortran, so no serialization is needed for the array structure.
# Array can simply be passed with with as.integer()
serialize_int_array <- function(arr, filename) {
  flat <- as.integer(arr)
  dims <- as.integer(dim(arr))
  ndim <- as.integer(length(dims))
  ascii <- utf8ToInt(filename)

  .Fortran("serialize_int_flat_r",
           arr = flat,
           dims = dims,
           ndim = ndim,
           filename_ascii = as.integer(ascii),
           fn_len = as.integer(length(ascii)))
}

# BASE R arrays are column-major just like fortran, so no serialization is needed for the array structure.
# Array can simply be passed with with as.double() to pass it in a flat format.
serialize_real_array <- function(arr, filename) {
  flat <- as.double(arr)

  dims <- if (is.null(dim(arr))) {
    as.integer(length(arr))  # 1D-Vector
  } else {
    as.integer(dim(arr))
  }

  ndim <- as.integer(length(dims))
  ascii <- utf8ToInt(filename)

  .Fortran("serialize_real_flat_r",
           arr = flat,
           dims = dims,
           ndim = ndim,
           filename_ascii = as.integer(ascii),
           fn_len = as.integer(length(ascii)))
}

# Serializes a character array to a file, encoding it as an integer matrix
# Each character is converted to its ASCII integer representation
# The matrix is then serialized with Fortran
serialize_char_array <- function(arr, filename) {
  stopifnot(is.character(arr))
  arr <- as.array(arr)
  dims <- dim(arr)
  if (is.null(dims)) dims <- length(arr)
  clen <- max(nchar(arr, type = "chars"))

  # encode to integer matrix
  # Chars can not be passed via .Fortran directly
  mat <- matrix(0L, nrow = clen, ncol = length(arr))
  for (i in seq_along(arr)) {
    chars <- utf8ToInt(substr(arr[i], 1, clen))
    mat[seq_along(chars), i] <- chars
  }

  invisible(.Fortran("serialize_char_flat_r",
    ascii_arr = as.integer(mat),
    dims = as.integer(dims),
    ndim = as.integer(length(dims)),
    clen = as.integer(clen),
    filename_ascii = utf8ToInt(filename),
    fn_len = nchar(filename)
  ))
}

# --- Testing for all functions ---
test_array_wrappers <- function(tmpdir = tempdir()) {
  cat("Start Wrapper-Tests...\n")
  fn <- function(name) file.path(tmpdir, name)

  cat("Integer tests...\n")
  cat("Dim 1: ")
  arr1 <- 1:10
  # Debug: check storage mode and length
  stopifnot(is.integer(arr1))
  stopifnot(length(arr1) == 10)
  serialize_int_array(arr1, fn("int1d.bin"))
  stopifnot(all(deserialize_int_array(fn("int1d.bin")) == arr1))
  cat("okay\n")

  cat("Dim 2: ")
  arr2 <- matrix(1:12, nrow=3, ncol=4)
  stopifnot(is.integer(arr2))
  stopifnot(length(arr2) == 12)
  serialize_int_array(arr2, fn("int2d.bin"))
  stopifnot(all(deserialize_int_array(fn("int2d.bin")) == arr2))
  cat("okay\n")

  cat("Dim 3: ")
  arr3 <- array(1:24, dim = c(2,3,4))
  serialize_int_array(arr3, fn("int3d.bin"))
  stopifnot(all(deserialize_int_array(fn("int3d.bin")) == arr3))
  cat("okay\n")

  cat("Dim 4: ")
  arr4 <- array(1:48, dim = c(2,3,4,2))
  serialize_int_array(arr4, fn("int4d.bin"))
  stopifnot(all(deserialize_int_array(fn("int4d.bin")) == arr4))
  cat("okay\n")

  cat("Dim 5: ")
  arr5 <- array(1:96, dim = c(2,3,4,2,2))
  serialize_int_array(arr5, fn("int5d.bin"))
  stopifnot(all(deserialize_int_array(fn("int5d.bin")) == arr5))
  cat("okay\n")

  # REAL
  cat("Real tests...\n")
  arr1r <- as.numeric(1:10) * 0.5
  serialize_real_array(arr1r, fn("real1d.bin"))
  cat("serialized real 1D arr\n")
  stopifnot(all(deserialize_real_array(fn("real1d.bin")) == arr1r))
  cat("okay\n")

  arr2r <- matrix(runif(12), nrow=3, ncol=4)
  serialize_real_array(arr2r, fn("real2d.bin"))
  stopifnot(all(deserialize_real_array(fn("real2d.bin")) == arr2r))
  cat("okay\n")

  arr3r <- array(runif(24), dim = c(2,3,4))
  serialize_real_array(arr3r, fn("real3d.bin"))
  stopifnot(all(deserialize_real_array(fn("real3d.bin")) == arr3r))
  cat("okay\n")

  arr4r <- array(runif(48), dim = c(2,3,4,2))
  serialize_real_array(arr4r, fn("real4d.bin"))
  stopifnot(all(deserialize_real_array(fn("real4d.bin")) == arr4r))
  cat("okay\n")

  arr5r <- array(runif(96), dim = c(2,3,4,2,2))
  serialize_real_array(arr5r, fn("real5d.bin"))
  stopifnot(all(deserialize_real_array(fn("real5d.bin")) == arr5r))
  cat("okay\n")
  cat("Real tests successful!\n")

  # CHAR
  cat("Character tests...\n")
  clen <- 8
  arr1c <- sprintf("%0*d", clen, 1:10)
  serialize_char_array(arr1c, fn("char1d.bin"))
  stopifnot(all(deserialize_char_array(fn("char1d.bin")) == arr1c))

  arr2c <- matrix(sprintf("%0*d", clen, 1:12), nrow=3, ncol=4)
  serialize_char_array(arr2c, fn("char2d.bin"))
  stopifnot(all(deserialize_char_array(fn("char2d.bin")) == arr2c))

  arr3c <- array(sprintf("%0*d", clen, 1:24), dim = c(2,3,4))
  serialize_char_array(arr3c, fn("char3d.bin"))
  stopifnot(all(deserialize_char_array(fn("char3d.bin")) == arr3c))

  arr4c <- array(sprintf("%0*d", clen, 1:48), dim = c(2,3,4,2))
  serialize_char_array(arr4c, fn("char4d.bin"))
  stopifnot(all(deserialize_char_array(fn("char4d.bin")) == arr4c))

  arr5c <- array(sprintf("%0*d", clen, 1:96), dim = c(2,3,4,2,2))
  serialize_char_array(arr5c, fn("char5d.bin"))
  stopifnot(all(deserialize_char_array(fn("char5d.bin")) == arr5c))

  # 1D-Array with different char lengths
  arr_ascii1 <- c("A", "G1", "GENE003", "BRCA1", "XYZ", "", "12345678", "SEQ")
  serialize_char_array(arr_ascii1, fn("char_ascii1d.bin"))
  stopifnot(all(deserialize_char_array(fn("char_ascii1d.bin")) == arr_ascii1))

  # 2D-Matrix with different length
  arr_ascii2 <- matrix(c("GENE1", "GENE22", "GENE333", "", "ID", "SEQ9999"), nrow = 2, byrow = TRUE)
  serialize_char_array(arr_ascii2, fn("char_ascii2d.bin"))
  stopifnot(all(deserialize_char_array(fn("char_ascii2d.bin")) == arr_ascii2))

  # 3D-Array with realistic data
  arr_ascii3 <- array(c("TP53", "BRCA1", "MT-ATP6", "CYTB", "ND1", "", "NRAS", "EGFR"), dim = c(2, 2, 2))
  serialize_char_array(arr_ascii3, fn("char_ascii3d.bin"))
  stopifnot(all(deserialize_char_array(fn("char_ascii3d.bin")) == arr_ascii3))

  cat("All ASCII string variation tests passed!\n")
}

# start tests
test_array_wrappers()