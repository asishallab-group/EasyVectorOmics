dyn.load("arrays/build/arrays.so")

deserialize_int_1d <- function(filename, n) {
  arr <- integer(n)
  ascii <- utf8ToInt(filename)
  res <- .Fortran("deserialize_int_1d_r",
                  arr = arr,
                  filename_ascii = as.integer(ascii),
                  fn_len = as.integer(length(ascii)),
                  PACKAGE = "arrays")
  res$arr
}

deserialize_int_2d <- function(filename, dim1, dim2) {
  arr <- integer(dim1 * dim2)
  res <- .Fortran("deserialize_int_2d_r",
                  arr = arr,
                  filename = as.character(filename),
                  PACKAGE = "arrays")
  array(res$arr, dim = c(dim1, dim2))
}

deserialize_int_3d <- function(filename, dim1, dim2, dim3) {
  arr <- integer(dim1 * dim2 * dim3)
  res <- .Fortran("deserialize_int_3d_r",
                  arr = arr,
                  filename = as.character(filename),
                  PACKAGE = "arrays")
  array(res$arr, dim = c(dim1, dim2, dim3))
}

deserialize_int_4d <- function(filename, dim1, dim2, dim3, dim4) {
  arr <- integer(dim1 * dim2 * dim3 * dim4)
  res <- .Fortran("deserialize_int_4d_r",
                  arr = arr,
                  filename = as.character(filename),
                  PACKAGE = "arrays")
  array(res$arr, dim = c(dim1, dim2, dim3, dim4))
}

deserialize_int_5d <- function(filename, dim1, dim2, dim3, dim4, dim5) {
  arr <- integer(dim1 * dim2 * dim3 * dim4 * dim5)
  res <- .Fortran("deserialize_int_5d_r",
                  arr = arr,
                  filename = as.character(filename),
                  PACKAGE = "arrays")
  array(res$arr, dim = c(dim1, dim2, dim3, dim4, dim5))
}

deserialize_real_1d <- function(filename, n) {
  arr <- double(n)
  res <- .Fortran("deserialize_real_1d_r",
                  arr = arr,
                  filename = as.character(filename),
                  PACKAGE = "arrays")
  res$arr
}

deserialize_real_2d <- function(filename, dim1, dim2) {
  arr <- double(dim1 * dim2)
  res <- .Fortran("deserialize_real_2d_r",
                  arr = arr,
                  filename = as.character(filename),
                  PACKAGE = "arrays")
  array(res$arr, dim = c(dim1, dim2))
}

deserialize_real_3d <- function(filename, dim1, dim2, dim3) {
  arr <- double(dim1 * dim2 * dim3)
  res <- .Fortran("deserialize_real_3d_r",
                  arr = arr,
                  filename = as.character(filename),
                  PACKAGE = "arrays")
  array(res$arr, dim = c(dim1, dim2, dim3))
}

deserialize_real_4d <- function(filename, dim1, dim2, dim3, dim4) {
  arr <- double(dim1 * dim2 * dim3 * dim4)
  res <- .Fortran("deserialize_real_4d_r",
                  arr = arr,
                  filename = as.character(filename),
                  PACKAGE = "arrays")
  array(res$arr, dim = c(dim1, dim2, dim3, dim4))
}

deserialize_real_5d <- function(filename, dim1, dim2, dim3, dim4, dim5) {
  arr <- double(dim1 * dim2 * dim3 * dim4 * dim5)
  res <- .Fortran("deserialize_real_5d_r",
                  arr = arr,
                  filename = as.character(filename),
                  PACKAGE = "arrays")
  array(res$arr, dim = c(dim1, dim2, dim3, dim4, dim5))
}

# Für char: Rückgabe als character-Vektor/Array, Annahme: maximale Länge clen muss bekannt sein
deserialize_char_1d <- function(filename, n, clen) {
  arr <- rep(strrep(" ", clen), n)
  res <- .Fortran("deserialize_char_1d_r",
                  arr = arr,
                  filename = as.character(filename),
                  PACKAGE = "arrays")
  trimws(res$arr)
}

deserialize_char_2d <- function(filename, dim1, dim2, clen) {
  arr <- rep(strrep(" ", clen), dim1 * dim2)
  res <- .Fortran("deserialize_char_2d_r",
                  arr = arr,
                  filename = as.character(filename),
                  PACKAGE = "arrays")
  matrix(trimws(res$arr), nrow = dim1, ncol = dim2)
}

deserialize_char_3d <- function(filename, dim1, dim2, dim3, clen) {
  arr <- rep(strrep(" ", clen), dim1 * dim2 * dim3)
  res <- .Fortran("deserialize_char_3d_r",
                  arr = arr,
                  filename = as.character(filename),
                  PACKAGE = "arrays")
  array(trimws(res$arr), dim = c(dim1, dim2, dim3))
}

deserialize_char_4d <- function(filename, dim1, dim2, dim3, dim4, clen) {
  arr <- rep(strrep(" ", clen), dim1 * dim2 * dim3 * dim4)
  res <- .Fortran("deserialize_char_4d_r",
                  arr = arr,
                  filename = as.character(filename),
                  PACKAGE = "arrays")
  array(trimws(res$arr), dim = c(dim1, dim2, dim3, dim4))
}

deserialize_char_5d <- function(filename, dim1, dim2, dim3, dim4, dim5, clen) {
  arr <- rep(strrep(" ", clen), dim1 * dim2 * dim3 * dim4 * dim5)
  res <- .Fortran("deserialize_char_5d_r",
                  arr = arr,
                  filename = as.character(filename),
                  PACKAGE = "arrays")
  array(trimws(res$arr), dim = c(dim1, dim2, dim3, dim4, dim5))
}

# --- Serialize Funktionen ---

serialize_int_1d <- function(arr, filename) {
    ascii <- utf8ToInt(filename)
    .Fortran("serialize_int_1d_r",
           arr = as.integer(arr),
           n1 = as.integer(length(arr)),
           filename_ascii = as.integer(ascii),
           fn_len = as.integer(length(ascii)),
           PACKAGE = "arrays")
}

serialize_int_2d <- function(arr, filename) {
  d <- dim(arr)
  .Fortran("serialize_int_2d_r",
           arr = as.integer(arr),
           n1 = as.integer(d[1]),
           n2 = as.integer(d[2]),
           filename = as.character(filename),
           PACKAGE = "arrays")
}

serialize_int_3d <- function(arr, filename) {
  d <- dim(arr)
  .Fortran("serialize_int_3d_r",
           arr = as.integer(arr),
           n1 = as.integer(d[1]),
           n2 = as.integer(d[2]),
           n3 = as.integer(d[3]),
           filename = as.character(filename),
           PACKAGE = "arrays")
}

serialize_int_4d <- function(arr, filename) {
  d <- dim(arr)
  .Fortran("serialize_int_4d_r",
           arr = as.integer(arr),
           n1 = as.integer(d[1]),
           n2 = as.integer(d[2]),
           n3 = as.integer(d[3]),
           n4 = as.integer(d[4]),
           filename = as.character(filename),
           PACKAGE = "arrays")
}

serialize_int_5d <- function(arr, filename) {
  d <- dim(arr)
  .Fortran("serialize_int_5d_r",
           arr = as.integer(arr),
           n1 = as.integer(d[1]),
           n2 = as.integer(d[2]),
           n3 = as.integer(d[3]),
           n4 = as.integer(d[4]),
           n5 = as.integer(d[5]),
           filename = as.character(filename),
           PACKAGE = "arrays")
}

serialize_real_1d <- function(arr, filename) {
  .Fortran("serialize_real_1d_r",
           arr = as.double(arr),
           n1 = as.integer(length(arr)),
           filename = as.character(filename),
           PACKAGE = "arrays")
}

serialize_real_2d <- function(arr, filename) {
  d <- dim(arr)
  .Fortran("serialize_real_2d_r",
           arr = as.double(arr),
           n1 = as.integer(d[1]),
           n2 = as.integer(d[2]),
           filename = as.character(filename),
           PACKAGE = "arrays")
}

serialize_real_3d <- function(arr, filename) {
  d <- dim(arr)
  .Fortran("serialize_real_3d_r",
           arr = as.double(arr),
           n1 = as.integer(d[1]),
           n2 = as.integer(d[2]),
           n3 = as.integer(d[3]),
           filename = as.character(filename),
           PACKAGE = "arrays")
}

serialize_real_4d <- function(arr, filename) {
  d <- dim(arr)
  .Fortran("serialize_real_4d_r",
           arr = as.double(arr),
           n1 = as.integer(d[1]),
           n2 = as.integer(d[2]),
           n3 = as.integer(d[3]),
           n4 = as.integer(d[4]),
           filename = as.character(filename),
           PACKAGE = "arrays")
}

serialize_real_5d <- function(arr, filename) {
  d <- dim(arr)
  .Fortran("serialize_real_5d_r",
           arr = as.double(arr),
           n1 = as.integer(d[1]),
           n2 = as.integer(d[2]),
           n3 = as.integer(d[3]),
           n4 = as.integer(d[4]),
           n5 = as.integer(d[5]),
           filename = as.character(filename),
           PACKAGE = "arrays")
}

serialize_char_1d <- function(arr, filename) {
  .Fortran("serialize_char_1d_r",
           arr = as.character(arr),
           n1 = as.integer(length(arr)),
           filename = as.character(filename),
           PACKAGE = "arrays")
}

serialize_char_2d <- function(arr, filename) {
  d <- dim(arr)
  .Fortran("serialize_char_2d_r",
           arr = as.character(arr),
           n1 = as.integer(d[1]),
           n2 = as.integer(d[2]),
           filename = as.character(filename),
           PACKAGE = "arrays")
}

serialize_char_3d <- function(arr, filename) {
  d <- dim(arr)
  .Fortran("serialize_char_3d_r",
           arr = as.character(arr),
           n1 = as.integer(d[1]),
           n2 = as.integer(d[2]),
           n3 = as.integer(d[3]),
           filename = as.character(filename),
           PACKAGE = "arrays")
}

serialize_char_4d <- function(arr, filename) {
  d <- dim(arr)
  .Fortran("serialize_char_4d_r",
           arr = as.character(arr),
           n1 = as.integer(d[1]),
           n2 = as.integer(d[2]),
           n3 = as.integer(d[3]),
           n4 = as.integer(d[4]),
           filename = as.character(filename),
           PACKAGE = "arrays")
}

serialize_char_5d <- function(arr, filename) {
  d <- dim(arr)
  .Fortran("serialize_char_5d_r",
           arr = as.character(arr),
           n1 = as.integer(d[1]),
           n2 = as.integer(d[2]),
           n3 = as.integer(d[3]),
           n4 = as.integer(d[4]),
           n5 = as.integer(d[5]),
           filename = as.character(filename),
           PACKAGE = "arrays")
}

# --- Testprogramm für alle Funktionen ---
test_array_wrappers <- function(tmpdir = tempdir()) {
  cat("Start Wrapper-Tests...\n")
  fn <- function(name) file.path(tmpdir, name)

  cat("Integer tests...\n")
  cat("Dim 1\n")
  arr1 <- 1:10
  # Debug: check storage mode and length
  stopifnot(is.integer(arr1))
  stopifnot(length(arr1) == 10)
  serialize_int_1d(arr1, fn("int1d.bin"))
  cat("serialized\n")
  cat("Serialized array: " , arr1)
  cat("\n")
  cat(deserialize_int_1d(fn("int1d.bin"), 10))
  stopifnot(all(deserialize_int_1d(fn("int1d.bin"), 10) == arr1))

  cat("Dim 2\n")
  arr2 <- matrix(1:12, nrow=3, ncol=4)
  stopifnot(is.integer(arr2))
  stopifnot(length(arr2) == 12)
  serialize_int_2d(arr2, fn("int2d.bin"))
  stopifnot(all(deserialize_int_2d(fn("int2d.bin"), 3, 4) == arr2))

  cat("Dim 3\n")
  arr3 <- array(1:24, dim = c(2,3,4))
  serialize_int_3d(arr3, fn("int3d.bin"))
  stopifnot(all(deserialize_int_3d(fn("int3d.bin"), 2,3,4) == arr3))

  cat("Dim 4\n")
  arr4 <- array(1:48, dim = c(2,3,4,2))
  serialize_int_4d(arr4, fn("int4d.bin"))
  stopifnot(all(deserialize_int_4d(fn("int4d.bin"), 2,3,4,2) == arr4))

  cat("Dim 5\n")
  arr5 <- array(1:96, dim = c(2,3,4,2,2))
  serialize_int_5d(arr5, fn("int5d.bin"))
  stopifnot(all(deserialize_int_5d(fn("int5d.bin"), 2,3,4,2,2) == arr5))

  # REAL
  arr1r <- as.numeric(1:10) * 0.5
  serialize_real_1d(arr1r, fn("real1d.bin"))
  stopifnot(all.equal(deserialize_real_1d(fn("real1d.bin"), 10), arr1r))

  arr2r <- matrix(runif(12), nrow=3, ncol=4)
  serialize_real_2d(arr2r, fn("real2d.bin"))
  stopifnot(all.equal(deserialize_real_2d(fn("real2d.bin"), 3, 4), arr2r))

  arr3r <- array(runif(24), dim = c(2,3,4))
  serialize_real_3d(arr3r, fn("real3d.bin"))
  stopifnot(all.equal(deserialize_real_3d(fn("real3d.bin"), 2,3,4), arr3r))

  arr4r <- array(runif(48), dim = c(2,3,4,2))
  serialize_real_4d(arr4r, fn("real4d.bin"))
  stopifnot(all.equal(deserialize_real_4d(fn("real4d.bin"), 2,3,4,2), arr4r))

  arr5r <- array(runif(96), dim = c(2,3,4,2,2))
  serialize_real_5d(arr5r, fn("real5d.bin"))
  stopifnot(all.equal(deserialize_real_5d(fn("real5d.bin"), 2,3,4,2,2), arr5r))

  # CHAR
  clen <- 8
  arr1c <- sprintf("%0*d", clen, 1:10)
  serialize_char_1d(arr1c, fn("char1d.bin"))
  stopifnot(all(deserialize_char_1d(fn("char1d.bin"), 10, clen) == arr1c))

  arr2c <- matrix(sprintf("%0*d", clen, 1:12), nrow=3, ncol=4)
  serialize_char_2d(arr2c, fn("char2d.bin"))
  stopifnot(all(deserialize_char_2d(fn("char2d.bin"), 3, 4, clen) == arr2c))

  arr3c <- array(sprintf("%0*d", clen, 1:24), dim = c(2,3,4))
  serialize_char_3d(arr3c, fn("char3d.bin"))
  stopifnot(all(deserialize_char_3d(fn("char3d.bin"), 2,3,4, clen) == arr3c))

  arr4c <- array(sprintf("%0*d", clen, 1:48), dim = c(2,3,4,2))
  serialize_char_4d(arr4c, fn("char4d.bin"))
  stopifnot(all(deserialize_char_4d(fn("char4d.bin"), 2,3,4,2, clen) == arr4c))

  arr5c <- array(sprintf("%0*d", clen, 1:96), dim = c(2,3,4,2,2))
  serialize_char_5d(arr5c, fn("char5d.bin"))
  stopifnot(all(deserialize_char_5d(fn("char5d.bin"), 2,3,4,2,2, clen) == arr5c))

  cat("Alle Wrapper-Tests erfolgreich!\n")
}

# Am Ende der Datei automatisch testen (optional, auskommentieren falls nicht gewünscht)
test_array_wrappers()