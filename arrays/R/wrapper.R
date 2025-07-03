dyn.load("arrays/build/array_r.so")

#' Serialize an integer array to file using Fortran routine
serialize_array <- function(arr, filename) {
    .Fortran("serialize_r",
            arr = arr,
            filename = as.character(filename))
    invisible(NULL)
}

#' Deserialize a real array from file using Fortran routine
deserialize_real_flat <- function(filename) {
    nmax <- 1e6 # adjust as needed
    flat <- double(nmax)
    dims <- integer(1)
    res <- .Fortran("deserialize_real_flat_r",
                    flat = flat,
                    dims = dims,
                    filename = as.character(filename))
    # Remove unused elements
    d <- res$dims[res$dims != 0]
    flat <- res$flat[seq_len(prod(d))]
    list(flat = flat, dims = d)
}

#' Deserialize an integer array from file using Fortran routine
deserialize_int_flat <- function(filename) {
    nmax <- 1e6 # adjust as needed
    flat <- integer(nmax)
    dims <- integer(1)
    res <- .Fortran("deserialize_int_flat_r",
                    flat = flat,
                    dims = dims,
                    filename = as.character(filename))
    d <- res$dims[res$dims != 0]
    flat <- res$flat[seq_len(prod(d))]
    list(flat = flat, dims = d)
}

#' Deserialize a character array from file using Fortran routine
deserialize_char_flat <- function(filename) {
  nmax <- 1e6 # adjust as needed
  clen <- integer(1)
  flat <- character(nmax)
  dims <- integer(1)
  res <- .Fortran("deserialize_char_flat_r",
                  flat = flat,
                  dims = dims,
                  clen = clen,
                  filename = as.character(filename))
  d <- res$dims[res$dims != 0]
  flat <- res$flat[seq_len(prod(d))]
  list(flat = flat, dims = d, clen = res$clen)
}

#' Test serialization and deserialization
test_array_serialization <- function() {
  # Integer test
  arr <- matrix(1L:12L, nrow = 3, ncol = 4)
  serialize_array(as.integer(arr), "test_int.bin")
  int_res <- deserialize_int_flat("test_int.bin")
  stopifnot(all(int_res$flat == as.integer(arr)))
  stopifnot(all(int_res$dims == dim(arr)))

  # Real test
  arr_real <- matrix(runif(12), nrow = 3, ncol = 4)
  # You would need a serialize_real_r wrapper for real arrays

  # Character test
  arr_char <- matrix(rep("A", 12), nrow = 3, ncol = 4)
  # You would need a serialize_char_r wrapper for character arrays

  cat("Integer serialization/deserialization test passed.\n")
}

test_array_serialization()