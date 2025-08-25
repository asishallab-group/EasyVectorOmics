dyn.load("Trees/build/libtrees.so") # Load the shared library

#' Build BST index (1D)
build_bst_index <- function(x) {
  if (!is.numeric(x)) {
    stop("Input x must be numeric")
  }
  
  n <- as.integer(length(x))
  ix <- integer(n)
  stack_left <- integer(n)
  stack_right <- integer(n)

  res <- .Fortran("build_bst_index_r",
                  x = as.double(x),
                  n = n,
                  ix = ix,
                  stack_left = stack_left,
                  stack_right = stack_right)

  res$ix
}

#' BST range query
bst_range_query <- function(x, ix, lo, hi) {
  if (length(x) != length(ix)) {
    stop("Length of x and ix must match")
  }
  
  n <- as.integer(length(x))
  out_ix <- integer(n)
  out_n <- integer(1)
  
  res <- .Fortran("bst_range_query_r",
                  x = as.double(x), 
                  ix = as.integer(ix), 
                  n = n, 
                  lo = as.double(lo), 
                  hi = as.double(hi), 
                  out_ix = out_ix, 
                  out_n = out_n)
  
  list(indices = res$out_ix[1:res$out_n], count = res$out_n)
}

#' Get sorted value from BST index
get_sorted_value <- function(x, ix, position) {
  if (position < 1 || position > length(ix)) {
    stop("Position must be between 1 and length(ix)")
  }
  x[ix[position]]
}

#' Build KD-Tree index (multidimensional)
build_kd_index <- function(X, dim_order = NULL) {
  if (!is.matrix(X)) {
    stop("Input X must be a matrix")
  }
  
  d <- as.integer(nrow(X))
  n <- as.integer(ncol(X))
  
  if (is.null(dim_order)) {
    dim_order <- 1:d  # Default: use dimensions in order
  }
  
  if (length(dim_order) != d) {
    stop("dim_order length must match number of dimensions")
  }
  
  kd_ix <- integer(n)
  work <- integer(n)
  subarray <- double(n)
  perm <- integer(n)
  stack_left <- integer(n)
  stack_right <- integer(n)
  
  res <- .Fortran("build_kd_index_r",
                  X = as.double(X), 
                  d = d, 
                  n = n, 
                  kd_ix = kd_ix, 
                  dim_order = as.integer(dim_order), 
                  work = work, 
                  subarray = subarray, 
                  perm = perm, 
                  stack_left = stack_left, 
                  stack_right = stack_right)
  
  res$kd_ix
}

#' Build Spherical KD-Tree index
build_spherical_kd <- function(V, dim_order = NULL) {
  if (!is.matrix(V)) {
    stop("Input V must be a matrix")
  }
  
  d <- as.integer(nrow(V))
  n <- as.integer(ncol(V))
  
  if (is.null(dim_order)) {
    dim_order <- 1:d  # Default: use dimensions in order
  }
  
  if (length(dim_order) != d) {
    stop("dim_order length must match number of dimensions")
  }
  
  sphere_ix <- integer(n)
  work <- integer(n)
  subarray <- double(n)
  perm <- integer(n)
  stack_left <- integer(n)
  stack_right <- integer(n)
  
  res <- .Fortran("build_spherical_kd_r",
                  V = as.double(V), 
                  d = d, 
                  n = n, 
                  sphere_ix = sphere_ix, 
                  dim_order = as.integer(dim_order), 
                  work = work, 
                  subarray = subarray, 
                  perm = perm, 
                  stack_left = stack_left, 
                  stack_right = stack_right)
  
  res$sphere_ix
}

#' Get point from KD-Tree index
get_kd_point <- function(X, kd_ix, position) {
  if (position < 1 || position > length(kd_ix)) {
    stop("Position must be between 1 and length(kd_ix)")
  }
  if (ncol(X) < max(kd_ix)) {
    stop("KD index contains invalid indices for matrix X")
  }
  X[, kd_ix[position]]
}