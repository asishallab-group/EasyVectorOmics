dyn.load("build/libtrees.so") # Load the shared library

#' Build BST index (1D)
build_bst_index <- function(x) {
  n <- as.integer(length(x))
  ix <- integer(n)
  stack_left <- integer(n)
  stack_right <- integer(n)

  cat("n =", n, "\n")
  cat("length(x) =", length(x), "\n")
  cat("length(ix) =", length(ix), "\n")
  cat("length(stack_left) =", length(stack_left), "\n")
  cat("length(stack_right) =", length(stack_right), "\n")

  res <- .Fortran("build_bst_index_r",
                  as.double(x),
                  as.integer(n),
                  ix,
                  stack_left,
                  stack_right)

  res[[3]]
}

#' BST range query
bst_range_query <- function(x, ix, lo, hi) {
  n <- as.integer(length(x))
  out_ix <- integer(n)
  out_n <- integer(n)
  res <- .Fortran("bst_range_query_r",
                  as.double(x), 
                  as.integer(ix), 
                  n, 
                  as.double(lo), 
                  as.double(hi), 
                  as.integer(out_ix), 
                  as.integer(out_n))
  count <- res[[7]][1]
  list(indices = res[[6]][seq_len(count)], count = count)
}

#' Build KD-Tree index (multidimensional)
build_kd_index <- function(X, dim_order) {
  d <- as.integer(nrow(X))
  n <- as.integer(ncol(X))
  kd_ix <- integer(n)
  work <- integer(n)
  subarray <- double(n)
  perm <- integer(n)
  stack_left <- integer(64)
  stack_right <- integer(64)
  res <- .Fortran("build_kd_index_r",
                  X = as.double(X), 
                  d = as.integer(d), 
                  n = as.integer(n), 
                  kd_ix = as.integer(kd_ix), 
                  dim_order = as.integer(dim_order), 
                  work = as.integer(work), 
                  subarray = as.double(subarray), 
                  perm = as.integer(perm), 
                  stack_left = as.integer(stack_left), 
                  stack_right = as.integer(stack_right))
  res[[4]]
}

#' Build Spherical KD-Tree index
build_spherical_kd <- function(V, dim_order) {
  d <- as.integer(nrow(V))
  n <- as.integer(ncol(V))
  sphere_ix <- integer(n)
  work <- integer(n)
  subarray <- double(n)
  perm <- integer(n)
  stack_left <- integer(64)
  stack_right <- integer(64)
  res <- .Fortran("build_spherical_kd_r",
                  as.double(V), d, n, sphere_ix, as.integer(dim_order), work, subarray, perm, stack_left, stack_right)
  res[[4]]
}

# --- Simple tests ---

cat("Testing BST index...\n")
set.seed(42)
x <- runif(10)
ix <- build_bst_index(x)
cat("Original x:", x, "\n")
cat("Sorted x[ix]:", x[ix], "\n")

cat("Testing BST range query...\n")
range_res <- bst_range_query(x, ix, 0.2, 0.8)
cat("Indices in range [0.2, 0.8]:", range_res$indices, "\n")
cat("Count:", range_res$count, "\n")

cat("Testing KD-Tree index...\n")
X <- matrix(runif(20), nrow=2)
dim_order <- c(1L, 2L)
kd_ix <- build_kd_index(X, dim_order)
cat("KD-Tree index:", kd_ix, "\n")

cat("Testing Spherical KD-Tree index...\n")
V <- matrix(rnorm(12), nrow=3)
V <- apply(V, 2, function(col) col / sqrt(sum(col^2)))
dim_order <- c(1L, 2L, 3L)
sphere_ix <- build_spherical_kd(V, dim_order)
cat("Spherical KD-Tree index:", sphere_ix, "\n")