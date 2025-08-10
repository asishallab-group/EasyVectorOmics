# Comprehensive R test suite for tissue versatility (mirrors Fortran unit tests)
# Load the shared library
dyn.load("build/libtensor-omics.so")

# Helper: call Fortran wrapper and return result as list
tv_call <- function(expr, select_vec, select_axes) {
  n_axes <- nrow(expr)
  n_vectors <- ncol(expr)
  n_selected_vectors <- sum(select_vec)
  n_selected_axes <- sum(select_axes)
  tissue_versatilities <- rep(0.0, n_selected_vectors)
  tissue_angles_deg <- rep(0.0, n_selected_vectors)
  .Fortran("compute_tissue_versatility_r",
           n_axes = as.integer(n_axes),
           n_vectors = as.integer(n_vectors),
           expression_vectors = as.double(expr),
           exp_vecs_selection_index = as.integer(select_vec),
           n_selected_vectors = as.integer(n_selected_vectors),
           axes_selection = as.integer(select_axes),
           n_selected_axes = as.integer(n_selected_axes),
           tissue_versatilities = as.double(tissue_versatilities),
           tissue_angles_deg = as.double(tissue_angles_deg))
}

# 1. Uniform expression (should yield TV=0)
test_uniform_expression <- function() {
  expr <- matrix(2, nrow=3, ncol=1)
  res <- tv_call(expr, c(1L), c(1L,1L,1L))
  stopifnot(abs(res$tissue_versatilities[1]) < 1e-12)
  stopifnot(abs(res$tissue_angles_deg[1]) < 1e-12)
  cat("test_uniform_expression passed\n")
}

# 2. Single axis expression (should yield TV=1)
test_single_axis_expression <- function() {
  expr <- matrix(c(0,0,5), nrow=3, ncol=1)
  res <- tv_call(expr, c(1L), c(1L,1L,1L))
  stopifnot(abs(res$tissue_versatilities[1] - 1) < 1e-12)
  stopifnot(res$tissue_angles_deg[1] > 0)
  cat("test_single_axis_expression passed\n")
}

# 3. Null vector (should yield TV=1, angle=90)
test_null_vector <- function() {
  expr <- matrix(0, nrow=3, ncol=1)
  res <- tv_call(expr, c(1L), c(1L,1L,1L))
  stopifnot(abs(res$tissue_versatilities[1] - 1) < 1e-12)
  stopifnot(abs(res$tissue_angles_deg[1] - 90) < 1e-12)
  cat("test_null_vector passed\n")
}

# 4. Partial axis selection (subspace)
test_partial_axis_selection <- function() {
  expr <- matrix(c(1,2,3), nrow=3, ncol=1)
  res <- tv_call(expr, c(1L), c(1L,0L,1L))
  stopifnot(res$tissue_versatilities[1] >= 0 && res$tissue_versatilities[1] <= 1)
  stopifnot(res$tissue_angles_deg[1] >= 0 && res$tissue_angles_deg[1] <= 90)
  cat("test_partial_axis_selection passed\n")
}

# 5. Mixed vectors (uniform, single axis, null)
test_mixed_vectors <- function() {
  expr <- matrix(c(1,1,1, 0,0,2, 0,0,0), nrow=3, ncol=3)
  res <- tv_call(expr, c(1L,1L,1L), c(1L,1L,1L))
  stopifnot(abs(res$tissue_versatilities[1]) < 1e-12)
  stopifnot(abs(res$tissue_versatilities[2] - 1) < 1e-12)
  stopifnot(abs(res$tissue_versatilities[3] - 1) < 1e-12)
  stopifnot(abs(res$tissue_angles_deg[1]) < 1e-12)
  stopifnot(res$tissue_angles_deg[2] > 0)
  stopifnot(abs(res$tissue_angles_deg[3] - 90) < 1e-12)
  cat("test_mixed_vectors passed\n")
}

# 6. Angle output in degrees for a known case (should be 45)
test_angle_degrees <- function() {
  expr <- matrix(c(1,0), nrow=2, ncol=1)
  res <- tv_call(expr, c(1L), c(1L,1L))
  stopifnot(abs(res$tissue_angles_deg[1] - 45) < 1e-12)
  cat("test_angle_degrees passed\n")
}

# 7. Multiple vectors selection
test_multiple_vectors_selection <- function() {
  expr <- matrix(c(1,1, 0,2, 0,0), nrow=2, ncol=3)
  res <- tv_call(expr, c(1L,0L,1L), c(1L,1L))
  stopifnot(abs(res$tissue_versatilities[1]) < 1e-12)
  stopifnot(abs(res$tissue_versatilities[2] - 1) < 1e-12)
  stopifnot(abs(res$tissue_angles_deg[1]) < 1e-12)
  stopifnot(abs(res$tissue_angles_deg[2] - 90) < 1e-12)
  cat("test_multiple_vectors_selection passed\n")
}

# 8. High-dimensional vectors (4D, 5D)
test_high_dimensional_vectors <- function() {
  expr4 <- matrix(1, nrow=4, ncol=1)
  expr5 <- matrix(2, nrow=5, ncol=1)
  res4 <- tv_call(expr4, c(1L), rep(1L,4))
  res5 <- tv_call(expr5, c(1L), rep(1L,5))
  stopifnot(abs(res4$tissue_versatilities[1]) < 1e-12)
  stopifnot(abs(res4$tissue_angles_deg[1]) < 1e-12)
  stopifnot(abs(res5$tissue_versatilities[1]) < 1e-12)
  stopifnot(abs(res5$tissue_angles_deg[1]) < 1e-12)
  cat("test_high_dimensional_vectors passed\n")
}

# 9. Randomized vectors and axes
test_randomized_vectors_axes <- function() {
  set.seed(42)
  n_axes <- 5
  n_vecs <- 4
  expr <- matrix(runif(n_axes * n_vecs), nrow=n_axes, ncol=n_vecs)
  res <- tv_call(expr, rep(1L,n_vecs), c(1L,0L,1L,0L,1L))
  stopifnot(all(res$tissue_versatilities >= 0 & res$tissue_versatilities <= 1))
  stopifnot(all(res$tissue_angles_deg >= 0 & res$tissue_angles_deg <= 90))
  cat("test_randomized_vectors_axes passed\n")
}

# 10. Numerical stability (very large/small values)
test_numerical_stability <- function() {
  expr <- matrix(c(1e-15,1e-15,1e-15, 1e15,1e15,1e15), nrow=3, ncol=2)
  res <- tv_call(expr, c(1L,1L), c(1L,1L,1L))
  stopifnot(abs(res$tissue_versatilities[1]) < 1e-12)
  stopifnot(abs(res$tissue_angles_deg[1]) < 1e-12)
  stopifnot(abs(res$tissue_versatilities[2]) < 1e-12)
  stopifnot(abs(res$tissue_angles_deg[2]) < 1e-12)
  cat("test_numerical_stability passed\n")
}

# 11. Invalid input: no axes selected (should return -1)
test_invalid_input_no_axes <- function() {
  expr <- matrix(c(1,2,3), nrow=3, ncol=1)
  res <- tv_call(expr, c(1L), c(0L,0L,0L))
  stopifnot(abs(res$tissue_versatilities[1] + 1) < 1e-12)
  cat("test_invalid_input_no_axes passed\n")
}

# 12. Multiple selection, partial axes
test_multiple_selection_partial_axes <- function() {
  expr <- matrix(c(1,2, 3,4, 5,6), nrow=2, ncol=3)
  res <- tv_call(expr, c(1L,0L,1L), c(1L,0L))
  stopifnot(length(res$tissue_versatilities) == 2)
  stopifnot(all(res$tissue_versatilities >= 0 & res$tissue_versatilities <= 1))
  cat("test_multiple_selection_partial_axes passed\n")
}

# Run all tests
cat("=================================================\n")
cat("    TISSUE VERSATILITY FULL R INTERFACE TESTS\n")
cat("=================================================\n\n")

test_uniform_expression()
test_single_axis_expression()
test_null_vector()
test_partial_axis_selection()
test_mixed_vectors()
test_angle_degrees()
test_multiple_vectors_selection()
test_high_dimensional_vectors()
test_randomized_vectors_axes()
test_numerical_stability()
test_invalid_input_no_axes()
test_multiple_selection_partial_axes()

cat("=================================================\n")
cat("             ALL TESTS COMPLETED\n")
cat("=================================================\n")
cat("If you see this message, all tissue versatility R interface tests passed! ✓\n")
