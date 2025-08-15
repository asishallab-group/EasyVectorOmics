# =====================
# Test cases for loess_smooth_2d_r Fortran wrapper
# =====================

dyn.load("build/libtensor-omics.so")

cat("\n[loess_smooth_2d_r] Case 1: Simple smoothing, no mask\n")
n_total <- as.integer(5)
n_target <- as.integer(2)
x_ref <- as.double(c(1, 2, 3, 4, 5))
y_ref <- matrix(as.double(c(10, 20, 30, 40, 50)), nrow = 1)
indices_used <- as.integer(1:5)
n_used <- as.integer(5)  # Using all indices
x_query <- as.double(c(2.5, 4.5))
kernel_sigma <- as.double(1.0)
kernel_cutoff <- as.double(2.0)
y_out <- matrix(0, nrow = 1, ncol = n_target)

res1 <- .Fortran("loess_smooth_2d_r",
  n_total = n_total,
  n_target = n_target,
  x_ref = x_ref,
  y_ref = y_ref,
  indices_used = indices_used,
  n_used = n_used,
  x_query = x_query,
  kernel_sigma = kernel_sigma,
  kernel_cutoff = kernel_cutoff,
  y_out = y_out
)
print(res1$y_out)
# Expected: y_out = smoothed values between y_ref points, using all points

cat("\n[loess_smooth_2d_r] Case 2: With mask (exclude some points)\n")
# Use only indices 1, 3, 5 (excluding points 2 and 4)
indices_used_filtered <- as.integer(c(1, 3, 5))
n_used_filtered <- as.integer(3)
y_out <- matrix(0, nrow = 1, ncol = n_target)
res2 <- .Fortran("loess_smooth_2d_r",
  n_total = n_total,
  n_target = n_target,
  x_ref = x_ref,
  y_ref = y_ref,
  indices_used = indices_used_filtered,
  n_used = n_used_filtered,
  x_query = x_query,
  kernel_sigma = kernel_sigma,
  kernel_cutoff = kernel_cutoff,
  y_out = y_out
)
print(res2$y_out)
# Expected: y_out = smoothed values using only x_ref/y_ref at positions 1, 3, 5

cat("\n[loess_smooth_2d_r] Case 3: Use only one point (should return that point's value)\n")
# Use only the first point
indices_used_single <- as.integer(c(1))
n_used_single <- as.integer(1)
y_out <- matrix(0, nrow = 1, ncol = n_target)
res3 <- .Fortran("loess_smooth_2d_r",
  n_total = n_total,
  n_target = n_target,
  x_ref = x_ref,
  y_ref = y_ref,
  indices_used = indices_used_single,
  n_used = n_used_single,
  x_query = x_query,
  kernel_sigma = kernel_sigma,
  kernel_cutoff = kernel_cutoff,
  y_out = y_out
)
print(res3$y_out)
# Expected: y_out = y_ref[1] for all queries (using only first point)


cat("\n[loess_smooth_2d_r] Case 4: x_query outside x_ref range\n")
# Expectation: y_out should extrapolate or fallback to nearest y_ref value
x_query <- as.double(c(-10, 100))
y_out <- matrix(0, nrow = 1, ncol = 2)
indices_used <- as.integer(1:5)
n_used <- as.integer(5)
res4 <- .Fortran("loess_smooth_2d_r",
  n_total = n_total,
  n_target = as.integer(2),
  x_ref = x_ref,
  y_ref = y_ref,
  indices_used = indices_used,
  n_used = n_used,
  x_query = x_query,
  kernel_sigma = kernel_sigma,
  kernel_cutoff = kernel_cutoff,
  y_out = y_out
)
print(res4$y_out)
# Expected: y_out close to y_ref[1] for -10, y_ref[5] for 100

cat("\n[loess_smooth_2d_r] Case 5: x_query == x_ref\n")
# Expectation: y_out should match y_ref at those points
x_query <- x_ref[1:2]
y_out <- matrix(0, nrow = 1, ncol = 2)
res5 <- .Fortran("loess_smooth_2d_r",
  n_total = n_total,
  n_target = as.integer(2),
  x_ref = x_ref,
  y_ref = y_ref,
  indices_used = indices_used,
  n_used = n_used,
  x_query = x_query,
  kernel_sigma = kernel_sigma,
  kernel_cutoff = kernel_cutoff,
  y_out = y_out
)
print(res5$y_out)
# Expected: y_out[1] == y_ref[1], y_out[2] == y_ref[2]

cat("\n[loess_smooth_2d_r] Case 6: y_ref contains NA\n")
# Expectation: Should handle NA gracefully, possibly returning NA or skipping in smoothing
y_ref_na <- y_ref
y_ref_na[1,3] <- NA_real_
y_out <- matrix(0, nrow = 1, ncol = n_target)
res6 <- tryCatch({
  .Fortran("loess_smooth_2d_r",
    n_total = n_total,
    n_target = n_target,
    x_ref = x_ref,
    y_ref = y_ref_na,
    indices_used = indices_used,
    n_used = n_used,
    x_query = x_query,
    kernel_sigma = kernel_sigma,
    kernel_cutoff = kernel_cutoff,
    y_out = y_out
  )
}, error = function(e) e)
print(res6$y_out)
# Expected: y_out with NA or interpolated values ignoring NA

cat("\n[loess_smooth_2d_r] Case 7: n_total = 1 (single point)\n")
# Expectation: y_out should be y_ref for all x_query
n_total1 <- as.integer(1)
x_ref1 <- as.double(42)
y_ref1 <- matrix(as.double(99), nrow = 1)
indices_used1 <- as.integer(1)
n_used1 <- as.integer(1)
x_query1 <- as.double(c(0, 42, 100))
y_out1 <- matrix(0, nrow = 1, ncol = 3)
res7 <- .Fortran("loess_smooth_2d_r",
  n_total = n_total1,
  n_target = as.integer(3),
  x_ref = x_ref1,
  y_ref = y_ref1,
  indices_used = indices_used1,
  n_used = n_used1,
  x_query = x_query1,
  kernel_sigma = kernel_sigma,
  kernel_cutoff = kernel_cutoff,
  y_out = y_out1
)
print(res7$y_out)
# Expected: y_out == 99 for all queries

cat("\n[loess_smooth_2d_r] Case 8: kernel_sigma = 0 (should return nearest y_ref)\n")
# Expectation: y_out should be y_ref at nearest x_ref for each x_query
kernel_sigma0 <- as.double(0)
indices_used <- as.integer(1:5)
n_used <- as.integer(5)
y_out <- matrix(0, nrow = 1, ncol = n_target)
res8 <- .Fortran("loess_smooth_2d_r",
  n_total = n_total,
  n_target = n_target,
  x_ref = x_ref,
  y_ref = y_ref,
  indices_used = indices_used,
  n_used = n_used,
  x_query = x_query,
  kernel_sigma = kernel_sigma0,
  kernel_cutoff = kernel_cutoff,
  y_out = y_out
)
print(res8$y_out)
# Expected: y_out == y_ref at nearest x_ref for each x_query

# =====================
# End of loess_smooth_2d_r tests
# =====================
