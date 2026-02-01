source("rcpp/tensoromics_functions.R")
source("rcpp/test_helpers.R")

TOL <- 1e-12

test_determine_shared_residual_range <- function() {

  # Helper
  approx_equal <- function(a, b) abs(a - b) < TOL

  # Test 1 — Basic correctness
  S1 <- array(c(
    1,2,3,4,
    5,6,-7,8,
    9,10,11,12,
    1,1,1,1
  ), dim = c(4, 2, 2))

  S2 <- array(c(
    2,-4,6,8,
    1,3,5,7,
    9,0,1,2
  ), dim = c(3, 2, 2))

  R <- tox_determine_shared_residual_range(S1, S2, 95)
  assertTrue(approx_equal(R, 10.65), "Test 1 failed: expected ~10.65")

  # Test 2 — Custom quantile
  R <- tox_determine_shared_residual_range(S1, S2, 50)
  assertTrue(approx_equal(R, 4.0), "Test 2 failed: expected ~4.0")

  # Test 3 — Quantile < 0 → error
  assertError(
    tox_determine_shared_residual_range(S1, S2, -1),
    "Test 3 failed: expected error for negative quantile"
  )

  # Test 4 — Quantile > 100 → error
  assertError(
    tox_determine_shared_residual_range(S1, S2, 150),
    "Test 4 failed: expected error for quantile > 100"
  )

  # Test 5 — NaNs ignored
  S1 <- array(c(
    NA_real_,2,3,
    4,5,6,
    -7,8,-11,
    9,10,12,
    1,1,1,1
  ), dim = c(4,2,2))
  S2 <- array(c(
    -1,2,3,
    4,5,6,
    -7,8,-11,
    9,10,NA_real_
  ), dim = c(3,2,2))

  R <- tox_determine_shared_residual_range(S1, S2, 95)
  assertTrue(approx_equal(R, 11.0), "Test 5 failed: expected ~11.0")

  # Test 6 — All zeros
  S1 <- array(0, dim = c(4,2,2))
  S2 <- array(0, dim = c(3,2,2))
  R <- tox_determine_shared_residual_range(S1, S2, 95)
  assertTrue(approx_equal(R, 0.0), "Test 6 failed: expected 0")

  # Test 7 — Single residual
  S1 <- array(3, dim = c(1, 1, 1))
  S2 <- array(-4, dim = c(1, 1, 1))
  R <- tox_determine_shared_residual_range(S1, S2, 95)
  assertTrue(approx_equal(R, 3.95), "Test 7 failed: expected ~3.95")
}

test_determine_shared_residual_range_expert <- function() {
  TOL <- 1e-12
  approx_equal <- function(a, b) abs(a - b) < TOL

  make_pool <- function(S1, S2) {
    pool <- abs(c(S1, S2))
     perm <- order(pool)
    list(pool = pool, perm = perm)
  }

  S1 <- array(c(
    1,2,3,4,
    5,6,-7,8,
    9,10,11,12,
    1,1,1,1
  ), dim = c(4, 2, 2))

  S2 <- array(c(
    2,-4,6,8,
    1,3,5,7,
    9,0,1,2
  ), dim = c(3, 2, 2))

  pp <- make_pool(S1, S2)

  # Test 1
  R <- tox_determine_shared_residual_range_expert(pp$pool, pp$perm, 95)
  assertTrue(approx_equal(R, 10.65), "Test 1 failed")

  # Test 2
  R <- tox_determine_shared_residual_range_expert(pp$pool, pp$perm, 50)
  assertTrue(approx_equal(R, 4.0), "Test 2 failed")

  # Test 3
  assertError(
    tox_determine_shared_residual_range_expert(pp$pool, pp$perm, -1),
    "Test 3 failed"
  )

  # Test 4
  assertError(
    tox_determine_shared_residual_range_expert(pp$pool, pp$perm, 150),
    "Test 4 failed"
  )

  # Test 5 — NaNs ignored
  S1 <- array(c(
    NA_real_,2,3,
    4,5,6,
    -7,8,-11,
    9,10,12,
    1,1,1,1
  ), dim = c(4,2,2))
  S2 <- array(c(
    -1,2,3,
    4,5,6,
    -7,8,-11,
    9,10,NA_real_
  ), dim = c(3,2,2))

  pp <- make_pool(S1, S2)
  R <- tox_determine_shared_residual_range_expert(pp$pool, pp$perm, 95)
  assertTrue(approx_equal(R, 11.0), "Test 5 failed")

  # Test 6 — All zeros
  S1 <- array(0, dim = c(4,2,2))
  S2 <- array(0, dim = c(3,2,2))
  pp <- make_pool(S1, S2)
  R <- tox_determine_shared_residual_range_expert(pp$pool, pp$perm, 95)
  assertTrue(approx_equal(R, 0.0), "Test 6 failed")

  # Test 7 — Single residual
  S1 <- array(3, dim = c(1, 1, 1))
  S2 <- array(-4, dim = c(1, 1, 1))
  pp <- make_pool(S1, S2)
  R <- tox_determine_shared_residual_range_expert(pp$pool, pp$perm, 95)
  assertTrue(approx_equal(R, 3.95), "Test 7 failed")
}

test_build_residual_histograms <- function() {
  TOL <- 1e-12
  R <- 2
  n_bins <- 4

  # Helper
  mat_equal <- function(a, b) all(a == b)
  mat_close <- function(a, b) all(abs(a - b) < TOL)

  # Test 1
  E <- array(c(
    -2, -0.5, 0.2, 1.7, 0.9, -1.2,
    0, 0, 0, 0, 0, 0,
    2.5, -3, 1.2, 0.4, -0.1, 0
  ), dim = c(3,2,3))

  out <- tox_build_residual_histograms(E, R, n_bins)

  expected_counts <- matrix(c(
    2,1,2,1,
    0,0,6,0,
    1,1,2,2
  ), 3, 4, byrow=TRUE)

  expected_pmf <- matrix(c(
    2/6,1/6,2/6,1/6,
    0,0,1,0,
    1/6,1/6,1/3,1/3
  ), 3, 4, byrow=TRUE)

  assertTrue(mat_equal(out$counts, expected_counts), "Test 1 counts mismatch")
  assertTrue(mat_close(out$pmf, expected_pmf), "Test 1 pmf mismatch")
  assertTrue(all(out$included_n_residuals == c(6,6,6)), "Test 1 included mismatch")

  # Test 2 — NaNs ignored
  E <- array(0, dim = c(3,2,3))
  E[1,1,1] <- NA_real_
  E[3,1,1] <- NA_real_
  E[2,1,2] <- NA_real_
  E[3,2,3] <- NA_real_

  out <- tox_build_residual_histograms(E, R, n_bins)

  expected_counts <- matrix(c(
    0,0,4,0,
    0,0,5,0,
    0,0,5,0
  ), 3, 4, byrow = TRUE)

  expected_pmf <- matrix(c(
    0,0,1,0,
    0,0,1,0,
    0,0,1,0
  ), 3, 4, byrow = TRUE)

  assertTrue(mat_equal(out$counts, expected_counts), "Test 2 counts mismatch")
  assertTrue(mat_close(out$pmf, expected_pmf), "Test 2 pmf mismatch")
  assertTrue(all(out$included_n_residuals == c(4,5,5)), "Test 2 included mismatch")

  # Test 3 — All NaN
  E <- array(NA_real_, dim = c(3,2,3))

  out <- tox_build_residual_histograms(E, R, n_bins)

  assertTrue(mat_equal(out$counts, matrix(0, 3, 4)), "Test 3 counts mismatch")
  assertTrue(mat_equal(out$pmf, matrix(0, 3, 4)), "Test 3 pmf mismatch")
  assertTrue(all(out$included_n_residuals == c(0,0,0)), "Test 3 included mismatch")

  # Test 4 — Boundary values
  E <- array(c(
    -2,-1,0,1,2,0,
    -2,-1,0,1,2,0,
    -2,-1,0,1,2,0
  ), dim = c(3,2,3))

  out <- tox_build_residual_histograms(E, R, n_bins)

  expected_counts <- matrix(c(
    1,1,2,2,
    1,1,2,2,
    1,1,2,2
  ), 3, 4, byrow = TRUE)

  expected_pmf <- expected_counts / 6

  assertTrue(mat_equal(out$counts, expected_counts), "Test 4 counts mismatch")
  assertTrue(mat_close(out$pmf, expected_pmf), "Test 4 pmf mismatch")
  assertTrue(all(out$included_n_residuals == c(6,6,6)), "Test 4 included mismatch")
}

test_compute_divergence_per_reference_point <- function() {
  TOL <- 1e-12
  approx_equal <- function(a, b) abs(a - b) < TOL

  # Test 1 — identical PMFs
  p <- matrix(c(
    0.1,0.2,0.3,0.4,
    0.25,0.25,0.25,0.25,
    1,0,0,0
  ), 3, 4, byrow = TRUE)

  q <- p

  jsd <- tox_compute_divergence_per_reference_point(p, q)
  assertTrue(all(abs(jsd) < TOL), "Test 1 failed")

  # Test 2 — disjoint PMFs
  p <- matrix(0, 3, 4)
  q <- matrix(0, 3, 4)
  p[1,1] <- 1
  q[1,2] <- 1

  jsd <- tox_compute_divergence_per_reference_point(p, q)
  assertTrue(approx_equal(jsd[1], log(2)), "Test 2 failed")
  assertTrue(jsd[2] == 0, "Test 2 row 2 failed")
  assertTrue(jsd[3] == 0, "Test 2 row 3 failed")

  # Test 3 — partial overlap
  p <- matrix(0, 3, 4)
  q <- matrix(0, 3, 4)
  p[1,] <- c(0.5,0.5,0,0)
  q[1,] <- c(0,1,0,0)

  jsd <- tox_compute_divergence_per_reference_point(p, q)

  expected <- 0.5 * (
    0.5 * log(2) +
    0.5 * log(2/3) +
    log(1/0.75)
  )

  assertTrue(approx_equal(jsd[1], expected), "Test 3 failed")

  # Test 4 — zero-probability bins
  p <- matrix(0, 3, 4)
  q <- matrix(0, 3, 4)
  p[1,1] <- 1
  q[1,3] <- 1

  jsd <- tox_compute_divergence_per_reference_point(p, q)
  assertTrue(approx_equal(jsd[1], log(2)), "Test 4 failed")

  # Test 5 — mixed patterns
  p <- matrix(
    c(
      0.2, 0.3, 0.5,
      0.0, 1.0, 0.0,
      0.0, 0.0, 0.25,
      0.25, 0.25, 0.25
    ),
    nrow = 3,
    ncol = 4
  )

  q <- matrix(
    c(
      0.2, 0.3, 0.5,
      0.0, 0.0, 1.0,
      0.0, 0.0, 0.25,
      0.25, 0.25, 0.25
    ),
    nrow = 3,
    ncol = 4
  )

  jsd <- tox_compute_divergence_per_reference_point(p, q)

  assertTrue(approx_equal(jsd[2], 0.5 * log(2)), "Test 5 row 2 failed")
  assertTrue(approx_equal(jsd[3], 0.5 * log(2)), "Test 5 row 3 failed")
}

test_compute_weighted_global_divergence <- function() {
  TOL <- 1e-12
  approx_equal <- function(a, b) abs(a - b) < TOL

  # Test 1 — uniform weights
  jsd <- c(0.1,0.2,0.3,0.4)
  n1 <- c(5L,5L,5L,5L)
  n2 <- c(5L,5L,5L,5L)
  out <- tox_compute_weighted_global_divergence(jsd, n1, n2)

  assertTrue(all(out$weights == 0.25), "Test 1 weights mismatch")
  assertTrue(approx_equal(out$global_js_divergence, 0.25), "Test 1 global JSD mismatch")

  # Test 2 — unequal sample counts
  jsd <- c(1,2,3,4)
  n1 <- c(10L,20L,30L,40L)
  n2 <- c(0L,10L,10L,10L)

  out <- tox_compute_weighted_global_divergence(jsd, n1, n2)

  expected_w <- c(10,30,40,50) / 130
  assertTrue(all(abs(out$weights - expected_w) < TOL), "Test 2 weights mismatch")

  expected <- sum(jsd * expected_w)
  assertTrue(approx_equal(out$global_js_divergence, expected), "Test 2 global JSD mismatch")

  # Test 3 — zero-sample neighborhoods
  jsd <- c(0.5,1.0,2.0,4.0)
  n1 <- c(0L,10L,0L,5L)
  n2 <- c(0L,0L,0L,5L)

  out <- tox_compute_weighted_global_divergence(jsd, n1, n2)

  expected_w <- c(0,0.5,0,0.5)
  assertTrue(all(abs(out$weights - expected_w) < TOL), "Test 3 weights mismatch")
  assertTrue(approx_equal(out$global_js_divergence, 2.5), "Test 3 global JSD mismatch")

  # Test 4 — all zero samples
  jsd <- c(1,2,3,4)
  n1 <- c(0L,0L,0L,0L)
  n2 <- c(0L,0L,0L,0L)

  out <- tox_compute_weighted_global_divergence(jsd, n1, n2)

  assertTrue(all(out$weights == 0), "Test 4 weights mismatch")
  assertTrue(out$global_js_divergence == 0, "Test 4 global JSD mismatch")

  # Test 5 — mixed jsd values, mixed sample counts
  jsd <- c(0.0, 0.5, 1.0, 2.0)
  n1  <- c(5L, 0L, 10L, 5L)
  n2  <- c(5L, 5L,  0L, 5L)

  out <- tox_compute_weighted_global_divergence(jsd, n1, n2)

  expected_w <- c(10, 5, 10, 10) / 35
  assertTrue(all(abs(out$weights - expected_w) < TOL),
             "Test 5 failed: weights mismatch")

  expected <- sum(jsd * expected_w)
  assertTrue(abs(out$global_js_divergence - expected) < TOL,
             "Test 5 failed: global JSD mismatch")
}

run_all_tests()