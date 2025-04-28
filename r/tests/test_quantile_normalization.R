library(testthat)

source("r/tensoromics_functions.R")

test_that("quantile_normalization preserves dimensions and names", {
  mat <- matrix(runif(6), nrow = 2)
  colnames(mat) <- c("sample1", "sample2", "sample3")
  rownames(mat) <- c("geneA", "geneB")

  result <- quantile_normalization(mat)

  expect_equal(dim(result), dim(mat))
  expect_equal(rownames(result), rownames(mat))
  expect_equal(colnames(result), colnames(mat))
})

test_that("quantile_normalization handles identical rows and maintains values", {
  mat <- matrix(rep(5, 6), nrow = 2)
  colnames(mat) <- c("sample1", "sample2", "sample3")
  rownames(mat) <- c("geneA", "geneB")

  result <- quantile_normalization(mat)

  expect_equal(dim(result), dim(mat))
  expect_true(all(is.finite(result)))

  # Build expected matrix WITH names
  expected <- matrix(5, nrow = 2, ncol = 3)
  colnames(expected) <- c("sample1", "sample2", "sample3")
  rownames(expected) <- c("geneA", "geneB")

  expect_equal(result, expected)
})


test_that("quantile_normalization does not introduce NAs and standardizes distribution across columns", {
  mat <- matrix(c(2, 0, 5, 3, 7, 1), nrow = 2)
  colnames(mat) <- c("sample1", "sample2", "sample3")
  rownames(mat) <- c("geneA", "geneB")

  result <- quantile_normalization(mat)

  expect_false(any(is.na(result)))
  expect_true(all(is.finite(result)))
  expect_equal(dim(result), dim(mat))

  # Now the important check: column distributions should be similar
  # Sort columns and they should be (approximately) equal
  sorted_cols <- apply(result, 2, sort)
  
  for (i in 2:ncol(sorted_cols)) {
    expect_equal(sorted_cols[, i], sorted_cols[, 1], tolerance = 1e-6)
  }
})
