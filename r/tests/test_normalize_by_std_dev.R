library(testthat)

source("r/tensoromics_functions.R")

test_that("normalize_by_std_dev normalizes values correctly", {
  mat <- matrix(c(2, 4, 6, 8), nrow = 2)
  colnames(mat) <- c("tissue1", "tissue2")
  rownames(mat) <- c("geneA", "geneB")
  
  result <- normalize_by_std_dev(mat)
  
  expect_equal(dim(result), dim(mat))
  expect_equal(rownames(result), rownames(mat))
  expect_equal(colnames(result), colnames(mat))
  
  # Manually calculate expected normalization
  std_dev <- apply(mat, 1, function(x) sqrt(mean(x^2)))
  expected <- sweep(mat, 1, std_dev, FUN = "/")
  
  expect_equal(as.numeric(result), as.numeric(expected))
})

test_that("normalize_by_std_dev handles constant rows (zero std dev -> normalized to 1)", {
  mat <- matrix(c(5, 5, 5, 5), nrow = 2)
  colnames(mat) <- c("tissue1", "tissue2")
  rownames(mat) <- c("geneA", "geneB")

  result <- normalize_by_std_dev(mat)

  expect_equal(dim(result), dim(mat))
  expect_false(any(is.na(result)))
  expect_true(all(is.finite(result)))

  # Now expected result is all 1s (5/5 = 1)
  expected <- matrix(1, nrow = 2, ncol = 2)
  expect_equal(as.numeric(result), as.numeric(expected))
})


test_that("normalize_by_std_dev normalizes large numbers properly", {
  mat <- matrix(c(1e6, 2e6, 1e6, 2e6), nrow = 2)
  colnames(mat) <- c("tissue1", "tissue2")
  rownames(mat) <- c("geneA", "geneB")

  result <- normalize_by_std_dev(mat)

  expect_equal(dim(result), dim(mat))
  expect_true(all(is.finite(result)))

  # Check that normalization rescales based on sqrt(mean(x^2))
  std_dev <- apply(mat, 1, function(x) sqrt(mean(x^2)))
  expected <- sweep(mat, 1, std_dev, FUN = "/")

  expect_equal(as.numeric(result), as.numeric(expected))
})
