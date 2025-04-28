library(testthat)

source("r/tensoromics_functions.R")

test_that("log2_transformation applies log2(x+1) properly to numeric values", {
  mat <- matrix(c(0, 3, 7, 15), nrow = 2)
  colnames(mat) <- c("sample1", "sample2")
  rownames(mat) <- c("geneA", "geneB")

  result <- log2_transformation(mat)

  expect_equal(dim(result), dim(mat))  # Dimension check
  expect_equal(rownames(result), rownames(mat))  # Row names preserved
  expect_equal(colnames(result), colnames(mat))  # Column names preserved
  
  expected_values <- log2(mat + 1)
  expect_equal(as.numeric(result), as.numeric(expected_values))  # Values
})

test_that("log2_transformation handles zeros correctly (log2(0+1) = 0)", {
  mat <- matrix(c(0, 0, 0, 0), nrow = 2)
  colnames(mat) <- c("sample1", "sample2")
  rownames(mat) <- c("geneA", "geneB")

  result <- log2_transformation(mat)

  expect_equal(dim(result), dim(mat))
  expect_true(all(result == 0))
  expect_equal(rownames(result), rownames(mat))
  expect_equal(colnames(result), colnames(mat))
})

test_that("log2_transformation preserves dimensions and names on random data", {
  mat <- matrix(runif(6), nrow = 2)
  colnames(mat) <- c("tissue1", "tissue2", "tissue3")
  rownames(mat) <- c("geneA", "geneB")

  result <- log2_transformation(mat)

  expect_equal(dim(result), dim(mat))
  expect_equal(rownames(result), rownames(mat))
  expect_equal(colnames(result), colnames(mat))
})
