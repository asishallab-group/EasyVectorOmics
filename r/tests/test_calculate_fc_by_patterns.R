library(testthat)

source("r/tensoromics_functions.R")

test_that("calculate_fc_by_patterns basic calculation and result correctness", {
  mat <- matrix(c(1, 2, 4, 8), nrow = 2)
  colnames(mat) <- c("muscle_dietM", "muscle_dietP")
  rownames(mat) <- c("geneA", "geneB")
  
  result <- calculate_fc_by_patterns(as.data.frame(mat),
                                      control_pattern = "dietM",
                                      condition_patterns = c("dietP"))

  expect_equal(dim(result), c(2, 1))

  # Build the expected data.frame
  expected_fc <- data.frame(
    muscle_dietP_logFC = c(3, 6),
    row.names = c("geneA", "geneB")
  )

  expect_equal(result, expected_fc)
})

test_that("calculate_fc_by_patterns handles missing matches with proper error", {
  mat <- matrix(c(1, 2, 4, 8), nrow = 2)
  colnames(mat) <- c("muscle_dietM", "brain_dietQ")
  rownames(mat) <- c("geneA", "geneB")

  expect_error(
    calculate_fc_by_patterns(as.data.frame(mat),
                              control_pattern = "dietM",
                              condition_patterns = c("dietZ")),
    "No valid control-condition pairs found!"
  )
})
