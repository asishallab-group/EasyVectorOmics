library(testthat)

source("r/tensoromics_functions.R")

test_that("calculate_tissue_averages averages replicates correctly across multiple tissues", {
  mat <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), nrow = 2)
  colnames(mat) <- c("tissue1_cont_1", "tissue1_cont_2", "tissue2_cont_1", "tissue2_cont_2", "tissue3_cont_1", "tissue3_cont_2")
  rownames(mat) <- c("geneA", "geneB")

  result <- calculate_tissue_averages(as.data.frame(mat))

  expect_equal(dim(result), c(2, 3))  # 3 tissues
  expect_true(all(colnames(result) %in% c("tissue1_cont", "tissue2_cont", "tissue3_cont")))
  
  # Check some values (mean of replicates)
  expect_equal(result["geneA", "tissue1_cont"], mean(c(1, 3)))
  expect_equal(result["geneA", "tissue2_cont"], mean(c(5, 7)))
  expect_equal(result["geneA", "tissue3_cont"], mean(c(9, 11)))
})


test_that("calculate_tissue_averages works with generic tissue replicates (3 tissues)", {
  mat <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), nrow = 2)
  colnames(mat) <- c("tissue1_1", "tissue1_2", "tissue2_1", "tissue2_2", "tissue3_1", "tissue3_2")
  rownames(mat) <- c("geneA", "geneB")

  result <- calculate_tissue_averages(as.data.frame(mat))

  expect_equal(dim(result), c(2, 3))  # 3 tissues
  expect_true(all(colnames(result) %in% c("tissue1", "tissue2", "tissue3")))
  
  # Check some values (mean of replicates)
  expect_equal(result["geneA", "tissue1"], mean(c(1, 3)))
  expect_equal(result["geneA", "tissue2"], mean(c(5, 7)))
  expect_equal(result["geneA", "tissue3"], mean(c(9, 11)))
})
