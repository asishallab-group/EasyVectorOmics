## Deterministic tests: quantile normalization
source("rcpp/tensoromics_functions.R")

test_preserves_dimensions <- function() {
	mat <- matrix(runif(2 * 3), nrow = 2, ncol = 3)
	res <- tox_quantile_normalization(mat)
	if (!all(dim(res) == dim(mat))) stop("quantile_normalization_test: dimensions not preserved")
	cat("quantile_normalization_test - preserves dimensions: PASS\n")
}

test_identical_rows <- function() {
	mat <- matrix(5, nrow = 2, ncol = 3)
	res <- tox_quantile_normalization(mat)
	if (!all.equal(as.numeric(res), as.numeric(mat), tolerance = 1e-12)) stop("quantile_normalization_test: identical rows not preserved")
	cat("quantile_normalization_test - identical rows: PASS\n")
}

test_single_row_behavior <- function() {
	mat <- matrix(c(1,2,3,4), nrow = 1)
	res <- tox_quantile_normalization(mat)
	# For a single row, quantile normalization should make every element equal to the row mean
	expected <- matrix(mean(mat[1, ]), nrow = 1, ncol = ncol(mat))
	if (!all.equal(as.numeric(res), as.numeric(expected), tolerance = 1e-12)) stop("quantile_normalization_test: single row expectation failed")
	cat("quantile_normalization_test - single row: PASS\n")
}

test_zero_matrix <- function() {
	mat <- matrix(0, nrow = 3, ncol = 3)
	res <- tox_quantile_normalization(mat)
	if (!all(res == 0)) stop("quantile_normalization_test: zero matrix not preserved")
	cat("quantile_normalization_test - zero matrix: PASS\n")
}

test_empty_input_error <- function() {
	mat <- matrix(numeric(0), nrow = 0, ncol = 0)
	ok <- FALSE
	tryCatch({
		tox_quantile_normalization(mat)
	}, error = function(e) {
		ok <<- grepl("Empty input arrays provided", e$message) || grepl("Empty input", e$message)
	})
	if (!ok) stop("quantile_normalization_test: empty input did not raise expected error")
	cat("quantile_normalization_test - empty input error: PASS\n")
}

# Run tests
test_preserves_dimensions()
test_identical_rows()
test_single_row_behavior()
test_zero_matrix()
test_empty_input_error()

cat("quantile_normalization_test: ALL PASSED\n")
