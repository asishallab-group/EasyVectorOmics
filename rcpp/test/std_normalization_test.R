## Deterministic tests: standard deviation normalization
source("rcpp/tensoromics_functions.R")

test_basic_normalization <- function() {
	# Small 2x2 example — compute expected by row-wise RMS (sqrt(mean(x^2)))
	mat <- matrix(c(2,6,4,8), nrow = 2, ncol = 2)
	# Each column is a sample; rows are genes
	expected <- mat
	for (i in seq_len(nrow(mat))) {
		rms <- sqrt(mean(mat[i, ]^2))
		expected[i, ] <- mat[i, ] / rms
	}

	res <- tox_normalize_by_std_dev(mat)
	if (!is.matrix(res)) stop("std_normalization_test: result is not a matrix")
	if (!all(dim(res) == dim(mat))) stop("std_normalization_test: dimension mismatch")
	if (!is.numeric(res)) stop("std_normalization_test: result not numeric")
	if (any(!is.finite(res))) stop("std_normalization_test: non-finite values found")

	if (!all.equal(as.numeric(res), as.numeric(expected), tolerance = 1e-12)) {
		stop("std_normalization_test: basic normalization numeric mismatch")
	}
	cat("std_normalization_test - basic: PASS\n")
}

test_constant_rows <- function() {
	mat <- matrix(5, nrow = 2, ncol = 3)
	res <- tox_normalize_by_std_dev(mat)
	# constant rows should normalize to 1.0 across each row
	if (!all(abs(res - 1.0) < 1e-12)) stop("std_normalization_test: constant rows did not normalize to 1")
	cat("std_normalization_test - constant rows: PASS\n")
}

test_zero_rows <- function() {
	mat <- matrix(0, nrow = 2, ncol = 4)
	res <- tox_normalize_by_std_dev(mat)
	if (!all(res == 0)) stop("std_normalization_test: zero rows did not remain zero")
	cat("std_normalization_test - zero rows: PASS\n")
}

test_single_row_col <- function() {
	# Single row
	mat1 <- matrix(c(2,4,6,8), nrow = 1)
	res1 <- tox_normalize_by_std_dev(mat1)
	rms1 <- sqrt(mean(mat1[1, ]^2))
	expected1 <- matrix(mat1[1, ] / rms1, nrow = 1)
	if (!all.equal(as.numeric(res1), as.numeric(expected1), tolerance = 1e-12)) stop("std_normalization_test: single row mismatch")

	# Single column
	mat2 <- matrix(c(2,4,6,8), nrow = 4, ncol = 1)
	res2 <- tox_normalize_by_std_dev(mat2)
	# For a single column, each row's RMS is absolute value; dividing yields sign-preserving 1 or 0
	if (!all(is.finite(res2))) stop("std_normalization_test: single column produced non-finite values")
	cat("std_normalization_test - single row/col: PASS\n")
}

test_empty_input_error <- function() {
	mat <- matrix(numeric(0), nrow = 0, ncol = 0)
	ok <- FALSE
	tryCatch({
		tox_normalize_by_std_dev(mat)
	}, error = function(e) {
		ok <<- grepl("Empty input arrays provided", e$message) || grepl("Empty input", e$message)
	})
	if (!ok) stop("std_normalization_test: empty input did not raise expected error")
	cat("std_normalization_test - empty input error: PASS\n")
}

# Run tests
test_basic_normalization()
test_constant_rows()
test_zero_rows()
test_single_row_col()
test_empty_input_error()

cat("std_normalization_test: ALL PASSED\n")
