## Deterministic tests: log2(x + 1) transformation
source("rcpp/tensoromics_functions.R")

test_basic_values <- function() {
	# Flattened example matching Fortran small-case: [0,3,7,15] -> [0,2,3,4]
	mat <- matrix(c(0,3,7,15), nrow = 2, ncol = 2)
	res <- tox_log2_transformation(mat)
	expected <- log2(mat + 1)
	if (!all.equal(as.numeric(res), as.numeric(expected), tolerance = 1e-12)) stop("log2_transformation_test: basic values mismatch")
	cat("log2_transformation_test - basic values: PASS\n")
}

test_zeros_handling <- function() {
	mat <- matrix(0, nrow = 2, ncol = 3)
	res <- tox_log2_transformation(mat)
	if (!all(res == 0)) stop("log2_transformation_test: zeros not preserved as zeros")
	cat("log2_transformation_test - zeros handling: PASS\n")
}

test_negative_handling <- function() {
	# Values between -1 and 0 are allowed (log2(x+1) defined)
	mat <- matrix(c(-0.5, -0.9, -0.99, -0.999), nrow = 2, ncol = 2)
	res <- tox_log2_transformation(mat)
	expected <- log2(mat + 1)
	if (!all.equal(as.numeric(res), as.numeric(expected), tolerance = 1e-12)) stop("log2_transformation_test: negative handling mismatch")
	cat("log2_transformation_test - negative handling: PASS\n")
}

test_empty_input_error <- function() {
	mat <- matrix(numeric(0), nrow = 0, ncol = 0)
	ok <- FALSE
	tryCatch({
		tox_log2_transformation(mat)
	}, error = function(e) {
		ok <<- grepl("Empty input arrays provided", e$message) || grepl("Empty input", e$message)
	})
	if (!ok) stop("log2_transformation_test: empty input did not raise expected error")
	cat("log2_transformation_test - empty input error: PASS\n")
}

# Run tests
test_basic_values()
test_zeros_handling()
test_negative_handling()
test_empty_input_error()

cat("log2_transformation_test: ALL PASSED\n")
