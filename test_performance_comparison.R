# Performance comparison: Rcpp vs .Fortran approaches

# Load Rcpp functions
source("rcpp/load_tensoromics.R")

# Store Rcpp function with different name to avoid conflicts
rcpp_euclidean_distance <- tox_euclidean_distance

# Load R/.Fortran functions (will overwrite tox_euclidean_distance)
source("r/tensoromics_functions.R")

# Store .Fortran function with different name  
fortran_euclidean_distance <- tox_euclidean_distance

# Test with vectors
n <- 1000000  # Larger test for more visible differences
vec1 <- rnorm(n)
vec2 <- rnorm(n)

cat("=== Performance Comparison: Rcpp vs .Fortran ===\n")
cat("Vector size:", length(vec1), "elements\n")
cat("Memory per vector:", round(as.numeric(object.size(vec1))/1024/1024, 2), "MB\n\n")

# Benchmark Rcpp (direct pointers)
cat("=== Rcpp Test (with direct pointers) ===\n")
time_rcpp <- system.time({
  result_rcpp <- rcpp_euclidean_distance(vec1, vec2)
})
cat("Rcpp Result:", result_rcpp, "\n")
print(time_rcpp)
cat("\n")

# Benchmark .Fortran (with copies)
cat("=== .Fortran Test (with data copies) ===\n")
time_fortran <- system.time({
  result_fortran <- fortran_euclidean_distance(vec1, vec2)
})
cat(".Fortran Result:", result_fortran, "\n")
print(time_fortran)
cat("\n")

# Compare results and performance
cat("=== Comparison ===\n")
cat("Results match:", all.equal(result_rcpp, result_fortran), "\n")
cat("Rcpp elapsed time:", time_rcpp[3], "seconds\n")
cat(".Fortran elapsed time:", time_fortran[3], "seconds\n")
speedup <- time_fortran[3] / time_rcpp[3]
cat("Speedup (Rcpp vs .Fortran):", round(speedup, 2), "x\n\n")

cat("=== Technical Details ===\n")
cat("Rcpp: v1.begin() -> direct pointer to data (no copying)\n") 
cat(".Fortran: as.double(v1) -> complete data copy before call\n")
cat("Memory overhead with .Fortran:", round(2 * as.numeric(object.size(vec1))/1024/1024, 2), "MB extra copies\n")
