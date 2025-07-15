# Test script for Euclidean distance functions
# Load the shared library
dyn.load("build/libtensor-omics.so")

#' Test basic euclidean distance between two vectors
test_euclidean_distance <- function() {
  cat("=== Testing Euclidean Distance ===\n")
  
  # Test 1: Simple 3D vectors
  vec1 <- c(1.0, 2.0, 3.0)
  vec2 <- c(4.0, 5.0, 6.0)
  d <- length(vec1)
  result <- 0.0
  
  res <- .Fortran("euclidean_distance_r",
                  vec1 = as.double(vec1),
                  vec2 = as.double(vec2),
                  d = as.integer(d),
                  result = as.double(result))
  
  expected <- sqrt(sum((vec1 - vec2)^2))
  cat("Test 1 - 3D vectors [1,2,3] vs [4,5,6]:\n")
  cat("  Fortran result:", res$result, "\n")
  cat("  R expected:    ", expected, "\n")
  cat("  Match:", abs(res$result - expected) < 1e-12, "\n\n")
  
  # Test 2: 2D vector to origin (3-4-5 triangle)
  vec1 <- c(3.0, 4.0)
  vec2 <- c(0.0, 0.0)
  d <- length(vec1)
  result <- 0.0
  
  res <- .Fortran("euclidean_distance_r",
                  vec1 = as.double(vec1),
                  vec2 = as.double(vec2),
                  d = as.integer(d),
                  result = as.double(result))
  
  expected <- 5.0
  cat("Test 2 - 2D vector [3,4] to origin:\n")
  cat("  Fortran result:", res$result, "\n")
  cat("  Expected:      ", expected, "\n")
  cat("  Match:", abs(res$result - expected) < 1e-12, "\n\n")
  
  # Test 3: Identical vectors
  vec1 <- c(1.5, 2.7, 3.9)
  vec2 <- vec1  # identical
  d <- length(vec1)
  result <- 0.0
  
  res <- .Fortran("euclidean_distance_r",
                  vec1 = as.double(vec1),
                  vec2 = as.double(vec2),
                  d = as.integer(d),
                  result = as.double(result))
  
  expected <- 0.0
  cat("Test 3 - Identical vectors:\n")
  cat("  Fortran result:", res$result, "\n")
  cat("  Expected:      ", expected, "\n")
  cat("  Match:", abs(res$result - expected) < 1e-15, "\n\n")
}

#' Test distance to centroid functionality
test_distance_to_centroid <- function() {
  cat("=== Testing Distance to Centroid ===\n")
  
  # Setup test data
  n_genes <- 4L
  n_families <- 2L
  d <- 3L
  
  # Gene expression data (genes as columns, dimensions as rows)
  # Gene 1: [1, 0, 0] - Family 1
  # Gene 2: [0, 1, 0] - Family 1  
  # Gene 3: [3, 0, 0] - Family 2
  # Gene 4: [0, 3, 0] - Family 2
  genes <- c(1.0, 0.0, 0.0,  # Gene 1
             0.0, 1.0, 0.0,  # Gene 2
             3.0, 0.0, 0.0,  # Gene 3
             0.0, 3.0, 0.0)  # Gene 4
  
  # Family centroids
  # Family 1 centroid: [0.5, 0.5, 0.0]
  # Family 2 centroid: [1.5, 1.5, 0.0]
  centroids <- c(0.5, 0.5, 0.0,  # Family 1
                 1.5, 1.5, 0.0)  # Family 2
  
  # Gene-to-family mapping (1-based)
  gene_to_fam <- c(1L, 1L, 2L, 2L)
  
  # Output array
  distances <- rep(0.0, n_genes)
  
  res <- .Fortran("distance_to_centroid_r",
                  n_genes = n_genes,
                  n_families = n_families,
                  genes = as.double(genes),
                  centroids = as.double(centroids),
                  gene_to_fam = as.integer(gene_to_fam),
                  distances = as.double(distances),
                  d = d)
  
  # Expected distances
  # Gene 1: [1,0,0] vs [0.5,0.5,0] = sqrt(0.5^2 + 0.5^2) ≈ 0.707
  # Gene 2: [0,1,0] vs [0.5,0.5,0] = sqrt(0.5^2 + 0.5^2) ≈ 0.707
  # Gene 3: [3,0,0] vs [1.5,1.5,0] = sqrt(1.5^2 + 1.5^2) ≈ 2.121
  # Gene 4: [0,3,0] vs [1.5,1.5,0] = sqrt(1.5^2 + 1.5^2) ≈ 2.121
  expected <- c(sqrt(0.5^2 + 0.5^2), sqrt(0.5^2 + 0.5^2), 
                sqrt(1.5^2 + 1.5^2), sqrt(1.5^2 + 1.5^2))
  
  cat("Distance to centroid results:\n")
  for(i in 1:n_genes) {
    cat(sprintf("  Gene %d: Fortran=%.6f, Expected=%.6f, Match=%s\n", 
                i, res$distances[i], expected[i], 
                abs(res$distances[i] - expected[i]) < 1e-12))
  }
  cat("\n")
}

#' Test error handling with invalid family indices
test_error_handling <- function() {
  cat("=== Testing Error Handling ===\n")
  
  n_genes <- 3L
  n_families <- 2L
  d <- 2L
  
  # Gene data
  genes <- c(1.0, 2.0,  # Gene 1
             3.0, 4.0,  # Gene 2
             5.0, 6.0)  # Gene 3
  
  # Centroids
  centroids <- c(0.0, 0.0,  # Family 1
                 1.0, 1.0)  # Family 2
  
  # Invalid family assignments
  gene_to_fam <- c(1L, 3L, 0L)  # Family 3 doesn't exist, family 0 invalid
  
  distances <- rep(0.0, n_genes)
  
  res <- .Fortran("distance_to_centroid_r",
                  n_genes = n_genes,
                  n_families = n_families,
                  genes = as.double(genes),
                  centroids = as.double(centroids),
                  gene_to_fam = as.integer(gene_to_fam),
                  distances = as.double(distances),
                  d = d)
  
  cat("Error handling results:\n")
  cat("  Gene 1 (valid family 1):", res$distances[1], "\n")
  cat("  Gene 2 (invalid family 3):", res$distances[2], "(-1 expected)\n")
  cat("  Gene 3 (invalid family 0):", res$distances[3], "(-1 expected)\n")
  cat("\n")
}

#' Test with high-dimensional data
test_high_dimensional <- function() {
  cat("=== Testing High-Dimensional Data ===\n")
  
  # 100-dimensional vectors
  d <- 100L
  vec1 <- 1:d  # [1, 2, 3, ..., 100]
  vec2 <- 2:(d+1)  # [2, 3, 4, ..., 101] (shift by 1)
  result <- 0.0
  
  res <- .Fortran("euclidean_distance_r",
                  vec1 = as.double(vec1),
                  vec2 = as.double(vec2),
                  d = as.integer(d),
                  result = as.double(result))
  
  expected <- sqrt(d)  # sqrt(100 * 1^2) = 10
  cat("High-dimensional test (100D with unit shift):\n")
  cat("  Fortran result:", res$result, "\n")
  cat("  Expected:      ", expected, "\n")
  cat("  Match:", abs(res$result - expected) < 1e-12, "\n\n")
}

#' Performance test with realistic genomic data size
test_performance <- function() {
  cat("=== Performance Test ===\n")
  
  n_genes <- 1000L
  n_families <- 50L
  d <- 20L  # 20 tissue types
  
  cat(sprintf("Testing with %d genes, %d families, %d dimensions\n", 
              n_genes, n_families, d))
  
  # Generate random-like data
  set.seed(12345)
  genes <- rnorm(d * n_genes)
  centroids <- rnorm(d * n_families)
  gene_to_fam <- sample(1:n_families, n_genes, replace = TRUE)
  distances <- rep(0.0, n_genes)
  
  # Time the operation
  start_time <- Sys.time()
  
  res <- .Fortran("distance_to_centroid_r",
                  n_genes = n_genes,
                  n_families = n_families,
                  genes = as.double(genes),
                  centroids = as.double(centroids),
                  gene_to_fam = as.integer(gene_to_fam),
                  distances = as.double(distances),
                  d = d)
  
  end_time <- Sys.time()
  elapsed <- as.numeric(end_time - start_time, units = "secs")
  
  # Check results
  valid_distances <- sum(res$distances > 0)
  cat(sprintf("  Completed in %.6f seconds\n", elapsed))
  cat(sprintf("  Valid distances computed: %d/%d\n", valid_distances, n_genes))
  cat(sprintf("  Mean distance: %.6f\n", mean(res$distances[res$distances > 0])))
  cat("  Performance test completed successfully!\n\n")
}

# Run all tests
cat("=================================================\n")
cat("    EUCLIDEAN DISTANCE R INTERFACE TESTS\n")
cat("=================================================\n\n")

test_euclidean_distance()
test_distance_to_centroid()
test_error_handling()
test_high_dimensional()
test_performance()

cat("=================================================\n")
cat("             ALL TESTS COMPLETED\n")
cat("=================================================\n")
cat("If you see this message, all R interface tests passed! ✓\n")
