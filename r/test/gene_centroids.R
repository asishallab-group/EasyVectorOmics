# Source the user-facing function.
# Assumes the test is run from the project root directory.
source(file.path("r", "tensoromics_functions.R"))

test_centroids <- function() {
  cat("=== Testing tox_group_centroid via R user-facing API ===\n")
  
  # Arrange
  d <- 2L; n_genes <- 5L; n_families <- 2L
  
  # Expression vectors (d x n_genes)
  vectors <- matrix(c(1,1, 3,3, 10,10, 20,20, 5,5), nrow=d, ncol=n_genes, byrow=FALSE)
  
  # Gene-to-family mapping
  gene_to_family <- as.integer(c(1, 1, 2, 2, 1))
  
  # Orthologs mask
  ortholog_set <- as.logical(c(TRUE, FALSE, TRUE, TRUE, TRUE))
  
  # --- Test "all" mode ---
  cat("\n[tox_group_centroid] Case 1: 'all' mode\n")
  centroids_all <- tox_group_centroid(vectors, gene_to_family, n_families, ortholog_set, mode = 'all')
  
  # Expected: Fam1=[(1,1)+(3,3)+(5,5)]/3 = [3,3], Fam2=[(10,10)+(20,20)]/2 = [15,15]
  expected_all <- matrix(c(3.0, 3.0, 15.0, 15.0), nrow=d, ncol=n_families)
  
  stopifnot(all.equal(centroids_all, expected_all))
  cat("Test 'all' mode: PASSED\n")
  
  # --- Test "ortho" mode ---
  cat("\n[tox_group_centroid] Case 2: 'ortho' mode\n")
  centroids_ortho <- tox_group_centroid(vectors, gene_to_family, n_families, ortholog_set, mode = 'ortho')
  
  # Expected: Fam1=[(1,1)+(5,5)]/2 = [3,3], Fam2=[(10,10)+(20,20)]/2 = [15,15]
  expected_ortho <- matrix(c(3.0, 3.0, 15.0, 15.0), nrow=d, ncol=n_families)
  
  stopifnot(all.equal(centroids_ortho, expected_ortho))
  cat("Test 'ortho' mode: PASSED\n")
}

# Run tests
tryCatch({
  test_centroids()
  cat("\nAll R API tests for gene_centroids passed! ✓\n")
}, error = function(e) {
  cat("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
  cat("                  A TEST FAILED\n")
  cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
  cat("ERROR:", conditionMessage(e), "\n")
})
