# Test script for the Fortran gene centroid library
lib_path <- file.path("build", "libtensor-omics.so")
if (!file.exists(lib_path)) {
  stop(paste("Shared library not found at:", lib_path, "\nPlease run './build.sh' first."))
}
dyn.load(lib_path)

test_centroids <- function() {
  cat("=== Testing Centroid Logic via R Interface ===\n")
  
  # Arrange
  d <- 2L; n_genes <- 5L; n_families <- 2L
  vectors <- matrix(c(1,1, 3,3, 10,10, 20,20, 5,5), nrow=d, ncol=n_genes, byrow=FALSE)
  gene_to_family <- as.integer(c(1, 1, 2, 2, 1))
  orthologs <- as.logical(c(TRUE, FALSE, TRUE, TRUE, TRUE))
  
  # --- Test "all" mode ---
  mode_all <- "all"
  mode_all_ascii <- as.integer(charToRaw(mode_all))
  centroid_matrix_all <- matrix(0.0, nrow=n_families, ncol=d)
  
  result_all <- .Fortran("group_centroid_r",
                         vectors = as.double(vectors),
                         d = d, n = n_genes,
                         gene_to_family_map = gene_to_family,
                         num_families = n_families,
                         centroid_matrix = centroid_matrix_all,
                         mode_ascii = mode_all_ascii,
                         mode_len = as.integer(nchar(mode_all)),
                         ortholog_set = orthologs,
                         selected_indices = integer(n_genes))
  
  expected_all_fam1 <- c(3.0, 3.0)
  expected_all_fam2 <- c(15.0, 15.0)
  stopifnot(all.equal(result_all$centroid_matrix[1,], expected_all_fam1))
  stopifnot(all.equal(result_all$centroid_matrix[2,], expected_all_fam2))
  cat("Test 'all' mode: PASSED\n")
  
  # --- Test "orthologs" mode ---
  mode_ortho <- "orthologs"
  mode_ortho_ascii <- as.integer(charToRaw(mode_ortho))
  centroid_matrix_ortho <- matrix(0.0, nrow=n_families, ncol=d)
  
  result_ortho <- .Fortran("group_centroid_r",
                           vectors = as.double(vectors),
                           d = d, n = n_genes,
                           gene_to_family_map = gene_to_family,
                           num_families = n_families,
                           centroid_matrix = centroid_matrix_ortho,
                           mode_ascii = mode_ortho_ascii,
                           mode_len = as.integer(nchar(mode_ortho)),
                           ortholog_set = orthologs,
                           selected_indices = integer(n_genes))

  expected_ortho_fam1 <- c(3.0, 3.0)
  stopifnot(all.equal(result_ortho$centroid_matrix[1,], expected_ortho_fam1))
  cat("Test 'orthologs' mode: PASSED\n")
}

# Run tests
test_centroids()
cat("\nAll R tests for gene_centroids passed! ✓\n")
