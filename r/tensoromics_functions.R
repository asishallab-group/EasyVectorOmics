# Load the shared library once when the package/script is loaded.
# Assumes the script is run from the project root directory.
.lib_path <- file.path("build", "libtensor-omics.so")
if (file.exists(.lib_path)) {
  dyn.load(.lib_path)
} else {
  warning(paste("Shared library not found at:", .lib_path, "\nRun './build.sh' from the project root."))
}


tox_group_centroid <- function(vectors, gene_to_family_map, num_families, ortholog_set, mode = 'all') {
  
  # 1) Validate inputs
  if (!is.matrix(vectors) || !is.numeric(vectors)) {
    stop("`vectors` must be a numeric matrix.")
  }
  d <- nrow(vectors)
  n_genes <- ncol(vectors)
  
  if (!is.integer(gene_to_family_map) || length(gene_to_family_map) != n_genes) {
    stop("`gene_to_family_map` must be an integer vector of length n_genes.")
  }
  if (!is.logical(ortholog_set) || length(ortholog_set) != n_genes) {
    stop("`ortholog_set` must be a logical vector of length n_genes.")
  }
  if (!mode %in% c('all', 'ortho')) {
    stop("`mode` must be either 'all' or 'ortho'.")
  }

  # 2) Prepare inputs/outputs for Fortran
  use_all_mode <- (mode == 'all')
  centroid_matrix_out <- matrix(0.0, nrow = d, ncol = num_families)
  selected_indices_ws <- integer(n_genes) # Workspace buffer

  # 3) Call Fortran
  result <- .Fortran("group_centroid_r",
                     vectors = as.double(vectors),
                     d = as.integer(d),
                     n = as.integer(n_genes),
                     gene_to_family_map = as.integer(gene_to_family_map),
                     num_families = as.integer(num_families),
                     centroid_matrix = centroid_matrix_out,
                     use_all_mode = as.logical(use_all_mode),
                     ortholog_set = as.logical(ortholog_set),
                     selected_indices = selected_indices_ws,
                     selected_indices_len = as.integer(n_genes))
  
  # 4) Return the populated output matrix
  return(result$centroid_matrix)
}
