# Performance comparison with real data: Rcpp vs .Fortran approaches
# Uses the same real gene expression data with ~10,000 genes

library(readr)

# Load Rcpp functions
source("rcpp/load_tensoromics.R")

# Store Rcpp function with different name to avoid conflicts
rcpp_distance_to_centroid <- tox_distance_to_centroid

# Load R/.Fortran functions (will overwrite tox_distance_to_centroid)
source("r/tensoromics_functions.R")

# Store .Fortran function with different name  
fortran_distance_to_centroid <- tox_distance_to_centroid

# Function to generate gene_to_family mapping (copied from real data example)
generate_gene_to_family_mapping <- function(orthogroups_file, centroids_file, gene_expression_file, 
                                           use_all_species = TRUE, target_species = NULL) {
  cat("Loading data files...\n")
  
  # Read files
  orthogroups <- read_tsv(orthogroups_file, show_col_types = FALSE)
  centroids <- read_tsv(centroids_file, show_col_types = FALSE)
  gene_expr <- read_tsv(gene_expression_file, show_col_types = FALSE)
  
  # Detect available species columns (all except "Orthogroup")
  species_columns <- colnames(orthogroups)[colnames(orthogroups) != "Orthogroup"]
  
  # Determine which columns to use for gene mapping
  if (!is.null(target_species)) {
    if (target_species %in% species_columns) {
      target_species_list <- target_species
    } else {
      stop("Error: species '", target_species, "' not available. Available species: ", 
           paste(species_columns, collapse = ", "))
    }
  } else if (use_all_species) {
    target_species_list <- species_columns
  }
  
  # Extract gene list from expression file
  gene_id_col <- names(gene_expr)[1]
  genes_in_expression <- gene_expr[[gene_id_col]]
  
  # Create gene_id to family index mapping
  gene_to_family <- rep(0, length(genes_in_expression))
  names(gene_to_family) <- genes_in_expression
  
  # Create orthogroup to numeric index mapping
  orthogroup_col <- names(centroids)[1]
  orthogroup_to_index <- setNames(1:nrow(centroids), centroids[[orthogroup_col]])
  
  # Generate gene-to-family mapping
  for (i in 1:nrow(orthogroups)) {
    orthogroup_id <- orthogroups$Orthogroup[i]
    
    if (!orthogroup_id %in% names(orthogroup_to_index)) {
      next
    }
    
    family_index <- orthogroup_to_index[orthogroup_id]
    
    for (species in target_species_list) {
      target_genes <- orthogroups[[species]][i]
      
      if (is.na(target_genes) || target_genes == "") {
        next
      }
      
      genes_list <- trimws(strsplit(target_genes, ",")[[1]])
      
      for (gene_id in genes_list) {
        if (gene_id %in% genes_in_expression) {
          if (gene_to_family[gene_id] == 0) {
            gene_to_family[gene_id] <- family_index
          }
        }
      }
    }
  }
  
  return(list(
    gene_to_family = as.integer(gene_to_family),
    gene_expression = gene_expr,
    centroids = centroids
  ))
}

# Main performance comparison function
run_real_data_performance_comparison <- function() {
  cat("=== PERFORMANCE COMPARISON WITH REAL DATA ===\n")
  cat("Testing distance_to_centroid with ~88,000 genes\n\n")
  
  # Define file paths
  base_path <- "material"
  orthogroups_file <- file.path(base_path, "filtered_families.tsv")
  centroids_file <- file.path(base_path, "centroids_orthologs_only.tsv")  
  gene_expression_file <- file.path(base_path, "normalization.tsv")
  
  # Verify files exist
  for (file in c(orthogroups_file, centroids_file, gene_expression_file)) {
    if (!file.exists(file)) {
      stop(paste("Error: File not found:", file))
    }
  }
  
  # Load and prepare data
  mapping_data <- generate_gene_to_family_mapping(
    orthogroups_file, 
    centroids_file, 
    gene_expression_file,
    use_all_species = TRUE
  )
  
  # Prepare matrices for Fortran (transpose for column-major order)
  gene_expr_matrix <- as.matrix(mapping_data$gene_expression[, -1])
  gene_expr_matrix <- t(gene_expr_matrix)
  
  centroids_matrix <- as.matrix(mapping_data$centroids[, -1])
  centroids_matrix <- t(centroids_matrix)
  
  # Convert matrices to vectors (column-major format for Fortran)
  genes_vector <- as.vector(gene_expr_matrix)
  centroids_vector <- as.vector(centroids_matrix)
  d <- nrow(gene_expr_matrix)
  
  # Display data info
  n_genes <- ncol(gene_expr_matrix)
  n_families <- ncol(centroids_matrix)
  genes_with_families <- sum(mapping_data$gene_to_family > 0)
  
  cat("Data loaded:\n")
  cat("- Genes:", n_genes, "\n")
  cat("- Families:", n_families, "\n") 
  cat("- Dimensions:", d, "\n")
  cat("- Genes with family assignment:", genes_with_families, "\n")
  cat("- Gene expression vector size:", round(object.size(genes_vector)/1024/1024, 2), "MB\n")
  cat("- Centroids vector size:", round(object.size(centroids_vector)/1024/1024, 2), "MB\n\n")
  
  # Benchmark Rcpp (direct pointers)
  cat("=== Rcpp Test (with direct pointers) ===\n")
  time_rcpp <- system.time({
    distances_rcpp <- rcpp_distance_to_centroid(genes_vector, centroids_vector, mapping_data$gene_to_family, d)
  })
  cat("Rcpp completed\n")
  print(time_rcpp)
  cat("\n")
  
  # Benchmark .Fortran (with copies)
  cat("=== .Fortran Test (with data copies) ===\n")
  time_fortran <- system.time({
    distances_fortran <- fortran_distance_to_centroid(genes_vector, centroids_vector, mapping_data$gene_to_family, d)
  })
  cat(".Fortran completed\n")
  print(time_fortran)
  cat("\n")
  
  # Compare results and performance
  cat("=== Comparison ===\n")
  cat("Results match:", all.equal(distances_rcpp, distances_fortran), "\n")
  cat("Rcpp elapsed time:", time_rcpp[3], "seconds\n")
  cat(".Fortran elapsed time:", time_fortran[3], "seconds\n")
  speedup <- time_fortran[3] / time_rcpp[3]
  cat("Speedup (Rcpp vs .Fortran):", round(speedup, 2), "x\n\n")
  
  # Calculate memory overhead
  total_data_size <- object.size(genes_vector) + object.size(centroids_vector) + object.size(mapping_data$gene_to_family)
  cat("=== Memory Analysis ===\n")
  cat("Total input data size:", round(as.numeric(total_data_size)/1024/1024, 2), "MB\n")
  cat("Memory overhead with .Fortran copies:", round(as.numeric(total_data_size)/1024/1024 * 2, 2), "MB\n")
  cat("Rcpp memory advantage: No data copying, direct pointer access\n\n")
  
  # Show sample results
  results_with_families <- which(mapping_data$gene_to_family > 0)
  if (length(results_with_families) > 0) {
    sample_indices <- head(results_with_families, 5)
    cat("=== Sample Results (first 5 genes with families) ===\n")
    cat("Gene Index | Family | Distance (Rcpp) | Distance (.Fortran)\n")
    for (i in sample_indices) {
      cat(sprintf("%10d | %6d | %15.6f | %18.6f\n", 
                  i, mapping_data$gene_to_family[i], distances_rcpp[i], distances_fortran[i]))
    }
  }
  
  return(list(
    rcpp_time = time_rcpp,
    fortran_time = time_fortran,
    speedup = speedup,
    distances_rcpp = distances_rcpp,
    distances_fortran = distances_fortran
  ))
}

# Run the comparison
if (interactive() || length(commandArgs(trailingOnly = TRUE)) == 0) {
  cat("Running performance comparison with real data...\n\n")
  results <- run_real_data_performance_comparison()
} else {
  cat("Script loaded. Use run_real_data_performance_comparison() to run the comparison.\n")
}
