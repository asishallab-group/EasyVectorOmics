# Example of using the Fortran module with real data
# This script demonstrates how to use distance_to_centroid with normalized TPM data,
# orthogroup centroids and gene-to-family mapping

library(readr)

# Function to generate gene_to_family mapping from Orthogroups.tsv file
generate_gene_to_family_mapping <- function(orthogroups_file, centroids_file, gene_expression_file, 
                                           use_all_species = TRUE, target_species = NULL) {
  cat("Loading data files...\n")
  
  # Read orthogroups file
  orthogroups <- read_tsv(orthogroups_file, show_col_types = FALSE)
  
  # Read centroids file
  centroids <- read_tsv(centroids_file, show_col_types = FALSE)
  
  # Read gene expression file
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
  gene_id_col <- names(gene_expr)[1]  # Use first column as gene_id
  genes_in_expression <- gene_expr[[gene_id_col]]
  
  # Create gene_id to family index mapping
  gene_to_family <- rep(0, length(genes_in_expression))  # 0 indicates no assignment
  names(gene_to_family) <- genes_in_expression
  
  # Create orthogroup to numeric index mapping
  orthogroup_col <- names(centroids)[1]  # Use first column as Orthogroup
  orthogroup_to_index <- setNames(1:nrow(centroids), centroids[[orthogroup_col]])
  
  # Generate gene-to-family mapping
  for (i in 1:nrow(orthogroups)) {
    orthogroup_id <- orthogroups$Orthogroup[i]
    
    # Check if this orthogroup has a centroid
    if (!orthogroup_id %in% names(orthogroup_to_index)) {
      next
    }
    
    family_index <- orthogroup_to_index[orthogroup_id]
    
    # Iterate over all target species
    for (species in target_species_list) {
      target_genes <- orthogroups[[species]][i]
      
      if (is.na(target_genes) || target_genes == "") {
        next
      }
      
      # Split genes by comma and clean spaces
      genes_list <- trimws(strsplit(target_genes, ",")[[1]])
      
      # Assign family index to each gene in expression file
      for (gene_id in genes_list) {
        if (gene_id %in% genes_in_expression) {
          # Only assign if not previously assigned (avoid overwriting)
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

# Function to call Fortran code
call_distance_to_centroid <- function(gene_expression_matrix, centroids_matrix, gene_to_family_vector) {
  # Verify dimensions
  n_genes <- ncol(gene_expression_matrix)
  n_families <- ncol(centroids_matrix)
  d <- nrow(gene_expression_matrix)
  
  # Verify dimension compatibility
  if (length(gene_to_family_vector) != n_genes) {
    stop("Error: gene_to_family length doesn't match number of genes")
  }
  
  if (nrow(centroids_matrix) != d) {
    stop("Error: centroid dimensions don't match gene expression")
  }
  
  # Prepare output vector
  distances <- numeric(n_genes)
  
  # Call Fortran function
  result <- .Fortran("distance_to_centroid_r",
                     n_genes = as.integer(n_genes),
                     n_families = as.integer(n_families),
                     genes = as.double(gene_expression_matrix),
                     centroids = as.double(centroids_matrix),
                     gene_to_fam = as.integer(gene_to_family_vector),
                     distances = as.double(distances),
                     d = as.integer(d))
  
  return(result$distances)
}

# Main function to run the complete example
run_real_data_example <- function() {
  cat("=== EUCLIDEAN DISTANCE CALCULATION WITH REAL DATA ===\n\n")
  
  # Define file paths
  base_path <- "results/sex_data/"
  orthogroups_file <- file.path(base_path, "filtered_families.tsv")
  centroids_file <- file.path(base_path, "centroids_orthologs_only_2.tsv")  
  gene_expression_file <- file.path(base_path, "normalization_2.tsv")
  
  # Verify files exist
  for (file in c(orthogroups_file, centroids_file, gene_expression_file)) {
    if (!file.exists(file)) {
      stop(paste("Error: File not found:", file))
    }
  }
  
  # Load Fortran library
  dyn.load("build/libtensor-omics.so")
  
  # Generate mapping and load data
  mapping_data <- generate_gene_to_family_mapping(
    orthogroups_file, 
    centroids_file, 
    gene_expression_file,
    use_all_species = TRUE
  )
  
  # Prepare matrices for Fortran (transpose for column-major order)
  gene_expr_matrix <- as.matrix(mapping_data$gene_expression[, -1])  # Exclude gene_id column
  gene_expr_matrix <- t(gene_expr_matrix)  # Transpose for Fortran column-major format
  
  centroids_matrix <- as.matrix(mapping_data$centroids[, -1])  # Exclude Orthogroup column
  centroids_matrix <- t(centroids_matrix)  # Transpose for Fortran column-major format
  
  # Call Fortran function
  distances <- call_distance_to_centroid(gene_expr_matrix, centroids_matrix, mapping_data$gene_to_family)
  
  # Create results dataframe
  results_df <- data.frame(
    gene_id = mapping_data$gene_expression$gene_id,
    family_index = mapping_data$gene_to_family,
    distance_to_centroid = distances,
    has_family = mapping_data$gene_to_family > 0,
    stringsAsFactors = FALSE
  )
  
  # Filter genes with assigned families
  results_with_families <- results_df[results_df$has_family, ]
  
  # Save results
  output_file <- file.path("results/distance_to_centroids_fortran.tsv")
  write_tsv(results_with_families, output_file)
  
  cat("Results saved to:", output_file, "\n")
  cat("Genes with families:", nrow(results_with_families), "\n")
  
  return(results_df)
}

# Run example if script is executed directly
if (interactive() || length(commandArgs(trailingOnly = TRUE)) == 0) {
  cat("Running example with real data...\n")
  results <- run_real_data_example()
} else {
  cat("Script loaded. Use run_real_data_example() to run the example.\n")
}
