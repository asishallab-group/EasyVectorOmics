# Load libraries
library(dplyr)
library(readr)

# Load logFC data
logfc_df <- read.table("/storage/EasyVectorOmics/FastQ_GSE125483_JK/kallisto_outputs/combined_output.tsv", header = TRUE, sep = "\t")

# Load orthologs data
orthologs_df <- read.table("/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/tree/majority/complete_tree_conserved_orthologs_2_formatted_filtered.tsv", header = TRUE, sep = "\t")

# Function to compute centroid per orthogroup 
calculate_orthogroup_centroids <- function(logfc_df, orthologs_df) {
  orthogroups <- unique(orthologs_df$Orthogroup)
#   orthogroups <- head(orthogroups, max_groups)  # Limit to first N orthogroups
  
  centroids_list <- list()
  
  for (og in orthogroups) {
    # cat("\nProcessing Orthogroup:", og, "\n")
    
    # Get genes in this orthogroup
    genes <- orthologs_df %>%
      filter(Orthogroup == og) %>%
      select(Gene1, Gene2) %>%
      unlist() %>%
      unique()
    
    # Get matching logFC values
    og_logfc <- logfc_df %>%
      filter(GeneID %in% genes)
    
    if (nrow(og_logfc) > 0) {
    #   print("Ortholog vectors:")
    #   print(og_logfc)
      
      # Compute centroid
      centroid <- colMeans(og_logfc[, -1, drop = FALSE], na.rm = TRUE)
    #   print("Centroid:")
    #   print(centroid)
      
      centroids_list[[og]] <- centroid
    } else {
      cat("No matching genes found in logFC data for", og, "\n")
    }
  }
  
  centroids_df <- do.call(rbind, centroids_list)
  centroids_df <- as.data.frame(centroids_df)
  centroids_df$Orthogroup <- rownames(centroids_df)
  rownames(centroids_df) <- NULL
  centroids_df <- centroids_df[, c("Orthogroup", setdiff(colnames(centroids_df), "Orthogroup"))]
  
  return(centroids_df)
}

centroids <- calculate_orthogroup_centroids(logfc_df, orthologs_df)

# Save result
write.table(centroids, file = "results/orthogroup_centroids_majorityrule.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
