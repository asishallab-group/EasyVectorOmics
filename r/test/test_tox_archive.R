source("r/tensoromics_functions_tox_data.R")

gene_ids <- c("Gene1", "Gene2", "Gene3")
n_genes <- 3
n_samples <- 1
expression_vectors <- c(c(1.0), c(2.0), c(3.0))
gene_to_fam <- c(1, 2, 3)
family_ids <- c("Fam1", "Fam2", "Fam3")
family_centroids <- c((1.0), (2.0), c(3.0))
shift_vectors <- c((0.0), c(0.0), c(0.0))

save_tox_data("test_archive_1_R.zip", gene_ids=gene_ids, gene_ids_name="gene_ids_v1.bin")
result <- read_tox_data("test_archive_1_R.zip", gene_ids = TRUE)
result2 <- read_tox_data("test_archive_1_R.zip", gene_ids = TRUE, family_ids=TRUE)

print(result2)