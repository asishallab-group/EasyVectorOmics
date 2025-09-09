source("r/tensoromics_functions_tox_data.R")

# ---- Example: replicate Fortran test logic in R ----

# Define your file lists (replace with your actual file paths)
files_6_replicates <- c(
  "material/kallisto_sex_data_Adipose.tsv",
  "material/kallisto_sex_data_Adrenal.tsv",
  "material/kallisto_sex_data_Colon.tsv",
  "material/kallisto_sex_data_Heart.tsv",
  "material/kallisto_sex_data_Liver.tsv",
  "material/kallisto_sex_data_Lung.tsv",
  "material/kallisto_sex_data_Muscle.tsv",
  "material/kallisto_sex_data_Skin.tsv",
  "material/kallisto_sex_data_Spleen.tsv",
  "material/kallisto_sex_data_Thyroid.tsv"
)
files_7_replicates <- c("material/kallisto_sex_data_Testis.tsv")
files_5_replicates <- c("material/kallisto_sex_data_Brain.tsv")
files_4_replicates <- c("material/kallisto_sex_data_Pituitary.tsv")

# Parameters
n_genes <- 88327
n_families <- 15512
gene_len <- 32
family_len <- 32
n_samples <- 10*6 + 7 + 5 + 4
value_cols_6 <- 2:7
value_cols_7 <- 2:8
value_cols_5 <- 2:6
value_cols_4 <- 2:5

# Allocate result matrices
kallisto_expr <- matrix(0, nrow = n_samples, ncol = n_genes)

# Read gene IDs from first file
cat("Reading gene IDs...\n")
res_gene <- read_gene_ids_from_file(files_6_replicates[1], n_genes, gene_len, n_header_rows = 1, gene_col = 1)
gene_ids <- res_gene$gene_ids
cat("ierr (gene ids):", res_gene$ierr, "\n")

# Read 6-replicate files
cat("Reading 6-replicate files...\n")
res_expr_6 <- read_expression_vectors(
  file_list = files_6_replicates,
  gene_ids = gene_ids,
  n_header_rows = 1,
  gene_col = 1,
  value_cols = value_cols_6,
  start_row = 1,
  delimiter = "\t",
  n_samples = 60
)
cat("ierr (6-replicates):", res_expr_6$ierr, "\n")
kallisto_expr[1:60, ] <- res_expr_6$expression_vectors

# Debug: Check dimensions and types before calling Fortran
cat("DEBUG: files_6_replicates length:", length(files_6_replicates), "\n")
cat("DEBUG: gene_ids length:", length(gene_ids), "\n")
cat("DEBUG: kallisto_expr dim:", dim(kallisto_expr), "\n")
cat("DEBUG: value_cols_6:", value_cols_6, "\n")
cat("DEBUG: n_samples for 6-replicates:", 60, "\n")
cat("DEBUG: gene_len:", gene_len, "\n")
cat("DEBUG: file_len:", max(nchar(files_6_replicates)), "\n")

# Read 7-replicate file
cat("Reading 7-replicate file...\n")
res_expr_7 <- read_expression_vectors(
  file_list = files_7_replicates,
  gene_ids = gene_ids,
  n_header_rows = 1,
  gene_col = 1,
  value_cols = value_cols_7,
  start_row = 1,
  delimiter = "\t",
  n_samples = 7
)
cat("ierr (7-replicates):", res_expr_7$ierr, "\n")
kallisto_expr[61:67, ] <- res_expr_7$expression_vectors

# Read 5-replicate file
cat("Reading 5-replicate file...\n")
res_expr_5 <- read_expression_vectors(
  file_list = files_5_replicates,
  gene_ids = gene_ids,
  n_header_rows = 1,
  gene_col = 1,
  value_cols = value_cols_5,
  start_row = 1,
  delimiter = "\t",
  n_samples = 5
)
cat("ierr (5-replicates):", res_expr_5$ierr, "\n")
kallisto_expr[68:72, ] <- res_expr_5$expression_vectors

# Read 4-replicate file
cat("Reading 4-replicate file...\n")
res_expr_4 <- read_expression_vectors(
  file_list = files_4_replicates,
  gene_ids = gene_ids,
  n_header_rows = 1,
  gene_col = 1,
  value_cols = value_cols_4,
  start_row = 1,
  delimiter = "\t",
  n_samples = 4
)
cat("ierr (4-replicates):", res_expr_4$ierr, "\n")
kallisto_expr[73:76, ] <- res_expr_4$expression_vectors

# Read family mapping
cat("Reading family file...\n")
res_family <- read_family_file("material/Orthogroups.tsv", gene_ids, n_families, family_len)
gene_family_ids <- res_family$family_ids
gene_to_fam <- res_family$gene_to_fam
cat("ierr (family):", res_family$ierr, "\n")
print(sample(gene_family_ids, 10))

# Filter out genes without family assignments
cat("Filtering unassigned genes...\n")
res_filter <- filter_unassigned_genes(gene_ids, kallisto_expr, gene_to_fam)
cat("ierr (filter):", res_filter$ierr, "\n")
gene_ids <- res_filter$gene_ids
kallisto_expr <- res_filter$expression_vectors
gene_to_fam <- res_filter$gene_to_fam
n_genes_kept <- res_filter$n_genes_kept
cat("n_genes_kept:", n_genes_kept, "\n")

# ...continue with further steps as needed...

