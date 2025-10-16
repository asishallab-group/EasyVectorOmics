source("r/tensoromics_functions_tox_data.R")

# ---- Example: replicate Fortran test logic in R ----

# Define your file lists (replace with your actual file paths)
file <- c("material/kallisto_sex_data_no_na.tsv")

# Parameters
n_genes <- 88327
n_families <- 15512
gene_len <- 32
family_len <- 32
n_samples <- 67
value_cols <- 2:68

# Allocate result matrices
kallisto_expr <- matrix(0, nrow = n_samples, ncol = n_genes)

# Read gene IDs from first file
cat("Reading gene IDs...\n")
res_gene <- read_gene_ids_from_file(file, n_genes, gene_len, n_header_rows = 1, gene_col = 1)
gene_ids <- res_gene$gene_ids
cat("ierr (gene ids):", res_gene$ierr, "\n")

# Read 6-replicate files
cat("Reading expression file...\n")
res_expr_6 <- read_expression_vectors(
  file_list = file,
  gene_ids = gene_ids,
  n_header_rows = 1,
  gene_col = 1,
  value_cols = value_cols,
  delimiter = "\t",
  n_samples = 67
)
cat("ierr :", res_expr_6$ierr, "\n")
kallisto_expr[1:67, ] <- res_expr_6$expression_vectors

# Debug: Check dimensions and types before calling Fortran
cat("DEBUG: gene_ids length:", length(gene_ids), "\n")
cat("DEBUG: kallisto_expr dim:", dim(kallisto_expr), "\n")
cat("DEBUG: value_cols:", value_cols, "\n")
cat("DEBUG: n_samples:", 67, "\n")
cat("DEBUG: gene_len:", gene_len, "\n")
cat("DEBUG: file_len:", max(nchar(file)), "\n")

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

# Data validation
cat("Starting data validation...\n")

# Set parameters for validation
n_genes <- n_genes_kept
n_samples <- nrow(kallisto_expr)

# Validate individual components
cat("Validating gene IDs uniqueness...\n")
validate_gene_ids_uniqueness(gene_ids, n_genes)

cat("Validating family IDs uniqueness...\n")
validate_family_ids_uniqueness(gene_family_ids, n_families)

cat("Validating gene-to-family mapping...\n")
validate_gene_to_family_mapping(gene_to_fam, n_genes, n_families)

cat("Validating expression data...\n")
validate_expression_data(kallisto_expr, n_genes, n_samples)

ortholog_set <- rep(TRUE, n_genes)
family_centroids <- tox_group_centroid(kallisto_expr, gene_to_fam, n_families, ortholog_set, mode='all')

cat(nrow(family_centroids), "x", ncol(family_centroids), "\n")

cat("Validating family centroids...\n")
validate_family_centroids(family_centroids, n_families, n_samples)

res <- tox_compute_shift_vector_field(kallisto_expr, family_centroids, gene_to_fam)
shift_vectors <- matrix(res$shift_vectors, nrow=2*n_samples, ncol=n_genes)

cat("Validating shift vectors...\n")
validate_shift_vectors(shift_vectors, kallisto_expr, family_centroids, gene_to_fam, n_samples)

# Comprehensive validation
cat("Performing comprehensive data validation...\n")
validate_all_data(
  n_genes, n_families, n_samples,
  gene_ids, gene_family_ids,
  gene_to_fam, kallisto_expr,
  family_centroids, shift_vectors
)

cat("All validations passed successfully!\n")

# Additional diagnostic output
cat("\n=== DATA SUMMARY ===\n")
cat("Samples:", n_samples, "\n")
cat("Genes:", n_genes, "\n")
cat("Families:", n_families, "\n")
cat("Expression matrix dimensions:", dim(kallisto_expr), "\n")
cat("Gene ID examples:", head(gene_ids), "\n")
cat("Family ID examples:", head(gene_family_ids), "\n")
cat("Gene-to-family mapping examples:", head(gene_to_fam), "\n")

cat("===Archive tests===\n")
save_tox_data(zip_filename="test_archive_1_R.zip", gene_ids=gene_ids, gene_ids_name="gene_ids_v1_R.bin",
                      expression_vectors=kallisto_expr, expression_vectors_name="expr_vecs_v1_R.bin", 
                      gene_to_fam=gene_to_fam, gene_to_fam_name="gene_to_fam_v1_R.bin",
                      family_ids=gene_family_ids, family_ids_name="family_ids_v1_R.bin", 
                      family_centroids=family_centroids, family_centroids_name="family_centroids_v1_R.bin",
                      shift_vectors=shift_vectors, shift_vectors_name="shift_vectors_v1_R.bin")

save_tox_data(zip_filename="test_archive_2_R.zip", gene_ids=gene_ids, gene_ids_name="gene_ids_v2_R.bin",
                      expression_vectors=kallisto_expr, expression_vectors_name="expr_vecs_v2_R.bin")

save_tox_data(zip_filename="test_archive_3_R.zip", gene_ids=gene_ids, gene_ids_name="gene_ids_v3_R.bin",
                      expression_vectors=kallisto_expr, expression_vectors_name="expr_vecs_v3_R.bin", gene_to_fam=gene_to_fam)

save_tox_data(zip_filename="test_archive_4_R.zip", family_centroids=family_centroids, family_centroids_name="family_centroids_v1_R.bin",
                      shift_vectors=shift_vectors, shift_vectors_name="shift_vectors_v1_R.bin")

save_tox_data(zip_filename="test_archive_5_R.zip")

result_1 <- read_tox_data(zip_filename="test_archive_1_R.zip", 
                                  gene_ids=TRUE,
                                  expression_vectors=TRUE,
                                  gene_to_fam=TRUE,
                                  family_ids=TRUE,
                                  family_centroids=TRUE,
                                  shift_vectors=TRUE)

result_2 <- read_tox_data(zip_filename="test_archive_1_R.zip",
                                  family_centroids=TRUE,
                                  shift_vectors=TRUE)

result_3 <- read_tox_data(zip_filename="test_archive_2_R.zip",
                                  gene_ids=TRUE,
                                  expression_vectors=TRUE,
                                  shift_vectors=TRUE)
tryCatch({
  result_py <- read_tox_data(zip_filename="test_archive_1_py.zip", 
                                  gene_ids=TRUE,
                                  expression_vectors=TRUE,
                                  gene_to_fam=TRUE,
                                  family_ids=TRUE,
                                  family_centroids=TRUE,
                                  shift_vectors=TRUE)
  cat("Successfully read from python archive")

  result_f <- read_tox_data(zip_filename="test_archive_1_f.zip", 
                                  gene_ids=TRUE,
                                  expression_vectors=TRUE,
                                  gene_to_fam=TRUE,
                                  family_ids=TRUE,
                                  family_centroids=TRUE,
                                  shift_vectors=TRUE)
  cat("Successfully read from fortran archive")
}, error = function(e) {
  cat("Cross-platform test skipped:", e$message, "\n")
})
cat("=== Testing create_zip_archive directly with dummy arrays ===\n")

# Erstelle Dummy-Daten
dummy_int_1d <- 1:10
dummy_int_2d <- matrix(1:12, nrow=3, ncol=4)
dummy_real_1d <- c(1.5, 2.5, 3.5, 4.5, 5.5)
dummy_real_2d <- matrix(rnorm(12), nrow=3, ncol=4)
dummy_char_1d <- c("apple", "banana", "cherry", "date", "elderberry")

# Serialisiere Dummy-Daten zu Dateien
cat("Serializing dummy arrays...\n")
tox_serialize_int_array(dummy_int_1d, "dummy_int_1d.bin")
tox_serialize_int_array(as.vector(dummy_int_2d), "dummy_int_2d.bin")  # Matrix zu Vektor
tox_serialize_real_array(dummy_real_1d, "dummy_real_1d.bin")
tox_serialize_real_array(dummy_real_2d, "dummy_real_2d.bin")
tox_serialize_char_array(dummy_char_1d, "dummy_char_1d.bin")

# Test 1: Erstelle Archiv mit allen Dummy-Dateien
cat("Test 1: Creating archive with all dummy files...\n")
create_zip_archive("test_dummy_all.zip",
                  keys = c("int_1d", "int_2d", "real_1d", "real_2d", "char_1d"),
                  filenames = c("dummy_int_1d.bin", "dummy_int_2d.bin", 
                                "dummy_real_1d.bin", "dummy_real_2d.bin", 
                                "dummy_char_1d.bin"))

# Test 2: Erstelle Archiv mit gemischten echten und Dummy-Daten
cat("Test 2: Creating archive with mixed real and dummy data...\n")
# Serialisiere einige echte Daten für diesen Test
tox_serialize_char_array(gene_ids[1:5], "sample_gene_ids.bin")
create_zip_archive("test_mixed_data.zip",
                  keys = c("sample_genes", "dummy_int", "dummy_real"),
                  filenames = c("sample_gene_ids.bin", "dummy_int_1d.bin", "dummy_real_1d.bin"))

# Test 3: Erstelle Archiv mit benutzerdefinierten Keys
cat("Test 3: Creating archive with custom keys...\n")
create_zip_archive("test_custom_keys.zip",
                  keys = c("experiment_config", "results_summary", "raw_data"),
                  filenames = c("dummy_char_1d.bin", "dummy_real_2d.bin", "dummy_int_2d.bin"))

# Test 4: Erstelle Archiv mit nur einer Datei
cat("Test 4: Creating archive with single file...\n")
create_zip_archive("test_single_file.zip",
                  keys = c("single_data"),
                  filenames = c("dummy_int_1d.bin"))

# Test 5: Teste Fehlerbehandlung mit ungleichen Arrays
cat("Test 5: Testing error handling with mismatched arrays...\n")
tryCatch({
  create_zip_archive("test_error.zip",
                    keys = c("key1", "key2"),
                    filenames = c("file1.bin"))  # Fehler: ungleiche Länge
}, error = function(e) {
  cat("Expected error caught:", e$message, "\n")
})

# Überprüfe ob Archive erstellt wurden
cat("Checking created archives...\n")
archives <- c("test_dummy_all.zip", "test_mixed_data.zip", 
              "test_custom_keys.zip", "test_single_file.zip")

for (archive in archives) {
  if (file.exists(archive)) {
    cat("✓ Archive created:", archive, "\n")
  } else {
    cat("✗ Archive missing:", archive, "\n")
  }
}

# Räume temporäre Dateien auf
cat("Cleaning up temporary files...\n")
temp_files <- c("dummy_int_1d.bin", "dummy_int_2d.bin", "dummy_real_1d.bin", 
                "dummy_real_2d.bin", "dummy_char_1d.bin", "sample_gene_ids.bin")

for (temp_file in temp_files) {
  if (file.exists(temp_file)) {
    file.remove(temp_file)
    cat("Removed:", temp_file, "\n")
  }
}

# 3. Lesen der Archive (bestehende Tests)
cat("=== Reading archives ===\n")
result_1 <- read_tox_data(zip_filename="test_archive_1_R.zip", 
                                  gene_ids=TRUE,
                                  expression_vectors=TRUE,
                                  gene_to_fam=TRUE,
                                  family_ids=TRUE,
                                  family_centroids=TRUE,
                                  shift_vectors=TRUE)

result_2 <- read_tox_data(zip_filename="test_archive_1_R.zip",
                                  family_centroids=TRUE,
                                  shift_vectors=TRUE)

result_3 <- read_tox_data(zip_filename="test_archive_2_R.zip",
                                  gene_ids=TRUE,
                                  expression_vectors=TRUE,
                                  shift_vectors=TRUE)

# 4. Teste Lesen der Dummy-Archive (falls read_tox_data generisch genug ist)
cat("=== Testing read of dummy archives ===\n")
tryCatch({
  # Versuche die Dummy-Archive zu lesen
  dummy_result <- read_tox_data(zip_filename="test_dummy_all.zip")
  cat("Successfully read dummy archive (structure only)\n")
}, error = function(e) {
  cat("Note: read_tox_data may not handle non-standard keys:", e$message, "\n")
})

# 5. Cross-platform Tests (bestehende Tests)
tryCatch({
  result_py <- read_tox_data(zip_filename="test_archive_1_py.zip", 
                                  gene_ids=TRUE,
                                  expression_vectors=TRUE,
                                  gene_to_fam=TRUE,
                                  family_ids=TRUE,
                                  family_centroids=TRUE,
                                  shift_vectors=TRUE)
  cat("Successfully read from python archive\n")

  result_f <- read_tox_data(zip_filename="test_archive_1_f.zip", 
                                  gene_ids=TRUE,
                                  expression_vectors=TRUE,
                                  gene_to_fam=TRUE,
                                  family_ids=TRUE,
                                  family_centroids=TRUE,
                                  shift_vectors=TRUE)
  cat("Successfully read from fortran archive\n")
}, error = function(e) {
  cat("Cross-platform test skipped:", e$message, "\n")
})

cat("=== All tests completed! ===\n")

# Zusätzliche Validierung: Prüfe Archiv-Inhalte
cat("=== Archive content validation ===\n")
library(utils)

# Liste alle erstellten Archive auf
archive_files <- list.files(pattern = "\\.zip$")
cat("Created archives:\n")
for (arch in archive_files) {
  file_info <- file.info(arch)
  cat(sprintf("  %s (%.2f KB)\n", arch, file_info$size/1024))
  
  # Versuche Archiv-Inhalt anzuzeigen (nur falls unzip verfügbar)
  tryCatch({
    content <- system(paste("unzip -l", arch), intern = TRUE)
    cat("    Contents:\n")
    for (line in tail(content, 5)) {  # Zeige letzte 5 Zeilen (ohne Header)
      if (grepl("\\.bin$|manifest\\.txt$", line)) {
        cat("     ", line, "\n")
      }
    }
  }, error = function(e) {
    # Ignoriere Fehler falls unzip nicht verfügbar
  })
}