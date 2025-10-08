import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from tensoromics_functions_tox_data import (
    read_expression_vectors,
    read_family_file,
    filter_unassigned_genes,
    read_gene_ids_from_file,
    validate_all_data,
    validate_data_structure,
    validate_expression_data,
    validate_family_centroids,
    validate_gene_ids_uniqueness,
    validate_family_ids_uniqueness,
    validate_gene_to_family_mapping,
    validate_shift_vectors,
    save_tox_data,
    read_tox_data
)
from tensoromics_functions import (
    tox_group_centroid,
    tox_compute_shift_vector_field
)

# ---- Example: replicate Fortran test logic in Python ----

# Define your file lists (replace with your actual file paths)
file = ["material/kallisto_sex_data_no_na.tsv"]
# Parameters
n_genes = 88327
n_families = 15512
gene_len = 32
family_len = 32
n_samples = 67
value_cols = np.arange(2, 68)

# Allocate result matrices
kallisto_expr = np.zeros((n_samples, n_genes), dtype=np.float64, order='F')

# Read gene IDs from first file
print("Reading gene IDs...")
gene_ids = read_gene_ids_from_file(file, n_genes, gene_len, n_header_rows=1, gene_col=1)
print("gene_ids sample:", gene_ids[:5])
print("gene_ids type:", type(gene_ids))

# Read 6-replicate files
print("Reading 6-replicate files...")
kallisto_expr = read_expression_vectors(
    file_list=file,
    gene_ids=gene_ids,  # Now accepts numpy array
    n_samples=67,
    n_header_rows=1,
    gene_col=1,
    value_cols=value_cols,
    delimiter="\t"
)
print("expression shape:", kallisto_expr.shape)

# Read family mapping
print("Reading family file...")
family_result = read_family_file("material/Orthogroups.tsv", gene_ids, family_len, n_families)
family_ids = family_result['family_ids']
gene_to_fam = family_result['gene_to_fam']
print("family_ids sample:", family_ids[:5])
print("family_ids type:", type(family_ids))
print("gene_to_fam sample:", gene_to_fam[:10])

# Filter out genes without family assignments
print("Filtering unassigned genes...")
filter_result = filter_unassigned_genes(gene_ids, gene_to_fam)
mask = filter_result['mask']
n_genes_kept = filter_result['n_genes_kept']
print("n_genes_kept:", n_genes_kept)
print("mask sample:", mask[:10])

print("First gene expression: ", kallisto_expr[:, 0])

# Filter arrays using mask
filtered_gene_ids = gene_ids[mask.astype(bool)]
filtered_kallisto_expr = kallisto_expr[:, mask.astype(bool)]
filtered_gene_to_fam = gene_to_fam[mask.astype(bool)]

print("filtered_gene_ids sample:", filtered_gene_ids[:5])
print("filtered_kallisto_expr shape:", filtered_kallisto_expr.shape)
print("filtered_gene_to_fam sample:", filtered_gene_to_fam[:10])

validate_gene_ids_uniqueness(gene_ids)
validate_family_ids_uniqueness(family_ids)
validate_expression_data(kallisto_expr, True)

ortholog_set = np.array([True for i in range(n_genes_kept)])

centroids = tox_group_centroid(filtered_kallisto_expr, filtered_gene_to_fam, n_families, ortholog_set)

print("Validating centroids...")
validate_family_centroids(centroids)
print("Centroids validated")

shift_vectors_result = tox_compute_shift_vector_field(filtered_kallisto_expr, centroids, filtered_gene_to_fam)
shift_vectors = shift_vectors_result["shift_vectors"]

print("Validating shift vectors...")
validate_shift_vectors(shift_vectors, filtered_kallisto_expr, centroids, filtered_gene_to_fam, 
                      n_genes_kept, n_samples, n_families)
print("Shift vectors validated")

validate_gene_to_family_mapping(filtered_gene_to_fam, n_families)

print("Number of genes kept after filtering:", n_genes_kept)
print("gene_family_size:", len(filtered_gene_to_fam))
print("filtered gene ids size:", len(filtered_gene_ids))

print("Validating data structure...")
validate_data_structure(n_genes_kept, n_families, n_samples, filtered_gene_ids, 
                        family_ids, filtered_gene_to_fam, filtered_kallisto_expr, centroids, shift_vectors)
print("Data structure validated")

print("Validating all data...")
validate_all_data(n_genes_kept, n_families, n_samples, filtered_gene_ids, family_ids, 
                  filtered_gene_to_fam, filtered_kallisto_expr, centroids, shift_vectors)
print("All data validated")

print("Testing save function...")
save_tox_data("test_archive_1_py.zip", gene_ids=filtered_gene_ids, gene_ids_name="gene_ids_v1.bin",
                      expression_vectors=filtered_kallisto_expr, expression_vectors_name="kallisto_data_v1.bin")

save_tox_data("test_archive_2_py.zip", family_centroids=centroids, family_centroids_name="centroids.bin")
try:
    save_tox_data("test_archive_3_py.zip", family_centroids=centroids, gene_ids=kallisto_expr)
except:
    print("Successfully threw exception")

save_tox_data("test_archive_4_py.zip", gene_ids=filtered_gene_ids, gene_ids_name="gene_ids_v1.bin",
                      expression_vectors=filtered_kallisto_expr, expression_vectors_name="kallisto_data_v1.bin",
                      gene_to_fam=filtered_gene_to_fam, gene_to_fam_name="gene_to_fam_v1.bin",
                      family_ids=family_ids, family_ids_name="family_ids_v1.bin",
                      family_centroids=centroids, family_centroids_name="family_centroids_v1.bin",
                      shift_vectors=shift_vectors, shift_vectors_name="shift_vectors_v1.bin")

result_1 = read_tox_data("test_archive_4_py.zip", load_gene_ids=True, load_expression_vectors=True, 
                                 load_gene_to_fam=True, load_family_ids=True, load_family_centroids=True, load_shift_vectors=True)
result_2 = read_tox_data("test_archive_4_py.zip", load_gene_ids=True, load_gene_to_fam=True)
try:
    result_3 = read_tox_data("test_archive_3_py.zip", load_gene_ids=True)
except:
    result_3 = None
    print("Successfully threw exception")
try:
    result_f = read_tox_data("test_archive_1_f.zip", True, True, True, True, True, True)
except:
    result_f = None 
    print("Could not read fortran archive. Does the file exist?")

try: 
    result_r = read_tox_data(zip_filename="test_archive_1_R.zip", gene_ids=True, expression_vectors=True, 
                                 gene_to_fam=True, family_ids=True, family_centroids=True, shift_vectors=True)
except:
    print("Could not read R archive. Does the file exist?")