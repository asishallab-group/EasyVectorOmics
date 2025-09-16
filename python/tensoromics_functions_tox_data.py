"""
TensorOmics Functions Module
Python wrapper functions for Fortran routines via C interface
"""

import numpy as np
import ctypes
import os

# Load library
dll_path = os.path.abspath("build/libtensor-omics.so")
ctypes.CDLL("libgomp.so.1", mode=ctypes.RTLD_GLOBAL)
lib = ctypes.CDLL(dll_path)

# Update the function signature for read_gene_ids_from_file_C
lib.read_gene_ids_from_file_C.argtypes = [
    ctypes.POINTER(ctypes.c_int),  # filename_ascii (as pointer)
    ctypes.c_int,                  # fn_len
    ctypes.POINTER(ctypes.c_int),  # gene_ids_ascii (as pointer)
    ctypes.c_int,                  # gene_ids_len
    ctypes.c_int,                  # n_genes
    ctypes.c_int,                  # n_header_rows
    ctypes.c_int,                  # gene_col
    ctypes.POINTER(ctypes.c_int)   # ierr
]
lib.read_gene_ids_from_file_C.restype = None

def read_gene_ids_from_file(filename, n_genes, gene_ids_len, n_header_rows, gene_col):
    # Convert filename to ASCII array
    fn_ascii = string_to_ascii_array(filename, len(filename))
    
    # Initialize gene_ids_ascii with spaces (ASCII 32)
    gene_ids_ascii = np.full((gene_ids_len, n_genes), 32, dtype=np.int32, order='F')
    
    ierr = ctypes.c_int()
    
    # Call C function
    lib.read_gene_ids_from_file_C(
        fn_ascii.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ctypes.c_int(len(filename)),
        gene_ids_ascii.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ctypes.c_int(gene_ids_len),
        ctypes.c_int(n_genes),
        ctypes.c_int(n_header_rows),
        ctypes.c_int(gene_col),
        ctypes.byref(ierr)
    )
    
    if ierr.value != 0:
        raise Exception(f"Error reading gene IDs: {ierr.value}")
    
    # Convert result back to strings
    gene_ids = []
    for i in range(n_genes):
        gene_ids.append(ascii_array_to_string(gene_ids_ascii[:, i]))
    
    return gene_ids

# read_expression_vectors_C
lib.read_expression_vectors_C.argtypes = [
    ctypes.POINTER(ctypes.c_int),  # file_list_ascii
    ctypes.c_int,                  # file_list_len
    ctypes.c_int,                  # n_files
    ctypes.POINTER(ctypes.c_int),  # gene_ids_ascii
    ctypes.c_int,                  # gene_ids_len
    ctypes.c_int,                  # n_genes
    ctypes.POINTER(ctypes.c_double), # expression_vectors_flat
    ctypes.c_int,                  # n_samples
    ctypes.c_int,                  # n_genes2
    ctypes.c_int,                  # n_header_rows
    ctypes.c_int,                  # gene_col
    ctypes.POINTER(ctypes.c_int),  # value_cols
    ctypes.c_int,                  # n_value_cols
    ctypes.c_int,                  # start_row
    ctypes.POINTER(ctypes.c_int),  # ierr
    ctypes.POINTER(ctypes.c_int),  # delimiter_ascii
    ctypes.c_int                   # dlen
]
lib.read_expression_vectors_C.restype = None

# read_family_file_C
lib.read_family_file_C.argtypes = [
    ctypes.POINTER(ctypes.c_int),  # filename_ascii
    ctypes.c_int,                  # fn_len
    ctypes.POINTER(ctypes.c_int),  # gene_ids_ascii
    ctypes.c_int,                  # gene_ids_len
    ctypes.c_int,                  # n_genes
    ctypes.POINTER(ctypes.c_int),  # family_ids_ascii
    ctypes.c_int,                  # family_ids_len
    ctypes.c_int,                  # n_families
    ctypes.POINTER(ctypes.c_int),  # gene_to_fam
    ctypes.POINTER(ctypes.c_int)   # ierr
]
lib.read_family_file_C.restype = None

# filter_unassigned_genes_C
lib.filter_unassigned_genes_C.argtypes = [
    ctypes.POINTER(ctypes.c_int),  # gene_ids_ascii
    ctypes.c_int,                  # gene_ids_len
    ctypes.c_int,                  # n_genes
    ctypes.POINTER(ctypes.c_double), # expression_vectors_flat
    ctypes.c_int,                  # n_samples
    ctypes.POINTER(ctypes.c_int),  # gene_to_fam
    ctypes.POINTER(ctypes.c_int),  # mask
    ctypes.POINTER(ctypes.c_int),  # n_genes_kept
    ctypes.POINTER(ctypes.c_int)   # ierr
]
lib.filter_unassigned_genes_C.restype = None

# --- Validation C function signatures ---

lib.validate_gene_to_family_mapping_C.argtypes = [
    ctypes.POINTER(ctypes.c_int), # gene_to_fam
    ctypes.c_int,                 # n_genes
    ctypes.c_int,                 # n_families
    ctypes.POINTER(ctypes.c_int)  # ierr
]
lib.validate_gene_to_family_mapping_C.restype = None

lib.validate_expression_data_C.argtypes = [
    ctypes.POINTER(ctypes.c_double), # expression_vectors
    ctypes.c_int,                    # n_genes
    ctypes.c_int,                    # n_samples
    ctypes.c_int,                    # check_non_negative
    ctypes.POINTER(ctypes.c_int)     # ierr
]
lib.validate_expression_data_C.restype = None

lib.validate_family_centroids_C.argtypes = [
    ctypes.POINTER(ctypes.c_double), # family_centroids
    ctypes.c_int,                    # n_families
    ctypes.c_int,                    # d
    ctypes.POINTER(ctypes.c_int)     # ierr
]
lib.validate_family_centroids_C.restype = None

lib.validate_shift_vectors_C.argtypes = [
    ctypes.POINTER(ctypes.c_double), # shift_vectors
    ctypes.POINTER(ctypes.c_double), # expression_vectors
    ctypes.POINTER(ctypes.c_double), # family_centroids
    ctypes.POINTER(ctypes.c_int),    # gene_to_fam
    ctypes.c_int,                    # d
    ctypes.c_int,                    # n_genes
    ctypes.c_int,                    # n_samples
    ctypes.c_int,                    # n_families
    ctypes.POINTER(ctypes.c_int)     # ierr
]
lib.validate_shift_vectors_C.restype = None

lib.validate_gene_ids_uniqueness_C.argtypes = [
    ctypes.POINTER(ctypes.c_int), # gene_ids_ascii
    ctypes.c_int,                 # gene_ids_len
    ctypes.c_int,                 # n_genes
    ctypes.POINTER(ctypes.c_int)  # ierr
]
lib.validate_gene_ids_uniqueness_C.restype = None

lib.validate_family_ids_uniqueness_C.argtypes = [
    ctypes.POINTER(ctypes.c_int), # family_ids_ascii
    ctypes.c_int,                 # fam_len
    ctypes.c_int,                 # n_families
    ctypes.POINTER(ctypes.c_int)  # ierr
]
lib.validate_family_ids_uniqueness_C.restype = None

lib.validate_data_structure_C.argtypes = [
    ctypes.c_int,                 # n_genes
    ctypes.c_int,                 # n_families
    ctypes.c_int,                 # n_samples
    ctypes.c_int,                 # d
    ctypes.POINTER(ctypes.c_int), # gene_ids_ascii
    ctypes.c_int,                 # gene_ids_len
    ctypes.POINTER(ctypes.c_int), # gene_family_ids_ascii
    ctypes.c_int,                 # fam_len
    ctypes.POINTER(ctypes.c_int), # gene_to_fam
    ctypes.POINTER(ctypes.c_double), # expression_vectors
    ctypes.POINTER(ctypes.c_double), # family_centroids
    ctypes.POINTER(ctypes.c_double), # shift_vectors
    ctypes.POINTER(ctypes.c_int)     # ierr
]
lib.validate_data_structure_C.restype = None

lib.validate_all_data_C.argtypes = [
    ctypes.c_int,                 # n_genes
    ctypes.c_int,                 # n_families
    ctypes.c_int,                 # n_samples
    ctypes.c_int,                 # d
    ctypes.POINTER(ctypes.c_int), # gene_ids_ascii
    ctypes.c_int,                 # gene_len
    ctypes.POINTER(ctypes.c_int), # gene_family_ids_ascii
    ctypes.c_int,                 # fam_len
    ctypes.POINTER(ctypes.c_int), # gene_to_fam
    ctypes.POINTER(ctypes.c_double), # expression_vectors
    ctypes.POINTER(ctypes.c_double), # family_centroids
    ctypes.POINTER(ctypes.c_double), # shift_vectors
    ctypes.POINTER(ctypes.c_int)     # ierr
]
lib.validate_all_data_C.restype = None

# Hilfsfunktionen für String-Konvertierung
def string_to_ascii_array(s, length):
    """Konvertiert einen String in ein ASCII-Integer-Array mit fester Länge"""
    ascii_array = np.zeros(length, dtype=np.int32, order='F')
    for i, char in enumerate(s):
        if i >= length:
            break
        ascii_array[i] = ord(char)
    return ascii_array

def ascii_array_to_string(ascii_array):
    """Konvertiert ein ASCII-Integer-Array zurück in einen String"""
    return ''.join(chr(x) for x in ascii_array if x != 0).strip()

# Funktion für read_expression_vectors_C
def read_expression_vectors(file_list, gene_ids, n_samples, n_header_rows, 
                           gene_col, value_cols, start_row, delimiter='\t'):
    # Konvertiere Eingaben in ASCII-Arrays
    file_list_ascii = np.zeros((max(len(f) for f in file_list), len(file_list)), 
                              dtype=np.int32, order='F')
    for i, file in enumerate(file_list):
        file_list_ascii[:, i] = string_to_ascii_array(file, file_list_ascii.shape[0])
    
    gene_ids_ascii = np.zeros((max(len(g) for g in gene_ids), len(gene_ids)), 
                             dtype=np.int32, order='F')
    for i, gene in enumerate(gene_ids):
        gene_ids_ascii[:, i] = string_to_ascii_array(gene, gene_ids_ascii.shape[0])
    
    # Vorbereiten der Ausgabearrays
    expression_vectors_flat = np.zeros(n_samples * len(gene_ids), dtype=np.float64, order='F')
    ierr = ctypes.c_int()
    delimiter_ascii = string_to_ascii_array(delimiter, 1)
    
    # Convert value_cols to ctypes array
    value_cols_ct = (ctypes.c_int * len(value_cols))(*value_cols)
    
    # C-Funktion aufrufen
    lib.read_expression_vectors_C(
        file_list_ascii.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ctypes.c_int(file_list_ascii.shape[0]),
        ctypes.c_int(len(file_list)),
        gene_ids_ascii.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ctypes.c_int(gene_ids_ascii.shape[0]),
        ctypes.c_int(len(gene_ids)),
        expression_vectors_flat.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_int(n_samples),
        ctypes.c_int(len(gene_ids)),
        ctypes.c_int(n_header_rows),
        ctypes.c_int(gene_col),
        value_cols_ct,
        ctypes.c_int(len(value_cols)),
        ctypes.c_int(start_row),
        ctypes.byref(ierr),
        delimiter_ascii.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ctypes.c_int(1)
    )
    
    if ierr.value != 0:
        raise Exception(f"Error reading expression vectors: {ierr.value}")
    
    # Konvertiere flaches Array zurück in 2D-Array
    expression_vectors = expression_vectors_flat.reshape((len(gene_ids), n_samples), order='F').T
    
    return expression_vectors

# Funktion für read_family_file_C
def read_family_file(filename, gene_ids, family_ids_len, n_families):
    # Konvertiere Eingaben in ASCII-Arrays
    fn_ascii = string_to_ascii_array(filename, len(filename))
    
    gene_ids_ascii = np.zeros((max(len(g) for g in gene_ids), len(gene_ids)), 
                             dtype=np.int32, order='F')
    for i, gene in enumerate(gene_ids):
        gene_ids_ascii[:, i] = string_to_ascii_array(gene, gene_ids_ascii.shape[0])
    
    # Vorbereiten der Ausgabearrays
    family_ids_ascii = np.zeros((family_ids_len, n_families), dtype=np.int32, order='F')
    gene_to_fam = np.zeros(len(gene_ids), dtype=np.int32, order='F')
    ierr = ctypes.c_int()
    
    # C-Funktion aufrufen
    lib.read_family_file_C(
        fn_ascii.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ctypes.c_int(len(filename)),
        gene_ids_ascii.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ctypes.c_int(gene_ids_ascii.shape[0]),
        ctypes.c_int(len(gene_ids)),
        family_ids_ascii.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ctypes.c_int(family_ids_len),
        ctypes.c_int(n_families),
        gene_to_fam.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ctypes.byref(ierr)
    )
    
    if ierr.value != 0:
        raise Exception(f"Error reading family file: {ierr.value}")
    
    # Konvertiere Ergebnis zurück zu Strings
    family_ids = []
    for i in range(n_families):
        family_ids.append(ascii_array_to_string(family_ids_ascii[:, i]))
    
    return family_ids, gene_to_fam

# Funktion für filter_unassigned_genes_C
def filter_unassigned_genes(gene_ids, expression_vectors, gene_to_fam):
    # Konvertiere Eingaben in ASCII-Arrays
    gene_ids_ascii = np.zeros((max(len(g) for g in gene_ids), len(gene_ids)), 
                             dtype=np.int32, order='F')
    for i, gene in enumerate(gene_ids):
        gene_ids_ascii[:, i] = string_to_ascii_array(gene, gene_ids_ascii.shape[0])
    
    # Vorbereiten der Ausgabearrays
    n_samples, n_genes = expression_vectors.shape
    expression_vectors_flat = expression_vectors.T.reshape(-1, order='F')
    mask = np.zeros(n_genes, dtype=ctypes.c_int())
    n_genes_kept = ctypes.c_int()
    ierr = ctypes.c_int()
    
    # C-Funktion aufrufen
    lib.filter_unassigned_genes_C(
        gene_ids_ascii.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ctypes.c_int(gene_ids_ascii.shape[0]),
        ctypes.c_int(n_genes),
        expression_vectors_flat.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_int(n_samples),
        gene_to_fam.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        mask.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ctypes.byref(n_genes_kept),
        ctypes.byref(ierr)
    )
    
    if ierr.value != 0:
        raise Exception(f"Error filtering unassigned genes: {ierr.value}")
    
    return mask, n_genes_kept.value

# --- Python wrappers for validation ---

def validate_gene_to_family_mapping(gene_to_fam, n_families):
    gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    n_genes = gene_to_fam.size
    ierr = ctypes.c_int()
    lib.validate_gene_to_family_mapping_C(
        gene_to_fam.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ctypes.c_int(n_genes),
        ctypes.c_int(n_families),
        ctypes.byref(ierr)
    )
    if ierr.value != 0:
        raise Exception(f"Validation failed: gene_to_fam (error {ierr.value})")

def validate_expression_data(expression_vectors, check_non_negative=True):
    arr = np.asfortranarray(expression_vectors, dtype=np.float64)
    n_genes, n_samples = arr.shape
    ierr = ctypes.c_int()
    lib.validate_expression_data_C(
        arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_int(n_genes),
        ctypes.c_int(n_samples),
        ctypes.c_int(1 if check_non_negative else 0),
        ctypes.byref(ierr)
    )
    if ierr.value != 0:
        raise Exception(f"Validation failed: expression_vectors (error {ierr.value})")

def validate_family_centroids(family_centroids):
    arr = np.asfortranarray(family_centroids, dtype=np.float64)
    d, n_families = arr.shape
    ierr = ctypes.c_int()
    lib.validate_family_centroids_C(
        arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_int(n_families),
        ctypes.c_int(d),
        ctypes.byref(ierr)
    )
    if ierr.value != 0:
        raise Exception(f"Validation failed: family_centroids (error {ierr.value})")

def validate_shift_vectors(shift_vectors, expression_vectors, family_centroids, gene_to_fam, d, n_genes, n_samples, n_families):
    shift_vectors = np.asfortranarray(shift_vectors, dtype=np.float64)
    expression_vectors = np.asfortranarray(expression_vectors, dtype=np.float64)
    family_centroids = np.asfortranarray(family_centroids, dtype=np.float64)
    gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    ierr = ctypes.c_int()
    lib.validate_shift_vectors_C(
        shift_vectors.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        expression_vectors.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        family_centroids.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        gene_to_fam.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ctypes.c_int(d),
        ctypes.c_int(n_genes),
        ctypes.c_int(n_samples),
        ctypes.c_int(n_families),
        ctypes.byref(ierr)
    )
    if ierr.value != 0:
        raise Exception(f"Validation failed: shift_vectors (error {ierr.value})")

def validate_gene_ids_uniqueness(gene_ids):
    gene_ids_len = max(len(g) for g in gene_ids)
    n_genes = len(gene_ids)
    gene_ids_ascii = np.zeros((gene_ids_len, n_genes), dtype=np.int32, order='F')
    for i, gene in enumerate(gene_ids):
        gene_ids_ascii[:, i] = string_to_ascii_array(gene, gene_ids_len)
    ierr = ctypes.c_int()
    lib.validate_gene_ids_uniqueness_C(
        gene_ids_ascii.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ctypes.c_int(gene_ids_len),
        ctypes.c_int(n_genes),
        ctypes.byref(ierr)
    )
    if ierr.value != 0:
        raise Exception(f"Validation failed: gene_ids uniqueness (error {ierr.value})")

def validate_family_ids_uniqueness(family_ids):
    fam_len = max(len(f) for f in family_ids)
    n_families = len(family_ids)
    family_ids_ascii = np.zeros((fam_len, n_families), dtype=np.int32, order='F')
    for i, fam in enumerate(family_ids):
        family_ids_ascii[:, i] = string_to_ascii_array(fam, fam_len)
    ierr = ctypes.c_int()
    lib.validate_family_ids_uniqueness_C(
        family_ids_ascii.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ctypes.c_int(fam_len),
        ctypes.c_int(n_families),
        ctypes.byref(ierr)
    )
    if ierr.value != 0:
        raise Exception(f"Validation failed: family_ids uniqueness (error {ierr.value})")

def validate_data_structure(n_genes, n_families, n_samples, d, gene_ids, gene_family_ids, gene_to_fam, expression_vectors, family_centroids, shift_vectors):
    gene_ids_len = max(len(g) for g in gene_ids)
    fam_len = max(len(f) for f in gene_family_ids)
    gene_ids_ascii = np.zeros((gene_ids_len, n_genes), dtype=np.int32, order='F')
    gene_family_ids_ascii = np.zeros((fam_len, n_genes), dtype=np.int32, order='F')
    for i in range(n_genes):
        gene_ids_ascii[:, i] = string_to_ascii_array(gene_ids[i], gene_ids_len)
    for i in range(n_families):
        gene_family_ids_ascii[:, i] = string_to_ascii_array(gene_family_ids[i], fam_len)
    gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    expression_vectors = np.asfortranarray(expression_vectors, dtype=np.float64)
    family_centroids = np.asfortranarray(family_centroids, dtype=np.float64)
    shift_vectors = np.asfortranarray(shift_vectors, dtype=np.float64)
    ierr = ctypes.c_int()
    lib.validate_data_structure_C(
        ctypes.c_int(n_genes),
        ctypes.c_int(n_families),
        ctypes.c_int(n_samples),
        ctypes.c_int(d),
        gene_ids_ascii.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ctypes.c_int(gene_ids_len),
        gene_family_ids_ascii.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ctypes.c_int(fam_len),
        gene_to_fam.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        expression_vectors.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        family_centroids.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        shift_vectors.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.byref(ierr)
    )
    if ierr.value != 0:
        raise Exception(f"Validation failed: data structure (error {ierr.value})")

def validate_all_data(n_genes, n_families, n_samples, d, gene_ids, gene_family_ids, gene_to_fam, expression_vectors, family_centroids, shift_vectors):
    gene_ids_len = max(len(g) for g in gene_ids)
    fam_len = max(len(f) for f in gene_family_ids)
    gene_ids_ascii = np.zeros((gene_ids_len, n_genes), dtype=np.int32, order='F')
    gene_family_ids_ascii = np.zeros((fam_len, n_families), dtype=np.int32, order='F')
    for i in range(n_genes):
        gene_ids_ascii[:, i] = string_to_ascii_array(gene_ids[i], gene_ids_len)
    for i in range(n_families):
        gene_family_ids_ascii[:, i] = string_to_ascii_array(gene_family_ids[i], fam_len)
    gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    expression_vectors = np.asfortranarray(expression_vectors, dtype=np.float64)
    family_centroids = np.asfortranarray(family_centroids, dtype=np.float64)
    shift_vectors = np.asfortranarray(shift_vectors, dtype=np.float64)
    ierr = ctypes.c_int()
    lib.validate_all_data_C(
        ctypes.c_int(n_genes),
        ctypes.c_int(n_families),
        ctypes.c_int(n_samples),
        ctypes.c_int(d),
        gene_ids_ascii.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ctypes.c_int(gene_ids_len),
        gene_family_ids_ascii.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ctypes.c_int(fam_len),
        gene_to_fam.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        expression_vectors.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        family_centroids.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        shift_vectors.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.byref(ierr)
    )
    if ierr.value != 0:
        raise Exception(f"Validation failed: all data (error {ierr.value})")