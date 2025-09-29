"""
TensorOmics Functions Module
Python wrapper functions for Fortran routines via C interface
"""

import numpy as np
import ctypes
import os
from tensoromics_functions import (
    tox_deserialize_char_nd,
    tox_serialize_char_nd,
    tox_serialize_int_nd,
    tox_deserialize_int_nd,
    tox_serialize_real_nd,
    tox_deserialize_real_nd
)

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
    ctypes.c_int,                  # n_header_rows
    ctypes.c_int,                  # gene_col
    ctypes.POINTER(ctypes.c_int),  # value_cols
    ctypes.c_int,                  # n_value_cols
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
    ctypes.c_int,                    # n_samples
    ctypes.POINTER(ctypes.c_int)     # ierr
]
lib.validate_family_centroids_C.restype = None

lib.validate_shift_vectors_C.argtypes = [
    ctypes.POINTER(ctypes.c_double), # shift_vectors
    ctypes.POINTER(ctypes.c_double), # expression_vectors
    ctypes.POINTER(ctypes.c_double), # family_centroids
    ctypes.POINTER(ctypes.c_int),    # gene_to_fam
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

# Helper functions for string conversion
def string_to_ascii_array(s, length):
    """Converts a string to a fixed-length ASCII integer array"""
    ascii_array = np.zeros(length, dtype=np.int32, order='F')
    for i, char in enumerate(s):
        if i >= length:
            break
        ascii_array[i] = ord(char)
    return ascii_array

def ascii_array_to_string(ascii_array):
    """Converts an ASCII integer array back to a string"""
    return ''.join(chr(x) for x in ascii_array if x != 0).strip()

def _ensure_string_array(arr):
    """Ensure input is a numpy array of strings"""
    if not isinstance(arr, np.ndarray):
        arr = np.array(arr, dtype=object)
    return arr

def _ensure_float_array(arr):
    """Ensure input is a numpy array of floats"""
    if not isinstance(arr, np.ndarray):
        arr = np.array(arr, dtype=np.float64)
    return arr

def _ensure_int_array(arr):
    """Ensure input is a numpy array of ints"""
    if not isinstance(arr, np.ndarray):
        arr = np.array(arr, dtype=np.int32)
    return arr

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
    
    # Convert result back to strings and return as numpy array
    gene_ids = []
    for i in range(n_genes):
        gene_ids.append(ascii_array_to_string(gene_ids_ascii[:, i]))
    
    return np.array(gene_ids, dtype='U')  # Return as numpy array

# Function for read_expression_vectors_C
def read_expression_vectors(file_list, gene_ids, n_samples, n_header_rows, 
                           gene_col, value_cols, delimiter='\t'):
    # Ensure inputs are numpy arrays
    gene_ids = _ensure_string_array(gene_ids)
    
    # Convert inputs to ASCII arrays
    file_list_ascii = np.zeros((max(len(f) for f in file_list), len(file_list)), 
                              dtype=np.int32, order='F')
    for i, file in enumerate(file_list):
        file_list_ascii[:, i] = string_to_ascii_array(file, file_list_ascii.shape[0])
    
    gene_ids_ascii = np.zeros((max(len(g) for g in gene_ids), len(gene_ids)), 
                             dtype=np.int32, order='F')
    for i, gene in enumerate(gene_ids):
        gene_ids_ascii[:, i] = string_to_ascii_array(gene, gene_ids_ascii.shape[0])
    
    # Prepare output arrays
    expression_vectors = np.zeros((n_samples, len(gene_ids)), dtype=np.float64, order='F')
    ierr = ctypes.c_int()
    delimiter_ascii = string_to_ascii_array(delimiter, 1)
    
    # Convert value_cols to ctypes array
    value_cols_ct = (ctypes.c_int * len(value_cols))(*value_cols)
    
    # Call C function
    lib.read_expression_vectors_C(
        file_list_ascii.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ctypes.c_int(file_list_ascii.shape[0]),
        ctypes.c_int(len(file_list)),
        gene_ids_ascii.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ctypes.c_int(gene_ids_ascii.shape[0]),
        ctypes.c_int(len(gene_ids)),
        expression_vectors.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_int(n_samples),
        ctypes.c_int(n_header_rows),
        ctypes.c_int(gene_col),
        value_cols_ct,
        ctypes.c_int(len(value_cols)),
        ctypes.byref(ierr),
        delimiter_ascii.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ctypes.c_int(1)
    )
    
    if ierr.value != 0:
        raise Exception(f"Error reading expression vectors: {ierr.value}")
        
    return expression_vectors

# Function for read_family_file_C
def read_family_file(filename, gene_ids, family_ids_len, n_families):
    # Ensure inputs are numpy arrays
    gene_ids = _ensure_string_array(gene_ids)
    
    # Convert inputs to ASCII arrays
    fn_ascii = string_to_ascii_array(filename, len(filename))
    
    gene_ids_ascii = np.zeros((max(len(g) for g in gene_ids), len(gene_ids)), 
                             dtype=np.int32, order='F')
    for i, gene in enumerate(gene_ids):
        gene_ids_ascii[:, i] = string_to_ascii_array(gene, gene_ids_ascii.shape[0])
    
    # Prepare output arrays
    family_ids_ascii = np.zeros((family_ids_len, n_families), dtype=np.int32, order='F')
    gene_to_fam = np.zeros(len(gene_ids), dtype=np.int32, order='F')
    ierr = ctypes.c_int()
    
    # Call C function
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
    
    # Convert result back to strings as numpy array
    family_ids = []
    for i in range(n_families):
        family_ids.append(ascii_array_to_string(family_ids_ascii[:, i]))
    
    return {
        'family_ids': np.array(family_ids, dtype='U'),
        'gene_to_fam': gene_to_fam
    }

# Function for filter_unassigned_genes_C
def filter_unassigned_genes(gene_ids, gene_to_fam):
    # Ensure inputs are numpy arrays
    gene_ids = _ensure_string_array(gene_ids)
    gene_to_fam = _ensure_int_array(gene_to_fam)
    
    # Convert inputs to ASCII arrays
    gene_ids_ascii = np.zeros((max(len(g) for g in gene_ids), len(gene_ids)), 
                             dtype=np.int32, order='F')
    for i, gene in enumerate(gene_ids):
        gene_ids_ascii[:, i] = string_to_ascii_array(gene, gene_ids_ascii.shape[0])
    
    # Prepare output arrays
    n_genes = len(gene_ids)
    mask = np.zeros(n_genes, dtype=ctypes.c_int())
    n_genes_kept = ctypes.c_int()
    ierr = ctypes.c_int()
    
    # Call C function
    lib.filter_unassigned_genes_C(
        gene_ids_ascii.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ctypes.c_int(gene_ids_ascii.shape[0]),
        ctypes.c_int(n_genes),
        gene_to_fam.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        mask.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ctypes.byref(n_genes_kept),
        ctypes.byref(ierr)
    )
    
    if ierr.value != 0:
        raise Exception(f"Error filtering unassigned genes: {ierr.value}")
    
    return {
        'mask': mask,
        'n_genes_kept': n_genes_kept.value
    }

# --- Python wrappers for validation ---

def validate_gene_to_family_mapping(gene_to_fam, n_families):
    gene_to_fam = _ensure_int_array(gene_to_fam)
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
    arr = _ensure_float_array(expression_vectors)
    arr = np.asfortranarray(arr, dtype=np.float64)
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
    arr = _ensure_float_array(family_centroids)
    arr = np.asfortranarray(arr, dtype=np.float64)
    n_samples, n_families = arr.shape
    ierr = ctypes.c_int()
    lib.validate_family_centroids_C(
        arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_int(n_families),
        ctypes.c_int(n_samples),
        ctypes.byref(ierr)
    )
    if ierr.value != 0:
        raise Exception(f"Validation failed: family_centroids (error {ierr.value})")

def validate_shift_vectors(shift_vectors, expression_vectors, family_centroids, gene_to_fam, n_genes, n_samples, n_families):
    shift_vectors = _ensure_float_array(shift_vectors)
    expression_vectors = _ensure_float_array(expression_vectors)
    family_centroids = _ensure_float_array(family_centroids)
    gene_to_fam = _ensure_int_array(gene_to_fam)
    
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
        ctypes.c_int(n_genes),
        ctypes.c_int(n_samples),
        ctypes.c_int(n_families),
        ctypes.byref(ierr)
    )
    if ierr.value != 0:
        raise Exception(f"Validation failed: shift_vectors (error {ierr.value})")

def validate_gene_ids_uniqueness(gene_ids):
    gene_ids = _ensure_string_array(gene_ids)
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
    family_ids = _ensure_string_array(family_ids)
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

def validate_data_structure(n_genes, n_families, n_samples, gene_ids, gene_family_ids, gene_to_fam, expression_vectors, family_centroids, shift_vectors):
    gene_ids = _ensure_string_array(gene_ids)
    gene_family_ids = _ensure_string_array(gene_family_ids)
    expression_vectors = _ensure_float_array(expression_vectors)
    family_centroids = _ensure_float_array(family_centroids)
    shift_vectors = _ensure_float_array(shift_vectors)
    gene_to_fam = _ensure_int_array(gene_to_fam)
    
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

def validate_all_data(n_genes, n_families, n_samples, gene_ids, gene_family_ids, gene_to_fam, expression_vectors, family_centroids, shift_vectors):
    gene_ids = _ensure_string_array(gene_ids)
    gene_family_ids = _ensure_string_array(gene_family_ids)
    expression_vectors = _ensure_float_array(expression_vectors)
    family_centroids = _ensure_float_array(family_centroids)
    shift_vectors = _ensure_float_array(shift_vectors)
    gene_to_fam = _ensure_int_array(gene_to_fam)
    
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
    
def tox_save_data_archive(zip_filename,
                          gene_ids=None, gene_ids_name=None,
                          expression_vectors=None, expression_vectors_name=None,
                          gene_to_fam=None, gene_to_fam_name=None,
                          family_ids=None, family_ids_name=None,
                          family_centroids=None, family_centroids_name=None,
                          shift_vectors=None, shift_vectors_name=None):
    import zipfile

    # Flags initialisieren
    gene_ids_array_valid = gene_ids_name_valid = False
    expression_vectors_array_valid = expression_vectors_name_valid = False
    gene_to_fam_array_valid = gene_to_fam_name_valid = False
    family_ids_array_valid = family_ids_name_valid = False
    family_centroids_array_valid = family_centroids_name_valid = False
    shift_vectors_array_valid = shift_vectors_name_valid = False

    if not isinstance(zip_filename, str):
        raise RuntimeError("Type mismatch: zip_filename must be a string")

    # Validierungen
    if gene_ids is not None:
        if len(gene_ids.shape) == 1:
            gene_ids_array_valid = True
        else:
            raise Exception(f"Gene IDs dimensions mismatch: Expected 1 but got {len(gene_ids.shape)}")
    if gene_ids_name is not None:
        if isinstance(gene_ids_name, str):
            gene_ids_name_valid = True
        else:
            raise Exception("Gene IDs name must be a string")
    if gene_ids_array_valid ^ gene_ids_name_valid:
        print("Gene IDs: Either provide array and filename or none. Skipping.")

    if expression_vectors is not None:
        if len(expression_vectors.shape) == 2:
            expression_vectors_array_valid = True
        else:
            raise Exception(f"Expression vectors dim mismatch: Expected 2 but got {len(expression_vectors.shape)}")
    if expression_vectors_name is not None:
        if isinstance(expression_vectors_name, str):
            expression_vectors_name_valid = True
        else:
            raise Exception("Expression vector name must be a string")
    if expression_vectors_array_valid ^ expression_vectors_name_valid:
        print("Expression vectors: Either provide array and filename or none. Skipping.")

    if gene_to_fam is not None:
        if len(gene_to_fam.shape) == 1:
            gene_to_fam_array_valid = True
        else:
            raise Exception(f"Gene to family dim mismatch: Expected 1 but got {len(gene_to_fam.shape)}")
    if gene_to_fam_name is not None:
        if isinstance(gene_to_fam_name, str):
            gene_to_fam_name_valid = True
        else:
            raise Exception("Gene to family name must be a string")
    if gene_to_fam_array_valid ^ gene_to_fam_name_valid:
        print("Gene to family: Either provide array and filename or none. Skipping.")

    if family_ids is not None:
        if len(family_ids.shape) == 1:
            family_ids_array_valid = True
        else:
            raise Exception(f"Family IDs dim mismatch: Expected 1 but got {len(family_ids.shape)}")
    if family_ids_name is not None:
        if isinstance(family_ids_name, str):
            family_ids_name_valid = True
        else:
            raise Exception("Family IDs name must be a string")
    if family_ids_array_valid ^ family_ids_name_valid:
        print("Family IDs: Either provide array and filename or none. Skipping.")

    if family_centroids is not None:
        if len(family_centroids.shape) == 2:
            family_centroids_array_valid = True
        else:
            raise Exception(f"Family centroids dim mismatch: Expected 2 but got {len(family_centroids.shape)}")
    if family_centroids_name is not None:
        if isinstance(family_centroids_name, str):
            family_centroids_name_valid = True
        else:
            raise Exception("Family centroids name must be a string")
    if family_centroids_array_valid ^ family_centroids_name_valid:
        print("Family centroids: Either provide array and filename or none. Skipping.")

    if shift_vectors is not None:
        if len(shift_vectors.shape) == 2:
            shift_vectors_array_valid = True
        else:
            raise Exception(f"Shift vectors dim mismatch: Expected 2 but got {len(shift_vectors.shape)}")
    if shift_vectors_name is not None:
        if isinstance(shift_vectors_name, str):
            shift_vectors_name_valid = True
        else:
            raise Exception("Shift vectors name must be a string")
    if shift_vectors_array_valid ^ shift_vectors_name_valid:
        print("Shift vectors: Either provide array and filename or none. Skipping.")

    # Manifest als Liste
    manifest_lines = []

    with zipfile.ZipFile(zip_filename, mode="x") as archive:
        if gene_ids_array_valid and gene_ids_name_valid:
            tox_serialize_char_nd(gene_ids, gene_ids_name)
            archive.write(gene_ids_name)
            print(f"Wrote gene IDs to {gene_ids_name}")
            os.remove(gene_ids_name)
            manifest_lines.append(f"gene_ids={gene_ids_name}")

        if expression_vectors_array_valid and expression_vectors_name_valid:
            tox_serialize_real_nd(expression_vectors, expression_vectors_name)
            archive.write(expression_vectors_name)
            print(f"Wrote expression vectors to {expression_vectors_name}")
            os.remove(expression_vectors_name)
            manifest_lines.append(f"expression={expression_vectors_name}")

        if gene_to_fam_array_valid and gene_to_fam_name_valid:
            tox_serialize_int_nd(gene_to_fam, gene_to_fam_name)
            archive.write(gene_to_fam_name)
            print(f"Wrote gene to family to {gene_to_fam_name}")
            os.remove(gene_to_fam_name)
            manifest_lines.append(f"gene_to_family={gene_to_fam_name}")

        if family_ids_array_valid and family_ids_name_valid:
            tox_serialize_char_nd(family_ids, family_ids_name)
            archive.write(family_ids_name)
            print(f"Wrote family IDs to {family_ids_name}")
            os.remove(family_ids_name)
            manifest_lines.append(f"family_ids={family_ids_name}")

        if family_centroids_array_valid and family_centroids_name_valid:
            tox_serialize_real_nd(family_centroids, family_centroids_name)
            archive.write(family_centroids_name)
            print(f"Wrote family centroids to {family_centroids_name}")
            os.remove(family_centroids_name)
            manifest_lines.append(f"family_centroids={family_centroids_name}")

        if shift_vectors_array_valid and shift_vectors_name_valid:
            tox_serialize_real_nd(shift_vectors, shift_vectors_name)
            archive.write(shift_vectors_name)
            print(f"Wrote shift vectors to {shift_vectors_name}")
            os.remove(shift_vectors_name)
            manifest_lines.append(f"shift_vectors={shift_vectors_name}")

        # Manifest schreiben
        archive.writestr("manifest.txt", "\n".join(manifest_lines) + "\n")

def tox_read_data_archive(zip_filename,
                          gene_ids=None,
                          expression_vectors=None,
                          gene_to_fam=None,
                          family_ids=None,
                          family_centroids=None,
                          shift_vectors=None):
    import zipfile, os
    if not isinstance(zip_filename, str) or zip_filename == "":
        raise RuntimeError("Zip name needs to be a non-empty string")

    result = {
        "gene_ids": None,
        "expression_vectors": None,
        "gene_to_fam": None,
        "family_ids": None,
        "family_centroids": None,
        "shift_vectors": None
    }

    with zipfile.ZipFile(zip_filename, mode="r") as archive:
        with archive.open("manifest.txt") as manifest:
            content = manifest.read().decode("utf-8").splitlines()

        gene_ids_filename         = content[0].split('=')[1]
        expression_vectors_filename = content[1].split('=')[1]
        gene_to_fam_filename      = content[2].split('=')[1]
        family_ids_filename       = content[3].split('=')[1]
        family_centroids_filename = content[4].split('=')[1]
        shift_vectors_filename    = content[5].split('=')[1]

        if gene_ids is not None and gene_ids_filename:
            archive.extract(gene_ids_filename)
            result["gene_ids"] = tox_deserialize_char_nd(gene_ids_filename)
            print(f"Gene ids extracted from {gene_ids_filename}")
            os.remove(gene_ids_filename)

        if expression_vectors is not None and expression_vectors_filename:
            archive.extract(expression_vectors_filename)
            result["expression_vectors"] = tox_deserialize_real_nd(expression_vectors_filename)
            print(f"Expression Vectors extracted from {expression_vectors_filename}")
            os.remove(expression_vectors_filename)

        if gene_to_fam is not None and gene_to_fam_filename:
            archive.extract(gene_to_fam_filename)
            result["gene_to_fam"] = tox_deserialize_int_nd(gene_to_fam_filename)
            print(f"Gene to family mapping extracted from {gene_to_fam_filename}")
            os.remove(gene_to_fam_filename)

        if family_ids is not None and family_ids_filename:
            archive.extract(family_ids_filename)
            result["family_ids"] = tox_deserialize_char_nd(family_ids_filename)
            print(f"Extracted family IDs from {family_ids_filename}")
            os.remove(family_ids_filename)

        if family_centroids is not None and family_centroids_filename:
            archive.extract(family_centroids_filename)
            result["family_centroids"] = tox_deserialize_real_nd(family_centroids_filename)
            print(f"Extracted family centroids from {family_centroids_filename}")
            os.remove(family_centroids_filename)

        if shift_vectors is not None and shift_vectors_filename:
            archive.extract(shift_vectors_filename)
            result["shift_vectors"] = tox_deserialize_real_nd(shift_vectors_filename)
            print(f"Extracted shift vectors from {shift_vectors_filename}")
            os.remove(shift_vectors_filename)

    return result