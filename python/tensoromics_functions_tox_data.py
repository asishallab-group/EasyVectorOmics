"""
TensorOmics Functions Module
Python wrapper functions for Fortran routines via C interface
"""
from error_handling import check_err_code
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

# Example for read_gene_ids_from_file_C
lib.read_gene_ids_from_file_C.argtypes = [
    ctypes.POINTER(ctypes.c_char),  # filename_raw
    ctypes.c_int,                  # fn_len
    ctypes.POINTER(ctypes.c_char),  # gene_ids_raw
    ctypes.c_int,                  # gene_ids_len
    ctypes.c_int,                  # n_genes
    ctypes.c_int,                  # n_header_rows
    ctypes.c_int,                  # gene_col
    ctypes.POINTER(ctypes.c_int)   # ierr
]
lib.read_gene_ids_from_file_C.restype = None

# read_expression_vectors_C
lib.read_expression_vectors_C.argtypes = [
    ctypes.POINTER(ctypes.c_char),  # file_list_raw
    ctypes.c_int,                  # file_list_len
    ctypes.c_int,                  # n_files
    ctypes.POINTER(ctypes.c_char),  # gene_ids_raw
    ctypes.c_int,                  # gene_ids_len
    ctypes.c_int,                  # n_genes
    ctypes.POINTER(ctypes.c_double), # expression_vectors_flat
    ctypes.c_int,                  # n_samples
    ctypes.c_int,                  # n_header_rows
    ctypes.c_int,                  # gene_col
    ctypes.POINTER(ctypes.c_int),  # value_cols
    ctypes.c_int,                  # n_value_cols
    ctypes.POINTER(ctypes.c_int),  # ierr
    ctypes.POINTER(ctypes.c_char),  # delimiter_raw
]
lib.read_expression_vectors_C.restype = None

# read_family_file_C
lib.read_family_file_C.argtypes = [
    ctypes.POINTER(ctypes.c_char),  # filename_raw
    ctypes.c_int,                  # fn_len
    ctypes.POINTER(ctypes.c_char),  # gene_ids_raw
    ctypes.c_int,                  # gene_ids_len
    ctypes.c_int,                  # n_genes
    ctypes.POINTER(ctypes.c_char),  # family_ids_raw
    ctypes.c_int,                  # family_ids_len
    ctypes.c_int,                  # n_families
    ctypes.POINTER(ctypes.c_int),  # gene_to_fam
    ctypes.POINTER(ctypes.c_int)   # ierr
]
lib.read_family_file_C.restype = None

# filter_unassigned_genes_C
lib.filter_unassigned_genes_C.argtypes = [
    ctypes.POINTER(ctypes.c_char),  # gene_ids_raw
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
    ctypes.POINTER(ctypes.c_char), # gene_ids_raw
    ctypes.c_int,                 # gene_ids_len
    ctypes.c_int,                 # n_genes
    ctypes.POINTER(ctypes.c_int)  # ierr
]
lib.validate_gene_ids_uniqueness_C.restype = None

lib.validate_family_ids_uniqueness_C.argtypes = [
    ctypes.POINTER(ctypes.c_char), # family_ids_raw
    ctypes.c_int,                 # fam_len
    ctypes.c_int,                 # n_families
    ctypes.POINTER(ctypes.c_int)  # ierr
]
lib.validate_family_ids_uniqueness_C.restype = None

lib.validate_data_structure_C.argtypes = [
    ctypes.c_int,                 # n_genes
    ctypes.c_int,                 # n_families
    ctypes.c_int,                 # n_samples
    ctypes.POINTER(ctypes.c_char), # gene_ids_raw
    ctypes.c_int,                 # gene_ids_len
    ctypes.POINTER(ctypes.c_char), # gene_family_ids_raw
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
    ctypes.POINTER(ctypes.c_char), # gene_ids_raw
    ctypes.c_int,                 # gene_len
    ctypes.POINTER(ctypes.c_char), # gene_family_ids_raw
    ctypes.c_int,                 # fam_len
    ctypes.POINTER(ctypes.c_int), # gene_to_fam
    ctypes.POINTER(ctypes.c_double), # expression_vectors
    ctypes.POINTER(ctypes.c_double), # family_centroids
    ctypes.POINTER(ctypes.c_double), # shift_vectors
    ctypes.POINTER(ctypes.c_int)     # ierr
]
lib.validate_all_data_C.restype = None

def string_to_c_char_array(s, length):
    """Convert string to c_char array with null termination"""
    if s is None:
        s = ""
    # Create array of c_char with specified length
    arr = (ctypes.c_char * length)()
    # Encode string and copy to array
    encoded = s.encode('ascii')
    for i in range(min(length, len(encoded))):
        arr[i] = encoded[i]
    # Ensure null termination if there's space
    if len(encoded) < length:
        arr[len(encoded)] = b'\x00'
    return arr

def c_char_array_to_string(c_array):
    """Convert c_char array back to string"""
    # Find null terminator or use full length
    bytes_list = []
    for i in range(len(c_array)):
        byte = c_array[i]
        if byte == b'\x00':
            break
        bytes_list.append(byte)
    return b''.join(bytes_list).decode('ascii').strip()

def strings_to_c_char_matrix(strings, max_length):
    """Convert list of strings to flat c_char array (Fortran-compatible, NumPy-wrapped)"""
    import numpy as np
    n_strings = len(strings)
    total_size = n_strings * max_length
    
    # Flaches ctypes-Array erzeugen
    matrix_type = ctypes.c_char * total_size
    matrix = matrix_type()

    # Initialisiere alles mit Nullbytes
    for i in range(total_size):
        matrix[i] = b'\x00'

    # Strings einfüllen (Fortran order)
    for i, s in enumerate(strings):
        encoded = s.encode('ascii')
        for j in range(min(max_length, len(encoded))):
            index = j + i * max_length
            matrix[index] = encoded[j:j+1]
        if len(encoded) < max_length:
            matrix[len(encoded) + i * max_length] = b'\x00'

    arr = np.ctypeslib.as_array(matrix)
    arr = arr.reshape((n_strings, max_length), order='F')

    return arr

def c_char_matrix_to_strings(matrix, max_length, n_strings):
    """Convert 2D c_char matrix or NumPy array back to list of strings (Fortran order)"""
    import numpy as np
    
    # Falls NumPy-Array: flach als 1D betrachten (Fortran order)
    if isinstance(matrix, np.ndarray):
        flat = matrix.ravel(order='F')
    else:
        flat = matrix

    strings = []
    for i in range(n_strings):
        chars = []
        for j in range(max_length):
            index = j + i * max_length
            if index >= len(flat):
                break

            char = flat[index]

            # 🔧 Fix: sicherstellen, dass 'char' ein einzelnes Byte ist
            if isinstance(char, np.ndarray):
                char = char.item()  # extrahiere den 0D-Wert
            
            # in Bytes umwandeln, falls es ein int ist
            if isinstance(char, (np.integer, int)):
                char = bytes([char])
            elif isinstance(char, (bytes, bytearray)):
                pass
            elif isinstance(char, str):
                char = char.encode('ascii')
            else:
                char = bytes([int(char)])

            # Null-Terminator prüfen
            if char == b'\x00':
                break

            chars.append(char)

        s = b''.join(chars).decode('ascii')
        strings.append(s)

    return strings


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
    # Ensure filename is a string (single file)
    if isinstance(filename, list):
        if len(filename) > 0:
            filename = filename[0]
        else:
            raise ValueError("filename cannot be an empty list for read_gene_ids_from_file")
    
    # Convert filename to c_char array
    fn_array = string_to_c_char_array(filename, len(filename))
    
    # Create output array for gene IDs - proper 2D array
    matrix_size = gene_ids_len * n_genes
    gene_ids_array = (ctypes.c_char * matrix_size)()
    
    ierr = ctypes.c_int()
    
    print(f"Debug: Creating gene_ids_array with size {matrix_size} = {gene_ids_len} * {n_genes}")
    
    # Call C function
    lib.read_gene_ids_from_file_C(
        ctypes.cast(ctypes.byref(fn_array), ctypes.POINTER(ctypes.c_char)),
        ctypes.c_int(len(filename)),
        ctypes.cast(ctypes.byref(gene_ids_array), ctypes.POINTER(ctypes.c_char)),
        ctypes.c_int(gene_ids_len),
        ctypes.c_int(n_genes),
        ctypes.c_int(n_header_rows),
        ctypes.c_int(gene_col),
        ctypes.byref(ierr)
    )
    
    check_err_code(ierr.value)
    
    # Convert result back to strings
    gene_ids = c_char_matrix_to_strings(gene_ids_array, gene_ids_len, n_genes)
    
    return np.array(gene_ids, dtype='U')

# Function for read_expression_vectors_C
def read_expression_vectors(file_list, gene_ids, n_samples, n_header_rows, 
                           gene_col, value_cols, delimiter='\t'):
    # Ensure file_list is a list (multiple files)
    if not isinstance(file_list, list):
        file_list = [file_list]  # Convert single file to list
    
    # Ensure inputs are numpy arrays
    gene_ids = _ensure_string_array(gene_ids)
    
    # Convert inputs to c_char matrices
    max_file_len = max(len(f) for f in file_list)
    file_list_matrix = strings_to_c_char_matrix(file_list, max_file_len)
    
    max_gene_len = max(len(g) for g in gene_ids)
    gene_ids_matrix = strings_to_c_char_matrix(gene_ids, max_gene_len)
    
    # Prepare output arrays
    expression_vectors = np.zeros((n_samples, len(gene_ids)), dtype=np.float64, order='F')
    ierr = ctypes.c_int()
    delimiter_array = string_to_c_char_array(delimiter, len(delimiter))
    
    # Convert value_cols to ctypes array
    value_cols_ct = (ctypes.c_int * len(value_cols))(*value_cols)
    
    # Call C function
    lib.read_expression_vectors_C(
        file_list_matrix.ctypes.data_as(ctypes.POINTER(ctypes.c_char)),
        ctypes.c_int(max_file_len),
        ctypes.c_int(len(file_list)),
        gene_ids_matrix.ctypes.data_as(ctypes.POINTER(ctypes.c_char)),
        ctypes.c_int(max_gene_len),
        ctypes.c_int(len(gene_ids)),
        expression_vectors.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.c_int(n_samples),
        ctypes.c_int(n_header_rows),
        ctypes.c_int(gene_col),
        value_cols_ct,
        ctypes.c_int(len(value_cols)),
        ctypes.byref(ierr),
        ctypes.cast(ctypes.byref(delimiter_array), ctypes.POINTER(ctypes.c_char))
    )
    
    check_err_code(ierr.value)
        
    return expression_vectors

# Function for read_family_file_C
def read_family_file(filename, gene_ids, family_ids_len, n_families):
    # Ensure filename is a string (single file)
    if isinstance(filename, list):
        if len(filename) > 0:
            filename = filename[0]  # Take first element if it's a list
        else:
            raise ValueError("filename cannot be an empty list for read_family_file")
    
    # Ensure inputs are numpy arrays
    gene_ids = _ensure_string_array(gene_ids)
    
    # Convert inputs to c_char arrays
    fn_array = string_to_c_char_array(filename, len(filename))
    
    max_gene_len = max(len(g) for g in gene_ids)
    gene_ids_matrix = strings_to_c_char_matrix(gene_ids, max_gene_len)
    
    # Prepare output arrays
    family_ids_matrix = strings_to_c_char_matrix([""] * n_families, family_ids_len)
    gene_to_fam = np.zeros(len(gene_ids), dtype=np.int32, order='F')
    ierr = ctypes.c_int()
    
    # Call C function
    lib.read_family_file_C(
        ctypes.cast(ctypes.byref(fn_array), ctypes.POINTER(ctypes.c_char)),
        ctypes.c_int(len(filename)),
        gene_ids_matrix.ctypes.data_as(ctypes.POINTER(ctypes.c_char)),
        ctypes.c_int(max_gene_len),
        ctypes.c_int(len(gene_ids)),
        family_ids_matrix.ctypes.data_as(ctypes.POINTER(ctypes.c_char)),
        ctypes.c_int(family_ids_len),
        ctypes.c_int(n_families),
        gene_to_fam.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ctypes.byref(ierr)
    )

    
    check_err_code(ierr.value)
    
    # Convert result back to strings as numpy array
    family_ids = c_char_matrix_to_strings(family_ids_matrix, family_ids_len, n_families)
    
    return {
        'family_ids': np.array(family_ids, dtype='U'),
        'gene_to_fam': gene_to_fam
    }

# Function for filter_unassigned_genes_C
def filter_unassigned_genes(gene_ids, gene_to_fam):
    # Ensure inputs are numpy arrays
    gene_ids = _ensure_string_array(gene_ids)
    gene_to_fam = _ensure_int_array(gene_to_fam)
    
    # Convert inputs to c_char matrix
    max_gene_len = max(len(g) for g in gene_ids)
    gene_ids_matrix = strings_to_c_char_matrix(gene_ids, max_gene_len)
    
    # Prepare output arrays
    n_genes = len(gene_ids)
    mask = np.zeros(n_genes, dtype=ctypes.c_int())
    n_genes_kept = ctypes.c_int()
    ierr = ctypes.c_int()
    
    # Call C function
    lib.filter_unassigned_genes_C(
        gene_ids_matrix.ctypes.data_as(ctypes.POINTER(ctypes.c_char)),
        ctypes.c_int(max_gene_len),
        ctypes.c_int(n_genes),
        gene_to_fam.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        mask.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        ctypes.byref(n_genes_kept),
        ctypes.byref(ierr)
    )
    
    check_err_code(ierr.value)
    
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
    check_err_code(ierr.value)

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
    check_err_code(ierr.value)

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
    check_err_code(ierr.value)

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
    check_err_code(ierr.value)

def validate_gene_ids_uniqueness(gene_ids):
    gene_ids = _ensure_string_array(gene_ids)
    gene_ids_len = max(len(g) for g in gene_ids)
    n_genes = len(gene_ids)
    gene_ids_raw = strings_to_c_char_matrix(gene_ids, gene_ids_len)
    ierr = ctypes.c_int()
    lib.validate_gene_ids_uniqueness_C(
        gene_ids_raw.ctypes.data_as(ctypes.POINTER(ctypes.c_char)),
        ctypes.c_int(gene_ids_len),
        ctypes.c_int(n_genes),
        ctypes.byref(ierr)
    )
    check_err_code(ierr.value)

def validate_family_ids_uniqueness(family_ids):
    family_ids = _ensure_string_array(family_ids)
    fam_len = max(len(f) for f in family_ids)
    n_families = len(family_ids)
    family_ids_raw = strings_to_c_char_matrix(family_ids, fam_len)
    ierr = ctypes.c_int()
    lib.validate_family_ids_uniqueness_C(
        family_ids_raw.ctypes.data_as(ctypes.POINTER(ctypes.c_char)),
        ctypes.c_int(fam_len),
        ctypes.c_int(n_families),
        ctypes.byref(ierr)
    )
    check_err_code(ierr.value)

def validate_data_structure(n_genes, n_families, n_samples, gene_ids, gene_family_ids, gene_to_fam, expression_vectors, family_centroids, shift_vectors):
    gene_ids = _ensure_string_array(gene_ids)
    gene_family_ids = _ensure_string_array(gene_family_ids)
    expression_vectors = _ensure_float_array(expression_vectors)
    family_centroids = _ensure_float_array(family_centroids)
    shift_vectors = _ensure_float_array(shift_vectors)
    gene_to_fam = _ensure_int_array(gene_to_fam)
    
    gene_ids_len = max(len(g) for g in gene_ids)
    fam_len = max(len(f) for f in gene_family_ids)
    gene_ids_raw = strings_to_c_char_matrix(gene_ids, gene_ids_len)
    gene_family_ids_raw = strings_to_c_char_matrix(gene_family_ids, fam_len)
    
    gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    expression_vectors = np.asfortranarray(expression_vectors, dtype=np.float64)
    family_centroids = np.asfortranarray(family_centroids, dtype=np.float64)
    shift_vectors = np.asfortranarray(shift_vectors, dtype=np.float64)
    ierr = ctypes.c_int()
    
    lib.validate_data_structure_C(
        ctypes.c_int(n_genes),
        ctypes.c_int(n_families),
        ctypes.c_int(n_samples),
        gene_ids_raw.ctypes.data_as(ctypes.POINTER(ctypes.c_char)),
        ctypes.c_int(gene_ids_len),
        gene_family_ids_raw.ctypes.data_as(ctypes.POINTER(ctypes.c_char)),
        ctypes.c_int(fam_len),
        gene_to_fam.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        expression_vectors.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        family_centroids.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        shift_vectors.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.byref(ierr)
    )
    check_err_code(ierr.value)

def validate_all_data(n_genes, n_families, n_samples, gene_ids, gene_family_ids, gene_to_fam, expression_vectors, family_centroids, shift_vectors):
    gene_ids = _ensure_string_array(gene_ids)
    gene_family_ids = _ensure_string_array(gene_family_ids)
    expression_vectors = _ensure_float_array(expression_vectors)
    family_centroids = _ensure_float_array(family_centroids)
    shift_vectors = _ensure_float_array(shift_vectors)
    gene_to_fam = _ensure_int_array(gene_to_fam)
    
    gene_ids_len = max(len(g) for g in gene_ids)
    fam_len = max(len(f) for f in gene_family_ids)
    gene_ids_raw = strings_to_c_char_matrix(gene_ids, gene_ids_len)
    gene_family_ids_raw = strings_to_c_char_matrix(gene_family_ids, fam_len)
    
    gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    expression_vectors = np.asfortranarray(expression_vectors, dtype=np.float64)
    family_centroids = np.asfortranarray(family_centroids, dtype=np.float64)
    shift_vectors = np.asfortranarray(shift_vectors, dtype=np.float64)
    ierr = ctypes.c_int()
    
    lib.validate_all_data_C(
        ctypes.c_int(n_genes),
        ctypes.c_int(n_families),
        ctypes.c_int(n_samples),
        gene_ids_raw.ctypes.data_as(ctypes.POINTER(ctypes.c_char)),
        ctypes.c_int(gene_ids_len),
        gene_family_ids_raw.ctypes.data_as(ctypes.POINTER(ctypes.c_char)),
        ctypes.c_int(fam_len),
        gene_to_fam.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        expression_vectors.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        family_centroids.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        shift_vectors.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ctypes.byref(ierr)
    )
    check_err_code(ierr.value)
    
def _prepare_string(s) -> tuple:
    """Prepare string for C interface (bytes + length)"""
    if s is None or s == "":
        # Return empty string with null terminator
        return b'\x00', 1
    # Include null terminator
    b = s.encode('utf-8') + b'\x00'
    return b, len(b)

def save_tox_data(zip_filename: str,
                 gene_ids = None,
                 expression_vectors = None,
                 gene_to_fam = None,
                 family_ids = None,
                 family_centroids = None,
                 shift_vectors = None,
                 gene_ids_name = None,
                 expression_vectors_name = None,
                 gene_to_fam_name = None,
                 family_ids_name = None,
                 family_centroids_name = None,
                 shift_vectors_name = None) -> None:
    """
    Save data to a zip archive - Python equivalent of R's save_tox_data
    """
    
    # First serialize the arrays to files (using your existing serialization functions)
    temp_files = []
    
    if gene_ids is not None and gene_ids_name:
        tox_serialize_char_nd(gene_ids, gene_ids_name)
        temp_files.append(gene_ids_name)
    
    if expression_vectors is not None and expression_vectors_name:
        tox_serialize_real_nd(expression_vectors, expression_vectors_name)
        temp_files.append(expression_vectors_name)
    
    if gene_to_fam is not None and gene_to_fam_name:
        tox_serialize_int_nd(gene_to_fam, gene_to_fam_name)
        temp_files.append(gene_to_fam_name)
    
    if family_ids is not None and family_ids_name:
        tox_serialize_char_nd(family_ids, family_ids_name)
        temp_files.append(family_ids_name)
    
    if family_centroids is not None and family_centroids_name:
        tox_serialize_real_nd(family_centroids, family_centroids_name)
        temp_files.append(family_centroids_name)
    
    if shift_vectors is not None and shift_vectors_name:
        tox_serialize_real_nd(shift_vectors, shift_vectors_name)
        temp_files.append(shift_vectors_name)
    
    # Prepare all strings with null terminators
    zip_b, zip_len = _prepare_string(zip_filename)
    gene_ids_b, gene_ids_len = _prepare_string(gene_ids_name)
    expr_b, expr_len = _prepare_string(expression_vectors_name)
    g2f_b, g2f_len = _prepare_string(gene_to_fam_name)
    fam_b, fam_len = _prepare_string(family_ids_name)
    centroids_b, centroids_len = _prepare_string(family_centroids_name)
    shift_b, shift_len = _prepare_string(shift_vectors_name)
    
    # Set up argument types for create_zip_archive_c - using c_int
    lib.create_zip_archive_c.argtypes = [
        ctypes.POINTER(ctypes.c_char), ctypes.c_int,  # zip_filename, zip_len
        ctypes.POINTER(ctypes.c_char), ctypes.c_int,  # gene_ids_file, gene_ids_file_len
        ctypes.POINTER(ctypes.c_char), ctypes.c_int,  # expression_file, expression_file_len
        ctypes.POINTER(ctypes.c_char), ctypes.c_int,  # gene_to_family_file, gene_to_family_file_len
        ctypes.POINTER(ctypes.c_char), ctypes.c_int,  # family_ids_file, family_ids_file_len
        ctypes.POINTER(ctypes.c_char), ctypes.c_int,  # family_centroids_file, family_centroids_file_len
        ctypes.POINTER(ctypes.c_char), ctypes.c_int,  # shift_vectors_file, shift_vectors_file_len
        ctypes.POINTER(ctypes.c_int)                  # ierr
    ]
    lib.create_zip_archive_c.restype = None
    
    ierr = ctypes.c_int()
    
    try:
        # Call the C-bound Fortran subroutine to create the zip archive
        lib.create_zip_archive_c(
            zip_b, ctypes.c_int(zip_len),
            gene_ids_b, ctypes.c_int(gene_ids_len),
            expr_b, ctypes.c_int(expr_len),
            g2f_b, ctypes.c_int(g2f_len),
            fam_b, ctypes.c_int(fam_len),
            centroids_b, ctypes.c_int(centroids_len),
            shift_b, ctypes.c_int(shift_len),
            ctypes.byref(ierr)
        )
        
        check_err_code(ierr.value)
        print(f"Successfully created archive: {zip_filename}")
        
    finally:
        # Clean up temporary files
        for temp_file in temp_files:
            if os.path.exists(temp_file):
                os.remove(temp_file)
                print(f"Removed temporary file: {temp_file}")

def read_tox_data(zip_filename: str,
                 load_gene_ids: bool = False,
                 load_expression_vectors: bool = False,
                 load_gene_to_fam: bool = False,
                 load_family_ids: bool = False,
                 load_family_centroids: bool = False,
                 load_shift_vectors: bool = False):
    """
    Read data from a zip archive - Python equivalent of R's read_tox_data
    """
    
    # Prepare the zip filename
    zip_b = string_to_c_char_array(zip_filename, len(zip_filename))
    zip_len = len(zip_filename)

    
    # Set up argument types for extract_zip_archive_c - using c_int
    lib.extract_zip_archive_c.argtypes = [
        ctypes.POINTER(ctypes.c_char), ctypes.c_int,  # zip_filename, filename_len
        ctypes.POINTER(ctypes.c_int)                  # ierr
    ]
    lib.extract_zip_archive_c.restype = None
    
    ierr = ctypes.c_int()
    
    # Call the C-bound Fortran subroutine to extract the zip archive
    lib.extract_zip_archive_c(
        ctypes.cast(ctypes.byref(zip_b), ctypes.POINTER(ctypes.c_char)),
        ctypes.c_int(zip_len),
        ctypes.byref(ierr)
    )

    check_err_code(ierr.value)
    
    print(f"Successfully extracted archive: {zip_filename}")
    
    result = {
        'gene_ids': None,
        'expression_vectors': None,
        'gene_to_fam': None,
        'family_ids': None,
        'family_centroids': None,
        'shift_vectors': None
    }
    
    # Read manifest file
    manifest_path = "manifest.txt"
    if not os.path.exists(manifest_path):
        raise FileNotFoundError("Manifest file not found in archive")
    
    # Parse manifest
    file_mapping = {}
    with open(manifest_path, 'r') as f:
        for line in f:
            parts = line.strip().split('=')
            if len(parts) == 2:
                file_mapping[parts[0]] = parts[1]
    
    print("Files found in archive:", file_mapping)
    
    # Load requested data using your existing deserialization functions
    try:
        if load_gene_ids and "gene_ids" in file_mapping:
            filename = file_mapping["gene_ids"]
            if os.path.exists(filename):
                result['gene_ids'] = tox_deserialize_char_nd(filename)
                print(f"Gene ids extracted from {filename}")
        
        if load_expression_vectors and "expression" in file_mapping:
            filename = file_mapping["expression"]
            if os.path.exists(filename):
                result['expression_vectors'] = tox_deserialize_real_nd(filename)
                print(f"Expression vectors extracted from {filename}")
        
        if load_gene_to_fam and "gene_to_family" in file_mapping:
            filename = file_mapping["gene_to_family"]
            if os.path.exists(filename):
                result['gene_to_fam'] = tox_deserialize_int_nd(filename)
                print(f"Gene to family mapping extracted from {filename}")
        
        if load_family_ids and "family_ids" in file_mapping:
            filename = file_mapping["family_ids"]
            if os.path.exists(filename):
                result['family_ids'] = tox_deserialize_char_nd(filename)
                print(f"Family IDs extracted from {filename}")
        
        if load_family_centroids and "family_centroids" in file_mapping:
            filename = file_mapping["family_centroids"]
            if os.path.exists(filename):
                result['family_centroids'] = tox_deserialize_real_nd(filename)
                print(f"Family centroids extracted from {filename}")
        
        if load_shift_vectors and "shift_vectors" in file_mapping:
            filename = file_mapping["shift_vectors"]
            if os.path.exists(filename):
                result['shift_vectors'] = tox_deserialize_real_nd(filename)
                print(f"Shift vectors extracted from {filename}")
                
    finally:
        # Cleanup extracted files
        files_to_remove = [manifest_path]
        for filename in file_mapping.values():
            if os.path.exists(filename):
                files_to_remove.append(filename)
        
        for file in files_to_remove:
            if os.path.exists(file):
                os.remove(file)
                print(f"Cleaned up: {file}")
    
    return result