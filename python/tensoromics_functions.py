import ctypes
import numpy as np
from error_handling import check_err_code

# Load compiled libraries
lib = ctypes.CDLL('build/libtensor-omics.so')  # For KD-Tree functions

# Configure BST argument types
lib.build_bst_index_C.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),  # values
    ctypes.c_int32,                                                  # num_values
    np.ctypeslib.ndpointer(dtype=np.int32),                          # sorted_indices (out)
    np.ctypeslib.ndpointer(dtype=np.int32),                          # stack_left
    np.ctypeslib.ndpointer(dtype=np.int32),                          # stack_right
    ctypes.POINTER(ctypes.c_int)                                     # ierr   
]

lib.bst_range_query_C.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64),  # values
    np.ctypeslib.ndpointer(dtype=np.int32),    # sorted_indices
    ctypes.c_int32,                            # num_values
    ctypes.c_double,                           # low
    ctypes.c_double,                           # high
    np.ctypeslib.ndpointer(dtype=np.int32),    # out_indices (out)
    ctypes.POINTER(ctypes.c_int32),            # number_matches (out)
    ctypes.POINTER(ctypes.c_int)
]

# Configure KD-Tree argument types
lib.build_kd_index_C.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64, flags='F_CONTIGUOUS'),  # X_flat (col-major)
    ctypes.c_int32,                                                  # d
    ctypes.c_int32,                                                  # n
    np.ctypeslib.ndpointer(dtype=np.int32),                          # kd_ix (out)
    np.ctypeslib.ndpointer(dtype=np.int32),                          # dim_order
    np.ctypeslib.ndpointer(dtype=np.int32),                          # work
    np.ctypeslib.ndpointer(dtype=np.float64),                        # subarray
    np.ctypeslib.ndpointer(dtype=np.int32),                          # perm
    np.ctypeslib.ndpointer(dtype=np.int32),                          # stack_left
    np.ctypeslib.ndpointer(dtype=np.int32),                          # stack_right
    ctypes.POINTER(ctypes.c_int)
]

# --- BST Functions ---
def build_bst_index(values):
    """
    Build a BST index for the given values.
    
    Parameters:
    values (np.array): 1D array of values to index
    
    Returns:
    np.array: BST indices (0-based for Python)
    """
    n = len(values)
    indices = np.empty(n, dtype=np.int32)
    stack_left = np.empty(n, dtype=np.int32)
    stack_right = np.empty(n, dtype=np.int32)
    ierr = ctypes.c_int()
    
    # Build BST index
    lib.build_bst_index_C(values, n, indices, stack_left, stack_right, ctypes.byref(ierr))
    check_err_code(ierr.value)
    
    # Convert from Fortran 1-based to Python 0-based indexing
    return indices - 1

def bst_range_query(values, indices, lower_bound, upper_bound):
    """
    Perform a range query on BST-indexed values.
    
    Parameters:
    values (np.array): Original values array
    indices (np.array): BST indices from build_bst_index
    lower_bound (float): Lower bound of range (inclusive)
    upper_bound (float): Upper bound of range (inclusive)
    
    Returns:
    tuple: (matching_indices, count) where matching_indices are 0-based Python indices
    """
    n = len(values)
    output_indices = np.empty(n, dtype=np.int32)
    match_count = ctypes.c_int32(0)
    ierr = ctypes.c_int()
    
    # Convert indices back to 1-based for Fortran
    indices_1based = indices + 1
    
    # Perform range query
    lib.bst_range_query_C(values, indices_1based, n, lower_bound, upper_bound, 
                         output_indices, ctypes.byref(match_count), ctypes.byref(ierr))
    check_err_code(ierr.value)
    
    # Convert from Fortran 1-based to Python 0-based indexing
    matching_indices = output_indices[:match_count.value] - 1
    
    return matching_indices, match_count.value

# --- KD-Tree Functions ---
def build_kd_index(points, dimension_order=None):
    """
    Build a KD-Tree index for the given points.
    
    Parameters:
    points (np.array): 2D array of points (d x n, Fortran order)
    dimension_order (np.array): Order of dimensions for splitting (1-based)
    
    Returns:
    np.array: KD-Tree indices (0-based for Python)
    """
    d, n = points.shape
    
    # Use default dimension order if not provided
    if dimension_order is None:
        dimension_order = np.arange(1, d + 1, dtype=np.int32)  # 1-based dimensions
    
    # Ensure Fortran order (column-major)
    if not points.flags.f_contiguous:
        points = np.asfortranarray(points)
    
    # Initialize arrays
    kd_indices = np.empty(n, dtype=np.int32)
    workspace = np.empty(n, dtype=np.int32)
    value_buffer = np.empty(n, dtype=np.float64)
    permutation = np.empty(n, dtype=np.int32)
    stack_left = np.empty(n, dtype=np.int32)
    stack_right = np.empty(n, dtype=np.int32)
    ierr = ctypes.c_int()
    
    # Build KD-Tree index
    lib.build_kd_index_C(points, d, n, kd_indices, dimension_order, workspace, 
                        value_buffer, permutation, stack_left, stack_right, ctypes.byref(ierr))
    check_err_code(ierr.value)
    
    # Convert from Fortran 1-based to Python 0-based indexing
    return kd_indices - 1

def build_spherical_kd(vectors, dimension_order=None):
    """
    Build a spherical KD-Tree index for the given unit vectors.
    
    Parameters:
    vectors (np.array): 2D array of unit vectors (d x n, Fortran order)
    dimension_order (np.array): Order of dimensions for splitting (1-based)
    
    Returns:
    np.array: Spherical KD-Tree indices (0-based for Python)
    """
    # For spherical KD-Tree, we use the same implementation as regular KD-Tree
    # but with a different name for clarity
    return build_kd_index(vectors, dimension_order)

