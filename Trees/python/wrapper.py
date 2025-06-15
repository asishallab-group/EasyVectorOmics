import ctypes
import numpy as np

# Shared Library laden
lib = ctypes.CDLL('./libtrees.so')

# BST Index Builder
lib.build_bst_index_C.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64),  # x
    ctypes.c_int32,                            # n
    np.ctypeslib.ndpointer(dtype=np.int32),    # ix (out)
    np.ctypeslib.ndpointer(dtype=np.int32),    # stack_left
    np.ctypeslib.ndpointer(dtype=np.int32)     # stack_right
]

def build_bst_index_C(x):
    """Python-Wrapper für build_bst_index_C"""
    n = len(x)
    ix = np.empty(n, dtype=np.int32)
    stack_left = np.empty(n, dtype=np.int32)
    stack_right = np.empty(n, dtype=np.int32)
    
    lib.build_bst_index_C(x, n, ix, stack_left, stack_right)
    return ix

# BST Range Query
lib.bst_range_query_C.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64),  # x
    np.ctypeslib.ndpointer(dtype=np.int32),    # ix
    ctypes.c_int32,                            # n
    ctypes.c_double,                           # lo
    ctypes.c_double,                           # hi
    np.ctypeslib.ndpointer(dtype=np.int32),    # out_ix (out)
    ctypes.POINTER(ctypes.c_int32)             # out_n (out)
]

def bst_range_query_C(x, ix, lo, hi):
    """Python-Wrapper für bst_range_query_C"""
    n = len(x)
    out_ix = np.empty(n, dtype=np.int32)
    out_n = ctypes.c_int32(0)
    
    lib.bst_range_query_C(x, ix, n, lo, hi, out_ix, ctypes.byref(out_n))
    return out_ix[:out_n.value]

# KD-Tree Builder
lib.build_kd_index_C.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64, flags='F_CONTIGUOUS'),  # X_flat (d x n)
    ctypes.c_int32,                                                  # d
    ctypes.c_int32,                                                  # n
    np.ctypeslib.ndpointer(dtype=np.int32),                          # kd_ix (out)
    np.ctypeslib.ndpointer(dtype=np.int32),                          # dim_order
    np.ctypeslib.ndpointer(dtype=np.int32),                          # work
    np.ctypeslib.ndpointer(dtype=np.float64),                        # subarray
    np.ctypeslib.ndpointer(dtype=np.int32),                          # perm
    np.ctypeslib.ndpointer(dtype=np.int32),                          # stack_left
    np.ctypeslib.ndpointer(dtype=np.int32)                           # stack_right
]

def build_kd_index_C(X, dim_order):
    """Python-Wrapper für build_kd_index_C"""
    d, n = X.shape
    kd_ix = np.empty(n, dtype=np.int32)
    work = np.empty(n, dtype=np.int32)
    subarray = np.empty(n, dtype=np.float64)
    perm = np.empty(n, dtype=np.int32)
    stack_left = np.empty(64, dtype=np.int32)
    stack_right = np.empty(64, dtype=np.int32)
    
    lib.build_kd_index_C(
        np.asarray(X, order='F'), d, n, kd_ix, dim_order,
        work, subarray, perm, stack_left, stack_right
    )
    return kd_ix

# Spherical KD-Tree Builder
lib.build_spherical_kd_C.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64, flags='F_CONTIGUOUS'),  # V_flat (d x n)
    ctypes.c_int32,                                                  # d
    ctypes.c_int32,                                                  # n
    np.ctypeslib.ndpointer(dtype=np.int32),                          # sphere_ix (out)
    np.ctypeslib.ndpointer(dtype=np.int32),                          # dim_order
    np.ctypeslib.ndpointer(dtype=np.int32),                          # work
    np.ctypeslib.ndpointer(dtype=np.float64),                        # subarray
    np.ctypeslib.ndpointer(dtype=np.int32),                          # perm
    np.ctypeslib.ndpointer(dtype=np.int32),                          # stack_left
    np.ctypeslib.ndpointer(dtype=np.int32)                           # stack_right
]

def build_spherical_kd_C(V, dim_order):
    """Python-Wrapper für build_spherical_kd_C"""
    d, n = V.shape
    sphere_ix = np.empty(n, dtype=np.int32)
    work = np.empty(n, dtype=np.int32)
    subarray = np.empty(n, dtype=np.float64)
    perm = np.empty(n, dtype=np.int32)
    stack_left = np.empty(64, dtype=np.int32)
    stack_right = np.empty(64, dtype=np.int32)
    
    lib.build_spherical_kd_C(
        np.asarray(V, order='F'), d, n, sphere_ix, dim_order,
        work, subarray, perm, stack_left, stack_right
    )
    return sphere_ix

# Test Cases for the Python Wrappers
if __name__ == "__main__":
    
    # --- Test BST ---
    print("Testing BST Index...")
    x = np.array([3.0, 1.0, 4.0, 1.5, 2.0], dtype=np.float64)
    ix = build_bst_index_C(x)
    print(f"Original: {x}\nIndices: {ix}\nSorted: {x[ix]}")
    
    # --- Test Range Query ---
    print("\nTesting BST Range Query [1.5, 3.0]...")
    matches = bst_range_query_C(x, ix, 1.5, 3.0)
    print(f"Matches: {matches} (Values: {x[matches]})")
    
    # --- Test KD-Tree ---
    print("\nTesting KD-Tree...")
    X = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], dtype=np.float64, order='F')
    kd_ix = build_kd_index_C(X, np.array([0, 1], dtype=np.int32))
    print(f"KD-Tree Indices: {kd_ix}")
    
    # --- Test Spherical KD-Tree ---
    print("\nTesting Spherical KD-Tree...")
    V = np.random.randn(3, 4)
    V = V / np.linalg.norm(V, axis=0)  # Normalize columns
    sphere_ix = build_spherical_kd_C(V, np.array([0, 1, 2], dtype=np.int32))
    print(f"Spherical Indices: {sphere_ix}")