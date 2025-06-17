import ctypes
import numpy as np

# Load compiled libraries
lib = ctypes.CDLL('Trees/build/libtrees.so')  # For KD-Tree functions

# --- BST Functions ---
# Configure BST argument types
lib.build_bst_index_C.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),  # x
    ctypes.c_int32,                                                  # n
    np.ctypeslib.ndpointer(dtype=np.int32),                          # ix (out)
    np.ctypeslib.ndpointer(dtype=np.int32),                          # stack_left
    np.ctypeslib.ndpointer(dtype=np.int32)                           # stack_right
]

lib.bst_range_query_C.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64),  # x
    np.ctypeslib.ndpointer(dtype=np.int32),    # ix
    ctypes.c_int32,                            # n
    ctypes.c_double,                           # lo
    ctypes.c_double,                           # hi
    np.ctypeslib.ndpointer(dtype=np.int32),    # out_ix (out)
    ctypes.POINTER(ctypes.c_int32)             # out_n (out)
]

# --- KD-Tree Functions ---
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
    np.ctypeslib.ndpointer(dtype=np.int32)                           # stack_right
]

# --- Test Cases ---
def test_bst():
    print("=== Testing BST ===")
    x = np.array([3.0, 1.0, 4.0, 2.0], dtype=np.float64)
    n = len(x)
    ix = np.empty(n, dtype=np.int32)
    stack_left = np.empty(n, dtype=np.int32)
    stack_right = np.empty(n, dtype=np.int32)
    
    # Build BST index
    lib.build_bst_index_C(x, n, ix, stack_left, stack_right)
    # Fortran indices are 1-based, Python is 0-based: subtract 1 for correct indexing
    ix_py = ix - 1
    print(f"BST indices (Fortran): {ix}")
    print(f"BST indices (Python 0-based): {ix_py}")
    print(f"Sorted values: {x[ix_py]}")

    # Range query
    out_ix = np.empty(n, dtype=np.int32)
    out_n = ctypes.c_int32(0)
    lib.bst_range_query_C(x, ix, n, 1.5, 3.5, out_ix, ctypes.byref(out_n))
    out_ix_py = out_ix[:out_n.value] - 1
    print(f"Range [1.5, 3.5] matches: {out_ix_py} (values: {x[out_ix_py]})")

def test_kdtree():
    print("\n=== Testing KD-Tree ===")
    X = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], dtype=np.float64, order='F')  # Column-major!
    d, n = X.shape
    dim_order = np.array([0, 1], dtype=np.int32)  # Split order: dim 0, then dim 1
    kd_ix = np.empty(n, dtype=np.int32)
    work = np.empty(n, dtype=np.int32)
    subarray = np.empty(n, dtype=np.float64)
    perm = np.empty(n, dtype=np.int32)
    stack_left = np.empty(n, dtype=np.int32)
    stack_right = np.empty(n, dtype=np.int32)
    
    # Build KD-Tree
    lib.build_kd_index_C(X, d, n, kd_ix, dim_order, work, subarray, perm, stack_left, stack_right)
    # Fortran indices are 1-based, Python is 0-based: subtract 1 for correct indexing
    kd_ix_py = kd_ix - 1
    print(f"KD-Tree indices (Fortran): {kd_ix}")
    print(f"KD-Tree indices (Python 0-based): {kd_ix_py}")
    print(f"Points in order:\n{X[:, kd_ix_py]}")

def test_bst_edge_cases():
    print("\n=== BST Edge Cases ===")
    # Empty array
    x = np.array([], dtype=np.float64)
    n = len(x)
    ix = np.empty(n, dtype=np.int32)
    stack_left = np.empty(n, dtype=np.int32)
    stack_right = np.empty(n, dtype=np.int32)
    try:
        lib.build_bst_index_C(x, n, ix, stack_left, stack_right)
        print("Empty array: No error (expected: no-op or error)")
    except Exception as e:
        print(f"Empty array: Exception caught as expected: {e}")

    # Single element
    x = np.array([42.0], dtype=np.float64)
    n = len(x)
    ix = np.empty(n, dtype=np.int32)
    stack_left = np.empty(n, dtype=np.int32)
    stack_right = np.empty(n, dtype=np.int32)
    lib.build_bst_index_C(x, n, ix, stack_left, stack_right)
    print(f"Single element BST indices: {ix}")

    # Wrong dtype (int instead of float)
    try:
        x_wrong = np.array([1, 2, 3], dtype=np.int32)
        ix_wrong = np.empty(3, dtype=np.int32)
        stack_left_wrong = np.empty(3, dtype=np.int32)
        stack_right_wrong = np.empty(3, dtype=np.int32)
        lib.build_bst_index_C(x_wrong, 3, ix_wrong, stack_left_wrong, stack_right_wrong)
        print("Wrong dtype: No error (unexpected, should fail)")
    except Exception as e:
        print(f"Wrong dtype: Exception caught as expected: {e}")

def test_kdtree_edge_cases():
    print("\n=== KD-Tree Edge Cases ===")
    # Empty matrix
    X = np.empty((2, 0), dtype=np.float64, order='F')
    d, n = X.shape
    kd_ix = np.empty(n, dtype=np.int32)
    dim_order = np.array([0, 1], dtype=np.int32)
    work = np.empty(n, dtype=np.int32)
    subarray = np.empty(n, dtype=np.float64)
    perm = np.empty(n, dtype=np.int32)
    stack_left = np.empty(n, dtype=np.int32)
    stack_right = np.empty(n, dtype=np.int32)
    try:
        lib.build_kd_index_C(X, d, n, kd_ix, dim_order, work, subarray, perm, stack_left, stack_right)
        print("Empty matrix: No error (expected: no-op or error)")
    except Exception as e:
        print(f"Empty matrix: Exception caught as expected: {e}")

    # Single point
    X = np.array([[1.0], [2.0]], dtype=np.float64, order='F')
    d, n = X.shape
    kd_ix = np.empty(n, dtype=np.int32)
    dim_order = np.array([0, 1], dtype=np.int32)
    work = np.empty(n, dtype=np.int32)
    subarray = np.empty(n, dtype=np.float64)
    perm = np.empty(n, dtype=np.int32)
    stack_left = np.empty(n, dtype=np.int32)
    stack_right = np.empty(n, dtype=np.int32)
    lib.build_kd_index_C(X, d, n, kd_ix, dim_order, work, subarray, perm, stack_left, stack_right)
    print(f"Single point KD-Tree indices: {kd_ix}")

    # Wrong dtype (float32 instead of float64)
    try:
        X_wrong = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32, order='F')
        kd_ix_wrong = np.empty(2, dtype=np.int32)
        dim_order_wrong = np.array([0, 1], dtype=np.int32)
        work_wrong = np.empty(2, dtype=np.int32)
        subarray_wrong = np.empty(2, dtype=np.float64)
        perm_wrong = np.empty(2, dtype=np.int32)
        stack_left_wrong = np.empty(2, dtype=np.int32)
        stack_right_wrong = np.empty(2, dtype=np.int32)
        lib.build_kd_index_C(X_wrong, 2, 2, kd_ix_wrong, dim_order_wrong, work_wrong, subarray_wrong, perm_wrong, stack_left_wrong, stack_right_wrong)
        print("Wrong dtype (float32): No error (unexpected, should fail)")
    except Exception as e:
        print(f"Wrong dtype (float32): Exception caught as expected: {e}")

    # Wrong shape (not Fortran contiguous)
    try:
        X_wrong_shape = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float64, order='C')
        kd_ix_wrong = np.empty(2, dtype=np.int32)
        dim_order_wrong = np.array([0, 1], dtype=np.int32)
        work_wrong = np.empty(2, dtype=np.int32)
        subarray_wrong = np.empty(2, dtype=np.float64)
        perm_wrong = np.empty(2, dtype=np.int32)
        stack_left_wrong = np.empty(2, dtype=np.int32)
        stack_right_wrong = np.empty(2, dtype=np.int32)
        lib.build_kd_index_C(X_wrong_shape, 2, 2, kd_ix_wrong, dim_order_wrong, work_wrong, subarray_wrong, perm_wrong, stack_left_wrong, stack_right_wrong)
        print("Wrong shape (C order): No error (unexpected, should fail)")
    except Exception as e:
        print(f"Wrong shape (C order): Exception caught as expected: {e}")

if __name__ == "__main__":
    test_bst()
    test_kdtree()
    test_bst_edge_cases()
    test_kdtree_edge_cases()