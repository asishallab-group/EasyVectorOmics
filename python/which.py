"""
'which' utility for Python, like in R/MATLAB.
Includes usage examples.
"""

import numpy as np

def which(cond):
    """
    Calls the Fortran which_c subroutine from libtensor-omics.so using ctypes.
    cond: np.ndarray of dtype int32 (0/1)
    Returns: idx_out (np.ndarray, int32, length n), m_out (int)
    """
    import ctypes
    import numpy as np
    n = cond.size
    idx_out = np.zeros(n, dtype=np.int32)
    m_max = n
    m_out = np.zeros(1, dtype=np.int32)
    # Load the shared library
    lib = ctypes.CDLL("build/libtensor-omics.so")
    # Prepare argument types
    which_c = lib.which_c
    which_c.argtypes = [
        ctypes.POINTER(ctypes.c_int32), # mask
        ctypes.c_int32,                 # n
        ctypes.POINTER(ctypes.c_int32), # idx_out
        ctypes.c_int32,                 # m_max
        ctypes.POINTER(ctypes.c_int32)  # m_out
    ]
    # Call Fortran subroutine
    which_c(
        cond.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
        ctypes.c_int32(n),
        idx_out.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
        ctypes.c_int32(m_max),
        m_out.ctypes.data_as(ctypes.POINTER(ctypes.c_int32))
    )
    return idx_out, int(m_out[0])

if __name__ == "__main__":
    # Test cases using the Fortran which_c subroutine via ctypes
    print("\n[which] Case 1: Simple mask (calls Fortran)")
    mask = np.array([1, 0, 1, 0, 0], dtype=np.int32)
    idx_out, m_out = which(mask)
    print(f"idx_out: {idx_out}")
    print(f"m_out: {m_out}")
    # Expected: idx_out = [1 3 0 0 0], m_out = 2

    print("\n[which] Case 2: All FALSE (calls Fortran)")
    mask = np.zeros(5, dtype=np.int32)
    idx_out, m_out = which(mask)
    print(f"idx_out: {idx_out}")
    print(f"m_out: {m_out}")
    # Expected: idx_out = [0 0 0 0 0], m_out = 0

    print("\n[which] Case 3: All TRUE (calls Fortran)")
    mask = np.ones(5, dtype=np.int32)
    idx_out, m_out = which(mask)
    print(f"idx_out: {idx_out}")
    print(f"m_out: {m_out}")
    # Expected: idx_out = [1 2 3 4 5], m_out = 5

    print("\n[which] Case 4: Empty mask (calls Fortran)")
    mask = np.zeros(0, dtype=np.int32)
    idx_out, m_out = which(mask)
    print(f"idx_out: {idx_out}")
    print(f"m_out: {m_out}")
    # Expected: idx_out = [], m_out = 0
