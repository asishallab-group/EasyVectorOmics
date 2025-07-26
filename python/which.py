"""
'which' utility for Python, like in R/MATLAB.
Includes usage examples.
"""

import numpy as np

def which(cond):
    """
    Emulates the Fortran which_c subroutine: accepts an int32 array (0/1), returns idx_out (length n, 1-based, zeros after m_out), m_out (int).
    cond: np.ndarray of dtype int32 (0/1)
    Returns: idx_out (np.ndarray, int32, length n), m_out (int)
    """
    n = cond.size
    idxs = np.where(cond != 0)[0]  # 0-based
    m_out = len(idxs)
    idx_out = np.zeros(n, dtype=np.int32)
    if m_out > 0:
        idx_out[:m_out] = idxs + 1  # Fortran 1-based
    return idx_out, m_out

if __name__ == "__main__":
    # Test cases matching the R Fortran wrapper tests
    print("\n[which] Case 1: Simple mask")
    mask = np.array([1, 0, 1, 0, 0], dtype=np.int32)
    idx_out, m_out = which(mask)
    print(idx_out)
    print(m_out)
    # Expected: idx_out = [1 3 0 0 0], m_out = 2

    print("\n[which] Case 2: All FALSE")
    mask = np.zeros(5, dtype=np.int32)
    idx_out, m_out = which(mask)
    print(idx_out)
    print(m_out)
    # Expected: idx_out = [0 0 0 0 0], m_out = 0

    print("\n[which] Case 3: All TRUE")
    mask = np.ones(5, dtype=np.int32)
    idx_out, m_out = which(mask)
    print(idx_out)
    print(m_out)
    # Expected: idx_out = [1 2 3 4 5], m_out = 5

    print("\n[which] Case 4: Empty mask")
    mask = np.zeros(0, dtype=np.int32)
    idx_out, m_out = which(mask)
    print(idx_out)
    print(m_out)
    # Expected: idx_out = [], m_out = 0
