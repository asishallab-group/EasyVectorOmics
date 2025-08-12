#!/usr/bin/env python3
"""
Test script for the gene centroid computation Fortran routines.
This script uses ctypes to call the C-compatible Fortran functions
and verifies their output against expected results.
"""

import ctypes
import numpy as np
import os
import platform

# --- 1. Library Loading ---

def _load_fortran_library():
    """
    Finds and loads the compiled Fortran shared library (.so, .dll, .dylib)
    created by the build.sh script.
    """
    # The build.sh script creates a symlink in the build directory
    lib_name = "libtensor-omics.so"
        
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    # The build.sh script places the library in the root of the build directory
    lib_path = os.path.join(project_root, "build", lib_name)

    if not os.path.exists(lib_path):
        raise OSError(f"Could not find library at '{lib_path}'. "
                      "Please build the Fortran project using the './build.sh' script.")

    # On some systems, dependent libraries like libgomp need to be preloaded.
    try:
        if platform.system() == "Linux":
            ctypes.CDLL("libgomp.so.1", mode=ctypes.RTLD_GLOBAL)
    except OSError:
        print("Warning: Could not preload libgomp.so.1. This may cause issues if OpenMP is used.")

    return ctypes.CDLL(lib_path)

# --- 2. Setup Function for the Fortran Procedure ---

def setup_group_centroid(lib):
    """Setup the group_centroid_c function interface."""
    func = lib.group_centroid_c
    func.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64), # vectors
        ctypes.c_int,                            # d
        ctypes.c_int,                            # n
        np.ctypeslib.ndpointer(dtype=np.int32), # gene_to_family_map
        ctypes.c_int,                            # num_families
        np.ctypeslib.ndpointer(dtype=np.float64), # centroid_matrix
        np.ctypeslib.ndpointer(dtype=np.int32), # mode_ascii
        ctypes.c_int,                            # mode_len
        np.ctypeslib.ndpointer(dtype=np.bool_),  # ortholog_set
        np.ctypeslib.ndpointer(dtype=np.int32), # selected_indices
        ctypes.c_int                             # selected_indices_len
    ]
    func.restype = None
    return func

# --- 3. Test Function ---

def test_centroid_logic():
    """Test the full centroid logic from Python."""
    print("=== Testing Centroid Logic ===")
    
    lib = _load_fortran_library()
    group_centroid_func = setup_group_centroid(lib)
    
    # Arrange: Common test data
    d, n_genes, n_families = 2, 5, 2
    
    # Expression vectors (d x n), Fortran-ordered
    vectors = np.array([
        [1.0, 3.0, 10.0, 20.0, 5.0],  # Dimension 1
        [1.0, 3.0, 10.0, 20.0, 5.0]   # Dimension 2
    ], dtype=np.float64, order='F')
    
    # Gene-to-family mapping (1-based for Fortran)
    gene_to_family = np.array([1, 1, 2, 2, 1], dtype=np.int32)
    
    # Orthologs mask
    orthologs = np.array([True, False, True, True, True], dtype=np.bool_)
    
    # Pre-allocated workspace buffer
    selected_indices = np.zeros(n_genes, dtype=np.int32)

    # --- Test "all" mode ---
    print("\n[group_centroid_c] Case 1: 'all' mode")
    mode_all_ascii = np.array([ord(c) for c in "all"], dtype=np.int32)
    centroids_out_all = np.zeros((n_families, d), dtype=np.float64, order='F')

    group_centroid_func(vectors, d, n_genes, gene_to_family, n_families, centroids_out_all,
                        mode_all_ascii, len("all"), orthologs, selected_indices, len(selected_indices))
    
    # Expected results for "all" mode
    # Fam 1 (genes 1,2,5): mean of [1,1], [3,3], [5,5] -> [3,3]
    # Fam 2 (genes 3,4):   mean of [10,10], [20,20] -> [15,15]
    expected_all_fam1 = np.array([3.0, 3.0])
    expected_all_fam2 = np.array([15.0, 15.0])
    
    np.testing.assert_allclose(centroids_out_all[0,:], expected_all_fam1)
    np.testing.assert_allclose(centroids_out_all[1,:], expected_all_fam2)
    print("Test 'all' mode: PASSED")

    # --- Test "orthologs" mode ---
    print("\n[group_centroid_c] Case 2: 'orthologs' mode")
    mode_ortho_ascii = np.array([ord(c) for c in "orthologs"], dtype=np.int32)
    centroids_out_ortho = np.zeros((n_families, d), dtype=np.float64, order='F')

    group_centroid_func(vectors, d, n_genes, gene_to_family, n_families, centroids_out_ortho,
                        mode_ortho_ascii, len("orthologs"), orthologs, selected_indices, len(selected_indices))

    # Expected results for "orthologs" mode
    # Fam 1 (genes 1,5 are orthologs): mean of [1,1], [5,5] -> [3,3]
    # Fam 2 (genes 3,4 are orthologs): mean of [10,10], [20,20] -> [15,15]
    expected_ortho_fam1 = np.array([3.0, 3.0])
    expected_ortho_fam2 = np.array([15.0, 15.0])
    
    np.testing.assert_allclose(centroids_out_ortho[0,:], expected_ortho_fam1)
    np.testing.assert_allclose(centroids_out_ortho[1,:], expected_ortho_fam2)
    print("Test 'orthologs' mode: PASSED")

# --- 4. Main Execution Block ---

def main():
    """Load the library and run all tests."""
    print("=================================================")
    print("      PYTHON WRAPPER TESTS FOR GENE CENTROIDS")
    print("=================================================")
    try:
        test_centroid_logic()
        print("\n=================================================")
        print("           ALL PYTHON TESTS COMPLETED")
        print("=================================================")
        print("If you see this message, all Python interface tests passed! ✓")

    except (ImportError, OSError, AssertionError) as e:
        print("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print("                  A TEST FAILED")
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print(f"ERROR: {e}")

if __name__ == "__main__":
    main()
