import ctypes
import numpy as np
import os
import platform

def _load_fortran_library():
    lib_name = "libtensor-omics.so"
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    lib_path = os.path.join(project_root, "build", lib_name)
    if not os.path.exists(lib_path):
        raise OSError(f"Could not find library at '{lib_path}'. Please run './build.sh'.")
    try:
        if platform.system() == "Linux":
            ctypes.CDLL("libgomp.so.1", mode=ctypes.RTLD_GLOBAL)
    except OSError:
        print("Warning: Could not preload libgomp.so.1.")
    return ctypes.CDLL(lib_path)

def setup_group_centroid(lib):
    """Setup group_centroid_c function"""
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
        # CORRECTED: The ortholog set is passed as an integer array (0 or 1)
        np.ctypeslib.ndpointer(dtype=np.int32),  # ortholog_set_int
        np.ctypeslib.ndpointer(dtype=np.int32), # selected_indices
        ctypes.c_int                             # selected_indices_len
    ]
    func.restype = None
    return func

def test_centroid_logic():
    """Test the full centroid logic from Python"""
    print("=== Testing Centroid Logic ===")
    
    lib = _load_fortran_library()
    group_centroid_func = setup_group_centroid(lib)
    
    # Arrange: Test data
    d, n_genes, n_families = 2, 5, 2
    vectors = np.array([
        [1.0, 3.0, 10.0, 20.0, 5.0],  # Dimension 1
        [1.0, 3.0, 10.0, 20.0, 5.0]   # Dimension 2
    ], dtype=np.float64, order='F')
    
    gene_to_family = np.array([1, 1, 2, 2, 1], dtype=np.int32)
    # CORRECTED: Create the ortholog mask as an integer array (1 for True, 0 for False)
    orthologs_int = np.array([1, 0, 1, 1, 1], dtype=np.int32)
    
    # Pre-allocated workspace buffer
    selected_indices = np.zeros(n_genes, dtype=np.int32)

    # --- Test "all" mode ---
    print("\n[group_centroid_c] Case 1: 'all' mode")
    mode_all_ascii = np.array([ord(c) for c in "all"], dtype=np.int32)
    centroids_out_all = np.zeros((n_families, d), dtype=np.float64, order='F')

    group_centroid_func(vectors, d, n_genes, gene_to_family, n_families, centroids_out_all,
                        mode_all_ascii, len("all"), orthologs_int, selected_indices, len(selected_indices))
    
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
                        mode_ortho_ascii, len("orthologs"), orthologs_int, selected_indices, len(selected_indices))

    expected_ortho_fam1 = np.array([3.0, 3.0])
    expected_ortho_fam2 = np.array([15.0, 15.0])
    
    np.testing.assert_allclose(centroids_out_ortho[0,:], expected_ortho_fam1)
    np.testing.assert_allclose(centroids_out_ortho[1,:], expected_ortho_fam2)
    print("Test 'orthologs' mode: PASSED")

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
