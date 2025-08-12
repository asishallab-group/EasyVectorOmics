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
        np.ctypeslib.ndpointer(dtype=np.bool_),  # ortholog_set
        np.ctypeslib.ndpointer(dtype=np.int32)  # selected_indices
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
    
    gene_to_family = np.array([1, 1, 2, 2, 1], dtype=np.int32) # Genes 1,2,5 in Fam 1; 3,4 in Fam 2
    orthologs = np.array([True, False, True, True, True], dtype=np.bool_)
    
    # --- Test "all" mode ---
    mode_all_ascii = np.array([ord(c) for c in "all"], dtype=np.int32)
    centroids_out_all = np.zeros((n_families, d), dtype=np.float64, order='F')
    selected_indices = np.zeros(n_genes, dtype=np.int32) # Workspace

    group_centroid_func(vectors, d, n_genes, gene_to_family, n_families, centroids_out_all,
                        mode_all_ascii, len("all"), orthologs, selected_indices)
    
    expected_all_fam1 = np.array([3.0, 3.0]) # mean of [1,1], [3,3], [5,5]
    expected_all_fam2 = np.array([15.0, 15.0]) # mean of [10,10], [20,20]
    
    np.testing.assert_allclose(centroids_out_all[0,:], expected_all_fam1)
    np.testing.assert_allclose(centroids_out_all[1,:], expected_all_fam2)
    print("Test 'all' mode: PASSED")

    # --- Test "orthologs" mode ---
    mode_ortho_ascii = np.array([ord(c) for c in "orthologs"], dtype=np.int32)
    centroids_out_ortho = np.zeros((n_families, d), dtype=np.float64, order='F')

    group_centroid_func(vectors, d, n_genes, gene_to_family, n_families, centroids_out_ortho,
                        mode_ortho_ascii, len("orthologs"), orthologs, selected_indices)

    expected_ortho_fam1 = np.array([3.0, 3.0]) # mean of [1,1], [5,5] (gene 2 is not ortholog)
    expected_ortho_fam2 = np.array([15.0, 15.0]) # mean of [10,10], [20,20]
    
    np.testing.assert_allclose(centroids_out_ortho[0,:], expected_ortho_fam1)
    np.testing.assert_allclose(centroids_out_ortho[1,:], expected_ortho_fam2)
    print("Test 'orthologs' mode: PASSED")

if __name__ == "__main__":
    test_centroid_logic()
    print("\nAll Python tests for gene_centroids passed! ✓")
