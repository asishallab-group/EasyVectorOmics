import numpy as np
import os
import sys

# Add the parent directory to the path to find the tensoromics_functions module
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from tensoromics_functions import tox_group_centroid

def test_centroid_logic():
    """Test the user-facing tox_group_centroid function."""
    print("=== Testing tox_group_centroid via user-facing API ===")
    
    # Arrange: Test data
    d, n_genes, n_families = 2, 5, 2
    vectors = np.array([
        [1.0, 3.0, 10.0, 20.0, 5.0],  # Dimension 1
        [1.0, 3.0, 10.0, 20.0, 5.0]   # Dimension 2
    ], dtype=np.float64, order='F')
    
    gene_to_family = np.array([1, 1, 2, 2, 1], dtype=np.int32)
    # CORRECTED: The ortholog set must be True for genes 1 and 5 (indices 0 and 4)
    # to produce the expected centroid for family 1.
    ortholog_set = np.array([True, False, True, True, True], dtype=np.bool_)
    
    # --- Test "all" mode ---
    print("\n[tox_group_centroid] Case 1: 'all' mode")
    centroids_all = tox_group_centroid(vectors, gene_to_family, n_families, ortholog_set, mode='all')
    
    expected_all = np.array([[3.0, 15.0], [3.0, 15.0]], order='F')
    
    np.testing.assert_allclose(centroids_all, expected_all)
    assert not centroids_all.flags.writeable
    print("Test 'all' mode: PASSED")

    # --- Test "ortho" mode ---
    print("\n[tox_group_centroid] Case 2: 'ortho' mode")
    centroids_ortho = tox_group_centroid(vectors, gene_to_family, n_families, ortholog_set, mode='ortho')
    
    expected_ortho = np.array([[3.0, 15.0], [3.0, 15.0]], order='F')

    np.testing.assert_allclose(centroids_ortho, expected_ortho)
    assert not centroids_ortho.flags.writeable
    print("Test 'ortho' mode: PASSED")

def main():
    """Run all tests."""
    print("=================================================")
    print("     PYTHON API TESTS FOR GENE CENTROIDS")
    print("=================================================")
    try:
        test_centroid_logic()
        print("\n=================================================")
        print("           ALL PYTHON TESTS COMPLETED")
        print("=================================================")
        print("If you see this message, all Python API tests passed! ✓")

    except Exception as e:
        print("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print("                  A TEST FAILED")
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print(f"ERROR: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
