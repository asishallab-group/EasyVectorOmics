import numpy as np
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from tensoromics_functions import tox_group_centroid

def test_centroid_logic():
    """Test the user-facing tox_group_centroid function."""
    print("=== Testing tox_group_centroid via user-facing API ===")
    
    # Test data 
    d, n_genes, n_families = 2, 5, 2
    vectors = np.array([
        [1.0, 9.0, 10.0, 20.0, 5.0],  # Dimension 1
        [1.0, 9.0, 10.0, 20.0, 5.0]   # Dimension 2
    ], dtype=np.float64, order='F')
    
    gene_to_family = np.array([1, 1, 2, 2, 1], dtype=np.int32)
    # Gene at index 1 is the only non-ortholog in Family 1
    ortholog_set = np.array([True, False, True, True, True], dtype=np.bool_)
    
    # --- Test "all" mode ---
    print("\n[tox_group_centroid] Case 1: 'all' mode")
    centroids_all = tox_group_centroid(vectors, gene_to_family, n_families, ortholog_set, mode='all')
    
    # Expected: (1+9+5)/3=5 for Fam1, (10+20)/2=15 for Fam2
    expected_all = np.array([[5.0, 15.0], [5.0, 15.0]], order='F')
    
    np.testing.assert_allclose(centroids_all, expected_all)
    assert not centroids_all.flags.writeable
    print("Test 'all' mode: PASSED")

    # --- Test "ortho" mode ---
    print("\n[tox_group_centroid] Case 2: 'ortho' mode")
    centroids_ortho = tox_group_centroid(vectors, gene_to_family, n_families, ortholog_set, mode='ortho')
    
    # Expected: (1+5)/2=3 for Fam1, (10+20)/2=15 for Fam2
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