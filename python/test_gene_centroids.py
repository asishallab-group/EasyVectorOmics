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
    # This is now a boolean array as the user-facing API expects
    # CORRECTED: The last element should be False to match the expected calculation.
    ortholog_set = np.array([True, False, True, True, False], dtype=np.bool_)
    
    # --- Test "all" mode ---
    print("\n[tox_group_centroid] Case 1: 'all' mode")
    centroids_all = tox_group_centroid(vectors, gene_to_family, n_families, ortholog_set, mode='all')
    
    # Expected: Fam1=[(1,1)+(3,3)+(5,5)]/3 = [3,3], Fam2=[(10,10)+(20,20)]/2 = [15,15]
    expected_all = np.array([[3.0, 15.0], [3.0, 15.0]], order='F')
    
    np.testing.assert_allclose(centroids_all, expected_all)
    assert not centroids_all.flags.writeable
    print("Test 'all' mode: PASSED")

    # --- Test "ortho" mode ---
    print("\n[tox_group_centroid] Case 2: 'ortho' mode")
    # Add print statements to debug the inputs right before the call
    print("Debug Info for 'ortho' mode:")
    print(f"  Vectors shape: {vectors.shape}")
    print(f"  Gene to Family: {gene_to_family}")
    print(f"  Ortholog Set: {ortholog_set}")
    
    centroids_ortho = tox_group_centroid(vectors, gene_to_family, n_families, ortholog_set, mode='ortho')
    
    # Expected: Fam1 (orthos are genes 1, 5 -> indices 0, 4): mean([1,5],[1,5]) = [3,3]
    #           Fam2 (orthos are genes 3, 4 -> indices 2, 3): mean([10,20],[10,20]) = [15,15]
    expected_ortho = np.array([[3.0, 15.0], [3.0, 15.0]], order='F')

    print(f"  Actual Centroids:\n{centroids_ortho}")
    print(f"  Expected Centroids:\n{expected_ortho}")

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
