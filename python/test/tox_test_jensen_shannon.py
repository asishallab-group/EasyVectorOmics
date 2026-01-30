# tox_test_jensen_shannon_test.py
import numpy as np
import sys
import os
from pathlib import Path

# Add parent directory to path to import tensoromics_functions
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensoromics_functions import tox_euclidean_distance, tox_distance_to_centroid

# Import your functions
from tensoromics_functions import (
    tox_compute_gene_means,
    tox_compute_residuals,
    tox_pool_means,
    tox_construct_neighborhoods
)

def test_compute_gene_means():
    """Test tox_compute_gene_means function"""
    print("[test_compute_gene_means] Testing compute_gene_means...")
    
    # Test 1: Basic 3x3 matrix
    expr = np.array([[10, 20, 30],
                     [12, 22, 32],
                     [14, 24, 34]], dtype=np.float64, order='F')
    
    means = tox_compute_gene_means(expr)
    
    print(f"  Input shape: {expr.shape}")
    print(f"  Output shape: {means.shape}")
    print(f"  Means: {means}")
    
    # Verify calculations
    expected = np.array([12.0, 22.0, 32.0])  # (10+12+14)/3, (20+22+24)/3, (30+32+34)/3
    np.testing.assert_array_almost_equal(means, expected, decimal=10)
    assert means.shape == (3,)
    
    # Test 2: With NaN values
    expr_with_nan = np.array([[10, 20, np.nan],
                              [12, np.nan, 32],
                              [14, 24, 34]], dtype=np.float64, order='F')
    
    means_with_nan = tox_compute_gene_means(expr_with_nan)
    
    print(f"  With NaN - Means: {means_with_nan}")
    
    # Gene 1: (10+12+14)/3 = 12
    # Gene 2: (20+24)/2 = 22 (NaN excluded)
    # Gene 3: (32+34)/2 = 33 (NaN excluded)
    expected_nan = np.array([12.0, 22.0, 33.0])
    np.testing.assert_array_almost_equal(means_with_nan, expected_nan, decimal=10)
    
    print("  ✓ test_compute_gene_means passed")

def test_compute_residuals():
    """Test tox_compute_residuals function"""
    print("\n[test_compute_residuals] Testing compute_residuals...")
    
    expr = np.array([[10, 20, 30],
                     [12, 22, 32],
                     [14, 24, 34]], dtype=np.float64, order='F')
    
    means = np.array([12.0, 22.0, 32.0], dtype=np.float64, order='F')
    
    residuals = tox_compute_residuals(expr, means)
    
    print(f"  Input expr shape: {expr.shape}")
    print(f"  Input means shape: {means.shape}")
    print(f"  Output residuals shape: {residuals.shape}")
    print(f"  Residuals:\n{residuals}")
    
    # Verify calculations: residuals = expr - means
    expected = np.array([[10-12, 20-22, 30-32],
                         [12-12, 22-22, 32-32],
                         [14-12, 24-22, 34-32]], dtype=np.float64)
    
    np.testing.assert_array_almost_equal(residuals, expected, decimal=10)
    assert residuals.shape == (3, 3)
    
    # Test with NaN values
    expr_nan = np.array([[10, np.nan, 30],
                         [12, 22, np.nan],
                         [np.nan, 24, 34]], dtype=np.float64, order='F')
    
    means_nan = np.array([11.0, 23.0, 32.0], dtype=np.float64, order='F')
    
    residuals_nan = tox_compute_residuals(expr_nan, means_nan)
    
    print(f"  With NaN - Residuals:\n{residuals_nan}")
    # Check that residuals at NaN positions are also NaN
    assert np.isnan(residuals_nan[0, 1])
    assert np.isnan(residuals_nan[1, 2])
    assert np.isnan(residuals_nan[2, 0])
    
    print("  ✓ test_compute_residuals passed")

def test_pool_means():
    """Test tox_pool_means function"""
    print("\n[test_pool_means] Testing pool_means...")
    
    # Test 1: Basic pooling
    mean_S1 = np.array([1.0, 3.0, 5.0, 7.0, 9.0], dtype=np.float64, order='F')
    mean_S2 = np.array([2.0, 4.0, 6.0, 8.0, 10.0], dtype=np.float64, order='F')
    n_points = 3
    
    result = tox_pool_means(mean_S1, mean_S2, n_points)
    
    print(f"  mean_S1: {mean_S1}")
    print(f"  mean_S2: {mean_S2}")
    print(f"  n_points: {n_points}")
    print(f"  N_pool: {result['N_pool']}")
    print(f"  x_star: {result['x_star']}")
    
    # Verify N_pool (no NaN values, so should be 10)
    assert result['N_pool'] == 10
    
    # Verify x_star has correct length
    assert len(result['x_star']) == n_points
    
    # Test 2: With NaN values
    mean_S1_nan = np.array([1.0, np.nan, 5.0, 7.0, 9.0], dtype=np.float64, order='F')
    mean_S2_nan = np.array([2.0, 4.0, np.nan, 8.0, 10.0], dtype=np.float64, order='F')
    
    result_nan = tox_pool_means(mean_S1_nan, mean_S2_nan, n_points)
    
    print(f"  With NaN - N_pool: {result_nan['N_pool']}")
    print(f"  With NaN - x_star: {result_nan['x_star']}")
    
    # Should exclude NaN values: 4 from S1 + 4 from S2 = 8
    assert result_nan['N_pool'] == 8
    
    print("  ✓ test_pool_means passed")

def test_construct_neighborhoods():
    """Test tox_construct_neighborhoods function"""
    print("\n[test_construct_neighborhoods] Testing construct_neighborhoods...")
    
    # Create test data
    n_points = 3
    x_star = np.array([10.0, 20.0, 30.0], dtype=np.float64, order='F')
    
    # 5 genes, 2 replicates
    mean_S = np.array([8.0, 12.0, 18.0, 22.0, 28.0], dtype=np.float64, order='F')
    resid_S = np.array([[1.0, -1.0, 2.0, -2.0, 3.0],
                        [-1.0, 1.0, -2.0, 2.0, -3.0]], 
                       dtype=np.float64, order='F')
    
    N_pool = 10
    neighborhood_size = 50  # Explicit size
    
    result = tox_construct_neighborhoods(x_star, mean_S, resid_S, N_pool, neighborhood_size)
    
    print(f"  x_star: {x_star}")
    print(f"  mean_S: {mean_S}")
    print(f"  resid_S shape: {resid_S.shape}")
    print(f"  N_pool: {N_pool}")
    print(f"  neighborhood_size: {neighborhood_size}")
    print(f"  k_x: {result['k_x']}")
    print(f"  neighborhood_residuals shape: {result['neighborhood_residuals'].shape}")
    print(f"  neighborhood_indices shape: {result['neighborhood_indices'].shape}")
    
    # Basic assertions
    assert result['k_x'] == neighborhood_size
    assert result['neighborhood_residuals'].shape[0] == n_points
    assert result['neighborhood_indices'].shape[0] == n_points
    assert result['neighborhood_indices'].shape[1] == 1000  # Maximum from Fortran
    
    # Check that indices are 1-based (from Fortran)
    # Should be >= 1 (valid) or -1 (invalid)
    indices = result['neighborhood_indices']
    valid_mask = indices != -1
    if np.any(valid_mask):
        assert np.all(indices[valid_mask] >= 1)
    
    # Test with automatic neighborhood size
    print("\n  Testing with automatic neighborhood size...")
    result_auto = tox_construct_neighborhoods(x_star, mean_S, resid_S, N_pool, -1)
    
    print(f"  Automatic k_x: {result_auto['k_x']}")
    assert result_auto['k_x'] >= 100  # Should be at least 100 with auto calculation
    
    print("  ✓ test_construct_neighborhoods passed")

def test_integration():
    """Test integration of multiple functions"""
    print("\n[test_integration] Testing function integration...")
    
    # Create test expression data
    expr_S1 = np.array([[10.0, 20.0, 30.0, 40.0],
                        [12.0, 22.0, 32.0, 42.0],
                        [14.0, 24.0, 34.0, 44.0]], 
                       dtype=np.float64, order='F')
    
    expr_S2 = np.array([[15.0, 25.0, 35.0, 45.0],
                        [17.0, 27.0, 37.0, 47.0],
                        [19.0, 29.0, 39.0, 49.0],
                        [21.0, 31.0, 41.0, 51.0]], 
                       dtype=np.float64, order='F')
    
    # Step 1: Compute means for both studies
    print("  Step 1: Computing gene means...")
    means_S1 = tox_compute_gene_means(expr_S1)
    means_S2 = tox_compute_gene_means(expr_S2)
    
    print(f"    means_S1 shape: {means_S1.shape}, values: {means_S1}")
    print(f"    means_S2 shape: {means_S2.shape}, values: {means_S2}")
    
    assert means_S1.shape == (4,)
    assert means_S2.shape == (4,)
    
    # Step 2: Compute residuals
    print("  Step 2: Computing residuals...")
    resid_S1 = tox_compute_residuals(expr_S1, means_S1)
    resid_S2 = tox_compute_residuals(expr_S2, means_S2)
    
    print(f"    resid_S1 shape: {resid_S1.shape}")
    print(f"    resid_S2 shape: {resid_S2.shape}")
    
    assert resid_S1.shape == expr_S1.shape
    assert resid_S2.shape == expr_S2.shape
    
    # Step 3: Pool means
    print("  Step 3: Pooling means...")
    n_points = 5
    pool_result = tox_pool_means(means_S1, means_S2, n_points)
    
    print(f"    N_pool: {pool_result['N_pool']}")
    print(f"    x_star shape: {pool_result['x_star'].shape}")
    print(f"    x_star: {pool_result['x_star']}")
    
    assert pool_result['N_pool'] == 8  # 4 + 4, no NaN values
    assert len(pool_result['x_star']) == n_points
    
    # Step 4: Construct neighborhoods
    print("  Step 4: Constructing neighborhoods...")
    neighborhoods = tox_construct_neighborhoods(
        pool_result['x_star'], means_S1, resid_S1, 
        pool_result['N_pool'], neighborhood_size=50
    )
    
    print(f"    k_x: {neighborhoods['k_x']}")
    print(f"    residuals shape: {neighborhoods['neighborhood_residuals'].shape}")
    print(f"    indices shape: {neighborhoods['neighborhood_indices'].shape}")
    
    assert neighborhoods['k_x'] == 50
    assert neighborhoods['neighborhood_residuals'].shape[0] == n_points
    assert neighborhoods['neighborhood_indices'].shape[0] == n_points
    
    print("  ✓ test_integration passed")

def test_error_handling():
    """Test error handling for invalid inputs"""
    print("\n[test_error_handling] Testing error handling...")
    
    # Test 1: Dimension mismatch in compute_residuals
    print("  Testing dimension mismatch...")
    expr = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.float64, order='F')
    means_wrong = np.array([1, 2], dtype=np.float64, order='F')  # Wrong length
    
    try:
        tox_compute_residuals(expr, means_wrong)
        print("  ERROR: Should have raised ValueError!")
        assert False, "Should have raised ValueError for dimension mismatch"
    except ValueError as e:
        print(f"  ✓ Correctly raised ValueError: {e}")
    
    # Test 2: Invalid n_points in pool_means
    print("  Testing invalid n_points...")
    mean_S1 = np.array([1, 2, 3], dtype=np.float64, order='F')
    mean_S2 = np.array([4, 5, 6], dtype=np.float64, order='F')
    
    try:
        tox_pool_means(mean_S1, mean_S2, n_points=0)
        print("  ERROR: Should have raised ValueError!")
        assert False, "Should have raised ValueError for n_points=0"
    except ValueError as e:
        print(f"  ✓ Correctly raised ValueError: {e}")
    
    # Test 3: Invalid neighborhood_size
    print("  Testing invalid neighborhood_size...")
    x_star = np.array([1, 2, 3], dtype=np.float64, order='F')
    mean_S = np.array([1, 2, 3, 4, 5], dtype=np.float64, order='F')
    resid_S = np.array([[1, 2, 3, 4, 5]], dtype=np.float64, order='F')
    
    try:
        tox_construct_neighborhoods(x_star, mean_S, resid_S, N_pool=10, neighborhood_size=0)
        print("  ERROR: Should have raised ValueError!")
        assert False, "Should have raised ValueError for neighborhood_size=0"
    except ValueError as e:
        print(f"  ✓ Correctly raised ValueError: {e}")
    
    print("  ✓ test_error_handling passed")

def run_all_tests():
    """Run all tests"""
    print("=" * 60)
    print("Running tests for tox_jensen_shannon_test Python wrappers")
    print("=" * 60)
    
    test_compute_gene_means()
    test_compute_residuals()
    test_pool_means()
    test_construct_neighborhoods()
    test_integration()
    test_error_handling()
    
    print("\n" + "=" * 60)
    print("All tests completed successfully! ✓")
    print("=" * 60)

if __name__ == "__main__":
    run_all_tests()