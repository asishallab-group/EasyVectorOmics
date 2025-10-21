import numpy as np
import sys
import os

# Add the path to your compiled shared library
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

# Import your wrapper functions
from tensoromics_functions import (
    calc_spike_thresholds,
    calc_spike_thresholds_alloc,
    calc_integrated_threshold,
    calc_integrated_threshold_alloc,
    detect_outliers_integrated,
    detect_outliers_spike
)

def create_sample_spike_data():
    """Create sample spike contributions data: 5 samples × 3 genes"""
    data = np.array([
        [1.0, 2.0, 3.0],  # Sample 1
        [2.0, 4.0, 6.0],  # Sample 2  
        [3.0, 6.0, 9.0],  # Sample 3
        [4.0, 8.0, 12.0], # Sample 4
        [5.0, 10.0, 15.0] # Sample 5
    ], dtype=np.float64, order='F')
    return data

def create_sample_permutation():
    """Create sample permutation array for spike data"""
    perm = np.array([
        [1, 1, 1],  # Gene 1 sorted indices
        [2, 2, 2],  # Gene 2 sorted indices
        [3, 3, 3],  # Gene 3 sorted indices
        [4, 4, 4],
        [5, 5, 5]
    ], dtype=np.int32, order='F')
    return perm

def create_sample_integrated_data():
    """Create sample integrated contributions data"""
    data = np.array([10.0, 20.0, 30.0, 40.0, 50.0], 
                   dtype=np.float64, order='F')
    return data

def create_sample_integrated_perm():
    """Create sample permutation for integrated data"""
    perm = np.array([1, 2, 3, 4, 5], dtype=np.int32, order='F')
    return perm

def test_calc_spike_thresholds_basic():
    """Test basic spike threshold calculation with pre-computed permutation"""
    print("Testing calc_spike_thresholds_basic...")
    
    spike_data = create_sample_spike_data()
    permutation = create_sample_permutation()
    
    thresholds = calc_spike_thresholds(spike_data, 80.0, permutation)
    
    # Validate output
    assert thresholds.shape == (3,), f"Expected shape (3,), got {thresholds.shape}"
    assert thresholds.dtype == np.float64, f"Expected float64, got {thresholds.dtype}"
    assert not thresholds.flags.writeable, "Output should be read-only"
    
    # Expected 80th percentiles
    expected = np.array([4.2, 8.4, 12.6], dtype=np.float64)
    np.testing.assert_array_almost_equal(thresholds, expected, decimal=10)
    
    print("calc_spike_thresholds_basic passed")

def test_calc_spike_thresholds_alloc_basic():
    """Test spike threshold calculation with internal allocation"""
    print("Testing calc_spike_thresholds_alloc_basic...")
    
    spike_data = create_sample_spike_data()
    thresholds = calc_spike_thresholds_alloc(spike_data, 80.0)
    
    # Validate output
    assert thresholds.shape == (3,), f"Expected shape (3,), got {thresholds.shape}"
    assert thresholds.dtype == np.float64, f"Expected float64, got {thresholds.dtype}"
    assert not thresholds.flags.writeable, "Output should be read-only"
    
    # Should get same results as performance version
    expected = np.array([4.2, 8.4, 12.6], dtype=np.float64)
    np.testing.assert_array_almost_equal(thresholds, expected, decimal=10)
    
    print("calc_spike_thresholds_alloc_basic passed")

def test_calc_integrated_threshold_basic():
    """Test basic integrated threshold calculation"""
    print("Testing calc_integrated_threshold_basic...")
    
    contributions = create_sample_integrated_data()
    permutation = create_sample_integrated_perm()
    
    threshold = calc_integrated_threshold(contributions, 80.0, permutation)
    
    # Validate output
    assert isinstance(threshold, float), f"Expected float, got {type(threshold)}"
    # 80th percentile of [10,20,30,40,50] is 42.0
    assert abs(threshold - 42.0) < 1e-10, f"Expected 42.0, got {threshold}"
    
    print("calc_integrated_threshold_basic passed")

def test_calc_integrated_threshold_alloc_basic():
    """Test integrated threshold calculation with internal allocation"""
    print("Testing calc_integrated_threshold_alloc_basic...")
    
    contributions = create_sample_integrated_data()
    threshold = calc_integrated_threshold_alloc(contributions, 80.0)
    
    # Validate output
    assert isinstance(threshold, float), f"Expected float, got {type(threshold)}"
    assert abs(threshold - 42.0) < 1e-10, f"Expected 42.0, got {threshold}"
    
    print("calc_integrated_threshold_alloc_basic passed")

def test_detect_outliers_integrated_basic():
    """Test integrated outlier detection"""
    print("Testing detect_outliers_integrated_basic...")
    
    contributions = create_sample_integrated_data()
    threshold = 35.0
    outliers = detect_outliers_integrated(contributions, threshold)
    
    # Validate output
    assert outliers.shape == (5,), f"Expected shape (5,), got {outliers.shape}"
    assert outliers.dtype == bool, f"Expected bool, got {outliers.dtype}"
    assert not outliers.flags.writeable, "Output should be read-only"
    
    # Values > 35.0 should be outliers: [10,20,30,40,50] -> [F,F,F,T,T]
    expected = np.array([False, False, False, True, True])
    np.testing.assert_array_equal(outliers, expected)
    
    print("detect_outliers_integrated_basic passed")

def test_detect_outliers_spike_basic():
    """Test spike outlier detection"""
    print("Testing detect_outliers_spike_basic...")
    
    spike_data = create_sample_spike_data()
    thresholds = np.array([3.5, 7.0, 10.5], dtype=np.float64, order='F')
    outliers = detect_outliers_spike(spike_data, thresholds)
    
    # Validate output
    assert outliers.shape == (5, 3), f"Expected shape (5, 3), got {outliers.shape}"
    assert outliers.dtype == bool, f"Expected bool, got {outliers.dtype}"
    assert not outliers.flags.writeable, "Output should be read-only"
    
    # Expected outliers pattern
    expected = np.array([
        [False, False, False],
        [False, False, False], 
        [False, False, False],
        [True,  True,  True],
        [True,  True,  True]
    ])
    np.testing.assert_array_equal(outliers, expected)
    
    print("detect_outliers_spike_basic passed")

def test_spike_thresholds_edge_cases():
    """Test spike thresholds with edge cases like duplicates and unsorted data"""
    print("Testing spike_thresholds_edge_cases...")
    
    # Data with duplicates and unsorted values
    spike_data = np.array([
        [5.0, 2.0],  # Sample 1
        [1.0, 2.0],  # Sample 2
        [5.0, 8.0],  # Sample 3  
        [3.0, 2.0],  # Sample 4
    ], dtype=np.float64, order='F')
    
    thresholds = calc_spike_thresholds_alloc(spike_data, 50.0)  # Median
    
    # Gene1: [1,3,5,5] median = 4.0, Gene2: [2,2,2,8] median = 2.0
    expected = np.array([4.0, 2.0], dtype=np.float64)
    np.testing.assert_array_almost_equal(thresholds, expected, decimal=10)
    
    print("spike_thresholds_edge_cases passed")

def test_integrated_threshold_edge_cases():
    """Test integrated thresholds with edge cases"""
    print("Testing integrated_threshold_edge_cases...")
    
    # Single element
    single_data = np.array([100.0], dtype=np.float64, order='F')
    threshold = calc_integrated_threshold_alloc(single_data, 50.0)
    assert abs(threshold - 100.0) < 1e-10, f"Expected 100.0, got {threshold}"
    
    # All identical values
    identical_data = np.array([5.0, 5.0, 5.0, 5.0], dtype=np.float64, order='F')
    threshold = calc_integrated_threshold_alloc(identical_data, 75.0)
    assert abs(threshold - 5.0) < 1e-10, f"Expected 5.0, got {threshold}"
    
    print("integrated_threshold_edge_cases passed")

def test_outlier_detection_edge_cases():
    """Test outlier detection with values exactly at threshold"""
    print("Testing outlier_detection_edge_cases...")
    
    # Values exactly at threshold should NOT be outliers
    contributions = np.array([5.0, 10.0, 10.0, 15.0], dtype=np.float64, order='F')
    threshold = 10.0
    
    outliers = detect_outliers_integrated(contributions, threshold)
    expected = np.array([False, False, False, True])  # Only > 10.0 is outlier
    np.testing.assert_array_equal(outliers, expected)
    
    print("outlier_detection_edge_cases passed")

def test_invalid_percentiles():
    """Test error handling for invalid percentile values"""
    print("Testing invalid_percentiles...")
    
    spike_data = create_sample_spike_data()
    
    try:
        calc_spike_thresholds_alloc(spike_data, -10.0)
        assert False, "Should have raised ValueError for negative percentile"
    except ValueError as e:
        assert "percentile_val must be between 0.0 and 100.0" in str(e)
    
    try:
        calc_spike_thresholds_alloc(spike_data, 150.0)
        assert False, "Should have raised ValueError for percentile > 100"
    except ValueError as e:
        assert "percentile_val must be between 0.0 and 100.0" in str(e)
    
    print("invalid_percentiles passed")

def test_dimension_mismatch_errors():
    """Test error handling for dimension mismatches"""
    print("Testing dimension_mismatch_errors...")
    
    spike_data = create_sample_spike_data()
    permutation = create_sample_permutation()
    
    # Wrong permutation dimensions
    wrong_perm = np.array([[1, 2], [3, 4]], dtype=np.int32, order='F')  # 2×2 vs 5×3
    try:
        calc_spike_thresholds(spike_data, 50.0, wrong_perm)
        assert False, "Should have raised ValueError for wrong permutation dimensions"
    except ValueError as e:
        assert "permutation shape" in str(e)
    
    # Wrong thresholds dimensions for spike outliers
    wrong_thresholds = np.array([1.0, 2.0], dtype=np.float64, order='F')  # 2 vs 3 genes
    try:
        detect_outliers_spike(spike_data, wrong_thresholds)
        assert False, "Should have raised ValueError for wrong thresholds dimensions"
    except ValueError as e:
        assert "thresholds size" in str(e)
    
    print("dimension_mismatch_errors passed")

def test_memory_layout_handling():
    """Test that functions handle different memory layouts correctly"""
    print("Testing memory_layout_handling...")
    
    # Create data in C order (row-major)
    spike_data_c = np.array([
        [1.0, 2.0, 3.0],
        [2.0, 4.0, 6.0], 
    ], dtype=np.float64, order='C')  # C order
    
    # Function should convert to Fortran order internally
    thresholds = calc_spike_thresholds_alloc(spike_data_c, 50.0)
    assert thresholds.shape == (3,), f"Expected shape (3,), got {thresholds.shape}"
    assert thresholds.dtype == np.float64, f"Expected float64, got {thresholds.dtype}"
    
    print("memory_layout_handling passed")

def test_data_type_conversion():
    """Test automatic data type conversion"""
    print("Testing data_type_conversion...")
    
    # Test with different input types that should be converted to float64
    spike_data_float32 = np.array([
        [1.0, 2.0],
        [3.0, 4.0]
    ], dtype=np.float32, order='F')
    
    thresholds = calc_spike_thresholds_alloc(spike_data_float32, 50.0)
    assert thresholds.dtype == np.float64, f"Expected float64, got {thresholds.dtype}"
    
    # Test with integer input
    spike_data_int = np.array([
        [1, 2],
        [3, 4]
    ], dtype=np.int32, order='F')
    
    thresholds = calc_spike_thresholds_alloc(spike_data_int, 50.0)
    assert thresholds.dtype == np.float64, f"Expected float64, got {thresholds.dtype}"
    
    print("data_type_conversion passed")

def test_large_scale_performance():
    """Test with larger datasets to ensure no memory issues"""
    print("Testing large_scale_performance...")
    
    # Use smaller sizes for unit tests, but test the pattern
    n_samples, n_genes = 100, 50
    large_data = np.random.rand(n_samples, n_genes).astype(np.float64, order='F')
    
    # These should complete without errors
    thresholds = calc_spike_thresholds_alloc(large_data, 95.0)
    assert thresholds.shape == (n_genes,), f"Expected shape ({n_genes},), got {thresholds.shape}"
    
    outliers = detect_outliers_spike(large_data, thresholds)
    assert outliers.shape == (n_samples, n_genes), f"Expected shape ({n_samples}, {n_genes}), got {outliers.shape}"
    
    print("large_scale_performance passed")

def test_consistency_between_versions():
    """Test that alloc and non-alloc versions give consistent results"""
    print("Testing consistency_between_versions...")
    
    spike_data = create_sample_spike_data()
    permutation = create_sample_permutation()
    
    # Test spike thresholds consistency
    thresholds_perf = calc_spike_thresholds(spike_data, 80.0, permutation)
    thresholds_alloc = calc_spike_thresholds_alloc(spike_data, 80.0)
    
    np.testing.assert_array_almost_equal(thresholds_perf, thresholds_alloc, decimal=10)
    
    # Test integrated thresholds consistency
    integrated_data = np.sum(spike_data, axis=1)
    integrated_perm = np.array([1, 2, 3, 4, 5], dtype=np.int32, order='F')
    
    threshold_perf = calc_integrated_threshold(integrated_data, 75.0, integrated_perm)
    threshold_alloc = calc_integrated_threshold_alloc(integrated_data, 75.0)
    
    assert abs(threshold_perf - threshold_alloc) < 1e-10, "Thresholds should be consistent"
    
    print("consistency_between_versions passed")

def run_all_tests():
    """Run all test functions"""
    print("Running comprehensive tests for outlier detection Python wrappers...")
    print("=" * 60)
    
    test_functions = [
        test_calc_spike_thresholds_basic,
        test_calc_spike_thresholds_alloc_basic,
        test_calc_integrated_threshold_basic,
        test_calc_integrated_threshold_alloc_basic,
        test_detect_outliers_integrated_basic,
        test_detect_outliers_spike_basic,
        test_spike_thresholds_edge_cases,
        test_integrated_threshold_edge_cases,
        test_outlier_detection_edge_cases,
        test_invalid_percentiles,
        test_dimension_mismatch_errors,
        test_memory_layout_handling,
        test_data_type_conversion,
        test_large_scale_performance,
        test_consistency_between_versions
    ]
    
    passed = 0
    failed = 0
    
    for test_func in test_functions:
        try:
            test_func()
            passed += 1
        except Exception as e:
            failed += 1
            print(f"✗ {test_func.__name__} failed: {e}")
            import traceback
            traceback.print_exc()
        print()
    
    print("=" * 60)
    print(f"Test Results: {passed} passed, {failed} failed")
    
    if failed == 0:
        print("All tests passed successfully!")
    else:
        print(f"{failed} tests failed")
    
    return failed == 0

if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)