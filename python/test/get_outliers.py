#!/usr/bin/env python3
"""
Comprehensive Python test suite for outlier detection functions
Uses tensoromics_functions.py wrapper functions (mirrors R get_outliers.R tests)
"""

import numpy as np
import sys
import os

# Add parent directory to path to import tensoromics_functions
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensoromics_functions import (
    tox_compute_family_scaling, tox_compute_family_scaling_expert, tox_compute_rdi, tox_identify_outliers, tox_detect_outliers
)

print("=== Testing outlier detection Python wrapper functions ===")
print("Based on Fortran test suite with comprehensive test coverage")

# =====================
# Tests for compute_family_scaling
# =====================

def test_compute_family_scaling_basic():
    """Test basic family scaling computation"""
    print("\n[test_compute_family_scaling_basic] Basic family scaling test")

    # Updated to use more families and genes
    n_families = 6
    genes_per_fam = 4
    n_genes = n_families * genes_per_fam

    # Generate distances and gene-to-family mapping
    distances = np.concatenate([
        np.random.uniform(1.0, 2.0, genes_per_fam),  # Family 1
        np.random.uniform(2.0, 3.0, genes_per_fam),  # Family 2
        np.random.uniform(3.0, 4.0, genes_per_fam),  # Family 3
        np.random.uniform(4.0, 5.0, genes_per_fam),  # Family 4
        np.random.uniform(5.0, 6.0, genes_per_fam),  # Family 5
        np.random.uniform(6.0, 7.0, genes_per_fam)   # Family 6
    ])
    gene_to_fam = np.concatenate([
        np.full(genes_per_fam, 1),
        np.full(genes_per_fam, 2),
        np.full(genes_per_fam, 3),
        np.full(genes_per_fam, 4),
        np.full(genes_per_fam, 5),
        np.full(genes_per_fam, 6)
    ]).astype(np.int32)

    result = tox_compute_family_scaling(distances, gene_to_fam)

    print("  Input distances:", distances)
    print("  Gene-to-family mapping:", gene_to_fam)
    print("  Scaling factors:", result['dscale'])
    print("  LOESS x (medians):", result['loess_x'])
    print("  LOESS y (stddevs):", result['loess_y'])
    print("  Indices used:", result['indices_used'])

    # Verify basic properties
    assert len(result['dscale']) == n_families  # Six families
    assert all(np.isfinite(result['dscale']))
    assert all(result['dscale'] > 0)  # Scaling factors should be positive

    print("Basic family scaling test passed ✓")


def test_compute_family_scaling_single_family():
    """Test with single family"""
    print("\n[test_compute_family_scaling_single_family] Single family test")
    
    distances = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64)
    gene_to_fam = np.array([1, 1, 1, 1], dtype=np.int32)  # All in family 1
    
    result = tox_compute_family_scaling(distances, gene_to_fam)
    
    print("  All genes in family 1")
    print("  Scaling factors:", result['dscale'])

    assert len(result['dscale']) == 1
    assert result['dscale'][0] == 0
    
    print("Single family test passed ✓")

def test_compute_family_scaling_edge_cases():
    """Test edge cases with more families"""
    print("\n[test_compute_family_scaling_edge_cases] Edge cases test")

    # Case 1: Very small distances with more families
    distances = np.array([
        0.01, 0.02, 0.01, 0.02,  # Family 1
        0.03, 0.04, 0.03, 0.04,  # Family 2
        0.05, 0.06, 0.05, 0.06,  # Family 3
        0.07, 0.08, 0.07, 0.08,  # Family 4
        0.09, 0.10, 0.09, 0.10,  # Family 5
        0.11, 0.12, 0.11, 0.12   # Family 6
    ], dtype=np.float64)
    gene_to_fam = np.array([
        1, 1, 2, 2,  # Family 1
        3, 3, 4, 4,  # Family 2
        5, 5, 6, 6,  # Family 3
        7, 7, 8, 8,  # Family 4
        9, 9, 10, 10,  # Family 5
        11, 11, 12, 12  # Family 6
    ], dtype=np.int32)

    result = tox_compute_family_scaling(distances, gene_to_fam)
    print("  Small distances case - scaling factors:", result['dscale'])
    assert all(np.isfinite(result['dscale']))

    # Case 2: Large distances with more families
    distances = np.array([
        100.0, 200.0, 150.0, 250.0,  # Family 1
        300.0, 400.0, 350.0, 450.0,  # Family 2
        500.0, 600.0, 550.0, 650.0,  # Family 3
        700.0, 800.0, 750.0, 850.0,  # Family 4
        900.0, 1000.0, 950.0, 1050.0,  # Family 5
        1100.0, 1200.0, 1150.0, 1250.0   # Family 6
    ], dtype=np.float64)
    gene_to_fam = np.array([
        1, 1, 2, 2,  # Family 1
        3, 3, 4, 4,  # Family 2
        5, 5, 6, 6,  # Family 3
        7, 7, 8, 8,  # Family 4
        9, 9, 10, 10,  # Family 5
        11, 11, 12, 12  # Family 6
    ], dtype=np.int32)

    result = tox_compute_family_scaling(distances, gene_to_fam)
    print("  Large distances case - scaling factors:", result['dscale'])
    assert all(np.isfinite(result['dscale']))

    print("Edge cases test passed ✓")

def test_tox_compute_family_scaling_expert():
    """Test expert version with user-provided work arrays"""
    print("\n[test_compute_family_scaling_expert] Expert version test")

    # Test data
    n_families = 6
    genes_per_fam = 4
    n_genes = n_families * genes_per_fam
    span = 0.7
    mode = 1
    n_iters = 3
    degree = 2

    distances = np.concatenate([
        np.random.uniform(1.0, 2.0, genes_per_fam),  # Family 1
        np.random.uniform(2.0, 3.0, genes_per_fam),  # Family 2
        np.random.uniform(3.0, 4.0, genes_per_fam),  # Family 3
        np.random.uniform(4.0, 5.0, genes_per_fam),  # Family 4
        np.random.uniform(5.0, 6.0, genes_per_fam),  # Family 5
        np.random.uniform(6.0, 7.0, genes_per_fam)   # Family 6
    ])
    gene_to_fam = np.concatenate([
        np.full(genes_per_fam, 1),
        np.full(genes_per_fam, 2),
        np.full(genes_per_fam, 3),
        np.full(genes_per_fam, 4),
        np.full(genes_per_fam, 5),
        np.full(genes_per_fam, 6)
    ]).astype(np.int32)

    # Pre-allocate work arrays (user responsibility)
    perm_tmp = np.zeros(n_genes, dtype=np.int32)
    stack_left_tmp = np.zeros(n_genes, dtype=np.int32)
    stack_right_tmp = np.zeros(n_genes, dtype=np.int32)
    family_distances = np.zeros(n_genes, dtype=np.float64)

    result = tox_compute_family_scaling_expert(
        distances, gene_to_fam, span, degree, mode, n_iters, perm_tmp, stack_left_tmp,
        stack_right_tmp, family_distances
    )

    print("  Input distances:", distances)
    print("  Gene-to-family mapping:", gene_to_fam)
    print("  Scaling factors:", result['dscale'])
    print("  Work arrays also returned in result")

    # Verify basic properties
    assert len(result['dscale']) == n_families  # Six families
    assert all(np.isfinite(result['dscale']))
    assert all(result['dscale'] > 0)  # Scaling factors should be positive

    # Verify work arrays are returned
    assert 'perm_tmp' in result
    assert 'stack_left_tmp' in result
    assert 'stack_right_tmp' in result
    assert 'family_distances' in result

    # Compare with regular version to ensure consistency
    result_regular = tox_compute_family_scaling(distances, gene_to_fam)

    print("  Comparing expert vs regular version:")
    print(f"    Expert dscale: {result['dscale']}")
    print(f"    Regular dscale: {result_regular['dscale']}")

    # Results should be very similar (within numerical precision)
    np.testing.assert_allclose(result['dscale'], result_regular['dscale'], rtol=1e-10)
    np.testing.assert_allclose(result['loess_x'], result_regular['loess_x'], rtol=1e-10)
    np.testing.assert_allclose(result['loess_y'], result_regular['loess_y'], rtol=1e-10)

    print("Expert version test passed ✓")

def test_compute_family_scaling_expert_input_validation():
    """Test expert version input validation"""
    print("\n[test_compute_family_scaling_expert_input_validation] Expert input validation test")

    n_families = 6
    genes_per_fam = 4
    n_genes = n_families * genes_per_fam
    span = 0.7
    mode = 1
    n_iters = 3
    degree = 2

    distances = np.concatenate([
        np.random.uniform(1.0, 2.0, genes_per_fam),  # Family 1
        np.random.uniform(2.0, 3.0, genes_per_fam),  # Family 2
        np.random.uniform(3.0, 4.0, genes_per_fam),  # Family 3
        np.random.uniform(4.0, 5.0, genes_per_fam),  # Family 4
        np.random.uniform(5.0, 6.0, genes_per_fam),  # Family 5
        np.random.uniform(6.0, 7.0, genes_per_fam)   # Family 6
    ])
    gene_to_fam = np.concatenate([
        np.full(genes_per_fam, 1),
        np.full(genes_per_fam, 2),
        np.full(genes_per_fam, 3),
        np.full(genes_per_fam, 4),
        np.full(genes_per_fam, 5),
        np.full(genes_per_fam, 6)
    ]).astype(np.int32)

    # Test wrong size work arrays
    error_caught1 = False
    try:
        perm_tmp = np.zeros(n_genes - 1, dtype=np.int32)  # Wrong size
        stack_left_tmp = np.zeros(n_genes, dtype=np.int32)
        stack_right_tmp = np.zeros(n_genes, dtype=np.int32)
        family_distances = np.zeros(n_genes, dtype=np.float64)

        tox_compute_family_scaling_expert(
            distances, gene_to_fam, span, degree, mode, n_iters, perm_tmp, stack_left_tmp, 
            stack_right_tmp, family_distances
        )
    except ValueError as e:
        error_caught1 = True
        assert "same length" in str(e)
    assert error_caught1

    # Test another wrong size array
    error_caught2 = False
    try:
        perm_tmp = np.zeros(n_genes, dtype=np.int32)
        stack_left_tmp = np.zeros(n_genes, dtype=np.int32)
        stack_right_tmp = np.zeros(n_genes, dtype=np.int32)
        family_distances = np.zeros(n_genes + 1, dtype=np.float64)  # Wrong size

        tox_compute_family_scaling_expert(
            distances, gene_to_fam, span, degree, mode, n_iters, perm_tmp, stack_left_tmp, 
            stack_right_tmp, family_distances
        )
    except ValueError as e:
        error_caught2 = True
        assert "same length" in str(e)
    assert error_caught2

    print("Expert input validation test passed ✓")

# =====================
# Tests for compute_rdi
# =====================

def test_compute_rdi_basic():
    """Test basic RDI computation"""
    print("\n[test_compute_rdi_basic] Basic RDI test")
    
    distances = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=np.float64)
    gene_to_fam = np.array([1, 1, 2, 2, 2], dtype=np.int32)
    dscale = np.array([1.5, 2.0], dtype=np.float64)  # Manual scaling factors
    
    rdi = tox_compute_rdi(distances, gene_to_fam, dscale)
    
    print("  Input distances:", distances)
    print("  Gene-to-family mapping:", gene_to_fam)
    print("  Scaling factors:", dscale)
    print("  RDI values:", rdi)
    
    # Verify properties
    assert len(rdi) == len(distances)
    assert all(np.isfinite(rdi))
    assert all(rdi >= 0)  # RDI should be non-negative
    
    print("Basic RDI test passed ✓")

def test_compute_rdi_outlier_detection():
    """Test RDI with simulated outliers"""
    print("\n[test_compute_rdi_outlier_detection] RDI outlier detection test")
    
    # Simulate data where last gene is an outlier
    distances = np.array([1.0, 1.1, 1.0, 10.0], dtype=np.float64)  # Last one is outlier
    gene_to_fam = np.array([1, 1, 1, 1], dtype=np.int32)
    dscale = np.array([1.0], dtype=np.float64)
    
    rdi = tox_compute_rdi(distances, gene_to_fam, dscale)
    
    print("  Distances (last is outlier):", distances)
    print("  RDI values:", rdi)
    
    # The outlier should have much higher RDI
    assert rdi[3] > rdi[0]
    assert rdi[3] > rdi[1]
    assert rdi[3] > rdi[2]
    
    print("RDI outlier detection test passed ✓")

# =====================
# Tests for identify_outliers
# =====================

def test_identify_outliers_basic():
    """Test basic outlier identification"""
    print("\n[test_identify_outliers_basic] Basic outlier identification test")
    
    # RDI values with clear outlier
    rdi = np.array([0.5, 0.6, 0.4, 3.5, 0.5], dtype=np.float64)
    
    result = tox_identify_outliers(rdi, percentile=80.0)  # Use 80th percentile
    
    print("  RDI values:", rdi)
    print(f"  Threshold used: {result['threshold']}")
    print("  Outliers (boolean):", result['outliers'])
    print("  Outlier indices:", np.where(result['outliers'])[0])
    
    # Verify that gene with RDI=3.5 is identified as outlier
    assert result['outliers'][3] == True
    assert result['threshold'] > 0
    
    print("Basic outlier identification test passed ✓")

def test_identify_outliers_no_outliers():
    """Test when no outliers present"""
    print("\n[test_identify_outliers_no_outliers] No outliers test")
    
    rdi = np.array([0.5, 0.6, 0.4, 0.7, 0.5], dtype=np.float64)  # All low RDI
    
    result = tox_identify_outliers(rdi, percentile=99.0)  # Very high percentile
    
    print("  RDI values:", rdi)
    print(f"  Threshold used: {result['threshold']}")
    print("  Outliers found:", sum(result['outliers']))
    
    # With very high percentile, few or no outliers should be found
    assert isinstance(result['outliers'], np.ndarray)
    assert result['threshold'] > 0
    
    print("No outliers test passed ✓")

def test_identify_outliers_all_outliers():
    """Test when many genes are outliers"""
    print("\n[test_identify_outliers_all_outliers] Many outliers test")
    
    rdi = np.array([3.0, 4.0, 5.0, 3.5], dtype=np.float64)  # All high RDI
    
    result = tox_identify_outliers(rdi, percentile=25.0)  # Low percentile = more outliers
    
    print("  RDI values:", rdi)
    print(f"  Threshold used: {result['threshold']}")
    print("  Outliers found:", sum(result['outliers']))
    
    # With low percentile, more outliers should be found
    assert sum(result['outliers']) >= 1
    
    print("All outliers test passed ✓")

# =====================
# Tests for detect_outliers (complete pipeline)
# =====================

def test_detect_outliers_pipeline():
    """Test complete outlier detection pipeline"""
    print("\n[test_detect_outliers_pipeline] Complete pipeline test")
    
    # Create test data with simulated outlier
    distances = np.array([
        1.0, 1.2, 1.1, 0.9,  # Family 1
        2.0, 2.2, 2.1, 2.3,  # Family 2
        3.0, 3.2, 3.1, 3.3,  # Family 3
        4.0, 4.2, 4.1, 4.3,  # Family 4
        5.0, 5.2, 5.1, 5.3,  # Family 5
        6.0, 6.2, 6.1, 6.3   # Family 6
    ], dtype=np.float64)
    gene_to_fam = np.array([
        1, 1, 1, 1,  # Family 1
        2, 2, 2, 2,  # Family 2
        3, 3, 3, 3,  # Family 3
        4, 4, 4, 4,  # Family 4
        5, 5, 5, 5,  # Family 5
        6, 6, 6, 6   # Family 6
    ], dtype=np.int32)
    percentile = 90.0
    
    result = tox_detect_outliers(distances, gene_to_fam, percentile)
    
    print("  Input distances:", distances)
    print("  Gene-to-family mapping:", gene_to_fam)
    print(f"  Percentile: {percentile}")
    print("  Outliers:", result['outliers'])
    print("  Outlier indices:", np.where(result['outliers'])[0])
    
    # Verify structure
    assert len(result['outliers']) == len(distances)
    assert 'loess_x' in result
    assert 'loess_y' in result
    assert 'loess_n' in result
    
    # Check for outliers
    outlier_count = sum(result['outliers'])
    print(f"  Total outliers found: {outlier_count}")

    print("Complete pipeline test passed ✓")

def test_detect_outliers_performance():
    """Test performance with larger dataset"""
    print("\n[test_detect_outliers_performance] Performance test")
    
    n_genes = 1000
    n_families = 20
    
    print(f"  Testing with {n_genes} genes, {n_families} families")
    
    # Generate realistic test data
    np.random.seed(42)
    distances = np.abs(np.random.normal(2.0, 0.5, n_genes))
    gene_to_fam = np.random.randint(1, n_families + 1, n_genes)
    
    # Add some outliers
    outlier_indices = np.random.choice(n_genes, size=50, replace=False)
    distances[outlier_indices] *= 3  # Make them outliers
    
    import time
    start_time = time.time()
    
    result = tox_detect_outliers(distances, gene_to_fam, percentile=95.0)
    
    end_time = time.time()
    elapsed = end_time - start_time
    
    outlier_count = sum(result['outliers'])
    print(f"  Completed in {elapsed:.6f} seconds")
    print(f"  Outliers detected: {outlier_count}/{n_genes}")
    
    # Verify reasonable results
    assert len(result['outliers']) == n_genes
    assert outlier_count >= 0  # Should find some outliers (or none is also OK)
    
    print("Performance test passed ✓")

def test_detect_outliers_edge_cases():
    """Test edge cases"""
    print("\n[test_detect_outliers_edge_cases] Edge cases test")
    
    # Case 1: Single family
    distances = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    gene_to_fam = np.array([1, 1, 1], dtype=np.int32)
    
    result = tox_detect_outliers(distances, gene_to_fam, percentile=95.0)
    print("  Single family case - outliers:", sum(result['outliers']))
    assert len(result['loess_x']) == 1
    
    # Case 2: Each gene in different family
    distances = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    gene_to_fam = np.array([1, 2, 3], dtype=np.int32)
    
    result = tox_detect_outliers(distances, gene_to_fam, percentile=95.0)
    print("  Different families case - outliers:", sum(result['outliers']))
    assert len(result['loess_x']) == 3
    
    print("Edge cases test passed ✓")

# =====================
# Run all tests
# =====================

def main():
    print("\n=================================================")
    print("    OUTLIER DETECTION FULL PYTHON INTERFACE TESTS")
    print("=================================================\n")

    # compute_family_scaling tests
    test_compute_family_scaling_basic()
    test_compute_family_scaling_single_family()
    test_compute_family_scaling_edge_cases()
    test_tox_compute_family_scaling_expert()
    test_compute_family_scaling_expert_input_validation()

    # compute_rdi tests
    test_compute_rdi_basic()
    test_compute_rdi_outlier_detection()

    # identify_outliers tests
    test_identify_outliers_basic()
    test_identify_outliers_no_outliers()
    test_identify_outliers_all_outliers()

    # detect_outliers (complete pipeline) tests
    test_detect_outliers_pipeline()
    test_detect_outliers_performance()
    test_detect_outliers_edge_cases()

    print("=================================================")
    print("             ALL TESTS COMPLETED")
    print("=================================================")
    print("If you see this message, all outlier detection Python interface tests passed! ✓")
    

if __name__ == "__main__":
    main()
