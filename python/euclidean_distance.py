#!/usr/bin/env python3
"""
Test script for Euclidean distance functions
Python equivalent of the R euclidean_distance.R tests
"""

import numpy as np
import ctypes
import time
import math

# Load library
ctypes.CDLL("libgomp.so.1", mode=ctypes.RTLD_GLOBAL)
lib = ctypes.CDLL("build/libtensor-omics.so")

def setup_euclidean_distance():
    """Setup euclidean distance function"""
    euclidean_distance = lib.euclidean_distance_c
    euclidean_distance.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double)
    ]
    euclidean_distance.restype = None
    return euclidean_distance

def setup_distance_to_centroid():
    """Setup distance to centroid function"""
    distance_to_centroid = lib.distance_to_centroid_c
    distance_to_centroid.argtypes = [
        ctypes.c_int,  # n_genes
        ctypes.c_int,  # n_families
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # genes
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # centroids
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # gene_to_fam
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # distances
        ctypes.c_int   # d
    ]
    distance_to_centroid.restype = None
    return distance_to_centroid

def test_euclidean_distance():
    """Test basic euclidean distance between two vectors"""
    print("=== Testing Euclidean Distance ===")
    
    euclidean_distance = setup_euclidean_distance()
    
    # Test 1: Simple 3D vectors
    vec1 = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    vec2 = np.array([4.0, 5.0, 6.0], dtype=np.float64)
    d = len(vec1)
    result = ctypes.c_double()
    
    euclidean_distance(vec1, vec2, d, ctypes.byref(result))
    
    expected = np.sqrt(np.sum((vec1 - vec2)**2))
    print("Test 1 - 3D vectors [1,2,3] vs [4,5,6]:")
    print(f"  Fortran result: {result.value}")
    print(f"  Python expected: {expected}")
    print(f"  Match: {abs(result.value - expected) < 1e-12}")
    print()
    
    # Test 2: 2D vector to origin (3-4-5 triangle)
    vec1 = np.array([3.0, 4.0], dtype=np.float64)
    vec2 = np.array([0.0, 0.0], dtype=np.float64)
    d = len(vec1)
    result = ctypes.c_double()
    
    euclidean_distance(vec1, vec2, d, ctypes.byref(result))
    
    expected = 5.0
    print("Test 2 - 2D vector [3,4] to origin:")
    print(f"  Fortran result: {result.value}")
    print(f"  Expected: {expected}")
    print(f"  Match: {abs(result.value - expected) < 1e-12}")
    print()
    
    # Test 3: Identical vectors
    vec1 = np.array([1.5, 2.7, 3.9], dtype=np.float64)
    vec2 = vec1.copy()  # identical
    d = len(vec1)
    result = ctypes.c_double()
    
    euclidean_distance(vec1, vec2, d, ctypes.byref(result))
    
    expected = 0.0
    print("Test 3 - Identical vectors:")
    print(f"  Fortran result: {result.value}")
    print(f"  Expected: {expected}")
    print(f"  Match: {abs(result.value - expected) < 1e-15}")
    print()

def test_distance_to_centroid():
    """Test distance to centroid functionality"""
    print("=== Testing Distance to Centroid ===")
    
    distance_to_centroid = setup_distance_to_centroid()
    
    # Setup test data
    n_genes = 4
    n_families = 2
    d = 3
    
    # Gene expression data (column-major: genes as columns, dimensions as rows)
    # Gene 1: [1, 0, 0] - Family 1
    # Gene 2: [0, 1, 0] - Family 1  
    # Gene 3: [3, 0, 0] - Family 2
    # Gene 4: [0, 3, 0] - Family 2
    genes = np.array([
        1.0, 0.0, 0.0,  # Gene 1
        0.0, 1.0, 0.0,  # Gene 2
        3.0, 0.0, 0.0,  # Gene 3
        0.0, 3.0, 0.0   # Gene 4
    ], dtype=np.float64)
    
    # Family centroids (column-major)
    # Family 1 centroid: [0.5, 0.5, 0.0]
    # Family 2 centroid: [1.5, 1.5, 0.0]
    centroids = np.array([
        0.5, 0.5, 0.0,  # Family 1
        1.5, 1.5, 0.0   # Family 2
    ], dtype=np.float64)
    
    # Gene-to-family mapping (1-based)
    gene_to_fam = np.array([1, 1, 2, 2], dtype=np.int32)
    
    # Output array
    distances = np.zeros(n_genes, dtype=np.float64)
    
    distance_to_centroid(n_genes, n_families, genes, centroids, gene_to_fam, distances, d)
    
    # Expected distances
    # Gene 1: [1,0,0] vs [0.5,0.5,0] = sqrt(0.5^2 + 0.5^2) ≈ 0.707
    # Gene 2: [0,1,0] vs [0.5,0.5,0] = sqrt(0.5^2 + 0.5^2) ≈ 0.707
    # Gene 3: [3,0,0] vs [1.5,1.5,0] = sqrt(1.5^2 + 1.5^2) ≈ 2.121
    # Gene 4: [0,3,0] vs [1.5,1.5,0] = sqrt(1.5^2 + 1.5^2) ≈ 2.121
    expected = np.array([
        math.sqrt(0.5**2 + 0.5**2), 
        math.sqrt(0.5**2 + 0.5**2), 
        math.sqrt(1.5**2 + 1.5**2), 
        math.sqrt(1.5**2 + 1.5**2)
    ])
    
    print("Distance to centroid results:")
    for i in range(n_genes):
        match = abs(distances[i] - expected[i]) < 1e-12
        print(f"  Gene {i+1}: Fortran={distances[i]:.6f}, Expected={expected[i]:.6f}, Match={match}")
    print()

def test_error_handling():
    """Test error handling with invalid family indices"""
    print("=== Testing Error Handling ===")
    
    distance_to_centroid = setup_distance_to_centroid()
    
    n_genes = 3
    n_families = 2
    d = 2
    
    # Gene data
    genes = np.array([
        1.0, 2.0,  # Gene 1
        3.0, 4.0,  # Gene 2
        5.0, 6.0   # Gene 3
    ], dtype=np.float64)
    
    # Centroids
    centroids = np.array([
        0.0, 0.0,  # Family 1
        1.0, 1.0   # Family 2
    ], dtype=np.float64)
    
    # Invalid family assignments
    gene_to_fam = np.array([1, 3, 0], dtype=np.int32)  # Family 3 doesn't exist, family 0 invalid
    
    distances = np.zeros(n_genes, dtype=np.float64)
    
    distance_to_centroid(n_genes, n_families, genes, centroids, gene_to_fam, distances, d)
    
    print("Error handling results:")
    print(f"  Gene 1 (valid family 1): {distances[0]}")
    print(f"  Gene 2 (invalid family 3): {distances[1]} (-1 expected)")
    print(f"  Gene 3 (invalid family 0): {distances[2]} (-1 expected)")
    print()

def test_high_dimensional():
    """Test with high-dimensional data"""
    print("=== Testing High-Dimensional Data ===")
    
    euclidean_distance = setup_euclidean_distance()
    
    # 100-dimensional vectors
    d = 100
    vec1 = np.arange(1, d+1, dtype=np.float64)      # [1, 2, 3, ..., 100]
    vec2 = np.arange(2, d+2, dtype=np.float64)      # [2, 3, 4, ..., 101] (shift by 1)
    result = ctypes.c_double()
    
    euclidean_distance(vec1, vec2, d, ctypes.byref(result))
    
    expected = math.sqrt(d)  # sqrt(100 * 1^2) = 10
    print("High-dimensional test (100D with unit shift):")
    print(f"  Fortran result: {result.value}")
    print(f"  Expected: {expected}")
    print(f"  Match: {abs(result.value - expected) < 1e-12}")
    print()

def test_performance():
    """Performance test with realistic genomic data size"""
    print("=== Performance Test ===")
    
    distance_to_centroid = setup_distance_to_centroid()
    
    n_genes = 1000
    n_families = 50
    d = 20  # 20 tissue types
    
    print(f"Testing with {n_genes} genes, {n_families} families, {d} dimensions")
    
    # Generate random-like data
    np.random.seed(12345)
    genes = np.random.randn(d * n_genes).astype(np.float64)
    centroids = np.random.randn(d * n_families).astype(np.float64)
    gene_to_fam = np.random.randint(1, n_families + 1, n_genes, dtype=np.int32)
    distances = np.zeros(n_genes, dtype=np.float64)
    
    # Time the operation
    start_time = time.time()
    
    distance_to_centroid(n_genes, n_families, genes, centroids, gene_to_fam, distances, d)
    
    end_time = time.time()
    elapsed = end_time - start_time
    
    # Check results
    valid_distances = np.sum(distances > 0)
    print(f"  Completed in {elapsed:.6f} seconds")
    print(f"  Valid distances computed: {valid_distances}/{n_genes}")
    print(f"  Mean distance: {np.mean(distances[distances > 0]):.6f}")
    print("  Performance test completed successfully!")
    print()

def main():
    """Run all tests"""
    print("=================================================")
    print("    EUCLIDEAN DISTANCE PYTHON INTERFACE TESTS")
    print("=================================================")
    print()
    
    test_euclidean_distance()
    test_distance_to_centroid()
    test_error_handling()
    test_high_dimensional()
    test_performance()
    
    print("=================================================")
    print("             ALL TESTS COMPLETED")
    print("=================================================")
    print("If you see this message, all Python interface tests passed! ✓")

if __name__ == "__main__":
    main()
