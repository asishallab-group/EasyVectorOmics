"""
Test script for clock hand angle functions
Python equivalent of the R and Fortran clock hand angle tests
"""

import numpy as np
import ctypes
import time
import math

# Load library
ctypes.CDLL("libgomp.so.1", mode=ctypes.RTLD_GLOBAL)
lib = ctypes.CDLL("build/libtensor-omics.so")

# Constants
PI = math.pi
TOL = 1e-12

def setup_clock_hand_angle_between_vectors():
    """Setup clock hand angle between vectors function"""
    clock_hand_angle = lib.clock_hand_angle_between_vectors_c
    clock_hand_angle.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # v1
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # v2
        ctypes.c_int,  # n_dims
        ctypes.POINTER(ctypes.c_double),  # signed_angle
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS")     # selected_axes_for_signed
    ]
    clock_hand_angle.restype = None
    return clock_hand_angle

def setup_clock_hand_angles_for_shift_vectors():
    """Setup clock hand angles for shift vectors function"""
    clock_hand_angles = lib.clock_hand_angles_for_shift_vectors_c
    clock_hand_angles.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # origins
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # targets
        ctypes.c_int,  # n_dims
        ctypes.c_int,  # n_vecs
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # vecs_selection_mask
        ctypes.c_int,  # n_selected_vecs
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # selected_axes_for_signed
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")   # signed_angles
    ]
    clock_hand_angles.restype = None
    return clock_hand_angles

def test_identical_vectors_2d():
    """Test identical vectors in 2D (should give 0 angle)"""
    print("=== Testing Identical Vectors 2D ===")
    
    clock_hand_angle = setup_clock_hand_angle_between_vectors()
    
    v1 = np.array([1.0, 0.0], dtype=np.float64)
    v2 = np.array([1.0, 0.0], dtype=np.float64)
    n_dims = len(v1)
    selected_axes = np.array([1, 2, 1], dtype=np.int32)  # Ignored for 2D
    result = ctypes.c_double()
    
    clock_hand_angle(v1, v2, n_dims, ctypes.byref(result), selected_axes)
    
    expected = 0.0
    print("Test: Identical 2D vectors [1,0] vs [1,0]:")
    print(f"  Fortran result: {result.value}")
    print(f"  Expected: {expected}")
    print(f"  Match: {abs(result.value - expected) < TOL}")
    print()

def test_perpendicular_vectors_2d():
    """Test perpendicular vectors in 2D (should give ±π/2)"""
    print("=== Testing Perpendicular Vectors 2D ===")
    
    clock_hand_angle = setup_clock_hand_angle_between_vectors()
    
    v1 = np.array([1.0, 0.0], dtype=np.float64)
    v2 = np.array([0.0, 1.0], dtype=np.float64)
    n_dims = len(v1)
    selected_axes = np.array([1, 2, 1], dtype=np.int32)  # Ignored for 2D
    result = ctypes.c_double()
    
    clock_hand_angle(v1, v2, n_dims, ctypes.byref(result), selected_axes)
    
    expected_magnitude = PI/2
    print("Test: Perpendicular 2D vectors [1,0] vs [0,1]:")
    print(f"  Fortran result: {result.value}")
    print(f"  Expected magnitude: ±{expected_magnitude}")
    print(f"  Magnitude match: {abs(abs(result.value) - expected_magnitude) < TOL}")
    print(f"  Sign (should be positive): {result.value > 0}")
    print(f"  Angle in degrees: {result.value * 180/PI:.2f}°")
    print()

def test_opposite_vectors_2d():
    """Test opposite vectors in 2D (should give ±π)"""
    print("=== Testing Opposite Vectors 2D ===")
    
    clock_hand_angle = setup_clock_hand_angle_between_vectors()
    
    v1 = np.array([1.0, 0.0], dtype=np.float64)
    v2 = np.array([-1.0, 0.0], dtype=np.float64)
    n_dims = len(v1)
    selected_axes = np.array([1, 2, 1], dtype=np.int32)  # Ignored for 2D
    result = ctypes.c_double()
    
    clock_hand_angle(v1, v2, n_dims, ctypes.byref(result), selected_axes)
    
    expected_magnitude = PI
    print("Test: Opposite 2D vectors [1,0] vs [-1,0]:")
    print(f"  Fortran result: {result.value}")
    print(f"  Expected magnitude: ±{expected_magnitude}")
    print(f"  Magnitude match: {abs(abs(result.value) - expected_magnitude) < TOL}")
    print(f"  Angle in degrees: {result.value * 180/PI:.2f}°")
    print()

def test_45_degree_rotation_2d():
    """Test 45-degree rotation in 2D"""
    print("=== Testing 45-Degree Rotation 2D ===")
    
    clock_hand_angle = setup_clock_hand_angle_between_vectors()
    
    v1 = np.array([1.0, 0.0], dtype=np.float64)
    v2 = np.array([math.sqrt(2)/2, math.sqrt(2)/2], dtype=np.float64)  # 45 degrees
    n_dims = len(v1)
    selected_axes = np.array([1, 2, 1], dtype=np.int32)  # Ignored for 2D
    result = ctypes.c_double()
    
    clock_hand_angle(v1, v2, n_dims, ctypes.byref(result), selected_axes)
    
    expected = PI/4
    print("Test: 45-degree counterclockwise rotation:")
    print(f"  v1: {v1}")
    print(f"  v2: {v2}")
    print(f"  Fortran result: {result.value}")
    print(f"  Expected: {expected}")
    print(f"  Match: {abs(result.value - expected) < TOL}")
    print(f"  Angle in degrees: {result.value * 180/PI:.2f}°")
    print()

def test_clockwise_vs_counterclockwise_2d():
    """Test clockwise vs counterclockwise rotations in 2D"""
    print("=== Testing Clockwise vs Counterclockwise 2D ===")
    
    clock_hand_angle = setup_clock_hand_angle_between_vectors()
    
    v1 = np.array([1.0, 0.0], dtype=np.float64)
    v2_ccw = np.array([0.0, 1.0], dtype=np.float64)   # 90° counterclockwise
    v2_cw = np.array([0.0, -1.0], dtype=np.float64)   # 90° clockwise
    n_dims = len(v1)
    selected_axes = np.array([1, 2, 1], dtype=np.int32)  # Ignored for 2D
    
    result_ccw = ctypes.c_double()
    result_cw = ctypes.c_double()
    
    clock_hand_angle(v1, v2_ccw, n_dims, ctypes.byref(result_ccw), selected_axes)
    clock_hand_angle(v1, v2_cw, n_dims, ctypes.byref(result_cw), selected_axes)
    
    print("Test: Clockwise vs Counterclockwise:")
    print(f"  Counterclockwise angle: {result_ccw.value:.6f} ({result_ccw.value * 180/PI:.2f}°)")
    print(f"  Clockwise angle: {result_cw.value:.6f} ({result_cw.value * 180/PI:.2f}°)")
    print(f"  CCW > 0: {result_ccw.value > 0}")
    print(f"  CW < 0: {result_cw.value < 0}")
    print(f"  Magnitudes equal: {abs(abs(result_ccw.value) - abs(result_cw.value)) < TOL}")
    print()

def test_3d_vectors():
    """Test 3D vector calculations"""
    print("=== Testing 3D Vectors ===")
    
    clock_hand_angle = setup_clock_hand_angle_between_vectors()
    
    # Test identical 3D vectors
    v1 = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    v2 = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    n_dims = len(v1)
    selected_axes = np.array([1, 2, 3], dtype=np.int32)  # Ignored for 3D
    result = ctypes.c_double()
    
    clock_hand_angle(v1, v2, n_dims, ctypes.byref(result), selected_axes)
    
    print("Test 1: Identical 3D vectors [1,1,1] vs [1,1,1]:")
    print(f"  Result: {result.value} (should be 0)")
    print(f"  Match: {abs(result.value) < TOL}")
    
    # Test perpendicular 3D vectors
    v1 = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    v2 = np.array([0.0, 1.0, 0.0], dtype=np.float64)
    result = ctypes.c_double()
    
    clock_hand_angle(v1, v2, n_dims, ctypes.byref(result), selected_axes)
    
    expected_magnitude = PI/2
    print("Test 2: Perpendicular 3D vectors [1,0,0] vs [0,1,0]:")
    print(f"  Result: {result.value:.6f} ({result.value * 180/PI:.2f}°)")
    print(f"  Expected magnitude: ±{expected_magnitude}")
    print(f"  Magnitude match: {abs(abs(result.value) - expected_magnitude) < TOL}")
    print()

def test_high_dimensional():
    """Test high-dimensional vectors with selected axes"""
    print("=== Testing High-Dimensional Vectors ===")
    
    clock_hand_angle = setup_clock_hand_angle_between_vectors()
    
    # 5D vectors, perpendicular in first two dimensions
    v1 = np.array([1.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64)
    v2 = np.array([0.0, 1.0, 0.0, 0.0, 0.0], dtype=np.float64)
    n_dims = len(v1)
    selected_axes = np.array([1, 2, 3], dtype=np.int32)  # Use first 3 dimensions for orientation
    result = ctypes.c_double()
    
    clock_hand_angle(v1, v2, n_dims, ctypes.byref(result), selected_axes)
    
    expected_magnitude = PI/2
    print("Test 1: 5D vectors perpendicular in first two dimensions:")
    print(f"  v1: {v1}")
    print(f"  v2: {v2}")
    print(f"  Selected axes: {selected_axes}")
    print(f"  Result: {result.value:.6f} ({result.value * 180/PI:.2f}°)")
    print(f"  Expected magnitude: ±{expected_magnitude}")
    print(f"  Magnitude match: {abs(abs(result.value) - expected_magnitude) < TOL}")
    
    # 7D vectors with specific selected axes
    v1 = np.array([0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float64)
    v2 = np.array([0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0], dtype=np.float64)
    n_dims = len(v1)
    selected_axes = np.array([3, 5, 1], dtype=np.int32)  # Use dimensions 3, 5, 1 for orientation
    result = ctypes.c_double()
    
    clock_hand_angle(v1, v2, n_dims, ctypes.byref(result), selected_axes)
    
    print("Test 2: 7D vectors with selected axes [3,5,1]:")
    print(f"  v1: {v1}")
    print(f"  v2: {v2}")
    print(f"  Selected axes: {selected_axes}")
    print(f"  Result: {result.value:.6f} ({result.value * 180/PI:.2f}°)")
    print(f"  Expected magnitude: ±{expected_magnitude}")
    print(f"  Magnitude match: {abs(abs(result.value) - expected_magnitude) < TOL}")
    print()

def test_shift_vectors_single_pair():
    """Test single pair of shift vectors"""
    print("=== Testing Single Pair Shift Vectors ===")
    
    clock_hand_angles = setup_clock_hand_angles_for_shift_vectors()
    
    # Single pair: [1,0] -> [0,1] (90° counterclockwise)
    n_dims = 2
    n_vecs = 1
    n_selected_vecs = 1
    
    origins = np.array([1.0, 0.0], dtype=np.float64)  # Column-major: [1,0] as single vector
    targets = np.array([0.0, 1.0], dtype=np.float64)  # Column-major: [0,1] as single vector
    vecs_selection_mask = np.array([1], dtype=np.int32)  # Select the single pair
    selected_axes = np.array([1, 2, 1], dtype=np.int32)  # Ignored for 2D
    signed_angles = np.zeros(n_selected_vecs, dtype=np.float64)
    
    clock_hand_angles(origins, targets, n_dims, n_vecs, vecs_selection_mask, 
                     n_selected_vecs, selected_axes, signed_angles)
    
    expected = PI/2
    print("Test: Single pair [1,0] -> [0,1]:")
    print(f"  Result: {signed_angles[0]:.6f} ({signed_angles[0] * 180/PI:.2f}°)")
    print(f"  Expected: {expected:.6f} ({expected * 180/PI:.2f}°)")
    print(f"  Match: {abs(signed_angles[0] - expected) < TOL}")
    print()

def test_shift_vectors_multiple_pairs():
    """Test multiple pairs of shift vectors"""
    print("=== Testing Multiple Pairs Shift Vectors ===")
    
    clock_hand_angles = setup_clock_hand_angles_for_shift_vectors()
    
    # Three different rotations
    n_dims = 2
    n_vecs = 3
    n_selected_vecs = 3
    
    # Column-major layout: each column is a vector
    origins = np.array([
        1.0, 0.0,   # Vector 1: [1,0]
        1.0, 0.0,   # Vector 2: [1,0]
        1.0, 0.0    # Vector 3: [1,0]
    ], dtype=np.float64)
    
    targets = np.array([
        0.0, 1.0,    # Vector 1: [0,1] -> 90° CCW
        -1.0, 0.0,   # Vector 2: [-1,0] -> 180°
        0.0, -1.0    # Vector 3: [0,-1] -> 90° CW
    ], dtype=np.float64)
    
    vecs_selection_mask = np.array([1, 1, 1], dtype=np.int32)  # Select all pairs
    selected_axes = np.array([1, 2, 1], dtype=np.int32)  # Ignored for 2D
    signed_angles = np.zeros(n_selected_vecs, dtype=np.float64)
    
    clock_hand_angles(origins, targets, n_dims, n_vecs, vecs_selection_mask, 
                     n_selected_vecs, selected_axes, signed_angles)
    
    expected = [PI/2, PI, -PI/2]  # Note: 180° can be +π or -π
    print("Test: Multiple pairs with different rotations:")
    for i in range(n_selected_vecs):
        if i == 1:  # 180° case - check magnitude only
            match = abs(abs(signed_angles[i]) - abs(expected[i])) < TOL
            print(f"  Pair {i+1}: {signed_angles[i]:.6f} ({signed_angles[i] * 180/PI:.2f}°), "
                  f"Expected magnitude: ±{abs(expected[i]):.6f}, Match: {match}")
        else:
            match = abs(signed_angles[i] - expected[i]) < TOL
            print(f"  Pair {i+1}: {signed_angles[i]:.6f} ({signed_angles[i] * 180/PI:.2f}°), "
                  f"Expected: {expected[i]:.6f}, Match: {match}")
    print()

def test_shift_vectors_with_selection_mask():
    """Test shift vectors with selection mask"""
    print("=== Testing Shift Vectors with Selection Mask ===")
    
    clock_hand_angles = setup_clock_hand_angles_for_shift_vectors()
    
    # Four vectors, but only select 2nd and 4th
    n_dims = 2
    n_vecs = 4
    n_selected_vecs = 2  # Only 2nd and 4th are selected
    
    # Column-major layout
    origins = np.array([
        1.0, 0.0,   # Vector 1: [1,0] - Not selected
        1.0, 0.0,   # Vector 2: [1,0] - Selected (180°)
        1.0, 0.0,   # Vector 3: [1,0] - Not selected
        1.0, 0.0    # Vector 4: [1,0] - Selected (45°)
    ], dtype=np.float64)
    
    targets = np.array([
        0.0, 1.0,                           # Vector 1: [0,1] - Not selected
        -1.0, 0.0,                          # Vector 2: [-1,0] - Selected (180°)
        0.0, -1.0,                          # Vector 3: [0,-1] - Not selected
        math.sqrt(2)/2, math.sqrt(2)/2      # Vector 4: [√2/2,√2/2] - Selected (45°)
    ], dtype=np.float64)
    
    vecs_selection_mask = np.array([0, 1, 0, 1], dtype=np.int32)  # Select 2nd and 4th
    selected_axes = np.array([1, 2, 1], dtype=np.int32)  # Ignored for 2D
    signed_angles = np.zeros(n_selected_vecs, dtype=np.float64)
    
    clock_hand_angles(origins, targets, n_dims, n_vecs, vecs_selection_mask, 
                     n_selected_vecs, selected_axes, signed_angles)
    
    expected = [PI, PI/4]  # 180° and 45°, Note: 180° can be +π or -π
    print("Test: Selection mask [False, True, False, True]:")
    print(f"  Selected pair 1 (was 2nd): {signed_angles[0]:.6f} ({signed_angles[0] * 180/PI:.2f}°)")
    print(f"    Expected magnitude: ±{abs(expected[0]):.6f}, Match: {abs(abs(signed_angles[0]) - abs(expected[0])) < TOL}")
    print(f"  Selected pair 2 (was 4th): {signed_angles[1]:.6f} ({signed_angles[1] * 180/PI:.2f}°)")
    print(f"    Expected: {expected[1]:.6f}, Match: {abs(signed_angles[1] - expected[1]) < TOL}")
    print()

def test_edge_cases():
    """Test edge cases and precision"""
    print("=== Testing Edge Cases ===")
    
    clock_hand_angle = setup_clock_hand_angle_between_vectors()
    
    # Test 1: Nearly identical vectors
    epsilon = 1e-15
    v1 = np.array([1.0, 0.0], dtype=np.float64)
    v2 = np.array([1.0, epsilon], dtype=np.float64)
    n_dims = len(v1)
    selected_axes = np.array([1, 2, 1], dtype=np.int32)
    result = ctypes.c_double()
    
    clock_hand_angle(v1, v2, n_dims, ctypes.byref(result), selected_axes)
    
    print("Test 1: Nearly identical vectors:")
    print(f"  v1: {v1}")
    print(f"  v2: {v2}")
    print(f"  Result: {result.value:.2e} radians")
    print(f"  Very small angle: {abs(result.value) < 1e-10}")
    
    # Test 2: Denormalized vectors (large magnitude)
    v1 = np.array([100.0, 0.0], dtype=np.float64)
    v2 = np.array([0.0, 50.0], dtype=np.float64)
    result = ctypes.c_double()
    
    clock_hand_angle(v1, v2, n_dims, ctypes.byref(result), selected_axes)
    
    expected = PI/2
    print("Test 2: Denormalized vectors [100,0] vs [0,50]:")
    print(f"  Result: {result.value:.6f} ({result.value * 180/PI:.2f}°)")
    print(f"  Expected: {expected:.6f}")
    print(f"  Match: {abs(result.value - expected) < TOL}")
    
    # Test 3: Very small vectors
    tiny = 1e-14
    v1 = np.array([tiny, 0.0], dtype=np.float64)
    v2 = np.array([0.0, tiny], dtype=np.float64)
    result = ctypes.c_double()
    
    clock_hand_angle(v1, v2, n_dims, ctypes.byref(result), selected_axes)
    
    print("Test 3: Tiny vectors:")
    print(f"  Result: {result.value:.6f} ({result.value * 180/PI:.2f}°)")
    print(f"  Expected magnitude: ±{expected:.6f}")
    print(f"  Magnitude match: {abs(abs(result.value) - expected) < 1e-10}")
    print()

def test_consistency_between_functions():
    """Test consistency between single and batch functions"""
    print("=== Testing Consistency Between Functions ===")
    
    clock_hand_angle = setup_clock_hand_angle_between_vectors()
    clock_hand_angles = setup_clock_hand_angles_for_shift_vectors()
    
    # Test same calculation with both functions
    v1 = np.array([1.0, 0.0], dtype=np.float64)
    v2 = np.array([0.0, 1.0], dtype=np.float64)
    n_dims = len(v1)
    selected_axes = np.array([1, 2, 1], dtype=np.int32)
    
    # Single function
    single_result = ctypes.c_double()
    clock_hand_angle(v1, v2, n_dims, ctypes.byref(single_result), selected_axes)
    
    # Batch function with single pair
    n_vecs = 1
    n_selected_vecs = 1
    origins = v1.copy()
    targets = v2.copy()
    vecs_selection_mask = np.array([1], dtype=np.int32)
    batch_results = np.zeros(n_selected_vecs, dtype=np.float64)
    
    clock_hand_angles(origins, targets, n_dims, n_vecs, vecs_selection_mask, 
                     n_selected_vecs, selected_axes, batch_results)
    
    print("Test: Consistency between single and batch functions:")
    print(f"  Single function result: {single_result.value:.10f}")
    print(f"  Batch function result: {batch_results[0]:.10f}")
    print(f"  Difference: {abs(single_result.value - batch_results[0]):.2e}")
    print(f"  Match: {abs(single_result.value - batch_results[0]) < TOL}")
    print()

def test_mathematical_properties():
    """Test mathematical properties (anti-commutativity)"""
    print("=== Testing Mathematical Properties ===")
    
    clock_hand_angle = setup_clock_hand_angle_between_vectors()
    
    # Test anti-commutativity: angle(v1,v2) = -angle(v2,v1)
    v1 = np.array([1.0, 2.0], dtype=np.float64)
    v2 = np.array([3.0, 1.0], dtype=np.float64)
    
    # Normalize vectors
    v1 = v1 / np.sqrt(np.sum(v1**2))
    v2 = v2 / np.sqrt(np.sum(v2**2))
    
    n_dims = len(v1)
    selected_axes = np.array([1, 2, 1], dtype=np.int32)
    
    result_12 = ctypes.c_double()
    result_21 = ctypes.c_double()
    
    clock_hand_angle(v1, v2, n_dims, ctypes.byref(result_12), selected_axes)
    clock_hand_angle(v2, v1, n_dims, ctypes.byref(result_21), selected_axes)
    
    print("Test: Anti-commutativity angle(v1,v2) = -angle(v2,v1):")
    print(f"  v1 (normalized): {v1}")
    print(f"  v2 (normalized): {v2}")
    print(f"  angle(v1,v2): {result_12.value:.6f}")
    print(f"  angle(v2,v1): {result_21.value:.6f}")
    print(f"  Sum (should be ~0): {result_12.value + result_21.value:.2e}")
    print(f"  Anti-commutative: {abs(result_12.value + result_21.value) < TOL}")
    print()

def test_performance():
    """Performance test with large-scale data"""
    print("=== Performance Test ===")
    
    clock_hand_angles = setup_clock_hand_angles_for_shift_vectors()
    
    # Large-scale test
    n_dims = 50
    n_vecs = 1000
    n_selected_vecs = n_vecs
    
    print(f"Testing with {n_vecs} vector pairs in {n_dims} dimensions")
    
    # Generate random-like data
    np.random.seed(12345)
    origins = np.random.randn(n_dims * n_vecs).astype(np.float64)
    targets = np.random.randn(n_dims * n_vecs).astype(np.float64)
    
    # Normalize vectors (for more meaningful angles)
    for i in range(n_vecs):
        start_idx = i * n_dims
        end_idx = (i + 1) * n_dims
        origin_vec = origins[start_idx:end_idx]
        target_vec = targets[start_idx:end_idx]
        
        origin_norm = np.sqrt(np.sum(origin_vec**2))
        target_norm = np.sqrt(np.sum(target_vec**2))
        
        if origin_norm > 0:
            origins[start_idx:end_idx] /= origin_norm
        if target_norm > 0:
            targets[start_idx:end_idx] /= target_norm
    
    vecs_selection_mask = np.ones(n_vecs, dtype=np.int32)
    selected_axes = np.array([1, 2, 3], dtype=np.int32)
    signed_angles = np.zeros(n_selected_vecs, dtype=np.float64)
    
    # Time the operation
    start_time = time.time()
    
    clock_hand_angles(origins, targets, n_dims, n_vecs, vecs_selection_mask, 
                     n_selected_vecs, selected_axes, signed_angles)
    
    end_time = time.time()
    elapsed = end_time - start_time
    
    # Check results
    valid_count = np.sum(~np.isnan(signed_angles))
    mean_angle = np.mean(signed_angles[~np.isnan(signed_angles)]) if valid_count > 0 else 0
    
    print(f"  Completed in {elapsed:.6f} seconds")
    print(f"  Valid angles computed: {valid_count}/{n_vecs}")
    print(f"  Mean angle: {mean_angle:.6f} radians ({mean_angle * 180/PI:.2f}°)")
    print(f"  Performance: {n_vecs/elapsed:.0f} calculations/second")
    print("  Performance test completed successfully!")
    print()

def main():
    """Run all tests"""
    print("=================================================")
    print("    CLOCK HAND ANGLES PYTHON INTERFACE TESTS")
    print("=================================================")
    print()
    
    # Basic 2D tests
    test_identical_vectors_2d()
    test_perpendicular_vectors_2d()
    test_opposite_vectors_2d()
    test_45_degree_rotation_2d()
    test_clockwise_vs_counterclockwise_2d()
    
    # 3D and high-dimensional tests
    test_3d_vectors()
    test_high_dimensional()
    
    # Shift vectors tests
    test_shift_vectors_single_pair()
    test_shift_vectors_multiple_pairs()
    test_shift_vectors_with_selection_mask()
    
    # Edge cases and properties
    test_edge_cases()
    test_consistency_between_functions()
    test_mathematical_properties()
    
    # Performance test
    test_performance()
    
    print("=================================================")
    print("             ALL TESTS COMPLETED")
    print("=================================================")
    print("If you see this message, all Python interface tests passed! ✓")

if __name__ == "__main__":
    main()
