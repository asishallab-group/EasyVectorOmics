#!/usr/bin/env python3
"""
Comprehensive Python test suite for tissue versatility (mirrors Fortran and R unit tests)
Uses the C wrapper for the Fortran routine via ctypes
"""

import numpy as np
import ctypes
import math
import os

# Load library
dll_path = os.path.abspath("build/libtensor-omics.so")
ctypes.CDLL("libgomp.so.1", mode=ctypes.RTLD_GLOBAL)
lib = ctypes.CDLL(dll_path)

import numpy as np
import ctypes

def setup_tissue_versatility():
    """Setup tissue versatility C wrapper using ndpointer for safety"""
    tv = lib.compute_tissue_versatility_c
    tv.argtypes = [
        ctypes.c_int,  # n_axes
        ctypes.c_int,  # n_vectors

        # Use ndpointer for arrays
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"),  # expression_vectors (Fortran order)
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # exp_vecs_selection_index
        ctypes.c_int,  # n_selected_vectors
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # axes_selection
        ctypes.c_int,  # n_selected_axes
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # tissue_versatilities
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # tissue_angles_deg
    ]
    tv.restype = None
    return tv

def tv_call(expr, select_vec, select_axes):
    """Call the C wrapper for tissue versatility, return dict of results"""
    n_axes = expr.shape[0]
    n_vectors = expr.shape[1]

    # Ensure arrays have correct dtype and memory layout
    expr_f = np.asfortranarray(expr, dtype=np.float64)
    select_vec = np.ascontiguousarray(select_vec, dtype=np.int32)
    select_axes = np.ascontiguousarray(select_axes, dtype=np.int32)

    n_selected_vectors = int(np.sum(select_vec))
    n_selected_axes = int(np.sum(select_axes))
    tissue_versatilities = np.empty(n_selected_vectors, dtype=np.float64)
    tissue_angles_deg = np.empty(n_selected_vectors, dtype=np.float64)

    # Call the C/Fortran wrapper
    tv = setup_tissue_versatility()
    tv(
        n_axes,
        n_vectors,
        expr_f,
        select_vec,
        n_selected_vectors,
        select_axes,
        n_selected_axes,
        tissue_versatilities,
        tissue_angles_deg
    )

    return {
        'tissue_versatilities': tissue_versatilities,
        'tissue_angles_deg': tissue_angles_deg
    }

# 1. Uniform expression (should yield TV=0)
def test_uniform_expression():
    expr = np.full((3, 1), 2.0)
    res = tv_call(expr, [1], [1, 1, 1])
    assert abs(res['tissue_versatilities'][0]) < 1e-12
    assert abs(res['tissue_angles_deg'][0]) < 1e-12
    print("test_uniform_expression passed")


# 2. Single axis expression (should yield TV=1)
def test_single_axis_expression():
    expr = np.array([[0],[0],[5]], dtype=np.float64)
    res = tv_call(expr, [1], [1,1,1])
    assert abs(res['tissue_versatilities'][0] - 1) < 1e-12
    assert res['tissue_angles_deg'][0] > 0
    print("test_single_axis_expression passed")

# 3. Null vector (should yield TV=1, angle=90)
def test_null_vector():
    expr = np.zeros((3,1), dtype=np.float64)
    res = tv_call(expr, [1], [1,1,1])
    assert abs(res['tissue_versatilities'][0] - 1) < 1e-12
    assert abs(res['tissue_angles_deg'][0] - 90) < 1e-12
    print("test_null_vector passed")

# 4. Partial axis selection (subspace)
def test_partial_axis_selection():
    expr = np.array([[1],[2],[3]], dtype=np.float64)
    res = tv_call(expr, [1], [1,0,1])
    assert 0 <= res['tissue_versatilities'][0] <= 1
    assert 0 <= res['tissue_angles_deg'][0] <= 90
    print("test_partial_axis_selection passed")

# 5. Mixed vectors (uniform, single axis, null)
def test_mixed_vectors():
    expr = np.array([[1,0,0],[1,0,0],[1,2,0]], dtype=np.float64)
    res = tv_call(expr, [1,1,1], [1,1,1])
    assert abs(res['tissue_versatilities'][0]) < 1e-12
    assert abs(res['tissue_versatilities'][1] - 1) < 1e-12
    assert abs(res['tissue_versatilities'][2] - 1) < 1e-12
    assert abs(res['tissue_angles_deg'][0]) < 1e-12
    assert res['tissue_angles_deg'][1] > 0
    assert abs(res['tissue_angles_deg'][2] - 90) < 1e-12
    print("test_mixed_vectors passed")

# 6. Angle output in degrees for a known case (should be 45)
def test_angle_degrees():
    expr = np.array([[1],[0]], dtype=np.float64)
    res = tv_call(expr, [1], [1,1])
    assert abs(res['tissue_angles_deg'][0] - 45) < 1e-12
    print("test_angle_degrees passed")

# 7. Multiple vectors selection
def test_multiple_vectors_selection():
    expr = np.array([[1,0,0],[1,2,0]], dtype=np.float64)
    res = tv_call(expr, [1,0,1], [1,1])
    assert abs(res['tissue_versatilities'][0]) < 1e-12
    assert abs(res['tissue_versatilities'][1] - 1) < 1e-12
    assert abs(res['tissue_angles_deg'][0]) < 1e-12
    assert abs(res['tissue_angles_deg'][1] - 90) < 1e-12
    print("test_multiple_vectors_selection passed")

# 8. High-dimensional vectors (4D, 5D)
def test_high_dimensional_vectors():
    expr4 = np.full((4,1), 1.0)
    expr5 = np.full((5,1), 2.0)
    res4 = tv_call(expr4, [1], [1,1,1,1])
    res5 = tv_call(expr5, [1], [1,1,1,1,1])
    assert abs(res4['tissue_versatilities'][0]) < 1e-12
    assert abs(res4['tissue_angles_deg'][0]) < 1e-12
    assert abs(res5['tissue_versatilities'][0]) < 1e-12
    assert abs(res5['tissue_angles_deg'][0]) < 1e-12
    print("test_high_dimensional_vectors passed")

# 9. Randomized vectors and axes
def test_randomized_vectors_axes():
    np.random.seed(42)
    n_axes = 5
    n_vecs = 4
    expr = np.random.rand(n_axes, n_vecs)
    res = tv_call(expr, [1]*n_vecs, [1,0,1,0,1])
    assert np.all((res['tissue_versatilities'] >= 0) & (res['tissue_versatilities'] <= 1))
    assert np.all((res['tissue_angles_deg'] >= 0) & (res['tissue_angles_deg'] <= 90))
    print("test_randomized_vectors_axes passed")

# 10. Numerical stability (very large/small values)
def test_numerical_stability():
    expr = np.array([[1e-15,1e15],[1e-15,1e15],[1e-15,1e15]], dtype=np.float64)
    res = tv_call(expr, [1,1], [1,1,1])
    assert abs(res['tissue_versatilities'][0]) < 1e-12
    assert abs(res['tissue_angles_deg'][0]) < 1e-12
    assert abs(res['tissue_versatilities'][1]) < 1e-12
    assert abs(res['tissue_angles_deg'][1]) < 1e-12
    print("test_numerical_stability passed")

# 11. Invalid input: no axes selected (should return -1)
def test_invalid_input_no_axes():
    expr = np.array([[1],[2],[3]], dtype=np.float64)
    res = tv_call(expr, [1], [0,0,0])
    assert abs(res['tissue_versatilities'][0] + 1) < 1e-12
    print("test_invalid_input_no_axes passed")

# 12. Multiple selection, partial axes
def test_multiple_selection_partial_axes():
    expr = np.array([[1,3,5],[2,4,6]], dtype=np.float64)
    res = tv_call(expr, [1,0,1], [1,0])
    assert len(res['tissue_versatilities']) == 2
    assert np.all((res['tissue_versatilities'] >= 0) & (res['tissue_versatilities'] <= 1))
    print("test_multiple_selection_partial_axes passed")


def main():
    print("=================================================")
    print("    TISSUE VERSATILITY FULL PYTHON INTERFACE TESTS")
    print("=================================================\n")
    test_uniform_expression()
    test_single_axis_expression()
    test_null_vector()
    test_partial_axis_selection()
    test_mixed_vectors()
    test_angle_degrees()
    test_multiple_vectors_selection()
    test_high_dimensional_vectors()
    test_randomized_vectors_axes()
    test_numerical_stability()
    test_invalid_input_no_axes()
    test_multiple_selection_partial_axes()
    print("=================================================")
    print("             ALL TESTS COMPLETED")
    print("=================================================")
    print("If you see this message, all tissue versatility Python interface tests passed! ✓")

if __name__ == "__main__":
    main()
