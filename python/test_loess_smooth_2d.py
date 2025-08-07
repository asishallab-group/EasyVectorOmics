"""
Unit tests for LOESS smoothing via C wrapper (Python)
Replicates the edge cases from R/Fortran tests.
"""

import numpy as np
import ctypes
import os
from ctypes import c_double

# Load library (adjust path if needed)
ctypes.CDLL("libgomp.so.1", mode=ctypes.RTLD_GLOBAL)
lib = ctypes.CDLL("build/libtensor-omics.so")

def setup_loess_smooth_2d():
    f = lib.loess_smooth_2d_c
    f.argtypes = [
        ctypes.c_int, 
        ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.int32, flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),
        ctypes.c_double, 
        ctypes.c_double,
        np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.int32, flags='C_CONTIGUOUS'),
    ]
    f.restype = None
    return f



def main():
    f = setup_loess_smooth_2d()
    # === Python tests equivalentes a los de R ===
    print("\n[loess_smooth_2d] Case 1: Simple smoothing, no mask (with Fortran 1-based padding)")
    n_total = 5
    n_target = 2
    x_ref = np.array([1,2,3,4,5], dtype=np.float64)
    y_ref = np.array([10,20,30,40,50], dtype=np.float64)
    indices_used = np.arange(1, 6, dtype=np.int32) 
    x_query = np.array([2.5, 4.5], dtype=np.float64)
    kernel_sigma = np.float64(1.0)
    kernel_cutoff = np.float64(2.0)


    y_out = np.zeros(n_target, dtype=np.float64)
    workspace_weights = np.zeros(n_total, dtype=np.float64)
    workspace_values = np.zeros(n_total, dtype=np.float64)
    mask_in = np.ones(n_total, dtype=np.int32)
    f(
        n_total, n_target,
        x_ref, y_ref, indices_used, x_query,
        ctypes.c_double(kernel_sigma), ctypes.c_double(kernel_cutoff),
        y_out, workspace_weights, workspace_values, mask_in
    )
    print("y_out:", y_out)

    print("\n[loess_smooth_2d] Case 2: With mask (exclude some points)")
    mask_in = np.array([1,0,1,0,1], dtype=np.int32)
    y_out = np.zeros(n_target, dtype=np.float64)
    f(n_total, n_target, x_ref, y_ref, indices_used, x_query, kernel_sigma, kernel_cutoff, y_out, workspace_weights, workspace_values, mask_in)
    print("y_out:", y_out)

    print("\n[loess_smooth_2d] Case 3: All mask FALSE (should fallback to y_ref)")
    mask_in = np.zeros(n_total, dtype=np.int32)
    y_out = np.zeros(n_target, dtype=np.float64)
    f(n_total, n_target, x_ref, y_ref, indices_used, x_query, kernel_sigma, kernel_cutoff, y_out, workspace_weights, workspace_values, mask_in)
    print("y_out:", y_out)

    print("\n[loess_smooth_2d] Case 4: x_query outside x_ref range")
    x_query = np.array([-10, 100], dtype=np.float64)
    mask_in = np.ones(n_total, dtype=np.int32)
    y_out = np.zeros(2, dtype=np.float64)
    f(n_total, 2, x_ref, y_ref, indices_used, x_query, kernel_sigma, kernel_cutoff, y_out, workspace_weights, workspace_values, mask_in)
    print("y_out:", y_out)

    print("\n[loess_smooth_2d] Case 5: x_query == x_ref")
    x_query = x_ref[:2].copy()
    y_out = np.zeros(2, dtype=np.float64)
    f(n_total, 2, x_ref, y_ref, indices_used, x_query, kernel_sigma, kernel_cutoff, y_out, workspace_weights, workspace_values, mask_in)
    print("y_out:", y_out)

    print("\n[loess_smooth_2d] Case 6: y_ref contains np.nan")
    y_ref_nan = y_ref.copy()
    y_ref_nan[2] = np.nan
    y_out = np.zeros(n_target, dtype=np.float64)
    try:
        f(n_total, n_target, x_ref, y_ref_nan, indices_used, np.array([2.5, 4.5], dtype=np.float64), kernel_sigma, kernel_cutoff, y_out, workspace_weights, workspace_values, mask_in)
        print("y_out:", y_out)
    except Exception as e:
        print("Exception (expected if nan not handled):", e)

    print("\n[loess_smooth_2d] Case 7: n_total = 1 (single point)")
    n_total1 = 1
    n_target1 = 3
    x_ref1 = np.array([42.0], dtype=np.float64)
    y_ref1 = np.array([99.0], dtype=np.float64)
    indices_used1 = np.array([1], dtype=np.int32)
    x_query1 = np.array([0, 42, 100], dtype=np.float64)
    y_out1 = np.zeros(n_target1, dtype=np.float64)
    workspace_weights1 = np.zeros(n_total1, dtype=np.float64)
    workspace_values1 = np.zeros(n_total1, dtype=np.float64)
    mask_in1 = np.ones(n_total1, dtype=np.int32)
    # For n_total1=1, no need to pad, but ensure correct dtype
    f(n_total1, n_target1, x_ref1, y_ref1, indices_used1, x_query1, kernel_sigma, kernel_cutoff, y_out1, workspace_weights1, workspace_values1, mask_in1)
    print("y_out:", y_out1)

    print("\n[loess_smooth_2d] Case 8: kernel_sigma = 0 (should return nearest y_ref)")
    kernel_sigma0 = 0.0
    y_out = np.zeros(n_target, dtype=np.float64)
    f(n_total, n_target, x_ref, y_ref, indices_used, np.array([2.5, 4.5], dtype=np.float64), kernel_sigma0, kernel_cutoff, y_out, workspace_weights, workspace_values, mask_in)
    print("y_out:", y_out)
    print("\nAll loess_smooth_2d tests completed.")

if __name__ == "__main__":
    main()
