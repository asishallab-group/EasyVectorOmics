#!/usr/bin/env python3
"""
Test script for outlier detection Fortran routines via C wrappers
Python equivalent of the R get_outliers.R tests
All comments and documentation in English.
"""

import numpy as np
import ctypes
import os
from ctypes import c_double, byref


# Load library
ctypes.CDLL("libgomp.so.1", mode=ctypes.RTLD_GLOBAL)
lib = ctypes.CDLL("build/libtensor-omics.so")

def np2c(arr, dtype):
    return arr.ctypes.data_as(dtype)

def print_header(msg):
    print("\n" + "="*60)
    print(msg)
    print("="*60)

def setup_compute_family_scaling():
    f = lib.compute_family_scaling_c
    f.argtypes = [
        ctypes.c_int,  # n_genes
        ctypes.c_int,  # n_families
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # distances (n_genes,)
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # gene_to_fam (n_genes,)
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # dscale (n_families,)
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # loess_x (n_families,)
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # loess_y (n_families,)
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # indices_used (n_families,)
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # perm_tmp (n_genes,)
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # stack_left_tmp (n_genes,)
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # stack_right_tmp (n_genes,)
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # workspace_weights (n_families,)
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # workspace_values (n_families,)
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # error_code (1,)
    ]
    f.restype = None
    return f

def setup_compute_rdi():
    f = lib.compute_rdi_c
    f.argtypes = [
        ctypes.c_int,  # n_genes
        ctypes.c_int,  # n_families
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # distances (n_genes,)
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # gene_to_fam (n_genes,)
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # dscale (n_families,)
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # rdi (n_genes,)
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # sorted_rdi (n_genes,)
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # perm (n_genes,)
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # stack_left (n_genes,)
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # stack_right (n_genes,)
    ]
    f.restype = None
    return f

def setup_identify_outliers():
    f = lib.identify_outliers_c
    f.argtypes = [
        ctypes.c_int,  # n_genes
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # rdi (n_genes,)
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # sorted_rdi (n_genes,)
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # is_outlier (n_genes,)
        ctypes.POINTER(ctypes.c_double),                                  # threshold (scalar)
        ctypes.c_double                                                  # percentile (scalar)
    ]
    f.restype = None
    return f

def setup_detect_outliers():
    f = lib.detect_outliers_c
    f.argtypes = [
        ctypes.c_int,  # n_genes
        ctypes.c_int,  # n_families
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # distances (n_genes,)
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # gene_to_fam (n_genes,)
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # work_array (n_genes,)
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # perm (n_genes,)
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # stack_left (n_genes,)
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # stack_right (n_genes,)
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # is_outlier (n_genes,)
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # loess_x (n_families,)
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # loess_y (n_families,)
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # loess_n (n_families,)
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # workspace_weights (n_families,)
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # workspace_values (n_families,)
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # error_code (1,)
        ctypes.c_double                                                  # percentile (scalar)
    ]
    f.restype = None
    return f

# Example: test for compute_family_scaling_c
def test_compute_family_scaling():
    print_header("[compute_family_scaling_c] Python LOESS-only test cases")
    n_genes = 5
    n_families = 2
    distances = np.array([1, 2, 3, 4, 5], dtype=np.float64)
    gene_to_fam = np.array([1, 1, 2, 2, 2], dtype=np.int32)
    dscale = np.zeros(n_families, dtype=np.float64)
    loess_x = np.zeros(n_families, dtype=np.float64)
    loess_y = np.zeros(n_families, dtype=np.float64)
    indices_used = np.zeros(n_families, dtype=np.int32)
    perm_tmp = np.zeros(n_genes, dtype=np.int32)
    stack_left_tmp = np.zeros(n_genes, dtype=np.int32)
    stack_right_tmp = np.zeros(n_genes, dtype=np.int32)
    workspace_weights = np.zeros(n_families, dtype=np.float64)
    workspace_values = np.zeros(n_families, dtype=np.float64)
    error_code = np.zeros(1, dtype=np.int32)
    f = setup_compute_family_scaling()
    # Case 1: LOESS-only, valid families
    print("\nCase 1: LOESS-only, valid families")
    f(n_genes, n_families, distances, gene_to_fam, dscale, loess_x, loess_y, indices_used, perm_tmp, stack_left_tmp, stack_right_tmp, workspace_weights, workspace_values, error_code)
    print('dscale:', dscale)
    print('error_code:', error_code[0])
    # Case 2: Invalid family index (should return error)
    print("\nCase 2: Invalid family index (should return error)")
    gene_to_fam_invalid = np.array([1, 3, 2, 2, 2], dtype=np.int32)
    f(n_genes, n_families, distances, gene_to_fam_invalid, dscale, loess_x, loess_y, indices_used, perm_tmp, stack_left_tmp, stack_right_tmp, workspace_weights, workspace_values, error_code)
    print('dscale:', dscale)
    print('error_code:', error_code[0])
    # Case 3: Family with only one gene (dscale=0)
    print("\nCase 3: Family with only one gene (dscale=0)")
    gene_to_fam_single = np.array([1, 2, 2, 2, 2], dtype=np.int32)
    f(n_genes, n_families, distances, gene_to_fam_single, dscale, loess_x, loess_y, indices_used, perm_tmp, stack_left_tmp, stack_right_tmp, workspace_weights, workspace_values, error_code)
    print('dscale:', dscale)
    print('error_code:', error_code[0])

def test_compute_rdi():
    print_header("[compute_rdi_c] Python test cases")
    n_families = 2
    distances = np.array([1, 2, 3, 4, 5], dtype=np.float64)
    gene_to_fam = np.array([1, 1, 2, 2, 2], dtype=np.int32)
    dscale = np.array([2, 4], dtype=np.float64)
    rdi = np.zeros(len(distances), dtype=np.float64)
    f = setup_compute_rdi()
    # [compute_rdi_c] Case 1: normal input
    print("\n[compute_rdi_c] Case 1: normal input")
    sorted_rdi = np.zeros(len(distances), dtype=np.float64)
    perm = np.arange(1, len(distances)+1, dtype=np.int32)
    stack_left = np.zeros(len(distances), dtype=np.int32)
    stack_right = np.zeros(len(distances), dtype=np.int32)
    f(len(distances), n_families, distances, gene_to_fam, dscale, rdi, sorted_rdi, perm, stack_left, stack_right)
    print("rdi:", rdi)
    print("sorted_rdi:", sorted_rdi)

    # [compute_rdi_c] Case 2: dscale with zeros
    print("\n[compute_rdi_c] Case 2: dscale with zeros")
    dscale_zeros = np.zeros(n_families, dtype=np.float64)
    rdi2 = np.zeros(len(distances), dtype=np.float64)
    sorted_rdi2 = np.zeros(len(distances), dtype=np.float64)
    perm2 = np.arange(1, len(distances)+1, dtype=np.int32)
    stack_left2 = np.zeros(len(distances), dtype=np.int32)
    stack_right2 = np.zeros(len(distances), dtype=np.int32)
    f(len(distances), n_families, distances, gene_to_fam, dscale_zeros, rdi2, sorted_rdi2, perm2, stack_left2, stack_right2)
    print("rdi:", rdi2)
    print("sorted_rdi:", sorted_rdi2)

    # [compute_rdi_c] Case 3: gene_to_fam out of range
    print("\n[compute_rdi_c] Case 3: gene_to_fam out of range")
    gene_to_fam_bad = np.array([1, 3, 2, 2, 2], dtype=np.int32)
    rdi3 = np.zeros(len(distances), dtype=np.float64)
    sorted_rdi3 = np.zeros(len(distances), dtype=np.float64)
    perm3 = np.arange(1, len(distances)+1, dtype=np.int32)
    stack_left3 = np.zeros(len(distances), dtype=np.int32)
    stack_right3 = np.zeros(len(distances), dtype=np.int32)
    f(len(distances), n_families, distances, gene_to_fam_bad, np.array([2, 4], dtype=np.float64), rdi3, sorted_rdi3, perm3, stack_left3, stack_right3)
    print("rdi:", rdi3)
    print("sorted_rdi:", sorted_rdi3)

    # [compute_rdi_c] Case 5: dscale shorter than max(gene_to_fam)
    print("\n[compute_rdi_c] Case 5: dscale shorter than max(gene_to_fam)")
    dscale_short = np.array([2], dtype=np.float64)
    rdi5 = np.zeros(len(distances), dtype=np.float64)
    sorted_rdi5 = np.zeros(len(distances), dtype=np.float64)
    perm5 = np.arange(1, len(distances)+1, dtype=np.int32)
    stack_left5 = np.zeros(len(distances), dtype=np.int32)
    stack_right5 = np.zeros(len(distances), dtype=np.int32)
    try:
        f(len(distances), n_families, distances, gene_to_fam, dscale_short, rdi5, sorted_rdi5, perm5, stack_left5, stack_right5)
        print("rdi:", rdi5)
        print("sorted_rdi:", sorted_rdi5)
    except Exception as e:
        print("Unexpected error:", e)

    # [compute_rdi_c] Case 6: vectors with incompatible length
    print("\n[compute_rdi_c] Case 6: vectors with incompatible length")
    distances_short = np.array([1, 2, 3], dtype=np.float64)
    gene_to_fam_short = np.array([1, 1, 2], dtype=np.int32)
    rdi6 = np.zeros(len(distances_short), dtype=np.float64)
    sorted_rdi6 = np.zeros(len(distances_short), dtype=np.float64)
    perm6 = np.arange(1, len(distances_short)+1, dtype=np.int32)
    stack_left6 = np.zeros(len(distances_short), dtype=np.int32)
    stack_right6 = np.zeros(len(distances_short), dtype=np.int32)
    try:
        f(len(distances_short), n_families, distances_short, gene_to_fam_short, np.array([2, 4], dtype=np.float64), rdi6, sorted_rdi6, perm6, stack_left6, stack_right6)
        print("rdi:", rdi6)
        print("sorted_rdi:", sorted_rdi6)
    except Exception as e:
        print("Unexpected error:", e)

def test_identify_outliers():
    print_header("[identify_outliers_c] Python test cases")
    f = setup_identify_outliers()
    # [identify_outliers_c] Case 1: Simple RDI, percentile 50
    # Expectation: Top 50% (highest 3) should be outliers (True), others False.
    print("\n[identify_outliers_c] Case 1: Simple RDI, percentile 50")
    rdi = np.array([0.3, 0.1, 0.5, 0.2, 0.4], dtype=np.float64)  # intentionally unsorted
    n = len(rdi)
    sorted_rdi = np.sort(rdi[rdi >= 0])
    if len(sorted_rdi) < n:
        sorted_rdi = np.concatenate([sorted_rdi, np.zeros(n - len(sorted_rdi))])
    is_outlier = np.zeros(n, dtype=np.int32)
    threshold = c_double(0.0)
    percentile = 50.0
    f(n, rdi, sorted_rdi, is_outlier, threshold, percentile)
    print("is_outlier:", is_outlier)
    print("threshold:", threshold.value)

    # [identify_outliers_c] Case 2: RDI with negative values (should be ignored)
    # Expectation: Negative RDI values are ignored for outlier detection; only positive values considered.
    print("\n[identify_outliers_c] Case 2: RDI with negative values (should be ignored)")
    rdi2 = np.array([0.3, -1, 0.5, 0.2, 0.4], dtype=np.float64)  # intentionally unsorted, with negative
    sorted_rdi2 = np.sort(rdi2[rdi2 >= 0])
    if len(sorted_rdi2) < n:
        sorted_rdi2 = np.concatenate([sorted_rdi2, np.zeros(n - len(sorted_rdi2))])
    is_outlier2 = np.zeros(n, dtype=np.int32)
    threshold2 = c_double(0.0)
    percentile2 = 80.0
    f(n, rdi2, sorted_rdi2, is_outlier2, threshold2, percentile2)
    print("is_outlier:", is_outlier2)
    print("threshold:", threshold2.value)

    # [identify_outliers_c] Case 3: All RDI zeros
    # Expectation: No outliers detected; all is_outlier should be False.
    print("\n[identify_outliers_c] Case 3: All RDI zeros")
    rdi3 = np.zeros(n, dtype=np.float64)
    sorted_rdi3 = np.sort(rdi3[rdi3 >= 0])
    if len(sorted_rdi3) < n:
        sorted_rdi3 = np.concatenate([sorted_rdi3, np.zeros(n - len(sorted_rdi3))])
    is_outlier3 = np.zeros(n, dtype=np.int32)
    threshold3 = c_double(0.0)
    percentile3 = 90.0
    f(n, rdi3, sorted_rdi3, is_outlier3, threshold3, percentile3)
    print("is_outlier:", is_outlier3)
    print("threshold:", threshold3.value)

    # [identify_outliers_c] Case 4: Percentile 0 (all outliers)
    # Expectation: All genes should be outliers (all True).
    print("\n[identify_outliers_c] Case 4: Percentile 0 (all outliers)")
    rdi4 = np.array([0.3, 0.1, 0.5, 0.2, 0.4], dtype=np.float64)  # intentionally unsorted
    sorted_rdi4 = np.sort(rdi4[rdi4 >= 0])
    if len(sorted_rdi4) < n:
        sorted_rdi4 = np.concatenate([sorted_rdi4, np.zeros(n - len(sorted_rdi4))])
    is_outlier4 = np.zeros(n, dtype=np.int32)
    threshold4 = c_double(0.0)
    percentile4 = 0.0
    f(n, rdi4, sorted_rdi4, is_outlier4, threshold4, percentile4)
    print("is_outlier:", is_outlier4)
    print("threshold:", threshold4.value)

    # [identify_outliers_c] Case 5: Percentile 100 (1 outlier)
    # Expectation: Only the highest RDI should be outlier (True), rest False.
    print("\n[identify_outliers_c] Case 5: Percentile 100 (1 outlier)")
    rdi5 = np.array([0.3, 0.1, 0.5, 0.2, 0.4], dtype=np.float64)  # intentionally unsorted
    sorted_rdi5 = np.sort(rdi5[rdi5 >= 0])
    if len(sorted_rdi5) < n:
        sorted_rdi5 = np.concatenate([sorted_rdi5, np.zeros(n - len(sorted_rdi5))])
    is_outlier5 = np.zeros(n, dtype=np.int32)
    threshold5 = c_double(0.0)
    percentile5 = 100.0
    f(n, rdi5, sorted_rdi5, is_outlier5, threshold5, percentile5)
    print("is_outlier:", is_outlier5)
    print("threshold:", threshold5.value)

    # [identify_outliers_c] Case 6: All RDI negative (should be ignored)
    # Expectation: All negative RDI values should be ignored for outlier detection; all is_outlier should be False.
    print("\n[identify_outliers_c] Case 6: All RDI negative (should be ignored)")
    rdi6 = np.array([-0.1, -0.2, -0.3, -0.4, -0.5], dtype=np.float64)
    sorted_rdi6 = np.sort(rdi6[rdi6 >= 0])
    if len(sorted_rdi6) < n:
        sorted_rdi6 = np.concatenate([sorted_rdi6, np.zeros(n - len(sorted_rdi6))])
    is_outlier6 = np.zeros(n, dtype=np.int32)
    threshold6 = c_double(0.0)
    percentile6 = 80.0
    f(n, rdi6, sorted_rdi6, is_outlier6, threshold6, percentile6)
    print("is_outlier:", is_outlier6)
    print("threshold:", threshold6.value)

def test_detect_outliers():
    print_header("[detect_outliers_c] Python LOESS-only test cases")
    f = setup_detect_outliers()
    # Case 1: LOESS-only, valid families
    n_genes = 6
    n_families = 2
    distances = np.array([1, 2, 3, 4, 5, 6], dtype=np.float64)
    gene_to_fam = np.array([1, 1, 2, 2, 2, 2], dtype=np.int32)
    work_array = np.zeros(n_genes, dtype=np.float64)
    perm = np.arange(1, n_genes+1, dtype=np.int32)
    stack_left = np.zeros(n_genes, dtype=np.int32)
    stack_right = np.zeros(n_genes, dtype=np.int32)
    is_outlier = np.zeros(n_genes, dtype=np.int32)
    loess_x = np.zeros(n_families, dtype=np.float64)
    loess_y = np.zeros(n_families, dtype=np.float64)
    loess_n = np.zeros(n_families, dtype=np.int32)
    workspace_weights = np.zeros(n_families, dtype=np.float64)
    workspace_values = np.zeros((1, n_families), dtype=np.float64)
    error_code = np.zeros(1, dtype=np.int32)
    percentile = 80.0
    f(n_genes, n_families, distances, gene_to_fam, work_array, perm, stack_left, stack_right, is_outlier, loess_x, loess_y, loess_n, workspace_weights, workspace_values, error_code, percentile)
    print("is_outlier:", is_outlier)
    print("error_code:", error_code[0])
    # Case 2: Invalid family index (should return error)
    print("\nCase 2: Invalid family index (should return error)")
    is_outlier[:] = 0
    gene_to_fam_invalid = np.array([1, 3, 2, 2, 2, 2], dtype=np.int32)
    workspace_values = np.zeros((1, n_families), dtype=np.float64)
    f(n_genes, n_families, distances, gene_to_fam_invalid, work_array, perm, stack_left, stack_right, is_outlier, loess_x, loess_y, loess_n, workspace_weights, workspace_values, error_code, percentile)
    print("is_outlier:", is_outlier)
    print("error_code:", error_code[0])

def main():
    test_compute_family_scaling()
    test_compute_rdi()
    test_identify_outliers()
    test_detect_outliers()
    print("\nAll outlier detection tests completed.")

if __name__ == "__main__":
    main()
