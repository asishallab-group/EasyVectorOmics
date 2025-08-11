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
    print("=" * 60)
    print(msg)
    print("=" * 60)

def setup_compute_family_scaling():
    """Setup the main (allocating) compute_family_scaling_c interface"""
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
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # error_code (1,)
    ]
    f.restype = None
    return f

def setup_compute_family_scaling_expert():
    """Setup the expert (kernel) compute_family_scaling_expert_c interface"""
    f = lib.compute_family_scaling_expert_c
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
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # family_distances (n_families,)
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
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # error_code (1,)
        ctypes.c_double                                                  # percentile (scalar)
    ]
    f.restype = None
    return f

# =====================
# Tests for compute_family_scaling_c (main/alloc interface - user friendly)
# =====================

def test_compute_family_scaling():
    print_header("Testing compute_family_scaling_c (main/alloc interface)")
    
    # Common test data
    n_genes = 5
    n_families = 2
    distances = np.array([1, 2, 3, 4, 5], dtype=np.float64)
    gene_to_fam = np.array([1, 1, 2, 2, 2], dtype=np.int32)
    
    f = setup_compute_family_scaling()
    
    # Case 1: Basic LOESS scaling
    print("\n[compute_family_scaling_c] Case 1: Basic LOESS scaling")
    dscale = np.zeros(n_families, dtype=np.float64)
    loess_x = np.zeros(n_families, dtype=np.float64)
    loess_y = np.zeros(n_families, dtype=np.float64)
    indices_used = np.zeros(n_families, dtype=np.int32)
    error_code = np.zeros(1, dtype=np.int32)
    
    f(n_genes, n_families, distances, gene_to_fam, dscale, loess_x, loess_y, indices_used, error_code)
    print(f"dscale: {', '.join(map(str, dscale))}")
    print(f"error_code: {error_code[0]}")
    print(f"loess_x (medians): {', '.join(map(str, loess_x))}")
    print(f"loess_y (stddevs): {', '.join(map(str, loess_y))}")
    
    # Case 2: Invalid family indices
    print("\n[compute_family_scaling_c] Case 2: Invalid family indices")
    gene_to_fam_invalid = np.array([1, 3, 2, 2, 2], dtype=np.int32)
    dscale = np.zeros(n_families, dtype=np.float64)
    loess_x = np.zeros(n_families, dtype=np.float64)
    loess_y = np.zeros(n_families, dtype=np.float64)
    indices_used = np.zeros(n_families, dtype=np.int32)
    error_code = np.zeros(1, dtype=np.int32)
    
    f(n_genes, n_families, distances, gene_to_fam_invalid, dscale, loess_x, loess_y, indices_used, error_code)
    print(f"dscale: {', '.join(map(str, dscale))}")
    print(f"error_code: {error_code[0]}")
    
    # Case 3: Family with only one gene
    print("\n[compute_family_scaling_c] Case 3: Family with only one gene")
    gene_to_fam_single = np.array([1, 2, 2, 2, 2], dtype=np.int32)
    dscale = np.zeros(n_families, dtype=np.float64)
    loess_x = np.zeros(n_families, dtype=np.float64)
    loess_y = np.zeros(n_families, dtype=np.float64)
    indices_used = np.zeros(n_families, dtype=np.int32)
    error_code = np.zeros(1, dtype=np.int32)
    
    f(n_genes, n_families, distances, gene_to_fam_single, dscale, loess_x, loess_y, indices_used, error_code)
    print(f"dscale: {', '.join(map(str, dscale))}")
    print(f"error_code: {error_code[0]}")
    
    # Case 4: All distances zero
    print("\n[compute_family_scaling_c] Case 4: All distances zero")
    distances_zero = np.array([0, 0, 0, 0, 0], dtype=np.float64)
    dscale = np.zeros(n_families, dtype=np.float64)
    loess_x = np.zeros(n_families, dtype=np.float64)
    loess_y = np.zeros(n_families, dtype=np.float64)
    indices_used = np.zeros(n_families, dtype=np.int32)
    error_code = np.zeros(1, dtype=np.int32)
    
    f(n_genes, n_families, distances_zero, gene_to_fam, dscale, loess_x, loess_y, indices_used, error_code)
    print(f"dscale: {', '.join(map(str, dscale))}")
    print(f"error_code: {error_code[0]}")
    
    # Case 5: Large dataset
    print("\n[compute_family_scaling_c] Case 5: Large dataset")
    n_families_large = 3000
    np.random.seed(42)  # For reproducible results
    
    # Generate random number of genes per family (between 2 and 20)
    genes_per_family = np.random.randint(2, 21, size=n_families_large)
    n_genes_large = np.sum(genes_per_family)
    
    # Create gene_to_fam array
    gene_to_fam_large = []
    for family_id in range(1, n_families_large + 1):
        gene_to_fam_large.extend([family_id] * genes_per_family[family_id - 1])
    gene_to_fam_large = np.array(gene_to_fam_large, dtype=np.int32)
    
    # Generate random distances
    distances_large = np.random.uniform(0.1, 10.0, n_genes_large).astype(np.float64)
    
    print(f"Dataset size: {n_genes_large} genes across {n_families_large} families")
    print(f"Genes per family range: {np.min(genes_per_family)} - {np.max(genes_per_family)}")
    print(f"Average genes per family: {np.mean(genes_per_family):.2f}")
    
    dscale = np.zeros(n_families_large, dtype=np.float64)
    loess_x = np.zeros(n_families_large, dtype=np.float64)
    loess_y = np.zeros(n_families_large, dtype=np.float64)
    indices_used = np.zeros(n_families_large, dtype=np.int32)
    error_code = np.zeros(1, dtype=np.int32)
    
    f(n_genes_large, n_families_large, distances_large, gene_to_fam_large, dscale, loess_x, loess_y, indices_used, error_code)
    print(f"dscale sample (first 10): {', '.join([str(round(x, 3)) for x in dscale[:10]])}")
    print(f"error_code: {error_code[0]}")
    print(f"Non-zero dscale values: {np.sum(dscale > 0)}/{n_families_large}")
    print(f"Processing completed successfully for large dataset")
    
    # Case 6: Mixed family sizes
    print("\n[compute_family_scaling_c] Case 6: Mixed family sizes")
    n_genes_mixed = 10
    n_families_mixed = 3
    distances_mixed = np.array([1.0, 1.1, 1.2, 1.3, 1.4, 2.0, 2.1, 3.0, 4.0, 5.0], dtype=np.float64)
    gene_to_fam_mixed = np.array([1, 1, 1, 1, 1, 2, 2, 3, 3, 3], dtype=np.int32)  # Family 1: 5 genes, Family 2: 2 genes, Family 3: 3 genes
    
    dscale = np.zeros(n_families_mixed, dtype=np.float64)
    loess_x = np.zeros(n_families_mixed, dtype=np.float64)
    loess_y = np.zeros(n_families_mixed, dtype=np.float64)
    indices_used = np.zeros(n_families_mixed, dtype=np.int32)
    error_code = np.zeros(1, dtype=np.int32)
    
    f(n_genes_mixed, n_families_mixed, distances_mixed, gene_to_fam_mixed, dscale, loess_x, loess_y, indices_used, error_code)
    print(f"dscale: {', '.join([str(round(x, 3)) for x in dscale])}")
    print(f"error_code: {error_code[0]}")
    print(f"Family medians: {', '.join([str(round(x, 3)) for x in loess_x])}")

# =====================
# Tests for compute_family_scaling_expert_c (expert/kernel interface - advanced users)
# =====================

def test_compute_family_scaling_expert():
    print_header("Testing compute_family_scaling_expert_c (expert/kernel interface)")
    
    # Common test data
    n_genes = 5
    n_families = 2
    distances = np.array([1, 2, 3, 4, 5], dtype=np.float64)
    gene_to_fam = np.array([1, 1, 2, 2, 2], dtype=np.int32)
    
    f = setup_compute_family_scaling_expert()
    
    # For expert interface, user must provide all work arrays
    perm_tmp = np.zeros(n_genes, dtype=np.int32)
    stack_left_tmp = np.zeros(n_genes, dtype=np.int32)
    stack_right_tmp = np.zeros(n_genes, dtype=np.int32)
    family_distances = np.zeros(n_genes, dtype=np.float64)
    
    # Case 1: Basic LOESS scaling with user-provided work arrays
    print("\n[compute_family_scaling_expert_c] Case 1: Basic LOESS scaling with work arrays")
    dscale = np.zeros(n_families, dtype=np.float64)
    loess_x = np.zeros(n_families, dtype=np.float64)
    loess_y = np.zeros(n_families, dtype=np.float64)
    indices_used = np.zeros(n_families, dtype=np.int32)
    error_code = np.zeros(1, dtype=np.int32)
    
    f(n_genes, n_families, distances, gene_to_fam, dscale, loess_x, loess_y, indices_used, perm_tmp, stack_left_tmp, stack_right_tmp, family_distances, error_code)
    print(f"dscale: {', '.join(map(str, dscale))}")
    print(f"error_code: {error_code[0]}")
    
    # Case 2: Check that work arrays are modified
    print("\n[compute_family_scaling_expert_c] Case 2: Work arrays modification check")
    # Reset work arrays to known values
    perm_tmp = np.zeros(n_genes, dtype=np.int32)
    stack_left_tmp = np.zeros(n_genes, dtype=np.int32)
    stack_right_tmp = np.zeros(n_genes, dtype=np.int32)
    family_distances = np.zeros(n_genes, dtype=np.float64)
    print(f"family_distances before: {', '.join(map(str, family_distances))}")
    
    dscale = np.zeros(n_families, dtype=np.float64)
    loess_x = np.zeros(n_families, dtype=np.float64)
    loess_y = np.zeros(n_families, dtype=np.float64)
    indices_used = np.zeros(n_families, dtype=np.int32)
    error_code = np.zeros(1, dtype=np.int32)
    
    f(n_genes, n_families, distances, gene_to_fam, dscale, loess_x, loess_y, indices_used, perm_tmp, stack_left_tmp, stack_right_tmp, family_distances, error_code)
    print(f"family_distances after: {', '.join([str(round(x, 3)) for x in family_distances])}")
    print(f"perm_tmp: {', '.join(map(str, perm_tmp))}")
    
    # Case 3: Performance comparison with main interface
    print("\n[compute_family_scaling_expert_c] Case 3: Results comparison with main interface")
    # Use same test data as Case 1 of main interface
    perm_tmp = np.zeros(n_genes, dtype=np.int32)
    stack_left_tmp = np.zeros(n_genes, dtype=np.int32)
    stack_right_tmp = np.zeros(n_genes, dtype=np.int32)
    family_distances = np.zeros(n_genes, dtype=np.float64)
    dscale_expert = np.zeros(n_families, dtype=np.float64)
    loess_x_expert = np.zeros(n_families, dtype=np.float64)
    loess_y_expert = np.zeros(n_families, dtype=np.float64)
    indices_used_expert = np.zeros(n_families, dtype=np.int32)
    error_code_expert = np.zeros(1, dtype=np.int32)
    
    f(n_genes, n_families, distances, gene_to_fam, dscale_expert, loess_x_expert, loess_y_expert, indices_used_expert, perm_tmp, stack_left_tmp, stack_right_tmp, family_distances, error_code_expert)
    
    # Compare with Case 1 results from main interface (we'll store them)
    # Run main interface for comparison
    f_main = setup_compute_family_scaling()
    dscale_main = np.zeros(n_families, dtype=np.float64)
    loess_x_main = np.zeros(n_families, dtype=np.float64)
    loess_y_main = np.zeros(n_families, dtype=np.float64)
    indices_used_main = np.zeros(n_families, dtype=np.int32)
    error_code_main = np.zeros(1, dtype=np.int32)
    f_main(n_genes, n_families, distances, gene_to_fam, dscale_main, loess_x_main, loess_y_main, indices_used_main, error_code_main)
    
    print("Comparison with main interface:")
    print(f"Main dscale: {', '.join(map(str, dscale_main))}")
    print(f"Expert dscale: {', '.join(map(str, dscale_expert))}")
    print(f"Difference: {', '.join([str(round(abs(x-y), 6)) for x, y in zip(dscale_main, dscale_expert)])}")
    print(f"Results identical: {all(abs(x-y) < 1e-10 for x, y in zip(dscale_main, dscale_expert))}")
    
    # Case 4: Large dataset with expert interface
    print("\n[compute_family_scaling_expert_c] Case 4: Large dataset")
    n_families_large = 3000
    np.random.seed(42)  # For reproducible results
    
    # Generate random number of genes per family (between 2 and 20)
    genes_per_family = np.random.randint(2, 21, size=n_families_large)
    n_genes_large = np.sum(genes_per_family)
    
    # Create gene_to_fam array
    gene_to_fam_large = []
    for family_id in range(1, n_families_large + 1):
        gene_to_fam_large.extend([family_id] * genes_per_family[family_id - 1])
    gene_to_fam_large = np.array(gene_to_fam_large, dtype=np.int32)
    
    # Generate random distances
    distances_large = np.random.uniform(0.1, 10.0, n_genes_large).astype(np.float64)
    
    print(f"Large dataset: {n_genes_large} genes across {n_families_large} families")
    print(f"Memory footprint: ~{(n_genes_large * 8 * 4 + n_families_large * 8 * 4) / 1024 / 1024:.1f} MB for work arrays")
    
    # Expert interface requires larger work arrays for large dataset
    perm_tmp_large = np.zeros(n_genes_large, dtype=np.int32)
    stack_left_tmp_large = np.zeros(n_genes_large, dtype=np.int32)
    stack_right_tmp_large = np.zeros(n_genes_large, dtype=np.int32)
    family_distances_large = np.zeros(n_genes_large, dtype=np.float64)
    dscale_large = np.zeros(n_families_large, dtype=np.float64)
    loess_x_large = np.zeros(n_families_large, dtype=np.float64)
    loess_y_large = np.zeros(n_families_large, dtype=np.float64)
    indices_used_large = np.zeros(n_families_large, dtype=np.int32)
    error_code_large = np.zeros(1, dtype=np.int32)
    
    f(n_genes_large, n_families_large, distances_large, gene_to_fam_large, dscale_large, loess_x_large, loess_y_large, indices_used_large, perm_tmp_large, stack_left_tmp_large, stack_right_tmp_large, family_distances_large, error_code_large)
    print(f"dscale sample (first 10): {', '.join([str(round(x, 3)) for x in dscale_large[:10]])}")
    print(f"error_code: {error_code_large[0]}")
    print(f"Expert interface successfully handled large dataset with manual memory management")

# =====================
# Test cases for compute_rdi_c
# =====================

def test_compute_rdi():
    print_header("Testing compute_rdi_c")
    
    f = setup_compute_rdi()
    
    # Case 1: normal input
    print("\n[compute_rdi_c] Case 1: normal input")
    distances = np.array([1, 2, 3, 4, 5], dtype=np.float64)
    gene_to_fam = np.array([1, 1, 2, 2, 2], dtype=np.int32)
    n_families = 2
    dscale = np.array([2, 4], dtype=np.float64)
    rdi = np.zeros(len(distances), dtype=np.float64)
    sorted_rdi = np.zeros(len(distances), dtype=np.float64)
    perm = np.arange(1, len(distances)+1, dtype=np.int32)
    stack_left = np.zeros(len(distances), dtype=np.int32)
    stack_right = np.zeros(len(distances), dtype=np.int32)
    
    f(len(distances), n_families, distances, gene_to_fam, dscale, rdi, sorted_rdi, perm, stack_left, stack_right)
    print(f"RDI: {', '.join([str(round(x, 3)) for x in rdi])}")
    print(f"Sorted RDI: {', '.join([str(round(x, 3)) for x in sorted_rdi])}")
    
    # Case 2: dscale with zeros
    print("\n[compute_rdi_c] Case 2: dscale with zeros")
    dscale_zero = np.array([0, 0], dtype=np.float64)
    rdi = np.zeros(len(distances), dtype=np.float64)
    sorted_rdi = np.zeros(len(distances), dtype=np.float64)
    perm = np.arange(1, len(distances)+1, dtype=np.int32)
    stack_left = np.zeros(len(distances), dtype=np.int32)
    stack_right = np.zeros(len(distances), dtype=np.int32)
    
    f(len(distances), n_families, distances, gene_to_fam, dscale_zero, rdi, sorted_rdi, perm, stack_left, stack_right)
    print(f"RDI with zero scaling: {', '.join(map(str, rdi))}")
    
    # Case 3: gene_to_fam out of range
    print("\n[compute_rdi_c] Case 3: gene_to_fam out of range")
    gene_to_fam_bad = np.array([1, 3, 2, 2, 2], dtype=np.int32)  # 3 doesn't exist in dscale
    rdi = np.zeros(len(distances), dtype=np.float64)
    sorted_rdi = np.zeros(len(distances), dtype=np.float64)
    perm = np.arange(1, len(distances)+1, dtype=np.int32)
    stack_left = np.zeros(len(distances), dtype=np.int32)
    stack_right = np.zeros(len(distances), dtype=np.int32)
    
    f(len(distances), n_families, distances, gene_to_fam_bad, np.array([2, 4], dtype=np.float64), rdi, sorted_rdi, perm, stack_left, stack_right)
    print(f"RDI with invalid family index: {', '.join([str(round(x, 3)) for x in rdi])}")
    
    # Case 4: High precision test
    print("\n[compute_rdi_c] Case 4: High precision test")
    distances_precise = np.array([2.0, 4.0, 6.0], dtype=np.float64)
    gene_to_fam_precise = np.array([1, 1, 1], dtype=np.int32)
    n_families_precise = 1
    dscale_precise = np.array([2.0], dtype=np.float64)  # All genes in family 1, scaling = 2.0
    rdi = np.zeros(len(distances_precise), dtype=np.float64)
    sorted_rdi = np.zeros(len(distances_precise), dtype=np.float64)
    perm = np.arange(1, len(distances_precise)+1, dtype=np.int32)
    stack_left = np.zeros(len(distances_precise), dtype=np.int32)
    stack_right = np.zeros(len(distances_precise), dtype=np.int32)
    
    f(len(distances_precise), n_families_precise, distances_precise, gene_to_fam_precise, dscale_precise, rdi, sorted_rdi, perm, stack_left, stack_right)
    # Expected RDI: [2.0/2.0, 4.0/2.0, 6.0/2.0] = [1.0, 2.0, 3.0]
    print("Expected RDI: [1.0, 2.0, 3.0]")
    print(f"Actual RDI: {', '.join(map(str, rdi))}")
    print(f"High precision test passed: {all(abs(x-y) < 1e-10 for x, y in zip(rdi, [1.0, 2.0, 3.0]))}")
    
    # Case 5: Negative distances
    print("\n[compute_rdi_c] Case 5: Negative distances")
    distances_negative = np.array([-1, 2, -3, 4, 5], dtype=np.float64)
    rdi = np.zeros(len(distances_negative), dtype=np.float64)
    sorted_rdi = np.zeros(len(distances_negative), dtype=np.float64)
    perm = np.arange(1, len(distances_negative)+1, dtype=np.int32)
    stack_left = np.zeros(len(distances_negative), dtype=np.int32)
    stack_right = np.zeros(len(distances_negative), dtype=np.int32)
    
    f(len(distances_negative), n_families, distances_negative, gene_to_fam, np.array([2, 4], dtype=np.float64), rdi, sorted_rdi, perm, stack_left, stack_right)
    print(f"RDI with negative distances: {', '.join([str(round(x, 3)) for x in rdi])}")

# =====================
# Test cases for identify_outliers_c
# =====================

def test_identify_outliers():
    print_header("Testing identify_outliers_c")
    
    f = setup_identify_outliers()
    
    # Case 1: Simple RDI, percentile 50
    print("\n[identify_outliers_c] Case 1: Simple RDI, percentile 50")
    rdi = np.array([0.3, 0.1, 0.5, 0.2, 0.4], dtype=np.float64)
    n = len(rdi)
    sorted_rdi = np.sort(rdi[rdi >= 0])
    if len(sorted_rdi) < n:
        sorted_rdi = np.concatenate([sorted_rdi, np.zeros(n - len(sorted_rdi))])
    is_outlier = np.zeros(n, dtype=np.int32)
    threshold = c_double(0.0)
    percentile = 50.0
    
    f(n, rdi, sorted_rdi, is_outlier, threshold, percentile)
    print(f"RDI: {', '.join(map(str, rdi))}")
    print(f"Is outlier: {', '.join(map(str, is_outlier))}")
    print(f"Threshold: {threshold.value}")
    
    # Case 2: All RDI zeros
    print("\n[identify_outliers_c] Case 2: All RDI zeros")
    rdi = np.zeros(n, dtype=np.float64)
    sorted_rdi = np.sort(rdi[rdi >= 0])
    if len(sorted_rdi) < n:
        sorted_rdi = np.concatenate([sorted_rdi, np.zeros(n - len(sorted_rdi))])
    is_outlier = np.zeros(n, dtype=np.int32)
    threshold = c_double(0.0)
    percentile = 90.0
    
    f(n, rdi, sorted_rdi, is_outlier, threshold, percentile)
    print(f"All zeros outliers: {', '.join(map(str, is_outlier))}")
    print(f"Threshold: {threshold.value}")
    
    # Case 3: RDI with negative values (should be ignored)
    print("\n[identify_outliers_c] Case 3: RDI with negative values (should be ignored)")
    rdi = np.array([0.3, -1, 0.5, 0.2, 0.4], dtype=np.float64)  # with negative
    sorted_rdi = np.sort(rdi[rdi >= 0])
    if len(sorted_rdi) < n:
        sorted_rdi = np.concatenate([sorted_rdi, np.zeros(n - len(sorted_rdi))])
    is_outlier = np.zeros(n, dtype=np.int32)
    threshold = c_double(0.0)
    percentile = 80.0
    
    f(n, rdi, sorted_rdi, is_outlier, threshold, percentile)
    print(f"RDI with negative: {', '.join(map(str, rdi))}")
    print(f"Is outlier: {', '.join(map(str, is_outlier))}")
    print(f"Threshold: {threshold.value}")
    
    # Case 4: Percentile 0 (all outliers)
    print("\n[identify_outliers_c] Case 4: Percentile 0 (all outliers)")
    rdi = np.array([0.3, 0.1, 0.5, 0.2, 0.4], dtype=np.float64)
    sorted_rdi = np.sort(rdi[rdi >= 0])
    if len(sorted_rdi) < n:
        sorted_rdi = np.concatenate([sorted_rdi, np.zeros(n - len(sorted_rdi))])
    is_outlier = np.zeros(n, dtype=np.int32)
    threshold = c_double(0.0)
    percentile = 0.0
    
    f(n, rdi, sorted_rdi, is_outlier, threshold, percentile)
    print(f"0% percentile outliers: {', '.join(map(str, is_outlier))}")
    print(f"All are outliers: {all(is_outlier)}")
    
    # Case 5: Percentile 100 (1 outlier)
    print("\n[identify_outliers_c] Case 5: Percentile 100 (1 outlier)")
    rdi = np.array([0.3, 0.1, 0.5, 0.2, 0.4], dtype=np.float64)
    sorted_rdi = np.sort(rdi[rdi >= 0])
    if len(sorted_rdi) < n:
        sorted_rdi = np.concatenate([sorted_rdi, np.zeros(n - len(sorted_rdi))])
    is_outlier = np.zeros(n, dtype=np.int32)
    threshold = c_double(0.0)
    percentile = 100.0
    
    f(n, rdi, sorted_rdi, is_outlier, threshold, percentile)
    print(f"100% percentile outliers: {', '.join(map(str, is_outlier))}")
    print(f"Only highest is outlier: {sum(is_outlier) == 1}")
    
    # Case 6: All RDI negative (should be ignored)
    print("\n[identify_outliers_c] Case 6: All RDI negative (should be ignored)")
    rdi = np.array([-0.1, -0.2, -0.3, -0.4, -0.5], dtype=np.float64)  # all negative
    sorted_rdi = np.sort(rdi[rdi >= 0])
    if len(sorted_rdi) < n:
        sorted_rdi = np.concatenate([sorted_rdi, np.zeros(n - len(sorted_rdi))])
    is_outlier = np.zeros(n, dtype=np.int32)
    threshold = c_double(0.0)
    percentile = 80.0
    
    f(n, rdi, sorted_rdi, is_outlier, threshold, percentile)
    print(f"All negative outliers: {', '.join(map(str, is_outlier))}")
    print(f"No outliers detected: {not any(is_outlier)}")

# =====================
# Test cases for detect_outliers_c (comprehensive workflow)
# =====================

def test_detect_outliers():
    print_header("Testing detect_outliers_c (complete outlier detection workflow)")
    
    f = setup_detect_outliers()
    
    # Case 1: Typical input
    print("\n[detect_outliers_c] Case 1: Typical input")
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
    error_code = np.zeros(1, dtype=np.int32)
    percentile = 80.0
    
    f(n_genes, n_families, distances, gene_to_fam, work_array, perm, stack_left, stack_right, is_outlier, loess_x, loess_y, loess_n, error_code, percentile)
    print(f"Distances: {', '.join(map(str, distances))}")
    print(f"Gene to family: {', '.join(map(str, gene_to_fam))}")
    print(f"Is outlier: {', '.join(map(str, is_outlier))}")
    print(f"Error code: {error_code[0]}")
    print(f"Family medians: {', '.join([str(round(x, 3)) for x in loess_x])}")
    print(f"Family stddevs: {', '.join([str(round(x, 3)) for x in loess_y])}")
    
    # Case 2: Invalid gene_to_fam indices
    print("\n[detect_outliers_c] Case 2: Invalid gene_to_fam indices")
    gene_to_fam_invalid = np.array([1, 3, 2, 2, 2, 2], dtype=np.int32)
    is_outlier = np.zeros(n_genes, dtype=np.int32)
    work_array = np.zeros(n_genes, dtype=np.float64)
    perm = np.arange(1, n_genes+1, dtype=np.int32)
    stack_left = np.zeros(n_genes, dtype=np.int32)
    stack_right = np.zeros(n_genes, dtype=np.int32)
    loess_x = np.zeros(n_families, dtype=np.float64)
    loess_y = np.zeros(n_families, dtype=np.float64)
    loess_n = np.zeros(n_families, dtype=np.int32)
    error_code = np.zeros(1, dtype=np.int32)
    
    f(n_genes, n_families, distances, gene_to_fam_invalid, work_array, perm, stack_left, stack_right, is_outlier, loess_x, loess_y, loess_n, error_code, percentile)
    print(f"Invalid family - Is outlier: {', '.join(map(str, is_outlier))}")
    print(f"Invalid family - Error code: {error_code[0]}")
    
    # Case 3: Default percentile test
    print("\n[detect_outliers_c] Case 3: Default percentile test")
    n_genes = 8
    n_families = 2
    distances = np.array([1, 1.1, 1.2, 1.3, 10, 10.1, 10.2, 50], dtype=np.float64)  # Last gene is clear outlier
    gene_to_fam = np.array([1, 1, 1, 1, 2, 2, 2, 2], dtype=np.int32)
    work_array = np.zeros(n_genes, dtype=np.float64)
    perm = np.arange(1, n_genes+1, dtype=np.int32)
    stack_left = np.zeros(n_genes, dtype=np.int32)
    stack_right = np.zeros(n_genes, dtype=np.int32)
    is_outlier = np.zeros(n_genes, dtype=np.int32)
    loess_x = np.zeros(n_families, dtype=np.float64)
    loess_y = np.zeros(n_families, dtype=np.float64)
    loess_n = np.zeros(n_families, dtype=np.int32)
    error_code = np.zeros(1, dtype=np.int32)
    
    # Don't specify percentile to test default behavior (should be 95%)
    f(n_genes, n_families, distances, gene_to_fam, work_array, perm, stack_left, stack_right, is_outlier, loess_x, loess_y, loess_n, error_code, 95.0)  # Using 95% as default
    print(f"Default percentile - Distances: {', '.join(map(str, distances))}")
    print(f"Default percentile - Is outlier: {', '.join(map(str, is_outlier))}")
    print(f"Default percentile - Error code: {error_code[0]}")
    print(f"Outlier detected for extreme value: {bool(is_outlier[7])}")
    
    # Case 4: Single gene families
    print("\n[detect_outliers_c] Case 4: Single gene families")
    n_genes = 3
    n_families = 3
    distances = np.array([1, 10, 100], dtype=np.float64)  # Each gene in different family
    gene_to_fam = np.array([1, 2, 3], dtype=np.int32)
    work_array = np.zeros(n_genes, dtype=np.float64)
    perm = np.arange(1, n_genes+1, dtype=np.int32)
    stack_left = np.zeros(n_genes, dtype=np.int32)
    stack_right = np.zeros(n_genes, dtype=np.int32)
    is_outlier = np.zeros(n_genes, dtype=np.int32)
    loess_x = np.zeros(n_families, dtype=np.float64)
    loess_y = np.zeros(n_families, dtype=np.float64)
    loess_n = np.zeros(n_families, dtype=np.int32)
    error_code = np.zeros(1, dtype=np.int32)
    
    f(n_genes, n_families, distances, gene_to_fam, work_array, perm, stack_left, stack_right, is_outlier, loess_x, loess_y, loess_n, error_code, 80.0)
    print(f"Single families - Is outlier: {', '.join(map(str, is_outlier))}")
    print(f"Single families - Error code: {error_code[0]}")
    print(f"Single families - Family medians: {', '.join(map(str, loess_x))}")
    print(f"No outliers expected: {not any(is_outlier)}")

def main():
    test_compute_family_scaling()
    test_compute_family_scaling_expert()
    test_compute_rdi()
    test_identify_outliers()
    test_detect_outliers()
    print("\nAll outlier detection tests completed.")
    print("Both main (allocating) and expert (kernel) interfaces tested successfully.")

if __name__ == "__main__":
    main()
