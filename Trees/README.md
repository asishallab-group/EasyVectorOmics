# Tree Structures for Tensor-Omics

This directory provides efficient Fortran implementations of tree-based data structures for high-performance scientific computing, focusing on **Binary Search Trees (BSTs)** and **k-d Trees (KD-Trees)**. These structures are designed for fast range queries, sorting, and multidimensional data indexing.

## Overview

- **Binary Search Tree (BST):**  
  Flat-index-based utilities for 1D range queries and sorted access over real-valued arrays.

- **k-d Tree (KD-Tree):**  
  Stack-based, non-recursive construction for multidimensional data, supporting both Cartesian and spherical data.

## Directory Structure

```
src/
├── binary_search/
│   ├── binary_search_tree.f90   # BST utilities for 1D queries
│   └── bst_test.f90             # Unit tests for BST
├── k-d_tree/
│   ├── k_d_tree.f90             # KD-Tree construction and utilities
│   └── test_kd.f90              # Unit and edge-case tests for KD-Tree
├── tox_sorting.f90              # Indirect sorting utilities (quicksort)
```

## Features

- **BST:**
  - Build sorted index for real arrays.
  - Efficient 1D range queries.
  - Access sorted values without data movement.

- **KD-Tree:**
  - Non-recursive, stack-based construction.
  - Handles arbitrary dimensions and edge cases.
  - Supports both Cartesian and spherical data.
  - Dimension order can be set by variance or custom.

- **Sorting:**
  - Indirect, stable quicksort for real, integer, and character arrays.

## Usage

### Binary Search Tree

```fortran
use binary_search_tree
! ...initialize x(:), n...
call build_bst_index(x, n, ix, stack_left, stack_right)
val = get_sorted_value(x, ix, i)
call bst_range_query(x, ix, n, lo, hi, out_ix, out_n)
```

### KD-Tree

```fortran
use kd_tree
! ...initialize X(d, n), dim_order(d)...
call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, stack_left, stack_right)
call build_spherical_kd(V, d, n, sphere_ix, dim_order, work, subarray, perm, stack_left, stack_right)
```

### Sorting Utilities

```fortran
use tox_sorting
call sort_array(array, perm, stack_left, stack_right)
```

## Testing

- Run `bst_test.f90` for BST validation.
- Run `test_kd.f90` for KD-Tree and edge-case tests.

## References

- See source files in `src/` for detailed documentation and comments.
- Sorting routines are in `src/tox_sorting.f90`.

---

For further details, consult the docstrings and comments in each module.

The test_kd currently gives the following output:
```
 2D Cartesian KD Tree index:           1           1           2           4           3           3
 FAIL:2D cartesiannot a permutation:           1           1           2           4           3           3
 3D Spherical KD Tree index:           1           7           3           7           1           1           8           8
 FAIL:3D sphericalnot a permutation:           1           7           3           7           1           1           8           8
 Edge n=0 passed.
 Edge n=1 passed.
 PASS:identical points
 FAIL:unit vectorsnot a permutation:           3           4           1           1
 PASS:high-d low-n
 FAIL:1D sortedindex out of bounds:          49
 All edge-case tests completed.
```
suggesting that the build_kd_index is not working as intended.

The BST_test returns the following:
```
BST index test PASSED: x(ix) is monotonic non-decreasing.
 get_sorted_value(x, ix, n/2) =   0.49634853005409241     
 bst_range_query: Number of values in [0.2, 0.8] =         610
 First value in range:   0.20035192454263662     
 Last value in range:   0.79847741414482898 
```
Therefore, the binary search tree should be build correctly