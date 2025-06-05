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
implicit none
integer, parameter :: n = 1000
real(8) :: x(n)
integer :: ix(n), stack_left(64), stack_right(64)
! Fill x(:) with your data
call build_bst_index(x, n, ix, stack_left, stack_right)
! Access sorted value:
print *, get_sorted_value(x, ix, 1)
! Range query:
integer :: out_ix(n), out_n
call bst_range_query(x, ix, n, 0.2d0, 0.8d0, out_ix, out_n)
print *, 'Number in range:', out_n
```

### KD-Tree

```fortran
use kd_tree
implicit none
integer, parameter :: d = 3, n = 10
real(8) :: X(d, n)
integer :: kd_ix(n), dim_order(d), work(n), perm(n), stack_left(64), stack_right(64)
real(8) :: subarray(n)
! Fill X(:, :) with your data
dim_order = [1, 2, 3] ! or set by variance
call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, stack_left, stack_right)
! To get the i-th point in KD order:
real(8) :: point(d)
call get_kd_point(X, kd_ix, i, point)
```

#### Spherical KD-Tree

```fortran
use kd_tree
implicit none
integer, parameter :: d = 3, n = 10
real(8) :: V(d, n)
integer :: sphere_ix(n), dim_order(d), work(n), perm(n), stack_left(64), stack_right(64)
real(8) :: subarray(n)
! Fill V(:, :) with unit vectors
dim_order = [1, 2, 3]
call build_spherical_kd(V, d, n, sphere_ix, dim_order, work, subarray, perm, stack_left, stack_right)
```

### Sorting Utilities

```fortran
use tox_sorting
implicit none
real(8) :: arr(10)
integer :: perm(10), stack_left(64), stack_right(64)
! Fill arr(:) with your data
call sort_array(arr, perm, stack_left, stack_right)
```

## Testing

### Running the Tests

- **BST:**  
  Compile and run `src/binary_search/bst_test.f90` to validate the binary search tree implementation.
- **KD-Tree:**  
  Compile and run `src/k-d_tree/test_kd.f90` to validate the KD-Tree implementation, including a wide range of edge cases and random data for robustness.

Example (from project root):

```sh
cd src/
./compile_files.sh
cd ../build
./test_kd
./bst_test
```

### What is Tested

- **BST tests:**  
  - Sorted index monotonicity
  - Sorted value access
  - Range queries

- **KD-Tree tests:**  
  - 2D/3D/4D/5D/10D, various `n`, including edge cases (`n=0`, `n=1`, all identical points, sorted, reversed, all zeros, random data)
  - Both Cartesian and spherical data
  - Permutation checks for all outputs

### Example Output

All tests currently pass:

```
2D Cartesian KD Tree index:           1           5           2           4           3           6
 PASS:2D cartesian
 3D Spherical KD Tree index:           3           5           4           1           8           2           6           7
 PASS:3D spherical
 Edge n=0 passed.
 Edge n=1 passed.
 PASS:identical points
 PASS:unit vectors
 PASS:high-d low-n
 PASS:1D sorted
 PASS:2D n=2
 PASS:1D n=1
 PASS:1D n=2
 PASS:3D n=100
 PASS:5D n=10
 PASS:2D reversed
 PASS:2D all zeros
 All edge-case tests completed.
```

## Notes

- When calling KD-Tree routines, ensure all workspace arrays are allocated to at least the number of points (`n`), and stack arrays are large enough for the recursion depth (typically `64` is sufficient for most practical cases).
- For best results, set `dim_order` to the order of dimensions by variance (see `sort_dim_order` in the test file for an example).

## References

- See source files in `src/` for detailed documentation and comments.
- Sorting routines are in `src/tox_sorting.f90`.

---

For further details, consult the docstrings and comments in each module.