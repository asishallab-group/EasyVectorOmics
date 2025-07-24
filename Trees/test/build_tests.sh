#!/bin/bash
# Compile tox_sorting and BST
gfortran -c ../src/tox_sorting.f90
gfortran -c ../src/binary_search_tree.f90
gfortran -c bst_test.f90
gfortran -o test_bst tox_sorting.o binary_search_tree.o bst_test.o

# Compile k-d tree sources and test (cartesian and spherical)
gfortran -c ../src/k_d_tree.f90
gfortran -c test_kd.f90
gfortran -o test_kd tox_sorting.o k_d_tree.o test_kd.o

rm *.o *.mod