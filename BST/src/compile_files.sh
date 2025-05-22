#!/bin/bash
gfortran -c tox_sorting.f90
gfortran -c binary_search/binary_search_tree.f90
gfortran -c binary_search/bst_test.f90
gfortran -o ../build/bst_test tox_sorting.o binary_search_tree.o bst_test.o

# Compile k-d tree sources and test
gfortran -c k-d_tree/k_d_tree.f90
gfortran -c k-d_tree/test_kd.f90
gfortran -o ../build/test_kd tox_sorting.o k_d_tree.o test_kd.o

mv *.o ../build
mv *.mod ../build
