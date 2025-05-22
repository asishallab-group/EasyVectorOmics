#!/bin/bash
gfortran -c tox_sorting.f90
gfortran -c binary_search_tree.f90
gfortran -c bst_test.f90
gfortran -o ../build/bst_test tox_sorting.o binary_search_tree.o bst_test.o
