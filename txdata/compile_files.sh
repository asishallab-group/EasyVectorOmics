#!/bin/bash
lfortran -c TensorInterface_mod.f90
lfortran -c TensorOmics_mod.f90
lfortran -o test_tox test_tox.f90 TensorInterface_mod.o TensorOmics_mod.o
