#!/bin/bash

# remove old binary files and shared libraries
rm -f *.bin *.so *.mod ../build/*

gfortran -O3 -shared -fPIC -o arrays.so serialize_int.f90 \
   serialize_real.f90 \
   serialize_char.f90 \
   serialize.f90 \
   int_deserialize_mod.f90 \
   real_deserialize_mod.f90 \
   char_deserialize_mod.f90 \
   array_utils.f90
   
# compile the test program
gfortran -O3 -o test_arrays test_arrays.f90 \
   serialize_int.f90 \
   serialize_real.f90 \
   serialize_char.f90 \
   int_deserialize_mod.f90 \
   real_deserialize_mod.f90 \
   char_deserialize_mod.f90 \
   serialize.f90 \
   array_utils.f90

# move libraries and modules to the build directory
mv *.so ../build
mv *.mod ../build