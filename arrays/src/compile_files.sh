#!/bin/bash

# remove old binary files and shared libraries
rm -f *.bin *.so *.mod ../build/*

# compile new files
gfortran -O3 -shared -fPIC -o int_deserialize.so int_deserialize_mod.f90
gfortran -O3 -shared -fPIC -o real_deserialize.so real_deserialize_mod.f90
gfortran -O3 -shared -fPIC -o char_deserialize.so char_deserialize_mod.f90
gfortran -O3 -shared -fPIC -o array_utils.so array_utils.f90
gfortran -O3 -shared -fPIC -o serialize.so serialize.f90

gfortran -O3 -shared -fPIC -o arrays.so serialize.f90 \
   int_deserialize_mod.f90 \
   real_deserialize_mod.f90 \
   char_deserialize_mod.f90 \
   array_utils.f90
   
# compile the test program
gfortran -O3 -o test_arrays test_arrays.f90 \
   int_deserialize_mod.f90 \
   real_deserialize_mod.f90 \
   char_deserialize_mod.f90 \
   serialize.f90 \
   array_utils.f90

# move libraries and modules to the build directory
mv *.so ../build
mv *.mod ../build