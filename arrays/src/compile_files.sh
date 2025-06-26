rm *.bin
gfortran -O3 -shared -fPIC -o  arrays.so arrays.f90
gfortran -O3 -o test_arrays test_arrays.f90 arrays.f90
mv *.so ../build
mv *.mod ../build