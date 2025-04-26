# Clean and rebuild
rm -f *.o *.mod *.so
lfortran -c TensorOmics_C_Module.f90
lfortran --shared -o libtensoromics_l.so TensorOmics_LFortran_mod.f90 TensorOmics_C_Module.o
gfortran -c TensorOmics_C_Module.f90
gfortran -shared -fPIC -o libtensoromics_g.so TensorOmics_GFortran_mod.f90 TensorOmics_C_Module.o