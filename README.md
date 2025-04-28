# Tensor Omics

Ubuntu clang version 14.0.0-1ubuntu1.1
R version 4.4.3 (2025-02-28) -- "Trophy Case"
LFortran version: 0.51.0
GNU Fortran (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0

gfortran -fPIC -shared -o build/libtensoromics_g.so src/tensoromics_normalization.f90
Rscript r/normalization.R 