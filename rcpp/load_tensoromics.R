# TensorOmics Rcpp Function Loader
#
# This script sets up the compilation environment and loads all TensorOmics
# Rcpp wrapper functions that interface with the underlying Fortran library.
#
# What this script does:
# 1. Loads the Rcpp library for C integration
# 2. Sets up the library path to find the compiled Fortran library (libtensor-omics.so)
# 3. Configures PKG_LIBS environment variable with:
#    - Runtime library path (-Wl,-rpath) for finding shared libraries at runtime
#    - Library search path (-L) for compilation linking
#    - Links to tensor-omics library (-ltensor-omics)
#    - Links to Fortran runtime (-lgfortran)
# 4. Compiles and loads the Rcpp functions from tensoromics_functions.cpp
# 5. Makes all tox_* functions available in the global environment
#
# Usage:
#   source("rcpp/load_tensoromics.R")
#
# After sourcing, all TensorOmics Rcpp functions will be available.

library(Rcpp)

# Get absolute path to build directory containing the compiled Fortran library
lib_path <- normalizePath("build")

# Set up compilation flags for linking with Fortran library
Sys.setenv(PKG_LIBS = paste0("-Wl,-rpath,", lib_path, " -L", lib_path, " -ltensor-omics -lgfortran"))

# Compile and load all TensorOmics Rcpp wrapper functions (includes error_handling.cpp)
sourceCpp("rcpp/tensoromics_functions.cpp", env = .GlobalEnv)

cat("✓ TensorOmics Rcpp functions loaded successfully\n")