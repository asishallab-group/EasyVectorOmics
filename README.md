
# TensorOmics - Normalization and Fold Change Calculation

This repository provides the core normalization and transformation methods for **TensorOmics**, focusing on highly efficient matrix operations using **Fortran** for backend computation and **R** as a convenient wrapper.

This scripts were tested using:
- Ubuntu clang version 14.0.0-1ubuntu1.1
- R version 4.4.3 (2025-02-28) -- "Trophy Case"
- LFortran version: 0.51.0
- GNU Fortran (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0


## 🚀 Project Structure

- **Fortran 90/95 Modules (`tensoromics_normalization.f90`)**  
  Provides fast implementations of:
  - Standard deviation normalization
  - Quantile normalization
  - Log2(x+1) transformation
  - Tissue replicate averaging
  - Log2 fold-change calculation
  
- **R Interface (`normalization.R`, `tensoromics_functions.R`)**  
  User-friendly R functions that call the compiled Fortran subroutines and return easily manipulated R data frames.

## ⚙️ Why use `gfortran` instead of `lfortran`?

Although **LFortran** (based on LLVM) was originally preferred due to its WebAssembly (WASM) compatibility goals, it is currently **unable to correctly compile shared libraries** (`.so`) for direct dynamic loading in R.

Specifically:
- LFortran does not yet support full linking for `.so` binaries needed by `.C()` in R.
- Some `source` bindings (e.g., `bind(C, name="...")`) are not fully functional with `lfortran`'s current stable releases.
- `gfortran` provides stable `.so` generation fully compatible with R's foreign function interface (FFI).

👉 Therefore, **`gfortran` was used to compile** the shared library for now.  
If in the future LFortran improves `.so` support, migration can be considered.

## 🛠️ How to Build

Compile the Fortran code using `gfortran`:

```bash
gfortran -fPIC -shared -o build/libtensoromics_g.so src/tensoromics_normalization.f90
```

This will produce a shared object (`.so`) file that R can dynamically load.

Check if the subroutines were compiled successfully with:
```bash
nm -g libtensoromics_g.so
```


## 📦 How to Use in R

First, load the shared library inside your R session:

```r
dyn.load("libtensoromics.so")
```

Now you can use the provided R wrapper functions:

```r
# Normalize by standard deviation
normalized_matrix <- normalize_by_std_dev(input_matrix)

# Perform quantile normalization
quantile_normalized <- quantile_normalization(input_matrix)

# Apply log2(x + 1) transformation
log_matrix <- log2_transformation(input_matrix)

# Average replicates by tissue groups
averaged_matrix <- calculate_tissue_averages(input_matrix)

# Calculate fold changes
# Use control_pattern and condition_patterns to indicate which columns are control and conditions.
fc_matrix <- calculate_fc_by_patterns(averaged_matrix, control_pattern = "dietM", condition_patterns = c("dietP"))
```

Each function transparently sends the data to Fortran and reconstructs the result in R with proper row and column names.

## 📋 Documentation

All R functions are documented using **Roxygen2-style comments**, meaning that:
- Parameters, outputs, and examples are clearly explained.
- Easy to generate `.Rd` files for package building if needed.

The Fortran code is also carefully commented inline to facilitate further development.


# 📈 Future Work

- Investigate LFortran or Flang as a compiler once `.so` creation is stable.
- WASM target: fully port TensorOmics preprocessing to run inside browsers via WebAssembly.


# 🧹 Repository Structure Summary

```
/src
  └── tensoromics_normalization.f90    # Fortran backend code
/methods
  └── normalization.R                  # R wrappers for normalization functions
  └── tensoromics_functions.R           # Additional R helper functions
/results
  └── (results outputs here)
```

