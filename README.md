
# TensorOmics - Normalization and Fold Change Calculation

This repository provides the core normalization and transformation methods for **TensorOmics**, focusing on highly efficient matrix operations using **Fortran** for backend computation and **R** as a convenient wrapper.

The scripts were tested using:
- Ubuntu clang version 14.0.0-1ubuntu1.1
- R version 4.4.3 (2025-02-28) -- "Trophy Case"
- LFortran version: 0.51.0
- GNU Fortran (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0


## 🚀 Project Structure

- **Fortran 90/95 Modules (`tox_normalization.f90`)**  
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
gfortran -g -fcheck=all -Wall -fPIC -shared -o build/tox_normalization.so src/tox_sorting.f90 src/tox_normalization.f90
```

This will produce a shared object (`.so`) file that R can dynamically load.

Check if the subroutines were compiled successfully with:
```bash
nm -g build/tox_normalization.so
```

- Unit tests: sorting
  ```bash
  gfortran -c src/tests/test_sorting.f90 -o build/test_sorting.o
  gfortran -o build/test_runner src/tests/main.f90 build/test_sorting.o src/tox_sorting.f90
  ./build/test_runner 
  ```

## ⚡ Alternative Compilation Options (for future WASM or portability)

### 1. Compile only to Object File (`.o`)

Instead of building a shared object, you can compile to an intermediate object file:

```bash
gfortran -c tox_normalization.f90
```

- Generates `tox_normalization.o`.
- Useful if you want to later link manually or transform the `.o` file (e.g., for WASM).

---

### 2. Direct WebAssembly Compilation (Experimental)

Using **LFortran** (experimental feature), you can **directly compile to WebAssembly**:

```bash
lfortran --emit-wasm tox_normalization.f90 -o tox_normalization.wasm
```

- This produces a `.wasm` module directly from Fortran.
- However, **current LFortran WebAssembly output is experimental** and may not yet support all features (like bindings with R or complex memory layouts).
- Migration to pure WebAssembly workflows is considered a future goal.


## 📦 How to Use in R

First, load the shared library inside your R session:

```r
dyn.load("build/tox_normalization.so")
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

The Fortran code is also carefully commented using FORD:


### ✅ FORD (Fortran Online Reference Documentation)

[FORD](https://github.com/Fortran-FOSS-Programmers/ford) is a documentation generator specifically designed for Fortran projects. It allows developers to create clean, structured, and navigable HTML documentation from source code using lightweight markup embedded in comments.

* Designed specifically for Fortran (unlike Doxygen which is general-purpose).
* Supports documentation of modules, subroutines, functions, derived types, and more.
* Uses `!!!` or `!>` comment syntax to annotate code.
* Ideal for scientific and engineering projects using modern Fortran.
* Easy to integrate into Git-based workflows.

Example usage:

```bash
ford ford.yml
```

This generates an HTML site you can explore in a browser (`doc/index.html` by default).

# 📈 Future Work

- [ ] Investigate full WebAssembly compilation and browser execution.
- [ ] Migrate to LFortran when .so and .wasm generation becomes fully stable.

# 🧹 Repository Structure Summary

```
/src
  └── tox_normalization.f90    # Fortran backend code
/methods
  └── normalization.R                  # R wrappers for normalization functions
  └── tensoromics_functions.R           # Additional R helper functions
/results
  └── (results outputs here)
```



# 📄 Input Data Format Specification

This project expects input gene expression matrices to follow specific naming conventions for proper normalization and processing.


## 📦 General Input Matrix Structure

- **Rows**: Genes (e.g., `geneA`, `geneB`, `FBpp0070000`, etc.).
- **Columns**: Samples corresponding to tissues and/or experimental groups.


## 📋 Column Naming Rules

| Feature | How it should be named |
|:--------|:-----------------------|
| Tissue name | Should appear as a prefix |
| Control vs condition groups | Distinguished using suffixes |
| Replicates | Indicated with final `_1`, `_2`, `_3`, etc. |

### ✅ Examples

| Column Name | Meaning |
|:------------|:--------|
| `muscle_control_1` | Muscle tissue, control group, replicate 1 |
| `muscle_control_2` | Muscle tissue, control group, replicate 2 |
| `muscle_treatmentA_1` | Muscle tissue, condition "treatmentA", replicate 1 |
| `brain_control` | Brain tissue, control (no replicate if only one) |
| `heart_dietM_1` | Heart tissue, dietM condition, replicate 1 |
| `heart_dietP_2` | Heart tissue, dietP condition, replicate 2 |


## 🧠 Important Points

- **Replicates**:  
  - Must be indicated with an underscore and number (`_1`, `_2`, etc.).
  - Replicates are automatically detected and averaged using `calculate_tissue_averages()`.

- **Conditions (control vs experimental)**:  
  - Control groups and condition groups must be distinguishable by a recognizable **pattern** (e.g., `"control"`, `"dietM"`, `"dietP"`, etc.).
  - Fold-change is computed based on **pattern matching** using functions like `calculate_fc_by_patterns()`.

- **Missing Replicates**:  
  - If no replicate number exists, it's assumed to be a single sample for that tissue/condition.

- **Naming Consistency**:  
  - Always keep naming consistent (same format across all columns).
  - No spaces or unusual characters; use underscores `_` to separate parts.


## ⚙️ Supported Input Scenarios

| Scenario | Supported? | Notes |
|:---------|:-----------|:------|
| Single sample per tissue | ✅ | Simply use the tissue name (e.g., `brain_control`) |
| Multiple replicates per tissue group | ✅ | Must use `_1`, `_2`, `_3` convention |
| Multiple conditions for same tissue | ✅ | Control and conditions detected based on user-provided patterns |
| No replicates and no conditions | ✅ | Treated as simple tissues for averaging or fold-change |


## 🚫 Not Supported

- Spaces inside column names (`"muscle control"` ❌).
- Different separator characters (`"muscle-control-1"` ❌ use `"muscle_control_1"` ✅).
- Completely inconsistent naming.
