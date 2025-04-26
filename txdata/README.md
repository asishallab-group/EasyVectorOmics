
---
Title: README for TensorOmics Data Structures
Author: Aaron Schroeder
Date: 25/04/2025
---

# 📦 TensorOmics Project Overview

This project defines data structures and interfaces for working with tensor-based gene expression data, particularly for multi-condition biological experiments. The `TensorOmics_Type` is designed to manage and process large-scale datasets, handling operations such as memory management, patch updates, and binary file serialization.

## 📁 File Structure

- **`TensorOmics_LFortran_mod.f90`**  
  Implements the `TensorOmics_Type` in **LFortran**. This version of the module is designed to be compiled with the LFortran compiler and provides basic tensor operations, including memory estimation and patch insertion.

- **`TensorOmics_GFortran_mod.f90`**  
  Implements the `TensorOmics_Type` for **GFortran**. This module is compiled with the GFortran compiler and provides additional functionality, including saving the tensor data to binary files.

- **`TensorOmics_C_Module.f90`**  
  A module designed to handle C-compatible data structures for easier integration with C-based systems. This is used to facilitate interoperability between the Fortran modules and the R wrapper.

- **`RWrapper.R`**  
  The R wrapper file that provides the interface for R to interact with the compiled Fortran libraries (`libtensoromics_l.so` and `libtensoromics_g.so`). The wrapper defines R functions that correspond to the Fortran subroutines.

---

## 🛠 Compilation Instructions

### 1. **Compiling with LFortran**

LFortran is used to compile the core tensor operations in `TensorOmics_LFortran_mod.f90`. This module is compiled into a shared library:

```bash
lfortran -shared -o libtensoromics_l.so TensorOmics_LFortran_mod.f90
```

### 2. **Compiling with GFortran**

The GFortran version (`TensorOmics_GFortran_mod.f90`) is compiled into a shared library for saving tensor data to binary files:

```bash
gfortran -shared -fPIC -o libtensoromics_g.so TensorOmics_GFortran_mod.f90
```

### 3. **Compilation with Both Compilers**

Both the LFortran and GFortran versions need to be compiled and linked separately. Make sure that `libtensoromics_l.so` and `libtensoromics_g.so` are located in the appropriate directory for your system.

### 4. **Using the R Wrapper**

The R wrapper provides a way to interact with the Fortran subroutines. You can load the compiled shared libraries and call the corresponding functions from R. 

The R wrapper uses `.Fortran()` to call the subroutines from the Fortran code directly. Here's an example of how to call the functions:

```r
# Load shared libraries
dyn.load("./libtensoromics_l.so")
dyn.load("./libtensoromics_g.so")

# Initialize tensoromics object
tom <- init_tensoromics(data_estimates)

# Calculate memory requirements
calculate_memory(tom)

# Update tensoromics with a patch
patch_matrix <- matrix(rnorm(100), nrow = tom$n_conditions)
update_tensoromics(tom, patch_matrix)

# Save tensoromics data
save_tensoromics(tom, "tensor_data.bin")
```

---

## 🧪 Features Implemented

- ✅ **Abstract Interface**: A base interface for tensor data structures, making it easy to extend for new tensor types.
- ✅ **Memory Estimation**: Calculates the memory required for a given tensor with specified conditions and genes.
- ✅ **Initialization and Zeroing**: Initializes the tensor data, setting all values to zero.
- ✅ **Patch Insertion**: Updates the tensor with patches of data, keeping track of insertion indices.
- ✅ **Binary Serialization**: Saves the tensor data to a binary file for future use.

---

## ⚠ Known Limitations and Issues

### **LFortran-Specific Issues**
- **Unformatted Binary I/O**: Writing `character(len=...)` variables through unformatted binary I/O may not work in some builds of LFortran. Workarounds are included for this.
- **Access Mode Restrictions**: LFortran currently does not support `access='stream'` in the `open()` function for writing raw binary data. This affects writing files without formatting overhead. The current implementation uses `access='sequential'`, which may cause issues with binary compatibility, especially for cross-language file parsing.
- **Limited Polymorphism**: LFortran does not support advanced polymorphism and dynamic allocation of class types, so these features were avoided in this implementation.
- **K-D Tree Support**: K-D trees and some other advanced data structures are not supported in LFortran.

### **Incompatibilities During Compilation**
- **C-Fortran Interoperability**: The original design attempted to expose the Fortran subroutines with C interop (`bind(c)`), but this caused several issues:
  - **Problems with `.C()` Calls**: The `.C()` interface in R was problematic due to incompatibilities with how LFortran and GFortran handle memory and data types.
  - **Workaround**: To resolve this, we will probably switch from using `.C()` to `.Fortran()`, which is a more natural interface for interacting with Fortran subroutines from R. A detailed explanaition is found at the end.

### **R Wrapper and Function Calls**
- The R wrapper will use `.Fortran()` for calling Fortran subroutines, avoiding the complexity and bugs associated with `.C()`.
- This change is necessary because `.C()` requires handling **C interop**, which leads to multiple issues, including problems with pointer passing, memory management, and data conversion between C and Fortran.

### **GFortran-Specific Issues**
- The GFortran code is used for saving the tensor data to a binary file, but due to differences in how GFortran and LFortran handle memory and type conversions, there was a need to split the functionality across two compilers. This ensures that both memory management and file I/O can be handled in the most efficient way possible for each compiler.

---

## 🔜 Open TODOs

- [ ] **Implement More Robust Error Handling**: The current `tox_update` and other subroutines need better error handling, especially for edge cases.
- [ ] **Validate Binary Save/Load Cycle**: Test the binary save/load functionality with real datasets to ensure cross-compatibility and data integrity.
- [ ] **Extend Test Suite**: Add tests for more patch insertion scenarios, edge cases, and error scenarios to ensure the code handles diverse data configurations.
- [ ] **Refactor String Metadata Handling**: Consider replacing unsafe string serialization methods (e.g., `character(len=*)`) with more robust alternatives that work across all compilers.
- [ ] **.Fortran() call**: Switch from using `.C()` to `.Fortran()`

---

## 🧑‍💻 Notes for Developers

This project is designed with forward compatibility in mind. All modules are implemented with clean modularity, allowing easy extensions or changes to the tensor data types.

The current implementation uses both **LFortran** and **GFortran** to handle different parts of the system. LFortran is used for its modern Fortran support and ease of use, while GFortran is employed for specific functionality like file saving due to its more mature support for such features.

The main challenge faced during development was **C-Fortran Interoperability**, which led to issues when using `.C()` in the R wrapper. To resolve these, the code was refactored to use `.Fortran()`, which simplifies interaction with Fortran code and avoids complex C interoperability.

### **Compilation Issues**
- LFortran and GFortran have some incompatible behavior, particularly with regard to pointer handling, memory management, and C interoperability. These issues necessitated the use of two separate compilers to handle different tasks.
- Some Fortran subroutines had to be modified to ensure they worked across both compilers, and additional care was taken to ensure compatibility between the shared libraries generated by the two compilers.

---

### **Why `.C()` Was Initially Used (and why `.Fortran()` is the better choice)**

The `.C()` function was originally used in an attempt to interface R with Fortran by taking advantage of C-compatible Fortran subroutines (using `bind(c)`), under the assumption that R's `.C()` function would be the easiest way to call Fortran subroutines.

However, as we didn’t need to interoperate with C and were only aiming for R-Fortran communication, it turned out that this approach was not optimal. Specifically:

1. **Unnecessary C Interoperability**: R’s `.C()` function is designed to interface with C libraries and C-style Fortran subroutines (`bind(c)`). However, since we were not using C, we don’t need to go through C bindings at all. Our only goal was to directly call Fortran code from R.
   
2. **Issues with `.C()` Functionality**: Using `.C()` with Fortran subroutines that were `bind(c)` introduced additional complications:
   - R required specific handling of pointers, arrays, and data types, leading to inconsistent behavior when passing more complex structures like `TensorOmics_Type`.
   - There were also compatibility issues between the data structures expected by `.C()` and those used in Fortran, which led to errors when calling the subroutines.

Given these challenges, it became clear that **`.Fortran()`** would be a more appropriate choice because:
- **R-Fortran Direct Communication**: `.Fortran()` is tailored to work directly with Fortran subroutines without needing C compatibility. This simplifies the interaction and avoids the extra complexity of C bindings.
- **Avoiding Pointer/Memory Issues**: Using `.Fortran()` allows R to work with Fortran arrays and types more naturally, which reduces the need for additional memory management and pointer handling.
