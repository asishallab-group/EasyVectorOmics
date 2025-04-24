---
Title: README for TensorOmics data structures
Author: Aaron Schroeder
Date: 24/04/2025
---

## 📦 Overview

This repository provides a minimal and extensible infrastructure for working with multi-dimensional gene expression data structures in Fortran. The project consists of a base interface for tensor-like data and a full implementation that supports expression data vectors, shift vectors, and metadata. Placeholder support for K-D trees is also included.


## 📁 File Structure

- **`TensorInterface_mod.f90`**  
  Defines the **abstract interface** `TensorInterface`, which includes the required method signatures:
  - `calculate_memory_requirements`
  - `init`
  - `update`
  - `save`

- **`TensorOmics_mod.f90`**  
  Implements the `TensorOmics_Type`, which **extends** the interface. This module includes:
  - Memory estimation and allocation
  - Data insertion and shift vector computation
  - Simple metadata handling
  - Save/load functionality
  - Placeholder K-D tree structures 
  
- **`test_tensoromics.f90`**  
  A basic program to test functionality, including:
  - Memory calculation and allocation
  - Data patch insertion
  - Save and build-tree operations


## ⚙️ Compilation Instructions

To compile the modules:

```bash
lfortran -c TensorInterface_mod.f90
lfortran -c TensorOmics_mod.f90
```

> **Note**: `test_tensoromics.f90` currently **does not compile successfully in LFortran** due to limited support for some advanced features (e.g., abstract types and polymorphism-related behaviors). You may compile this test program using **GFortran** instead:

```bash
gfortran -o test_tensoromics TensorInterface_mod.f90 TensorOmics_mod.f90 test_tensoromics.f90
```

## ❗ Known Issues & Workarounds

### ❌ LFortran Limitations
- **Abstract types and derived-type polymorphism** are not yet fully supported in LFortran.
- **K-D Trees** are not supported in LFortran.
- Serialization of derived types with `write(unit) derived_type` may not work reliably.

### ✅ Workaround
- Compile the test program using **GFortran**, which has more complete Fortran support.


## 🔄 Cross-Compiler Compatibility

- It is technically possible to **mix GFortran and LFortran compiled modules**, but:
  - Interface compatibility must be ensured (data layout, types, ABI).
  - Prefer simple I/O formats like ASCII or binary with explicit layout.
  - Use a common linker and ensure runtime libraries are available (e.g., `-llfortran_runtime` for LFortran objects).

## 📌 Future ToDos
- Bug fixes
- More test cases
- Getting the test program to compile
- Extend functions

## ✅ Status Summary

| Component             | LFortran | GFortran |
|----------------------|----------|----------|
| TensorInterface_mod  | ✅       | ✅       |
| TensorOmics_mod      | ✅       | ✅       |
| test_tensoromics     | ❌       | ✅       |
