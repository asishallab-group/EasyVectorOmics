---
Title: README for TensorOmics Data Structures
Author: Aaron Schroeder
Date: 24/04/2025
---

# 📦 TensorOmics Project Overview

This project defines data structures and interfaces for working with tensor-based gene expression data, particularly for multi-condition biological experiments. It adheres to the requirements set in Issue #2.

## 📁 File Structure

- **`TensorInterface_mod.f90`**  
  Defines the abstract `TensorInterface` type and its deferred procedures, serving as a contract for compatible tensor structures.

- **`TensorOmics_mod.f90`**  
  Implements the `TensorOmics_Type`, a concrete extension of `TensorInterface`, with features for memory management, patch updates, and binary saving.

- **`test_tox.f90`**  
  A minimal working test program that initializes a `TensorOmics_Type`, inserts a small data patch, and saves the result to a binary file.

---

## 🛠 Compilation Instructions

Use the following commands to compile with **LFortran**:

```bash
lfortran -c TensorInterface_mod.f90
lfortran -c TensorOmics_mod.f90
lfortran -o test_tox test_tox.f90 TensorInterface_mod.o TensorOmics_mod.o
./test_tox
```

If you encounter any linker or runtime issues, use:

```bash
lfortran --show-backtrace
```

---

## 🧪 Features Implemented

- ✅ Abstract interface for tensor data modules  
- ✅ Memory requirement estimation for varying gene counts  
- ✅ Initialization and zeroing of all arrays  
- ✅ Patch insertion with tracking of insertion indices  
- ✅ Binary file serialization (including metadata)

---

## ⚠ Known Limitations (LFortran-specific)

- Writing `character(len=...) :: meta_str(:,:)` via unformatted binary IO may not work in some LFortran builds. This section can be commented out or adapted if needed.
- Assumed-length strings (`character(len=*)`) in arguments may need fixed-length alternatives for full compatibility.
- Advanced polymorphism and dynamic allocation of class types are intentionally avoided to ensure support.
- K-D Trees are not supported in LFortran

## Known Issues:
- `open(newunit=unit, file=filename, form='unformatted', status='replace')` in `TensorOmics_mod.f90` line 93 is used, because LFortran currently does not support `access='stream'` in `open()` statements, which is required for writing raw binary data without formatting overhead. As a workaround, the code uses standard unformatted access (`access='sequential'`), which may insert record markers and affect binary compatibility with other readers. This limitation may impact downstream binary parsing or cross-language interoperability.
- Writing binary curently ends in corrupted files. Checking the `data.txdata` with `hexdump -c data.txdata` provides an overview of the file, which apparently includes massive amounts of space symbols and NULL values.
---

## 🔜 Open TODOs

- [ ] Implement more robust error handling in `tox_update`
- [ ] Validate binary save/load cycle with real datasets
- [ ] Extend test suite with multiple patches and edge cases
- [ ] Optionally replace string metadata serialization with a safer method


## 🧑‍💻 Notes for Developers

All modules are implemented with forward compatibility and clean modularity in mind. The interface design allows additional tensor-based types to be created and plugged in easily by extending `TensorInterface`.

The current implementation is tailored for compatibility with **LFortran**, the compiler required by the team. It avoids unsupported features and includes workarounds for known bugs or gaps in compiler support.