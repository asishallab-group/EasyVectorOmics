# Array Serialization Library for Fortran

## Overview

This library provides a comprehensive set of Fortran modules for serializing and deserializing arrays of various types (integer, real, character) and dimensions (1D to 5D) to/from a custom binary file format. The implementation includes support for variable-length character strings and efficient handling of multi-dimensional arrays.

## Key Features

- Support for integer (int32), real (real64), and character arrays
- Handling of arrays from 1 to 5 dimensions
- Custom binary format with metadata including:
  - Magic number for file identification
  - Type code (1=int32, 2=real64, 3=character)
  - Array dimensionality and shape
  - Character length information for string arrays
- Memory-efficient serialization/deserialization
- Comprehensive error checking
- C and R language bindings for interoperability

## File Format Specification

The binary file format consists of:
1. Magic number (4 bytes, 'FA20' in hex)
2. Type code (4 bytes)
3. Number of dimensions (4 bytes)
4. Dimension sizes (4 bytes each)
5. For character arrays: maximum string length (4 bytes)
6. Array data:
   - For character arrays: each element preceded by its actual length (4 bytes)
   - Other types: contiguous binary data

## Modules

### Core Modules
- `serialize_mod.f90`: Main interface module for serialization, provides generic interface
- `serialize_int.f90`: Integer array serialization
- `serialize_real.f90`: Real array serialization
- `serialize_char.f90`: Character array serialization
- `deserialize_int.f90`: Integer array deserialization
- `deserialize_real.f90`: Real array deserialization
- `deserialize_char.f90`: Character array deserialization
- `array_utils.f90`: Utility functions for array metadata

### Language Bindings
- C bindings for all serialization/deserialization functions
- R-compatible flat array interfaces

## Build Instructions

1. Ensure GNU Fortran compiler is installed
2. Run the build script:
   ```bash
   ./compile_files.sh
   ```

This will generate:
- Shared library (`arrays.so`)
- Test executable (`test_arrays`)
- Module files (`.mod`)

## Usage Examples

### Fortran
```fortran
use serialize_mod
integer(int32), allocatable :: arr(:,:)
real(real64), allocatable :: rarr(:,:,:)
character(len=5), allocatable :: carr(:,:)

! Serialization
call serialize(arr, "array.bin")
call serialize(rarr, "real_array.bin")
call serialize(carr, "char_array.bin")

! Deserialization
call deserialize_int(arr, "array.bin")
call deserialize_real(rarr, "real_array.bin")
call deserialize_char(carr, "char_array.bin")
```

### C
All functions are exposed to C using the `_C` functions in each module. 

### R
Since R can only pass/accepts flat arrays, the serialization and deserialization is not needed for multidimensional arrays. Both R and fortran use the same memory layout for arrays, therefore they can be passed without problems.

Chars and pointers can also not be passed cia `.Fortran()`, therefore the char array is converted to ASCII and then passed as integer and in fortran reconverted. TO allow reading the data, the array needs to be preallocated by R. To allow efficient memory assignment, the metadata is read first and the array is then created based on size and dimensions read.
## Testing

The test program `test_arrays` validates all functionality including:
- Basic serialization/deserialization for all types and dimensions
- Edge cases (empty arrays, 1×1 arrays)
- Variable-length character strings
- Multi-dimensional array handling