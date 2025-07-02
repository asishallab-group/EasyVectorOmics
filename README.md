
# Tensor Omics

See the Tensor_Omics.tex file for details.

# TOX Project Structure

This repository contains the source code, methods, snippets and tests for the **Tensor Omics (TOX)** project.

## Folder Overview

```
/build
  └── ...       # Compiled Fortran binaries and intermediate build files

/doc
  └── ...       # Documentation generated automatically using FORD (Fortran documentation tool)

/python
  └── ...       # Python scripts that execute pipeline logic and invoke subroutines

/r
  └── ...       # R scripts that execute pipeline logic and invoke subroutines

/snippets
└── ...         #Code templates or reusable short logic blocks 

/src
  └── ...       # Fortran backend

/test
  └── ...       # Fortran testing

/helper
  └── helper_c_wrapper.py       #  helper script to generate c wrappers subroutines

build.sh        # Compile and generate shared libraries
ford.yml        # Generates documentation
fpm.toml        # Defines compilation options
test_runner.sh  # Compile and generate unit test

```

## Notes

* **`/build`** is used to store shared libraries, compiled and binary files resulting from Fortran compilation. It keeps the repo clean by separating source and compiled code.
* **`/doc`** contains the auto-generated documentation, which is built using [FORD](https://github.com/Fortran-FOSS-Programmers/ford) from annotated Fortran source files.
* **`/python`** includes python scripts that coordinate analysis workflows
* **`/r`** includes R scripts that coordinate analysis workflows
* **`/snippets/`** includes frequently used or testable units of logic reused across development stages. 
  - Snippets should be easy to create and use. The goal is to give the user access to the subroutine names along with their respective arguments, and nothing more. Example:
  ```
    {
        "Call to subroutine_name": {
          "prefix": "tox|f42:subroutine_name",
          "body": [
            "! Brief explanation on what the subroutine does",
            "call subroutine_name(arg1, arg2, arg3, arg4)"
          ],
          "description": "Insert a call to subroutine_name with brief explanation"
        }
    }
  ```
* **`/src`** contains performance-critical Fortran code. These are compiled during the build process.
  - All `.f90` files should include `precompiler_constants.f90`
* **`/test`** contains the unit tests for the Fortran subroutines.

  * The file `asserts.f90` must exist and can be modified if additional assert functions are needed.
  * There must be a central program called `run_tests.f90` which contains all the test calls defined in the modules.
  * Each subroutine's tests should be placed in independent modules (one file per tested subroutine). All of them must include a function `run_all_tests()` ... **TO BE COMPLETED**
  * All test modules must be named `mod_<subroutine_name>.f90` to ensure they are compiled before `run_tests.f90`. Otherwise, compilation errors may occur.
  * Full tests can be executed using `fpm test` or the script `test_runner.sh` ... **TO BE COMPLETED**

* **`/helper`** this folder will not be included in the final version of TOX. For now, it serves to help us create the C wrapper for the subroutines more quickly and easily. See details in `helper/readme.md`.

---

### FORD (Fortran Online Reference Documentation)

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


### Tensor Omics Snippets

Organize and place the snippets inside the appropriate snippet folders according to their functionality:

- Use the `f42:` prefix for F42-compliant infrastructure.
- Use the `tox:` prefix for application-specific Tensor Omics subroutines.

See `snippets/readme.md` for details.

---


### Compilation

The `build.sh` script will compile all the files located in the `src/` directory.

It creates a directory for the compiled objects under `/build/<compiler>/`, and the resulting shared library will be named `libtensor-omics.so`.

This `.so` file is the one that must be loaded from R or Python.

Every time the code is compiled, a new `/build/<compiler>/` directory is created. To simplify access, the script creates a symbolic link to the latest compiled shared library so that R and Python can always load the same file consistently.

Usage:

→ Uses the `gfortran` compiler without performance optimizations.

```bash
./build.sh
```

→ Uses the `gfortran` compiler with maximum performance flags.

```bash
./build.sh --max-performance
```

→ Uses the `ifx` compiler with maximum performance flags.

```bash
./build.sh --max-performance FC=ifx
```

---

Feel free to extend this README with additional information.
