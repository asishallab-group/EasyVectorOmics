
# Tensor Omics Snippets

Organize and place the snippets inside the appropriate snippet folders according to their functionality:

- Use the `f42:` prefix for F42-compliant infrastructure.
- Use the `tox:` prefix for application-specific Tensor Omics routines.

---
## How to Add Snippets in VSCode

To install a snippet in Visual Studio Code:

1. Open the command palette (`Ctrl+Shift+P` or `Cmd+Shift+P`).
2. Search for:  
   `Preferences: Configure Snippets`
3. Choose a global file (e.g. `global.code-snippets`) or the language-specific file for Fortran (e.g. `fortran.json`).
4. Paste the snippet block inside the JSON file.

---

## Snippets

### `f42:do-parallel`

Start typing `f42:do-parallel` and let the snippet autocomplete.

It includes tab stops for easy customization:
- `${1:i}` → loop index variable
- `${2:N}` → name of the loop limit variable
- `${3:10}` → initial value assigned to `N`
- Loop body → the cursor jumps there after filling in the loop setup

The snippet expands into a unified parallel loop structure that automatically selects the appropriate backend based on compiler flags:
- Coarray
- OpenMP
- Serial fallback

This allows backend-independent parallel execution using preprocessor directives.

---

#### Compilation Examples

Use one of the following compilation commands depending on the desired backend:

##### Serial (default)

```bash
gfortran -cpp src/test_parallel.f90 -o build/test_parallel
```

##### OpenMP

```bash
gfortran -cpp -DUSE_OPENMP -fopenmp src/test_parallel.f90 -o build/test_parallel
```

##### Coarray (single-image)

```bash
gfortran -cpp -DUSE_COARRAY -fcoarray=single src/test_parallel.f90 -o build/test_parallel
```

> Ensure that you're using `-cpp` if your source files have the `.f90` extension. For `.F90`, preprocessing is usually automatic.

---

Feel free to extend this README with additional snippets or compilation configurations.
