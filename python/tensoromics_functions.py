import ctypes
import numpy as np
import os
import platform

# ---------------------------------------------------------------------
# Load Fortran library
# ---------------------------------------------------------------------

def _load_fortran_library():
    """Loads the shared Fortran library."""
    lib_name = "libtensor-omics.so"
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    lib_path = os.path.join(project_root, "build", lib_name)
    
    if not os.path.exists(lib_path):
        raise OSError(f"Could not find library at '{lib_path}'. Please run './build.sh'.")
    
    try:
        if platform.system() == "Linux":
            ctypes.CDLL("libgomp.so.1", mode=ctypes.RTLD_GLOBAL)
    except OSError:
        print("Warning: Could not preload libgomp.so.1.")
        
    return ctypes.CDLL(lib_path)

# ---------------------------------------------------------------------
# Wrapper function setup
# ---------------------------------------------------------------------

def _setup_function_signatures(lib):
    """
    Sets up the argument and return types for the Fortran C-interface functions.
    """
    func = lib.group_centroid_c
    func.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"), # vectors
        ctypes.c_int,                                                  # d
        ctypes.c_int,                                                  # n
        np.ctypeslib.ndpointer(dtype=np.int32, flags="F_CONTIGUOUS"),  # gene_to_family_map
        ctypes.c_int,                                                  # num_families
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"), # centroid_matrix (out)
        ctypes.c_int,                                                  # use_all_mode (as int)
        np.ctypeslib.ndpointer(dtype=np.int32, flags="F_CONTIGUOUS"),   # ortholog_set (as int array)
        np.ctypeslib.ndpointer(dtype=np.int32, flags="F_CONTIGUOUS"),  # selected_indices
        ctypes.c_int                                                   # selected_indices_len
    ]
    func.restype = None
    return func

# --- Module-level setup execution ---
lib = _load_fortran_library()
_group_centroid_c = _setup_function_signatures(lib)

# ---------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------

def _readonly(*arrays: np.ndarray) -> None:
    """Mark all given NumPy arrays as read-only."""
    for a in arrays:
        if isinstance(a, np.ndarray):
            a.flags.writeable = False

# ---------------------------------------------------------------------
# User-facing API function
# ---------------------------------------------------------------------

def tox_group_centroid(vectors, gene_to_family_map, num_families, ortholog_set, mode='all'):
    """
    Computes expression centroids for groups of genes.

    Parameters
    ----------
    vectors : np.ndarray
        A 2D NumPy array (d x n_genes) of gene expression vectors.
    gene_to_family_map : np.ndarray
        A 1D NumPy array of length n_genes, mapping each gene to a family ID.
    num_families : int
        The total number of unique families.
    ortholog_set : np.ndarray
        A 1D boolean NumPy array of length n_genes, indicating ortholog membership.
    mode : str, optional
        The calculation mode. 'all' (default) or 'ortho'.

    Returns
    -------
    np.ndarray
        A read-only (d x num_families) NumPy array containing the computed centroids.
    """
    # 1) Validate inputs
    if not isinstance(vectors, np.ndarray) or vectors.ndim != 2:
        raise ValueError("`vectors` must be a 2D NumPy array.")
    d, n_genes = vectors.shape
    if not isinstance(gene_to_family_map, np.ndarray) or gene_to_family_map.size != n_genes:
        raise ValueError("`gene_to_family_map` must be a 1D NumPy array of size n_genes.")
    if not isinstance(ortholog_set, np.ndarray) or ortholog_set.size != n_genes:
        raise ValueError("`ortholog_set` must be a 1D NumPy array of size n_genes.")
    if mode not in ['all', 'ortho']:
        raise ValueError("`mode` must be either 'all' or 'ortho'.")

    # 2) Prepare input/output buffers
    vecs_f = np.asarray(vectors, dtype=np.float64, order="F")
    g2f_map_f = np.asarray(gene_to_family_map, dtype=np.int32, order="F")
    
    # Convert boolean inputs to integers for the C interface
    use_all_mode_int = 1 if mode == 'all' else 0
    ortho_set_int_f = np.asarray(ortholog_set, dtype=np.int32, order="F")
    
    centroids_out = np.zeros((d, num_families), dtype=np.float64, order="F")
    selected_indices = np.zeros(n_genes, dtype=np.int32, order="F")

    # 3) Call Fortran
    _group_centroid_c(
        vecs_f,
        d,
        n_genes,
        g2f_map_f,
        num_families,
        centroids_out,
        use_all_mode_int,
        ortho_set_int_f,
        selected_indices,
        n_genes
    )

    # 4) Mark output as read-only and return
    _readonly(centroids_out)
    return centroids_out