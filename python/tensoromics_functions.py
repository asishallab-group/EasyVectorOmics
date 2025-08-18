import ctypes
import numpy as np
import os
import platform

# ---------------------------------------------------------------------
# Load Fortran library
# ---------------------------------------------------------------------

def _load_fortran_library():
    """Dynamically loads the shared Fortran library."""
    lib_name = "libtensor-omics.so"
    # Assume this file is in python/, so go up one level to the project root
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    lib_path = os.path.join(project_root, "build", lib_name)
    
    if not os.path.exists(lib_path):
        raise OSError(f"Could not find library at '{lib_path}'. Please run './build.sh'.")
    
    try:
        # Preload libgomp for OpenMP support on Linux if available
        if platform.system() == "Linux":
            ctypes.CDLL("libgomp.so.1", mode=ctypes.RTLD_GLOBAL)
    except OSError:
        print("Warning: Could not preload libgomp.so.1.")
        
    return ctypes.CDLL(lib_path)

lib = _load_fortran_library()

# ---------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------

def _readonly(*arrays: np.ndarray) -> None:
    """Mark all given NumPy arrays as read-only."""
    for a in arrays:
        if isinstance(a, np.ndarray):
            a.flags.writeable = False

# ---------------------------------------------------------------------
# Wrapper function setup
# ---------------------------------------------------------------------

# Setup for group_centroid_c, done once on module load.
_group_centroid_c = lib.group_centroid_c
_group_centroid_c.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"), # vectors
    ctypes.c_int,                                                  # d
    ctypes.c_int,                                                  # n
    np.ctypeslib.ndpointer(dtype=np.int32, flags="F_CONTIGUOUS"),  # gene_to_family_map
    ctypes.c_int,                                                  # num_families
    np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"), # centroid_matrix (out)
    ctypes.c_bool,                                                 # use_all_mode
    np.ctypeslib.ndpointer(dtype=np.bool_, flags="F_CONTIGUOUS"),   # ortholog_set
]
_group_centroid_c.restype = None

# ---------------------------------------------------------------------
# User-facing API function
# ---------------------------------------------------------------------

def tox_group_centroid(vectors, gene_to_family_map, num_families, ortholog_set, mode='all'):
    """
    Computes expression centroids for groups of genes.

    Parameters
    ----------
    vectors : np.ndarray
        A 2D NumPy array (d x n_genes) of gene expression vectors,
        preferably in Fortran order.
    gene_to_family_map : np.ndarray
        A 1D NumPy array of length n_genes, mapping each gene to a family ID.
    num_families : int
        The total number of unique families.
    ortholog_set : np.ndarray
        A 1D boolean NumPy array of length n_genes, indicating ortholog membership.
    mode : str, optional
        The calculation mode. 'all' (default) uses all genes in a family.
        'ortho' uses only genes where `ortholog_set` is True.

    Returns
    -------
    np.ndarray
        A read-only (d x num_families) NumPy array containing the computed centroids.
        NOTE: Returned NumPy arrays are read-only for safety.
        If you need to modify them (e.g., for plotting), use `.copy()`.
    """
    # 1) Validate inputs
    if not isinstance(vectors, np.ndarray) or vectors.ndim != 2:
        raise ValueError("`vectors` must be a 2D NumPy array.")
    if not isinstance(gene_to_family_map, np.ndarray) or gene_to_family_map.ndim != 1:
        raise ValueError("`gene_to_family_map` must be a 1D NumPy array.")
    if not isinstance(ortholog_set, np.ndarray) or ortholog_set.ndim != 1:
        raise ValueError("`ortholog_set` must be a 1D NumPy array.")
    
    d, n_genes = vectors.shape
    if gene_to_family_map.size != n_genes or ortholog_set.size != n_genes:
        raise ValueError("Input array dimensions are inconsistent.")
        
    if mode not in ['all', 'ortho']:
        raise ValueError("`mode` must be either 'all' or 'ortho'.")

    # 2) Prepare input/output buffers
    vecs_f = np.asarray(vectors, dtype=np.float64, order="F")
    g2f_map_f = np.asarray(gene_to_family_map, dtype=np.int32, order="F")
    ortho_set_f = np.asarray(ortholog_set, dtype=np.bool_, order="F")
    
    use_all_mode = (mode == 'all')
    
    centroids_out = np.zeros((d, num_families), dtype=np.float64, order="F")

    # 3) Call Fortran (argtypes already set)
    _group_centroid_c(
        vecs_f,
        d,
        n_genes,
        g2f_map_f,
        num_families,
        centroids_out,
        use_all_mode,
        ortho_set_f
    )

    # 4) Mark output as read-only and return
    _readonly(centroids_out)
    return centroids_out
