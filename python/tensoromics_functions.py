"""
TensorOmics Functions Module
Python wrapper functions for Fortran routines via C interface
"""

import numpy as np
import ctypes
import os

# Load library
dll_path = os.path.abspath("build/libtensor-omics.so")
ctypes.CDLL("libgomp.so.1", mode=ctypes.RTLD_GLOBAL)
lib = ctypes.CDLL(dll_path)


def tox_errors(ierr):
    """
    Check error code and raise informative exception if needed
    
    Args:
        ierr: Error code from Fortran routine
    
    Raises:
        RuntimeError: If error code indicates failure
    """
    error_messages = {
        0: None,  # Success
        101: "File could not be opened",
        102: "Could not read magic number",
        103: "Could not read array type code", 
        104: "Could not read array dimension number",
        105: "Could not read array dimensions",
        106: "Could not read character length",
        107: "Could not read array data",
        200: "Invalid file format (magic number mismatch)",
        202: "No axes selected (empty input)",
        5002: "File not open or unit not connected",
        9999: "Unknown error"
    }
    
    msg = error_messages.get(ierr, f"Unknown Fortran error code: {ierr}")
    
    if msg is not None:
        raise RuntimeError(msg)

def _readonly(*arrays: np.ndarray) -> None:
    """Mark all given NumPy arrays as read-only."""
    for a in arrays:
        if isinstance(a, np.ndarray):
            a.flags.writeable = False
            # NOTE: Returned NumPy arrays are read-only for safety.
            # If you need to modify them (e.g., for plotting), use `.copy()`.


def tox_clock_hand_angle_between_vectors(v1, v2, selected_axes_for_signed):
    """
    Calculate clock hand angle between two vectors
    
    Args:
        v1: First vector (numpy array)
        v2: Second vector (numpy array) 
        selected_axes_for_signed: Boolean array indicating which axes to use for signed angle
    
    Returns:
        float: Signed angle between vectors in degrees
    """
    # Input validation and conversion
    v1 = np.ascontiguousarray(v1, dtype=np.float64)  # 1D vectors can stay C-contiguous
    v2 = np.ascontiguousarray(v2, dtype=np.float64)  # 1D vectors can stay C-contiguous
    selected_axes_for_signed = np.ascontiguousarray(selected_axes_for_signed, dtype=np.int32)
    
    n_dims = len(v1)
    
    if len(v2) != n_dims:
        raise ValueError("v1 and v2 must have same length")
    if len(selected_axes_for_signed) != n_dims:
        raise ValueError("selected_axes_for_signed must have same length as vectors")
    
    # Prepare output
    signed_angle = ctypes.c_double(0.0)
    
    # Setup C wrapper
    clock_hand_angle = lib.clock_hand_angle_between_vectors_c
    clock_hand_angle.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # v1
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # v2
        ctypes.c_int,  # n_dims
        ctypes.POINTER(ctypes.c_double),  # signed_angle
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS")     # selected_axes_for_signed
    ]
    clock_hand_angle.restype = None
    
    # Call Fortran routine
    clock_hand_angle(v1, v2, n_dims, ctypes.byref(signed_angle), selected_axes_for_signed)
    
    return signed_angle.value

def tox_clock_hand_angles_for_shift_vectors(origins, targets, vecs_selection_mask, selected_axes_for_signed):
    """
    Calculate clock hand angles for shift vectors
    
    Args:
        origins: Origin vectors (n_dims x n_vecs)
        targets: Target vectors (n_dims x n_vecs)
        vecs_selection_mask: Boolean array indicating which vectors to process
        selected_axes_for_signed: Boolean array indicating which axes to use for signed angle
    
    Returns:
        numpy.ndarray: Signed angles for selected vectors in degrees
    """
    # Input validation and conversion - PRESERVE FORTRAN ORDER
    origins = np.asfortranarray(origins, dtype=np.float64)  # Use Fortran order
    targets = np.asfortranarray(targets, dtype=np.float64)  # Use Fortran order
    vecs_selection_mask = np.ascontiguousarray(vecs_selection_mask, dtype=np.int32)
    selected_axes_for_signed = np.ascontiguousarray(selected_axes_for_signed, dtype=np.int32)
    
    n_dims, n_vecs = origins.shape
    
    if targets.shape != (n_dims, n_vecs):
        raise ValueError("origins and targets must have same shape")
    if len(vecs_selection_mask) != n_vecs:
        raise ValueError("vecs_selection_mask must match number of vectors")
    if len(selected_axes_for_signed) != n_dims:
        raise ValueError("selected_axes_for_signed must match number of dimensions")
    
    n_selected_vecs = int(np.sum(vecs_selection_mask))
    
    # Prepare output
    signed_angles = np.zeros(n_selected_vecs, dtype=np.float64)
    
    # Setup C wrapper - CHANGE FLAGS FOR FORTRAN ORDER
    clock_hand_angles = lib.clock_hand_angles_for_shift_vectors_c
    clock_hand_angles.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"),  # origins (Fortran order)
        np.ctypeslib.ndpointer(dtype=np.float64, flags="F_CONTIGUOUS"),  # targets (Fortran order)
        ctypes.c_int,  # n_dims
        ctypes.c_int,  # n_vecs
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # vecs_selection_mask
        ctypes.c_int,  # n_selected_vecs
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # selected_axes_for_signed
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")   # signed_angles
    ]
    clock_hand_angles.restype = None
    
    # Call Fortran routine
    clock_hand_angles(origins, targets, n_dims, n_vecs, vecs_selection_mask, 
                     n_selected_vecs, selected_axes_for_signed, signed_angles)
    
    # Mark output as read-only
    _readonly(signed_angles)
    
    return signed_angles

def relative_axes_changes_from_shift_vector(shift_vector):
    """
    Compute relative axis contributions from a shift vector (RAP space).
    Args:
        shift_vector (array-like): Input vector (1D)
    Returns:
        np.ndarray: Relative axis contributions (sum to 1)
    """
    vec = np.ascontiguousarray(shift_vector, dtype=np.float64)
    n_dims = len(vec)
    contrib = np.zeros(n_dims, dtype=np.float64)
    # Setup C wrapper
    relative_axes_changes = lib.relative_axes_changes_from_shift_vector_c
    relative_axes_changes.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # shift_vector
        ctypes.c_int,  # n_dims
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")   # contrib
    ]
    relative_axes_changes.restype = None
    relative_axes_changes(vec, n_dims, contrib)
    _readonly(contrib)
    return contrib

def relative_axes_expression_from_expression_vector(expression_vector):
    """
    Compute relative axis contributions from an expression vector (RAP space).
    Args:
        expression_vector (array-like): Input vector (1D)
    Returns:
        np.ndarray: Relative axis contributions (sum to 1)
    """
    vec = np.ascontiguousarray(expression_vector, dtype=np.float64)
    n_dims = len(vec)
    contrib = np.zeros(n_dims, dtype=np.float64)
    # Setup C wrapper
    relative_axes_changes = lib.relative_axes_expression_from_expression_vector_c
    relative_axes_changes.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # expression_vector
        ctypes.c_int,  # n_dims
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")   # contrib
    ]
    relative_axes_changes.restype = None
    relative_axes_changes(vec, n_dims, contrib)
    _readonly(contrib)
    return contrib