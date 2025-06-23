import numpy as np
import ctypes
import pytest

# Cargar la librería compartida
lib = ctypes.CDLL("build/tox_normalization.so")
normalize = lib.normalize_by_std_dev_c
normalize.argtypes = [
    ctypes.c_int, ctypes.c_int,
    np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
]
normalize.restype = None

def call_normalize(mat):
    n_genes, n_tissues = mat.shape
    input_flat = np.asfortranarray(mat).ravel(order='F')
    output_flat = np.zeros_like(input_flat)
    normalize(n_genes, n_tissues, input_flat, output_flat)
    return output_flat.reshape((n_genes, n_tissues), order='F')

def test_normalizes_values_correctly():
    mat = np.array([[2, 6], [4, 8]], dtype=np.float64, order='F')
    result = call_normalize(mat)
    # Manual normalization
    std_dev = np.sqrt(np.mean(mat**2, axis=1, keepdims=True))
    expected = mat / std_dev
    assert result.shape == mat.shape
    np.testing.assert_allclose(result, expected, rtol=1e-12, atol=1e-12)

def test_handles_constant_rows():
    mat = np.array([[5, 5], [5, 5]], dtype=np.float64, order='F')
    result = call_normalize(mat)
    expected = np.ones_like(mat)
    assert result.shape == mat.shape
    assert np.all(np.isfinite(result))
    np.testing.assert_allclose(result, expected, rtol=1e-12, atol=1e-12)

def test_normalizes_large_numbers():
    mat = np.array([[1e6, 1e6], [2e6, 2e6]], dtype=np.float64, order='F')
    result = call_normalize(mat)
    std_dev = np.sqrt(np.mean(mat**2, axis=1, keepdims=True))
    expected = mat / std_dev
    assert result.shape == mat.shape
    assert np.all(np.isfinite(result))
    np.testing.assert_allclose(result, expected, rtol=1e-12, atol=1e-12)


# Load the shared library
quantile_norm = lib.quantile_normalization_c
quantile_norm.argtypes = [
    ctypes.c_int, ctypes.c_int,
    np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
    np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
    ctypes.c_int
]
quantile_norm.restype = None

def call_quantile_normalization(mat):
    n_genes, n_tissues = mat.shape
    input_flat = np.asfortranarray(mat).ravel(order='F')
    output_flat = np.zeros_like(input_flat)
    temp_col = np.zeros(n_genes, dtype=np.float64, order='F')
    rank_means = np.zeros(n_genes, dtype=np.float64, order='F')
    perm = np.zeros(n_genes, dtype=np.int32, order='F')
    max_stack = max(2 * n_genes, 2)
    stack_left = np.zeros(max_stack, dtype=np.int32, order='F')
    stack_right = np.zeros(max_stack, dtype=np.int32, order='F')
    quantile_norm(
        n_genes, n_tissues,
        input_flat, output_flat,
        temp_col, rank_means, perm, stack_left, stack_right, max_stack
    )
    return output_flat.reshape((n_genes, n_tissues), order='F')

def test_preserves_dimensions_and_names():
    mat = np.random.rand(2, 3)
    result = call_quantile_normalization(mat)
    assert result.shape == mat.shape

def test_handles_identical_rows_and_maintains_values():
    mat = np.full((2, 3), 5.0)
    result = call_quantile_normalization(mat)
    expected = np.full((2, 3), 5.0)
    assert result.shape == mat.shape
    assert np.all(np.isfinite(result))
    np.testing.assert_allclose(result, expected, rtol=1e-12, atol=1e-12)

def test_no_nans_and_standardizes_column_distributions():
    mat = np.array([[2, 3, 7], [0, 5, 1]], dtype=np.float64, order='F')
    result = call_quantile_normalization(mat)
    assert not np.any(np.isnan(result))
    assert np.all(np.isfinite(result))
    assert result.shape == mat.shape

    # Check that sorted columns are (approximately) equal
    sorted_cols = np.sort(result, axis=0)
    for i in range(1, sorted_cols.shape[1]):
        np.testing.assert_allclose(sorted_cols[:, i], sorted_cols[:, 0], rtol=1e-6, atol=1e-6)

if __name__ == "__main__":
    pytest.main([__file__])