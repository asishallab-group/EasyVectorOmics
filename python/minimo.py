import numpy as np
import ctypes

ctypes.CDLL("libgomp.so.1", mode=ctypes.RTLD_GLOBAL)
lib = ctypes.CDLL("build/libtensor-omics.so")

f = lib.loess_smooth_2d_ctest
f.argtypes = [
    ctypes.c_int, ctypes.c_int,
    np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),
    ctypes.c_double, ctypes.c_double,
    np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, flags='C_CONTIGUOUS'),
]
f.restype = None

n_total = 5
n_target = 2
x_ref = np.array([1,2,3,4,5], dtype=np.float64)
y_ref = np.array([10,20,30,40,50], dtype=np.float64)
indices_used = np.arange(1,6, dtype=np.int32)
x_query = np.array([2.5, 4.5], dtype=np.float64)
kernel_sigma = 1.0
kernel_cutoff = 2.0
y_out = np.zeros(n_target, dtype=np.float64)
workspace_weights = np.zeros(n_total, dtype=np.float64)
workspace_values = np.zeros(n_total, dtype=np.float64)
mask_in = np.ones(n_total, dtype=np.int32)

f(
    n_total, n_target,
    x_ref, y_ref, indices_used, x_query,
    kernel_sigma, kernel_cutoff,  # <-- plain floats!
    y_out, workspace_weights, workspace_values, mask_in
)
print('y_out:', y_out)