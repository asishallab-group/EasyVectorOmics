import ctypes
import numpy as np
import os
from error_handling import check_err_code

# load fortran library
lib_path = os.path.join(os.path.dirname(__file__), "../build/libtensor-omics.so")
arrays_lib = ctypes.CDLL(lib_path)

# helper to convert a filename to ASCII chars to transfer it as integer to fortran
def _filename_to_ascii_array(filename):
    ascii_arr = np.array([ord(c) for c in filename], dtype=np.int32)
    return ascii_arr, np.int32(len(ascii_arr))

# Function to convert array of strings to integer ASCII array
def _string_array_to_ascii_matrix(strings: np.ndarray) -> tuple[np.ndarray, int]:
    """
    Transforms a char array to integer
    """
    if not isinstance(strings, np.ndarray) or strings.dtype.kind != 'U':
        raise ValueError("Input must be a numpy array of strings (dtype='U')")
    
    flat = strings.flatten(order='F')  # important: Fortran-order
    clen = max(len(s) for s in flat)
    total = flat.size

    mat = np.zeros((clen, total), dtype=np.int32, order='F')
    for i, s in enumerate(flat):
        codes = [ord(c) for c in s]
        mat[:len(codes), i] = codes

    return mat, clen

# Helper function to read dimensions of integer/real array
def tox_get_array_metadata(filename, max_dims=5, with_clen=False):
    """
    Reads dimensions (and optionally character length) of an array file.
    with_clen=True -> returns (dims, clen)
    with_clen=False -> returns dims only
    """
    filename_ascii, fn_len = _filename_to_ascii_array(filename)

    dims_out = np.zeros(max_dims, dtype=np.int32)
    ndims = ctypes.c_int()
    ierr = ctypes.c_int()
    clen = ctypes.c_int()  # immer vorhanden

    # Gemeinsamer Prototyp
    arrays_lib.get_array_metadata_C.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"), # filename_ascii
        ctypes.c_int,                                                         # fn_len
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"), # dims_out
        ctypes.POINTER(ctypes.c_int),                                         # ndims
        ctypes.POINTER(ctypes.c_int),                                         # ierr
        ctypes.POINTER(ctypes.c_int)                                          # clen
    ]
    arrays_lib.get_array_metadata_C.restype = None

    # Aufruf
    arrays_lib.get_array_metadata_C(
        filename_ascii,
        fn_len,
        dims_out,
        ctypes.byref(ndims),
        ctypes.byref(ierr),
        ctypes.byref(clen)
    )

    check_err_code(ierr)

    if with_clen:
        return dims_out[:ndims.value], clen.value
    else:
        return dims_out[:ndims.value]

# serilization of an n-dimensional integer array
def tox_serialize_int_nd(arr: np.ndarray, filename: str):
    """
    Serializes an n-dimensional integer32 array to a binary file
    """
    if not isinstance(arr, np.ndarray) or arr.dtype != np.int32:
        raise ValueError("arr must be a numpy array of int32")

    # Make sure layout is fortran compatible
    arr_f = np.asfortranarray(arr)

    # dimensions
    dims = np.array(arr.shape, dtype=np.int32)
    ndim = arr.ndim
    ierr = ctypes.c_int()

    # flat array to pass to fortran
    flat = arr_f.ravel(order='F')

    # prepare filename
    filename_ascii, fn_len = _filename_to_ascii_array(filename)

    arrays_lib.serialize_int_nd_C.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # arr
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # dims
        ctypes.c_int,  # ndim
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # filename_ascii
        ctypes.c_int,  # fn_len
        ctypes.POINTER(ctypes.c_int) 
    ]
    arrays_lib.serialize_int_nd_C.restype = None

    # call function
    arrays_lib.serialize_int_nd_C(
        flat,
        dims,
        ndim,
        filename_ascii,
        fn_len,
        ctypes.byref(ierr)
    )

    check_err_code(ierr)

# Deserialize an n dimensional integer array
def tox_deserialize_int_nd(filename):
    """
    Deserializes an n-dimensional int32-Array.
    """
    # read size of the array
    dims = tox_get_array_metadata(filename)
    print(f"Deserializing array with dimensions: {dims}")
    # create array with the proper size
    total_size = np.prod(dims)
    arr = np.zeros(total_size, dtype=np.int32, order='F')  # gets a 1D array
    ascii_arr, fn_len = _filename_to_ascii_array(filename)
    ierr = ctypes.c_int()

    arrays_lib.deserialize_int_C.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="F_CONTIGUOUS"),  # arr
        ctypes.c_int,                                                          # total size
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # filename_ascii
        ctypes.c_int,                                                          # fn_len
        ctypes.POINTER(ctypes.c_int)                                           # ierr
    ]
    arrays_lib.deserialize_int_C.restype = None

    arrays_lib.deserialize_int_C(arr, total_size, ascii_arr, fn_len, ctypes.byref(ierr))
    check_err_code(ierr)
    return arr.reshape(dims, order='F')  # Reshape to original dimensions

def tox_serialize_real_nd(arr: np.ndarray, filename: str):
    """
    Serializes an n-dimensional real64 array to a binary file
    """
    if not isinstance(arr, np.ndarray) or arr.dtype != np.float64:
        raise ValueError("arr must be a numpy array of float64")

    # make sure layout is fortran compatible
    arr_f = np.asfortranarray(arr)

    # dimensions
    dims = np.array(arr.shape, dtype=np.int32)
    ndim = arr.ndim
    ierr = ctypes.c_int()

    # flat array with fortran order
    flat = arr_f.ravel(order='F')

    # ASCII-Filename preparation
    filename_ascii, fn_len = _filename_to_ascii_array(filename)

    # declare args
    arrays_lib.serialize_real_nd_C.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"), # arr
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # dims
        ctypes.c_int,  # ndim
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # filename_ascii
        ctypes.c_int,  # fn_len
        ctypes.POINTER(ctypes.c_int)  # ierr
    ]
    arrays_lib.serialize_real_nd_C.restype = None

    # call function
    arrays_lib.serialize_real_nd_C(
        flat,
        dims,
        ndim,
        filename_ascii,
        fn_len,
        ctypes.byref(ierr)
    )
    check_err_code(ierr)

def tox_deserialize_real_nd(filename):
    """
    Deserializes an n-dimensional array of type real64
    """
    #read dimensions
    dims = tox_get_array_metadata(filename)
    print(f"Deserializing array with dimensions: {dims}")
    # create array with correct size
    total_size = np.prod(dims)
    arr = np.zeros(total_size, dtype=np.float64, order='F')  # accept flat array
    ascii_arr, fn_len = _filename_to_ascii_array(filename)
    ierr = ctypes.c_int()

    arrays_lib.deserialize_real_C.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"),  # arr
        ctypes.c_int,                                                          # total size
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # filename_ascii
        ctypes.c_int,                                                           # fn_len
        ctypes.POINTER(ctypes.c_int)                                           # ierr
    ]
    arrays_lib.deserialize_real_C.restype = None

    arrays_lib.deserialize_real_C(arr, total_size, ascii_arr, fn_len, ctypes.byref(ierr))
    check_err_code(ierr)
    return arr.reshape(dims, order='F')  # Reshape

def tox_serialize_char_nd(arr: np.ndarray, filename: str):
    """
    Serializes an n-dimensional character array to a binary file
    """
    if not isinstance(arr, np.ndarray) or arr.dtype.kind != 'U':
        raise ValueError("arr must be a numpy array of strings (dtype='U')")

    dims = np.array(arr.shape, dtype=np.int32)
    ndim = arr.ndim
    ierr = ctypes.c_int()

    ascii_mat, clen = _string_array_to_ascii_matrix(arr)
    ascii_ptr = np.ascontiguousarray(ascii_mat.ravel(order='F'))

    filename_ascii, fn_len = _filename_to_ascii_array(filename)

    arrays_lib.serialize_char_flat_C.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),  # ascii_ptr
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),  # dims
        ctypes.c_int,                                                          # ndim
        ctypes.c_int,                                                          # clen
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),  # filename_ascii
        ctypes.c_int,                                                           # fn_len
        ctypes.POINTER(ctypes.c_int)                                           # ierr
    ]
    arrays_lib.serialize_char_flat_C.restype = None

    arrays_lib.serialize_char_flat_C(
        ascii_ptr,
        dims,
        ndim,
        clen,
        filename_ascii,
        fn_len,
        ctypes.byref(ierr)
    )
    check_err_code(ierr)


def tox_deserialize_char_nd(filename: str, ndim_max=5):
    """
    Deserializes an n-dimensional character array from a file
    """
    #convert filename to ASCII array
    filename_ascii, fn_len = _filename_to_ascii_array(filename)

    # Get metadata
    dims, clen = tox_get_array_metadata(filename, max_dims=ndim_max, with_clen=True)
    ierr = ctypes.c_int()

    print(f"Deserializing char array with dimensions: {dims}, clen: {clen}")
    total_array_size = np.prod(dims)
    # SCII array of size clen x total
    ascii_arr = np.asfortranarray(np.zeros((clen, total_array_size), dtype=np.int32))

    # declare arguments
    arrays_lib.deserialize_char_flat_C.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags="F_CONTIGUOUS"), # ascii_arr
        ctypes.c_int,           # clen
        ctypes.c_int,           # total
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # filename_ascii
        ctypes.c_int,            # fn_len
        ctypes.POINTER(ctypes.c_int)                                           # ierr
    ]
    arrays_lib.deserialize_char_flat_C.restype = None

    arrays_lib.deserialize_char_flat_C(
        ascii_arr,
        clen,
        total_array_size,
        filename_ascii,
        fn_len,
        ctypes.byref(ierr)
    )
    check_err_code(ierr)

    shape = tuple(dims)
    total_array_size = np.prod(shape)

    ascii_arr_2d = ascii_arr.reshape((clen, total_array_size), order='F')

    # ASCII → String-Array
    chars = np.array([
    ''.join(chr(c) for c in ascii_arr_2d[:clen, i] if c > 0)
        for i in range(total_array_size)
    ], dtype=f'U{clen}')

    chars = chars.reshape(shape, order='F')

    return chars