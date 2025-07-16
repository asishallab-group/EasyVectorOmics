import ctypes
import numpy as np
import os

# Lade die Fortran-Bibliothek
lib_path = os.path.join(os.path.dirname(__file__), "../build/arrays.so")
arrays_lib = ctypes.CDLL(lib_path)

def _filename_to_ascii_array(filename):
    ascii_arr = np.array([ord(c) for c in filename], dtype=np.int32)
    return ascii_arr, np.int32(len(ascii_arr))

def _string_array_to_ascii_matrix(strings: np.ndarray) -> tuple[np.ndarray, int]:
    """
    Wandelt ein n-dimensionales Array von Strings in eine ASCII-Matrix (clen, total_size)
    """
    if not isinstance(strings, np.ndarray) or strings.dtype.kind != 'U':
        raise ValueError("Input must be a numpy array of strings (dtype='U')")
    
    flat = strings.flatten(order='F')  # wichtig: Fortran-order
    clen = max(len(s) for s in flat)
    total = flat.size

    mat = np.zeros((clen, total), dtype=np.int32, order='F')
    for i, s in enumerate(flat):
        codes = [ord(c) for c in s]
        mat[:len(codes), i] = codes

    return mat, clen

def get_array_dims(filename, max_dims=5):
    """
    Liest die Dimensionen eines gespeicherten Arrays über get_array_dims().
    """
    ascii_arr, fn_len = _filename_to_ascii_array(filename)
    dims_out = np.zeros(max_dims, dtype=np.int32, order='F')
    ndims = np.zeros(1, dtype=np.int32, order='F')

    arrays_lib.get_array_dims_C.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # filename_ascii
        ctypes.c_int,                                                          # fn_len
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # dims_out
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS")   # ndims
    ]
    arrays_lib.get_array_dims_C.restype = None

    arrays_lib.get_array_dims_C(ascii_arr, fn_len, dims_out, ndims)

    return dims_out[:ndims[0]]

def get_char_array_metadata(filename: str, max_ndims: int = 10):
    """
    Ruft die Metadaten einer char-Array-Datei ab: Dimensionen, Anzahl Dimensionen, Typcode, clen.
    Erwartet, dass die Fortran-Subroutine get_array_metadata_chars_C eingebunden ist.
    """
    # ASCII-Filename vorbereiten
    filename_ascii, fn_len = _filename_to_ascii_array(filename)

    # Ausgabe-Puffer vorbereiten
    dims_out = np.zeros(max_ndims, dtype=np.int32)
    ndims = ctypes.c_int()
    type_code = ctypes.c_int()
    clen = ctypes.c_int()

    # Argument-Typen deklarieren
    arrays_lib.get_array_metadata_chars_C.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # filename_ascii
        ctypes.c_int,                                                         # fn_len
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # dims_out
        ctypes.c_int,                                                         # dims_len
        ctypes.POINTER(ctypes.c_int),                                         # ndims
        ctypes.POINTER(ctypes.c_int),                                         # type_code_out
        ctypes.POINTER(ctypes.c_int),                                         # clen_out
    ]
    arrays_lib.get_array_metadata_chars_C.restype = None

    # Aufruf
    arrays_lib.get_array_metadata_chars_C(
        filename_ascii,
        fn_len,
        dims_out,
        max_ndims,
        ctypes.byref(ndims),
        ctypes.byref(type_code),
        ctypes.byref(clen)
    )

    # Rückgabe: nur relevante Dimensionen
    return dims_out[:ndims.value], clen.value

# Python-Wrapper
def serialize_int_nd(arr: np.ndarray, filename: str):
    """
    Serialisiert ein n-dimensionales int32-Array in eine Binärdatei mit Header.
    """
    if not isinstance(arr, np.ndarray) or arr.dtype != np.int32:
        raise ValueError("arr must be a numpy array of int32")

    # Sicherstellen: Fortran-kompatibles Layout
    arr_f = np.asfortranarray(arr)

    # Dimensionen
    dims = np.array(arr.shape, dtype=np.int32)
    ndim = arr.ndim

    # Flaches Array (Fortran-richtig)
    flat = arr_f.ravel(order='F')

    # ASCII-Filename vorbereiten
    filename_ascii, fn_len = _filename_to_ascii_array(filename)

    # Fortran-Funktion deklarieren
    arrays_lib.serialize_int_nd_C.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # arr
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # dims
        ctypes.c_int,  # ndim
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # filename_ascii
        ctypes.c_int  # fn_len
    ]
    arrays_lib.serialize_int_nd_C.restype = None

    # Aufruf
    arrays_lib.serialize_int_nd_C(
        flat,
        dims,
        ndim,
        filename_ascii,
        fn_len
    )


def deserialize_int_nd(filename):
    """
    Hauptfunktion: Deserialisiert ein beliebig-dimensionales int32-Array.
    """
    dims = get_array_dims(filename)
    print(f"Deserializing array with dimensions: {dims}")
    total_size = np.prod(dims)
    arr = np.zeros(total_size, dtype=np.int32, order='F')  # 1D flach, später reshapen
    ascii_arr, fn_len = _filename_to_ascii_array(filename)

    arrays_lib.deserialize_int_C.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # arr
        ctypes.c_int,                                                          # total size
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # filename_ascii
        ctypes.c_int                                                           # fn_len
    ]
    arrays_lib.deserialize_int_C.restype = None

    arrays_lib.deserialize_int_C(arr, total_size, ascii_arr, fn_len)
    return arr.reshape(dims, order='F')  # Reshape in die Originaldimensionen

def serialize_real_nd(arr: np.ndarray, filename: str):
    """
    Serialisiert ein n-dimensionales float64-Array in eine Binärdatei mit Header.
    """
    if not isinstance(arr, np.ndarray) or arr.dtype != np.float64:
        raise ValueError("arr must be a numpy array of float64")

    # Sicherstellen: Fortran-kompatibles Layout
    arr_f = np.asfortranarray(arr)

    # Dimensionen
    dims = np.array(arr.shape, dtype=np.int32)
    ndim = arr.ndim

    # Flaches Array (Fortran-richtig)
    flat = arr_f.ravel(order='F')

    # ASCII-Filename vorbereiten
    filename_ascii, fn_len = _filename_to_ascii_array(filename)

    # Fortran-Funktion deklarieren
    arrays_lib.serialize_real_nd_C.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"), # arr
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # dims
        ctypes.c_int,  # ndim
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # filename_ascii
        ctypes.c_int  # fn_len
    ]
    arrays_lib.serialize_real_nd_C.restype = None

    # Aufruf
    arrays_lib.serialize_real_nd_C(
        flat,
        dims,
        ndim,
        filename_ascii,
        fn_len
    )


def deserialize_real_nd(filename):
    """
    Hauptfunktion: Deserialisiert ein beliebig-dimensionales int32-Array.
    """
    dims = get_array_dims(filename)
    print(f"Deserializing array with dimensions: {dims}")
    total_size = np.prod(dims)
    arr = np.zeros(total_size, dtype=np.float64, order='F')  # 1D flach, später reshapen
    ascii_arr, fn_len = _filename_to_ascii_array(filename)

    arrays_lib.deserialize_real_C.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"),  # arr
        ctypes.c_int,                                                          # total size
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # filename_ascii
        ctypes.c_int                                                           # fn_len
    ]
    arrays_lib.deserialize_real_C.restype = None

    arrays_lib.deserialize_real_C(arr, total_size, ascii_arr, fn_len)
    return arr.reshape(dims, order='F')  # Reshape in die Originaldimensionen

def serialize_char_nd(arr: np.ndarray, filename: str):
    """
    Serialisiert ein n-dimensionales String-Array in eine Binärdatei (Fortran-kompatibel).
    """
    if not isinstance(arr, np.ndarray) or arr.dtype.kind != 'U':
        raise ValueError("arr must be a numpy array of strings (dtype='U')")

    dims = np.array(arr.shape, dtype=np.int32)
    ndim = arr.ndim

    ascii_mat, clen = _string_array_to_ascii_matrix(arr)
    ascii_ptr = np.ascontiguousarray(ascii_mat.ravel(order='F'))

    filename_ascii, fn_len = _filename_to_ascii_array(filename)

    arrays_lib.serialize_char_flat_C.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),  # ascii_ptr
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),  # dims
        ctypes.c_int,                                                          # ndim
        ctypes.c_int,                                                          # clen
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),  # filename_ascii
        ctypes.c_int                                                           # fn_len
    ]
    arrays_lib.serialize_char_flat_C.restype = None

    arrays_lib.serialize_char_flat_C(
        ascii_ptr,
        dims,
        ndim,
        clen,
        filename_ascii,
        fn_len
    )


def deserialize_char_nd(filename: str, ndim_max=5):
    """
    Deserialisiert ein n-dimensionales Array von Strings aus Binärformat.
    """
    #convert filename to ASCII array
    filename_ascii, fn_len = _filename_to_ascii_array(filename)

    dims_out = np.zeros(ndim_max, dtype=np.int32)

    # Puffergrößen (zunächst nur Dateigröße abschätzen, hier z. B. max. 1000 Strings)
    dims, clen = get_char_array_metadata(filename, ndim_max)

    print(f"Deserializing char array with dimensions: {dims}, clen: {clen}")
    total = np.prod(dims)
    # ASCII-Puffer vorbereiten (clen x total)
    ascii_arr = np.asfortranarray(np.zeros((clen, total), dtype=np.int32))

    # Rückgabewerte: byref-Variablen für ndim/clen_out
    ndim_out = ctypes.c_int()
    clen_out = ctypes.c_int()

    # Fortran-Funktion definieren
    arrays_lib.deserialize_char_flat_C.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags="F_CONTIGUOUS"), # ascii_arr
        ctypes.c_int,           # clen
        ctypes.c_int,           # total
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # dims_out
        ctypes.POINTER(ctypes.c_int),  # ndim_out
        ctypes.POINTER(ctypes.c_int),  # clen_out
        np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # filename_ascii
        ctypes.c_int,           # fn_len
        ctypes.c_int            # ndim_actual
    ]
    arrays_lib.deserialize_char_flat_C.restype = None

    arrays_lib.deserialize_char_flat_C(
        ascii_arr,
        clen,
        total,
        dims_out,
        ctypes.byref(ndim_out),
        ctypes.byref(clen_out),
        filename_ascii,
        fn_len,
        ndim_max
    )

    ndim = ndim_out.value
    clen = clen_out.value
    shape = tuple(dims_out[:ndim])
    total = np.prod(shape)

    # Wichtig: 2D wiederherstellen
    ascii_arr_2d = ascii_arr.reshape((clen, total), order='F')

    # ASCII → String-Array
    chars = np.array([
    ''.join(chr(c) for c in ascii_arr_2d[:clen, i] if c > 0)
        for i in range(total)
    ], dtype=f'U{clen}')

    chars = chars.reshape(shape, order='F')

    return chars


def int_test():
    array = np.array([1, 2, 3, 4, 5], dtype=np.int32, order='F')
    filename = "test_int_1d.bin"
    serialize_int_nd(array, filename)
    print(f"Serialized array to {filename}")
    res = deserialize_int_nd(filename)
    assert np.array_equal(res, array), "Deserialized array does not match original"
    print(res)

    array_2d = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.int32, order='F')
    filename_2d = "test_int_2d.bin"
    serialize_int_nd(array_2d, filename_2d)
    print(f"Serialized 2D array to {filename_2d}")
    res_2d = deserialize_int_nd(filename_2d)
    print(res_2d)
    assert np.array_equal(res_2d, array_2d), "Deserialized 2D array does not match original"

    array_3d = np.array([[[1, 2], [3, 4]], [[5, 6], [7, 8]]], dtype=np.int32, order='F')
    filename_3d = "test_int_3d.bin"
    serialize_int_nd(array_3d, filename_3d)
    print(f"Serialized 3D array to {filename_3d}")
    res_3d = deserialize_int_nd(filename_3d)
    print(res_3d)
    assert np.array_equal(res_3d, array_3d), "Deserialized 3D array does not match original"
    
def real_test():
    array = np.array([1.5, 2.3, 3.2, 4.0, 5.0], dtype=np.float64, order='F')
    filename = "test_real_1d.bin"
    serialize_real_nd(array, filename)
    print(f"Serialized array to {filename}")
    res = deserialize_real_nd(filename)
    assert np.array_equal(res, array), "Deserialized array does not match original"
    print(res)

    array_2d = np.array([[1.0, 2.0, 3.0], [4.0, 5.7, 6.0]], dtype=np.float64, order='F')
    filename_2d = "test_real_2d.bin"
    serialize_real_nd(array_2d, filename_2d)
    print(f"Serialized 2D array to {filename_2d}")
    res_2d = deserialize_real_nd(filename_2d)
    print(res_2d)
    assert np.array_equal(res_2d, array_2d), "Deserialized 2D array does not match original"

    array_3d = np.array([[[1.0, 2.0], [3.3, 4.0]], [[5.0, 6.8], [7.0, 8.0]]], dtype=np.float64, order='F')
    filename_3d = "test_real_3d.bin"
    serialize_real_nd(array_3d, filename_3d)
    print(f"Serialized 3D array to {filename_3d}")
    res_3d = deserialize_real_nd(filename_3d)
    print(res_3d)
    assert np.array_equal(res_3d, array_3d), "Deserialized 3D array does not match original"

    empty_array = np.array([], dtype=np.float64, order='F')
    empty_filename = "test_real_empty.bin"
    serialize_real_nd(empty_array, empty_filename)
    print(f"Serialized empty array to {empty_filename}")
    res_empty = deserialize_real_nd(empty_filename)
    assert res_empty.size == 0, "Deserialized empty array should be empty"

def char_test():
    array = np.array(["hello", "world"], dtype='U5', order='F')
    filename = "test_char_1d.bin"
    serialize_char_nd(array, filename)
    print(f"Serialized array to {filename}")
    res = deserialize_char_nd(filename)
    print(res)
    assert np.array_equal(res, array), "Deserialized array does not match original"
    

    array_2d = np.array([["foo", "bar"], ["baz", "qux"]], dtype='U5', order='F')
    filename_2d = "test_char_2d.bin"
    serialize_char_nd(array_2d, filename_2d)
    print(f"Serialized 2D array to {filename_2d}")
    res_2d = deserialize_char_nd(filename_2d)
    print(res_2d)
    assert np.array_equal(res_2d, array_2d), "Deserialized 2D array does not match original"

    array_3d = np.array([[["abb", "bbbbbbb"], ["cfs", "d"]], [["e", ""], ["g", "h"]]], dtype='U5', order='F')
    filename_3d = "test_char_3d.bin"
    serialize_char_nd(array_3d, filename_3d)
    print(f"Serialized 3D array to {filename_3d}")
    res_3d = deserialize_char_nd(filename_3d)
    print(res_3d)
    assert np.array_equal(res_3d, array_3d), "Deserialized 3D array does not match original"

print("============ INT TEST ============")
int_test()
print("============ REAL TEST ============")
real_test()
print("============ CHAR TEST ============")
char_test()