import ctypes
import numpy as np
import os

# Lade die Fortran-Bibliothek
lib_path = os.path.join(os.path.dirname(__file__), "../arrays.so")
arrays_lib = ctypes.CDLL(lib_path)

# Hilfsfunktion: Übergibt einen Dateinamen als C-String
def _to_c_char_p(filename):
    return ctypes.c_char_p(filename.encode('utf-8'))

def deserialize_int_1d(filename):
    """
    Wrapper für die Fortran-Funktion deserialize_int_1d_C.
    :param filename: Dateiname als String
    :return: 1D numpy array (int32)
    """
    # Annahme: Die Fortran-Funktion allokiert das Array und gibt einen Zeiger zurück
    # Signatur: void deserialize_int_1d_C(int **arr, char *filename)
    # Wir müssen einen Zeiger auf einen Zeiger übergeben
    arr_ptr = ctypes.POINTER(ctypes.c_int)()
    arrays_lib.deserialize_int_1d_C.argtypes = [ctypes.POINTER(ctypes.POINTER(ctypes.c_int)), ctypes.c_char_p]
    arrays_lib.deserialize_int_1d_C.restype = None
    arrays_lib.deserialize_int_1d_C(ctypes.byref(arr_ptr), _to_c_char_p(filename))
    # Die Größe muss separat bestimmt werden (z.B. durch Metadaten oder Rückgabewert)
    # Hier als Platzhalter: Rückgabe als raw pointer (ohne Größe)
    return arr_ptr  # In der Praxis: numpy.frombuffer(..., dtype=np.int32, count=...)

def deserialize_int_2d(filename):
    """
    Wrapper für die Fortran-Funktion deserialize_int_2d_C.
    :param filename: Dateiname als String
    :return: 2D numpy array (int32)
    """
    arr_ptr = ctypes.POINTER(ctypes.c_int)()
    arrays_lib.deserialize_int_2d_C.argtypes = [ctypes.POINTER(ctypes.POINTER(ctypes.c_int)), ctypes.c_char_p]
    arrays_lib.deserialize_int_2d_C.restype = None
    arrays_lib.deserialize_int_2d_C(ctypes.byref(arr_ptr), _to_c_char_p(filename))
    return arr_ptr

def deserialize_int_3d(filename):
    """
    Wrapper für die Fortran-Funktion deserialize_int_3d_C.
    :param filename: Dateiname als String
    :return: 3D numpy array (int32)
    """
    arr_ptr = ctypes.POINTER(ctypes.c_int)()
    arrays_lib.deserialize_int_3d_C.argtypes = [ctypes.POINTER(ctypes.POINTER(ctypes.c_int)), ctypes.c_char_p]
    arrays_lib.deserialize_int_3d_C.restype = None
    arrays_lib.deserialize_int_3d_C(ctypes.byref(arr_ptr), _to_c_char_p(filename))
    return arr_ptr

def deserialize_int_4d(filename):
    """
    Wrapper für die Fortran-Funktion deserialize_int_4d_C.
    :param filename: Dateiname als String
    :return: 4D numpy array (int32)
    """
    arr_ptr = ctypes.POINTER(ctypes.c_int)()
    arrays_lib.deserialize_int_4d_C.argtypes = [ctypes.POINTER(ctypes.POINTER(ctypes.c_int)), ctypes.c_char_p]
    arrays_lib.deserialize_int_4d_C.restype = None
    arrays_lib.deserialize_int_4d_C(ctypes.byref(arr_ptr), _to_c_char_p(filename))
    return arr_ptr

def deserialize_int_5d(filename):
    """
    Wrapper für die Fortran-Funktion deserialize_int_5d_C.
    :param filename: Dateiname als String
    :return: 5D numpy array (int32)
    """
    arr_ptr = ctypes.POINTER(ctypes.c_int)()
    arrays_lib.deserialize_int_5d_C.argtypes = [ctypes.POINTER(ctypes.POINTER(ctypes.c_int)), ctypes.c_char_p]
    arrays_lib.deserialize_int_5d_C.restype = None
    arrays_lib.deserialize_int_5d_C(ctypes.byref(arr_ptr), _to_c_char_p(filename))
    return arr_ptr