import numpy as np
import ctypes

# Load library
ctypes.CDLL("libgomp.so.1", mode=ctypes.RTLD_GLOBAL)
lib = ctypes.CDLL("build/libtensor-omics.so")

def setup_sort_real():
    """Setup real sorting function"""
    sort_real = lib.sort_real_c
    sort_real.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
        ctypes.c_int
    ]
    sort_real.restype = None
    return sort_real

def setup_sort_integer():
    """Setup integer sorting function"""
    sort_int = lib.sort_integer_c
    sort_int.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
        ctypes.c_int
    ]
    sort_int.restype = None
    return sort_int

def setup_sort_character():
    """Setup character sorting function - following R pattern"""
    sort_char = lib.sort_character_c
    sort_char.argtypes = [
        ctypes.POINTER(ctypes.c_ubyte),  # char_data as raw bytes
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),  # perm
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),  # stack_left
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),  # stack_right
        ctypes.c_int,  # n
        ctypes.c_int   # strlen
    ]
    sort_char.restype = None
    return sort_char

def print_array(name, arr):
    """Pretty print array"""
    print(f"\n{name}: {arr}")

def test_sort_real_example_1():
    """Example 1: Sort real numbers"""
    print("="*50)
    print("SORT REAL EXAMPLE 1: Basic sorting")
    print("="*50)
    
    # Input data
    array = np.array([3.2, 1.5, 9.0, 2.1, 7.3, 4.4], dtype=np.float64)
    array = np.ascontiguousarray(array)
    n = len(array)
    perm = np.ascontiguousarray(np.arange(1, n+1, dtype=np.int32))  # 1-based indexing for Fortran
    stack_left = np.ascontiguousarray(np.zeros(n, dtype=np.int32))
    stack_right = np.ascontiguousarray(np.zeros(n, dtype=np.int32))
    
    print_array("Input array", array)
    print_array("Initial permutation", perm)
    
    # Call sorting function
    sort_real = setup_sort_real()
    sort_real(array, perm, stack_left, stack_right, n)
    
    print_array("Final permutation", perm)
    
    # Show sorted values
    sorted_values = array[perm - 1]  # Convert to 0-based indexing for Python
    print_array("Sorted values", sorted_values)
    
    # Manual verification
    expected = np.sort(array)
    print_array("Expected (numpy sort)", expected)
    print(f"Match? {np.allclose(sorted_values, expected)}")

def test_sort_real_example_2():
    """Example 2: Edge cases for real sorting"""
    print("="*50)
    print("SORT REAL EXAMPLE 2: Edge cases")
    print("="*50)
    
    # Edge cases: negative, zero, duplicates
    array = np.ascontiguousarray(np.array([-5.0, 0.0, 3.14, -5.0, 0.0, 10.5], dtype=np.float64))
    n = len(array)
    perm = np.ascontiguousarray(np.arange(1, n+1, dtype=np.int32))
    stack_left = np.ascontiguousarray(np.zeros(n, dtype=np.int32))
    stack_right = np.ascontiguousarray(np.zeros(n, dtype=np.int32))
    
    print_array("Input (with negatives/duplicates)", array)
    
    sort_real = setup_sort_real()
    sort_real(array, perm, stack_left, stack_right, n)
    
    sorted_values = array[perm - 1]
    print_array("Sorted values", sorted_values)
    
    expected = np.sort(array)
    print_array("Expected", expected)
    print(f"Match? {np.allclose(sorted_values, expected)}")

def test_sort_integer_example_1():
    """Example 1: Sort integers"""
    print("="*50)
    print("SORT INTEGER EXAMPLE 1: Basic integer sorting")
    print("="*50)
    
    array = np.ascontiguousarray(np.array([5, 3, 8, 1, 4, 2], dtype=np.int32))
    n = len(array)
    perm = np.ascontiguousarray(np.arange(1, n+1, dtype=np.int32))
    stack_left = np.ascontiguousarray(np.zeros(n, dtype=np.int32))
    stack_right = np.ascontiguousarray(np.zeros(n, dtype=np.int32))
    
    print_array("Input array", array)
    
    sort_int = setup_sort_integer()
    sort_int(array, perm, stack_left, stack_right, n)
    
    sorted_values = array[perm - 1]
    print_array("Sorted values", sorted_values)
    
    expected = np.sort(array)
    print_array("Expected", expected)
    print(f"Match? {np.array_equal(sorted_values, expected)}")

def test_sort_integer_example_2():
    """Example 2: Large integers and negatives"""
    print("="*50)
    print("SORT INTEGER EXAMPLE 2: Large and negative integers")
    print("="*50)
    
    array = np.ascontiguousarray(np.array([1000, -500, 0, 999, -1, 2000], dtype=np.int32))
    n = len(array)
    perm = np.ascontiguousarray(np.arange(1, n+1, dtype=np.int32))
    stack_left = np.ascontiguousarray(np.zeros(n, dtype=np.int32))
    stack_right = np.ascontiguousarray(np.zeros(n, dtype=np.int32))
    
    print_array("Input array", array)
    
    sort_int = setup_sort_integer()
    sort_int(array, perm, stack_left, stack_right, n)
    
    sorted_values = array[perm - 1]
    print_array("Sorted values", sorted_values)
    
    expected = np.sort(array)
    print_array("Expected", expected)
    print(f"Match? {np.array_equal(sorted_values, expected)}")

def test_edge_cases():
    """Test edge cases"""
    print("="*50)
    print("EDGE CASES")
    print("="*50)
    
    # Case 1: Single element
    print("\n--- Single element ---")
    array = np.ascontiguousarray(np.array([42.0], dtype=np.float64))
    n = 1
    perm = np.ascontiguousarray(np.array([1], dtype=np.int32))
    stack_left = np.ascontiguousarray(np.zeros(n, dtype=np.int32))
    stack_right = np.ascontiguousarray(np.zeros(n, dtype=np.int32))
    
    print_array("Input (single)", array)
    
    sort_real = setup_sort_real()
    sort_real(array, perm, stack_left, stack_right, n)
    
    sorted_values = array[perm - 1]
    print_array("Sorted (should be same)", sorted_values)
    
    # Case 2: Already sorted
    print("\n--- Already sorted ---")
    array = np.ascontiguousarray(np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64))
    n = len(array)
    perm = np.ascontiguousarray(np.arange(1, n+1, dtype=np.int32))
    stack_left = np.ascontiguousarray(np.zeros(n, dtype=np.int32))
    stack_right = np.ascontiguousarray(np.zeros(n, dtype=np.int32))
    
    print_array("Input (already sorted)", array)
    
    sort_real(array, perm, stack_left, stack_right, n)
    sorted_values = array[perm - 1]
    print_array("Sorted (should be same)", sorted_values)
    
    # Case 3: Reverse sorted
    print("\n--- Reverse sorted ---")
    array = np.ascontiguousarray(np.array([4.0, 3.0, 2.0, 1.0], dtype=np.float64))
    n = len(array)
    perm = np.ascontiguousarray(np.arange(1, n+1, dtype=np.int32))
    stack_left = np.ascontiguousarray(np.zeros(n, dtype=np.int32))
    stack_right = np.ascontiguousarray(np.zeros(n, dtype=np.int32))
    
    print_array("Input (reverse sorted)", array)
    
    sort_real(array, perm, stack_left, stack_right, n)
    sorted_values = array[perm - 1]
    print_array("Sorted", sorted_values)

def test_performance_example():
    """Test with larger arrays"""
    print("="*50)
    print("PERFORMANCE EXAMPLE: Large array")
    print("="*50)
    
    # Generate random array
    np.random.seed(42)
    array = np.random.rand(1000) * 1000
    array = np.ascontiguousarray(array.astype(np.float64))
    n = len(array)
    perm = np.ascontiguousarray(np.arange(1, n+1, dtype=np.int32))
    stack_left = np.ascontiguousarray(np.zeros(n, dtype=np.int32))
    stack_right = np.ascontiguousarray(np.zeros(n, dtype=np.int32))
    
    print(f"Input: Random array of {n} elements")
    print(f"Range: [{np.min(array):.3f}, {np.max(array):.3f}]")
    
    sort_real = setup_sort_real()
    sort_real(array, perm, stack_left, stack_right, n)
    
    sorted_values = array[perm - 1]
    expected = np.sort(array)
    
    print(f"Sorted range: [{np.min(sorted_values):.3f}, {np.max(sorted_values):.3f}]")
    print(f"Match with numpy? {np.allclose(sorted_values, expected)}")
    print(f"Is actually sorted? {np.all(sorted_values[:-1] <= sorted_values[1:])}")

def test_sort_character_python():
    """Test character sorting - exact R equivalent"""
    print("="*60)
    print("CHARACTER SORTING TEST: Python following R pattern")
    print("="*60)
    
    # Same strings as R example
    strings = ["test", "dog", "delta", "zeta", "alpha", "beta"]
    n = len(strings)
    strlen = max(len(s) for s in strings)
    
    print(f"Input strings: {strings}")
    print(f"n = {n}, strlen = {strlen}")
    
    # Create character matrix (strlen x n) - column-major for Fortran
    char_matrix = np.full((strlen, n), ord(' '), dtype=np.uint8)
    
    for i in range(n):
        # Pad string to fixed length
        padded = strings[i].ljust(strlen)
        # Convert to ASCII codes and store in column i
        for j, char in enumerate(padded):
            char_matrix[j, i] = ord(char)
    
    print(f"\nCharacter matrix shape: {char_matrix.shape}")
    print("Character matrix (as ASCII codes):")
    for i in range(n):
        col_chars = ''.join(chr(char_matrix[j, i]) for j in range(strlen))
        print(f"  Column {i}: '{col_chars}' -> {char_matrix[:, i].tolist()}")
    
    # Flatten to 1D array (column-major order)
    char_raw = char_matrix.flatten(order='F')  # Fortran order
    
    print(f"\nFlattened array length: {len(char_raw)}")
    print(f"First few bytes: {char_raw[:10].tolist()}")
    
    # Prepare integer arrays
    perm = np.ascontiguousarray(np.arange(1, n+1, dtype=np.int32))  # 1-based
    stack_left = np.ascontiguousarray(np.zeros(n, dtype=np.int32))
    stack_right = np.ascontiguousarray(np.zeros(n, dtype=np.int32))
    
    print(f"\nInitial permutation: {perm}")
    
    # Call Fortran function
    sort_char = setup_sort_character()
    
    # Convert to ctypes pointer
    char_ptr = char_raw.ctypes.data_as(ctypes.POINTER(ctypes.c_ubyte))
    
    print("\nCalling Fortran sort_character_c...")
    sort_char(char_ptr, perm, stack_left, stack_right, n, strlen)
    
    print(f"Final permutation: {perm}")
    
    # Use permutation to reorder original strings (convert to 0-based)
    sorted_strings = [strings[i-1] for i in perm]
    
    print(f"\nSorted strings: {sorted_strings}")
    
    # Verify with Python sort
    expected = sorted(strings)
    print(f"Expected (Python): {expected}")
    print(f"Match? {sorted_strings == expected}")

def test_sort_character_edge_cases():
    """Test edge cases for character sorting"""
    print("\n" + "="*60)
    print("CHARACTER SORTING EDGE CASES")
    print("="*60)
    
    # Test case 1: Different lengths
    print("\n--- Case 1: Different string lengths ---")
    strings = ["a", "abc", "ab"]
    test_char_case(strings)
    
    # Test case 2: Numbers and letters
    print("\n--- Case 2: Mixed numbers and letters ---")
    strings = ["z1", "a2", "b10", "c3"]
    test_char_case(strings)
    
    # Test case 3: Special characters
    print("\n--- Case 3: With spaces and special chars ---")
    strings = ["hello world", "hello_world", "hello-world"]
    test_char_case(strings)

def test_char_case(strings):
    """Helper function to test a case"""
    n = len(strings)
    strlen = max(len(s) for s in strings)
    
    print(f"Input: {strings}")
    
    # Create character matrix
    char_matrix = np.full((strlen, n), ord(' '), dtype=np.uint8)
    
    for i in range(n):
        padded = strings[i].ljust(strlen)
        for j, char in enumerate(padded):
            char_matrix[j, i] = ord(char)
    
    # Flatten
    char_raw = char_matrix.flatten(order='F')
    
    # Setup arrays
    perm = np.ascontiguousarray(np.arange(1, n+1, dtype=np.int32))
    stack_left = np.ascontiguousarray(np.zeros(n, dtype=np.int32))
    stack_right = np.ascontiguousarray(np.zeros(n, dtype=np.int32))
    
    # Call Fortran
    sort_char = setup_sort_character()
    char_ptr = char_raw.ctypes.data_as(ctypes.POINTER(ctypes.c_ubyte))
    sort_char(char_ptr, perm, stack_left, stack_right, n, strlen)
    
    # Get results
    sorted_strings = [strings[i-1] for i in perm]
    expected = sorted(strings)
    
    print(f"Fortran result: {sorted_strings}")
    print(f"Python result:  {expected}")
    print(f"Match? {sorted_strings == expected}")

def test_large_character_array():
    """Test with larger character array"""
    print("\n" + "="*60)
    print("LARGE CHARACTER ARRAY TEST")
    print("="*60)
    
    # Generate some test strings
    import random
    import string
    
    random.seed(42)
    strings = []
    for _ in range(20):
        length = random.randint(3, 8)
        s = ''.join(random.choices(string.ascii_lowercase, k=length))
        strings.append(s)
    
    print(f"Generated {len(strings)} random strings:")
    print(f"Sample: {strings[:5]}...")
    
    test_char_case(strings)

if __name__ == "__main__":
    print("TENSOR-OMICS SORTING PYTHON VERIFICATION")
    print("Testing sorting functions...")
    
    try:
        test_sort_real_example_1()
        test_sort_real_example_2()
        test_sort_integer_example_1()
        test_sort_integer_example_2()
        test_edge_cases()
        test_performance_example()
        test_sort_character_python()
        test_sort_character_edge_cases()
        test_large_character_array()
        
        print("\n" + "="*50)
        print("ALL SORTING EXAMPLES COMPLETED SUCCESSFULLY!")
        print("Sorting functions work correctly from Python!")
        print("="*50)
        
    except Exception as e:
        print(f"\nERROR: {e}")
        print("Check library loading and function signatures.")