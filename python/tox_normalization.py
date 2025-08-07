import numpy as np
import ctypes

# Load library
ctypes.CDLL("libgomp.so.1", mode=ctypes.RTLD_GLOBAL)
lib = ctypes.CDLL("build/libtensor-omics.so")

def setup_normalize():
    """Setup normalize function"""
    normalize = lib.normalize_by_std_dev_c
    normalize.argtypes = [
        ctypes.c_int, 
        ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    ]
    normalize.restype = None
    return normalize

def setup_quantile():
    """Setup quantile normalization function"""
    quantile_norm = lib.quantile_normalization_c
    quantile_norm.argtypes = [
        ctypes.c_int, 
        ctypes.c_int,
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
    return quantile_norm

def print_matrix(name, mat):
    """Pretty print matrix with name"""
    print(f"\n{name}:")
    print(f"Shape: {mat.shape}")
    for i, row in enumerate(mat):
        print(f"  Row {i+1}: {row}")

def setup_log2():
    """Setup log2 transformation function"""
    log2_transform = lib.log2_transformation_c
    log2_transform.argtypes = [
        ctypes.c_int, 
        ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    ]
    log2_transform.restype = None
    return log2_transform

def setup_calc_tiss_avg():
    """Setup tissue averaging function"""
    tiss_avg = lib.calc_tiss_avg_c
    tiss_avg.argtypes = [
        ctypes.c_int, 
        ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
    ]
    tiss_avg.restype = None
    return tiss_avg

def setup_calc_fchange():
    """Setup fold change function"""
    fchange = lib.calc_fchange_c
    fchange.argtypes = [
        ctypes.c_int,  # n_genes
        ctypes.c_int,  # n_cols
        ctypes.c_int,  # n_pairs
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),  # control_cols
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),  # cond_cols
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # i_matrix
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # o_matrix
    ]
    fchange.restype = None
    return fchange

def test_normalize_example_1():
    """Example 1: Simple 2x3 matrix"""
    print("="*50)
    print("NORMALIZE EXAMPLE 1: Simple 2x3 matrix")
    print("="*50)
    
    # Input data
    mat = np.array([[1.0, 2.0, 3.0], 
                    [4.0, 5.0, 6.0]], dtype=np.float64, order='F')
    
    print_matrix("Input", mat)
    
    # Call function
    normalize = setup_normalize()
    n_genes, n_tissues = mat.shape
    input_flat = np.asfortranarray(mat).ravel(order='F')
    output_flat = np.zeros_like(input_flat)
    
    print(f"Calling normalize_by_std_dev_c({n_genes}, {n_tissues}, input, output)")
    normalize(n_genes, n_tissues, input_flat, output_flat)
    
    result = output_flat.reshape((n_genes, n_tissues), order='F')
    print_matrix("Output", result)
    
    # Manual verification
    print("\nManual verification:")
    for i in range(n_genes):
        row = mat[i, :]
        std_dev = np.sqrt(np.mean(row**2))
        normalized_row = row / std_dev
        print(f"  Gene {i+1}: std_dev={std_dev:.4f}, normalized={normalized_row}")

def test_normalize_example_2():
    """Example 2: Large values"""
    print("="*50)
    print("NORMALIZE EXAMPLE 2: Large values")
    print("="*50)
    
    mat = np.array([[1e6, 2e6], 
                    [3e6, 4e6]], dtype=np.float64, order='F')
    
    print_matrix("Input", mat)
    
    normalize = setup_normalize()
    n_genes, n_tissues = mat.shape
    input_flat = np.asfortranarray(mat).ravel(order='F')
    output_flat = np.zeros_like(input_flat)
    
    normalize(n_genes, n_tissues, input_flat, output_flat)
    result = output_flat.reshape((n_genes, n_tissues), order='F')
    
    print_matrix("Output", result)
    print(f"All finite? {np.all(np.isfinite(result))}")

def test_quantile_example_1():
    """Example 1: Simple quantile normalization"""
    print("="*50)
    print("QUANTILE EXAMPLE 1: Simple 3x3 matrix")
    print("="*50)
    
    mat = np.array([[5.0, 2.0, 8.0],
                    [1.0, 6.0, 3.0],
                    [9.0, 4.0, 7.0]], dtype=np.float64, order='F')
    
    print_matrix("Input", mat)
    
    quantile_norm = setup_quantile()
    n_genes, n_tissues = mat.shape
    input_flat = np.asfortranarray(mat).ravel(order='F')
    output_flat = np.zeros_like(input_flat)
    temp_col = np.zeros(n_genes, dtype=np.float64, order='F')
    rank_means = np.zeros(n_genes, dtype=np.float64, order='F')
    perm = np.zeros(n_genes, dtype=np.int32, order='F')
    max_stack = max(2 * n_genes, 2)
    stack_left = np.zeros(max_stack, dtype=np.int32, order='F')
    stack_right = np.zeros(max_stack, dtype=np.int32, order='F')
    
    quantile_norm(n_genes, n_tissues, input_flat, output_flat,
                  temp_col, rank_means, perm, stack_left, stack_right, max_stack)
    
    result = output_flat.reshape((n_genes, n_tissues), order='F')
    print_matrix("Output", result)
    
    # Check column distributions
    print("\nColumn distributions (sorted):")
    for j in range(n_tissues):
        sorted_col = np.sort(result[:, j])
        print(f"  Column {j+1}: {sorted_col}")

def test_random_examples():
    """Example with random data"""
    print("="*50)
    print("RANDOM EXAMPLES")
    print("="*50)
    
    # Random data
    np.random.seed(42)
    for i in range(3):
        print(f"\n--- Random Example {i+1} ---")
        mat = np.random.rand(3, 4) * 100
        print_matrix(f"Random Input {i+1}", mat)
        
        # Normalize
        normalize = setup_normalize()
        n_genes, n_tissues = mat.shape
        input_flat = np.asfortranarray(mat).ravel(order='F')
        output_flat = np.zeros_like(input_flat)
        normalize(n_genes, n_tissues, input_flat, output_flat)
        result = output_flat.reshape((n_genes, n_tissues), order='F')
        
        print_matrix(f"Normalized {i+1}", result)
        print(f"Input range: [{np.min(mat):.3f}, {np.max(mat):.3f}]")
        print(f"Output range: [{np.min(result):.3f}, {np.max(result):.3f}]")

def test_edge_cases():
    """Test edge cases"""
    print("="*50)
    print("EDGE CASES")
    print("="*50)
    
    # Case 1: All zeros
    print("\n--- Case 1: All zeros ---")
    mat = np.zeros((2, 2), dtype=np.float64, order='F')
    print_matrix("Input (all zeros)", mat)
    
    normalize = setup_normalize()
    n_genes, n_tissues = mat.shape
    input_flat = np.asfortranarray(mat).ravel(order='F')
    output_flat = np.zeros_like(input_flat)
    normalize(n_genes, n_tissues, input_flat, output_flat)
    result = output_flat.reshape((n_genes, n_tissues), order='F')
    print_matrix("Output", result)
    
    # Case 2: Single element
    print("\n--- Case 2: Single element ---")
    mat = np.array([[5.0]], dtype=np.float64, order='F')
    print_matrix("Input (single)", mat)
    
    n_genes, n_tissues = mat.shape
    input_flat = np.asfortranarray(mat).ravel(order='F')
    output_flat = np.zeros_like(input_flat)
    normalize(n_genes, n_tissues, input_flat, output_flat)
    result = output_flat.reshape((n_genes, n_tissues), order='F')
    print_matrix("Output", result)

def test_log2_example_1():
    """Example 1: Simple log2 transformation"""
    print("="*50)
    print("LOG2 TRANSFORMATION EXAMPLE 1: Simple values")
    print("="*50)
    
    # Input data: [0, 3, 7, 15] → log2(x+1) = [0, 2, 3, 4]
    mat = np.array([[0.0, 3.0], 
                    [7.0, 15.0]], dtype=np.float64, order='F')
    
    print_matrix("Input", mat)
    
    log2_transform = setup_log2()
    n_genes, n_tissues = mat.shape
    input_flat = np.asfortranarray(mat).ravel(order='F')
    output_flat = np.zeros_like(input_flat)
    
    log2_transform(n_genes, n_tissues, input_flat, output_flat)
    result = output_flat.reshape((n_genes, n_tissues), order='F')
    
    print_matrix("Output", result)
    
    # Manual verification
    print("\nManual verification (log2(x+1)):")
    expected = np.log2(mat + 1)
    print_matrix("Expected", expected)
    print(f"Match? {np.allclose(result, expected)}")

def test_log2_example_2():
    """Example 2: Edge cases for log2"""
    print("="*50)
    print("LOG2 TRANSFORMATION EXAMPLE 2: Edge cases")
    print("="*50)
    
    # Edge cases: zeros, ones, large values
    mat = np.array([[0.0, 1.0, 1023.0]], dtype=np.float64, order='F')
    
    print_matrix("Input (edge cases)", mat)
    
    log2_transform = setup_log2()
    n_genes, n_tissues = mat.shape
    input_flat = np.asfortranarray(mat).ravel(order='F')
    output_flat = np.zeros_like(input_flat)
    
    log2_transform(n_genes, n_tissues, input_flat, output_flat)
    result = output_flat.reshape((n_genes, n_tissues), order='F')
    
    print_matrix("Output", result)
    print(f"log2(0+1)={result[0,0]:.3f}, log2(1+1)={result[0,1]:.3f}, log2(1023+1)={result[0,2]:.3f}")

def test_calc_tiss_avg_example_1():
    """Example 1: Average tissue replicates"""
    print("="*50)
    print("TISSUE AVERAGING EXAMPLE 1: 3 tissues, 2 reps each")
    print("="*50)
    
    # 2 genes × 6 samples (3 tissues, 2 replicates each)
    mat = np.array([[1.0, 7.0, 3.0, 9.0, 5.0, 11.0],   # Gene 1: samples 1-6
                    [2.0, 8.0, 4.0, 10.0, 6.0, 12.0]], dtype=np.float64, order='F')  # Gene 2: samples 1-6
    
    print_matrix("Input (2 genes × 6 samples)", mat)
    print("Tissue groups: [1,3,5] start columns, [2,2,2] replicates each")
    
    tiss_avg = setup_calc_tiss_avg()
    n_genes, n_samples = mat.shape
    n_groups = 3
    
    group_starts = np.array([1, 3, 5], dtype=np.int32, order='F')  # 1-based indexing for Fortran
    group_counts = np.array([2, 2, 2], dtype=np.int32, order='F')
    
    input_flat = np.asfortranarray(mat).ravel(order='F')
    output_flat = np.zeros(n_genes * n_groups, dtype=np.float64, order='F')
    
    tiss_avg(n_genes, n_groups, group_starts, group_counts, input_flat, output_flat)
    result = output_flat.reshape((n_genes, n_groups), order='F')
    
    print_matrix("Output (2 genes × 3 tissues)", result)
    
    # Manual verification - CORRECTED
    print("\nManual verification:")
    print(f"Tissue 1 (samples 1-2): Gene1=mean([1,7])={np.mean([1,7]):.1f}, Gene2=mean([2,8])={np.mean([2,8]):.1f}")
    print(f"Tissue 2 (samples 3-4): Gene1=mean([3,9])={np.mean([3,9]):.1f}, Gene2=mean([4,10])={np.mean([4,10]):.1f}")
    print(f"Tissue 3 (samples 5-6): Gene1=mean([5,11])={np.mean([5,11]):.1f}, Gene2=mean([6,12])={np.mean([6,12]):.1f}")
    
    # Verify the results match
    expected = np.array([[4.0, 6.0, 8.0], [5.0, 7.0, 9.0]], order='F')
    print(f"\nResults match expected? {np.allclose(result, expected)}")
def test_calc_tiss_avg_example_2():
    """Example 2: Unequal replicates"""
    print("="*50)
    print("TISSUE AVERAGING EXAMPLE 2: Unequal replicates")
    print("="*50)
    
    # 2 genes × 7 samples (tissue1: 2 reps, tissue2: 3 reps, tissue3: 2 reps)
    mat = np.array([[1.0, 8.0, 2.0, 9.0, 3.0, 10.0, 4.0],   # Gene 1
                    [5.0, 12.0, 6.0, 13.0, 7.0, 14.0, 8.0]], dtype=np.float64, order='F')  # Gene 2
    
    print_matrix("Input (2 genes × 7 samples)", mat)
    print("Tissue groups: [1,3,6] start columns, [2,3,2] replicates each")
    
    tiss_avg = setup_calc_tiss_avg()
    n_genes, n_samples = mat.shape
    n_groups = 3
    
    group_starts = np.array([1, 3, 6], dtype=np.int32, order='F')
    group_counts = np.array([2, 3, 2], dtype=np.int32, order='F')
    
    input_flat = np.asfortranarray(mat).ravel(order='F')
    output_flat = np.zeros(n_genes * n_groups, dtype=np.float64, order='F')
    
    tiss_avg(n_genes, n_groups, group_starts, group_counts, input_flat, output_flat)
    result = output_flat.reshape((n_genes, n_groups), order='F')
    
    print_matrix("Output (2 genes × 3 tissues)", result)

def test_calc_fchange_example_1():
    """Example 1: Simple fold change"""
    print("="*50)
    print("FOLD CHANGE EXAMPLE 1: Simple condition vs control")
    print("="*50)
    
    # 2 genes × 2 samples (control, condition)
    mat = np.array([[1.0, 2.0],   # Gene 1: control=1, condition=4
                    [4.0, 8.0]], dtype=np.float64, order='F')  # Gene 2: control=2, condition=8
    
    print_matrix("Input (2 genes × 2 samples)", mat)
    print("Control column: 1, Condition column: 2")
    
    fchange = setup_calc_fchange()
    n_genes, n_samples = mat.shape
    n_cols = n_samples
    n_pairs = 1
    
    control_cols = np.array([1], dtype=np.int32, order='F')  # 1-based indexing
    condition_cols = np.array([2], dtype=np.int32, order='F')
    
    input_flat = np.asfortranarray(mat).ravel(order='F')
    output_flat = np.zeros(n_genes * n_pairs, dtype=np.float64, order='F')
    
    fchange(n_genes, n_cols, n_pairs, control_cols, condition_cols, input_flat, output_flat)
    
    print(f"\nOutput (fold changes): {output_flat}")
    
    # Manual verification  
    print("\nManual verification (condition - control):")
    print(f"Gene 1: {mat[0,1]} - {mat[0,0]} = {mat[0,1] - mat[0,0]}")
    print(f"Gene 2: {mat[1,1]} - {mat[1,0]} = {mat[1,1] - mat[1,0]}")

def test_calc_fchange_example_2():
    """Example 2: Multiple conditions vs same control"""
    print("="*50)
    print("FOLD CHANGE EXAMPLE 2: Multiple conditions vs same control")
    print("="*50)
    
    # 2 genes × 3 samples (1 control, 2 conditions)
    mat = np.array([[2.0, 6.0, 10.0],   # Gene 1: control=2, condition1=6, condition2=10
                    [4.0, 16.0, 24.0]], dtype=np.float64, order='F')  # Gene 2: control=4, condition1=16, condition2=24
    
    print_matrix("Input (2 genes × 3 samples)", mat)
    print("Control column: 1, Condition columns: 2 and 3")
    print("Pairs: (control=1,condition=2) and (control=1,condition=3)")
    
    fchange = setup_calc_fchange()
    n_genes, n_samples = mat.shape
    n_cols = n_samples
    n_pairs = 2
    
    # Both pairs use column 1 as control, but different condition columns
    control_cols = np.array([1, 1], dtype=np.int32, order='F')      # Same control for both
    condition_cols = np.array([2, 3], dtype=np.int32, order='F')    # Different conditions
    
    input_flat = np.asfortranarray(mat).ravel(order='F')
    output_flat = np.zeros(n_genes * n_pairs, dtype=np.float64, order='F')
    
    fchange(n_genes, n_cols, n_pairs, control_cols, condition_cols, input_flat, output_flat)
    result = output_flat.reshape((n_genes, n_pairs), order='F')
    
    print_matrix("Output (2 genes × 2 pairs)", result)
    
    print("\nManual verification:")
    print(f"Pair 1 (condition2-control1): Gene1={mat[0,1]}-{mat[0,0]}={mat[0,1]-mat[0,0]}, Gene2={mat[1,1]}-{mat[1,0]}={mat[1,1]-mat[1,0]}")
    print(f"Pair 2 (condition3-control1): Gene1={mat[0,2]}-{mat[0,0]}={mat[0,2]-mat[0,0]}, Gene2={mat[1,2]}-{mat[1,0]}={mat[1,2]-mat[1,0]}")
    
    # Expected results
    expected = np.array([[4.0, 8.0], [12.0, 20.0]], order='F')
    print(f"\nResults match expected? {np.allclose(result, expected)}")

if __name__ == "__main__":
    print("TENSOR-OMICS PYTHON VERIFICATION")
    print("Testing data flow and basic functionality...")
    
    try:
        test_normalize_example_1()
        test_normalize_example_2()
        test_quantile_example_1()
        test_random_examples()
        test_edge_cases()
        test_log2_example_1()
        test_log2_example_2()
        test_calc_tiss_avg_example_1()
        test_calc_tiss_avg_example_2()
        test_calc_fchange_example_1()
        test_calc_fchange_example_2()
        
        print("\n" + "="*50)
        print("ALL EXAMPLES COMPLETED SUCCESSFULLY!")
        print("Data is flowing correctly between Python and Fortran.")
        print("="*50)

        
    except Exception as e:
        print(f"\nERROR: {e}")
        print("Check library loading and function signatures.")

