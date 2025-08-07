#!/usr/bin/env python3
"""
Example of using the Fortran module with real data
This script demonstrates how to use distance_to_centroid with normalized TPM data,
orthogroup centroids and gene-to-family mapping
"""

import numpy as np
import pandas as pd
import ctypes
import os
from pathlib import Path

# Load library
ctypes.CDLL("libgomp.so.1", mode=ctypes.RTLD_GLOBAL)
lib = ctypes.CDLL("build/libtensor-omics.so")

def setup_euclidean_distance():
    """Setup euclidean distance function"""
    euclidean_distance = lib.euclidean_distance_c
    euclidean_distance.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_double)
    ]
    euclidean_distance.restype = None
    return euclidean_distance

def setup_distance_to_centroid():
    """Setup distance to centroid function"""
    distance_to_centroid = lib.distance_to_centroid_c
    distance_to_centroid.argtypes = [
        ctypes.c_int,  # n_genes
        ctypes.c_int,  # n_families
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # genes
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # centroids
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),    # gene_to_fam
        np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),  # distances
        ctypes.c_int   # d
    ]
    distance_to_centroid.restype = None
    return distance_to_centroid

def generate_gene_to_family_mapping(orthogroups_file, centroids_file, gene_expression_file, 
                                   use_all_species=True, target_species=None):
    """
    Generate gene_to_family mapping from Orthogroups.tsv file
    
    Parameters:
    -----------
    orthogroups_file : str
        Path to orthogroups TSV file
    centroids_file : str
        Path to centroids TSV file  
    gene_expression_file : str
        Path to gene expression TSV file
    use_all_species : bool
        Whether to use all species columns for mapping
    target_species : str or None
        Specific species to target (overrides use_all_species)
    
    Returns:
    --------
    dict : Dictionary with gene_to_family mapping, gene_expression, and centroids
    """
    print("Loading data files...")
    
    # Read files
    orthogroups = pd.read_csv(orthogroups_file, sep='\t')
    centroids = pd.read_csv(centroids_file, sep='\t')
    gene_expr = pd.read_csv(gene_expression_file, sep='\t')
    
    # Detect available species columns (all except "Orthogroup")
    species_columns = [col for col in orthogroups.columns if col != "Orthogroup"]
    
    # Determine which columns to use for gene mapping
    if target_species is not None:
        if target_species in species_columns:
            target_species_list = [target_species]
        else:
            raise ValueError(f"Error: species '{target_species}' not available. Available species: {', '.join(species_columns)}")
    elif use_all_species:
        target_species_list = species_columns

    
    # Extract gene list from expression file
    gene_id_col = gene_expr.columns[0]  # Use first column as gene_id
    genes_in_expression = gene_expr[gene_id_col].tolist()
    
    # Create gene_id to family index mapping
    gene_to_family = np.zeros(len(genes_in_expression), dtype=np.int32)  # 0 indicates no assignment
    gene_to_index = {gene: idx for idx, gene in enumerate(genes_in_expression)}
    
    # Create orthogroup to numeric index mapping
    orthogroup_col = centroids.columns[0]  # Use first column as Orthogroup
    orthogroup_to_index = {ortho: idx+1 for idx, ortho in enumerate(centroids[orthogroup_col])}  # 1-based for Fortran
    
    # Generate gene-to-family mapping
    for i, row in orthogroups.iterrows():
        orthogroup_id = row['Orthogroup']
        
        # Check if this orthogroup has a centroid
        if orthogroup_id not in orthogroup_to_index:
            continue
            
        family_index = orthogroup_to_index[orthogroup_id]
        
        # Iterate over all target species
        for species in target_species_list:
            if species not in row or pd.isna(row[species]) or row[species] == "":
                continue
                
            # Split genes by comma and clean spaces
            genes_list = [gene.strip() for gene in str(row[species]).split(',')]
            
            # Assign family index to each gene in expression file
            for gene_id in genes_list:
                if gene_id in gene_to_index:
                    gene_idx = gene_to_index[gene_id]
                    # Only assign if not previously assigned (avoid overwriting)
                    if gene_to_family[gene_idx] == 0:
                        gene_to_family[gene_idx] = family_index
    
    return {
        'gene_to_family': gene_to_family,
        'gene_expression': gene_expr,
        'centroids': centroids
    }

def call_distance_to_centroid(gene_expression_matrix, centroids_matrix, gene_to_family_vector):
    """
    Call Fortran distance_to_centroid function
    
    Parameters:
    -----------
    gene_expression_matrix : np.ndarray
        Gene expression matrix (d x n_genes) - column-major order
    centroids_matrix : np.ndarray
        Centroids matrix (d x n_families) - column-major order
    gene_to_family_vector : np.ndarray
        Gene-to-family mapping vector (1-based indices)
    
    Returns:
    --------
    np.ndarray : Distance values for each gene
    """
    # Verify dimensions
    d, n_genes = gene_expression_matrix.shape
    d_centroids, n_families = centroids_matrix.shape
    
    # Verify dimension compatibility
    if len(gene_to_family_vector) != n_genes:
        raise ValueError("Error: gene_to_family length doesn't match number of genes")
    
    if d_centroids != d:
        raise ValueError("Error: centroid dimensions don't match gene expression")
    
    # Prepare arrays for Fortran (C-contiguous)
    genes_flat = np.ascontiguousarray(gene_expression_matrix.ravel(order='F'))
    centroids_flat = np.ascontiguousarray(centroids_matrix.ravel(order='F'))
    gene_to_fam = np.ascontiguousarray(gene_to_family_vector, dtype=np.int32)
    distances = np.zeros(n_genes, dtype=np.float64)
    
    # Call Fortran function
    distance_to_centroid = setup_distance_to_centroid()
    distance_to_centroid(
        n_genes,
        n_families,
        genes_flat,
        centroids_flat,
        gene_to_fam,
        distances,
        d
    )
    
    return distances

def run_real_data_example():
    """
    Main function to run the complete example
    """
    print("=== EUCLIDEAN DISTANCE CALCULATION WITH REAL DATA ===\n")
    
    # Define file paths
    base_path = Path("material/")
    orthogroups_file = base_path / "filtered_families.tsv"
    centroids_file = base_path / "centroids_orthologs_only.tsv"
    gene_expression_file = base_path / "normalization.tsv"
    
    # Verify files exist
    for file_path in [orthogroups_file, centroids_file, gene_expression_file]:
        if not file_path.exists():
            raise FileNotFoundError(f"Error: File not found: {file_path}")
    
    # Generate mapping and load data
    mapping_data = generate_gene_to_family_mapping(
        str(orthogroups_file),
        str(centroids_file),
        str(gene_expression_file),
        use_all_species=True
    )
    
    # Prepare matrices for Fortran (column-major order)
    gene_expr_matrix = mapping_data['gene_expression'].iloc[:, 1:].values  # Exclude gene_id column
    gene_expr_matrix = np.asfortranarray(gene_expr_matrix.T)  # Transpose to d x n_genes
    
    centroids_matrix = mapping_data['centroids'].iloc[:, 1:].values  # Exclude Orthogroup column
    centroids_matrix = np.asfortranarray(centroids_matrix.T)  # Transpose to d x n_families
    
    print(f"Gene expression matrix shape: {gene_expr_matrix.shape}")
    print(f"Centroids matrix shape: {centroids_matrix.shape}")
    print(f"Gene-to-family vector length: {len(mapping_data['gene_to_family'])}")
    
    # Call Fortran function
    distances = call_distance_to_centroid(
        gene_expr_matrix, 
        centroids_matrix, 
        mapping_data['gene_to_family']
    )
    
    # Create results dataframe
    results_df = pd.DataFrame({
        'gene_id': mapping_data['gene_expression'].iloc[:, 0],
        'family_index': mapping_data['gene_to_family'],
        'distance_to_centroid': distances,
        'has_family': mapping_data['gene_to_family'] > 0
    })
    
    # Filter genes with assigned families
    results_with_families = results_df[results_df['has_family']].copy()
    
    # Save results
    output_file = "results/distance_to_centroids_python.tsv"
    results_with_families.to_csv(output_file, sep='\t', index=False)
    
    print(f"Results saved to: {output_file}")
    print(f"Genes with families: {len(results_with_families)}")
    
    return results_df

def test_basic_euclidean():
    """Test basic euclidean distance calculation"""
    print("=== TESTING BASIC EUCLIDEAN DISTANCE ===\n")
    
    # Test vectors
    vec1 = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    vec2 = np.array([4.0, 5.0, 6.0], dtype=np.float64)
    
    # Call Fortran function
    euclidean_distance = setup_euclidean_distance()
    result = ctypes.c_double()
    
    euclidean_distance(
        np.ascontiguousarray(vec1),
        np.ascontiguousarray(vec2),
        len(vec1),
        ctypes.byref(result)
    )
    
    # Expected result
    expected = np.sqrt(np.sum((vec1 - vec2)**2))
    
    print(f"Vector 1: {vec1}")
    print(f"Vector 2: {vec2}")
    print(f"Fortran result: {result.value:.6f}")
    print(f"Python expected: {expected:.6f}")
    print(f"Match: {abs(result.value - expected) < 1e-12}")
    
    return abs(result.value - expected) < 1e-12

if __name__ == "__main__":
    print("TENSOR-OMICS PYTHON EUCLIDEAN DISTANCE EXAMPLE")
    print("=" * 60)
    
    try:
        # Test basic functionality first
        if test_basic_euclidean():
            print("\n✓ Basic euclidean distance test passed")
        else:
            print("\n✗ Basic euclidean distance test failed")
            exit(1)
        
        # Run real data example
        print("\n" + "=" * 60)
        results = run_real_data_example()
        
        print("\n" + "=" * 60)
        print("EXAMPLE COMPLETED SUCCESSFULLY!")
        print(f"Processed {len(results)} genes total")
        print(f"Found {len(results[results['has_family']])} genes with family assignments")
        print("=" * 60)
        
    except Exception as e:
        print(f"\nERROR: {e}")
        print("Check file paths and library loading.")
        exit(1)
