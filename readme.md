# EasyVectorOmics Analysis Pipeline

## Workflow

1. **Data Normalization**  
   Script: `normalizer.py`  
   Input: A TSV file with gene expression values in the format:  
   GeneID | Tissue_X | Tissue_Y | Tissue_Z | ...  
   Output: A normalized TPM expression matrix

2. **Protein Similarity Calculation**  
   Script: `harmonic_mean_calc.py`  
   Input: The combined BLAST output file from OrthoFinder2 substitude back with original protein ids.
   Output: A BLAST file with an added column containing the harmonic mean similarity between protein pairs.

3. **Gene Classification via Phylogeny**  
   Scripts: `process_tree.py`, `tree_rest.py`  
   Inputs for `process_tree.py`:  
   - OrthoFinder Orthologues folder  
   - `harmonic_mean.pkl`  
   - `Duplications.tsv`  
   - Resolved gene trees folder  
   - `Orthogroups.tsv`  

   `process_tree.py` analyzes resolved gene trees to classify gene relationships into conserved orthologs, inparalogs, outparalogs, special inparalogs, and special outparalogs.  
   `tree_rest.py` handles orthogroups not represented in the resolved trees and assigns classifications using remaining input data.

4. **Tandem Gene Classification for Synteny Analysis**  
   Script: `TandemsPaper.py`  
   Input:  
   - All-vs-all BLASTn results  
   - Combined species GTF file   
   Output: A TSV file listing tandem gene clusters, their species, and the best representative gene for each cluster.

5. **Neighborhood-Based Synteny Classification**  
   Script: `neighborhood.py`  
   Input:  
   - Species-specific GTF files  
   - Tandems file from `TandemsPaper.py`  
   - All-vs-all BLASTn results  
   Output: A TSV file listing only source-to-source gene relationships.

6. **Gene Relationship Completion Based on Synteny**  
   Script: `classify_genes.py`  
   Input:  
   - `Orthogroups.tsv` from OrthoFinder  
   - Results from the neighborhood synteny classification  
   Output: TSV files classifying remaining genes into conserved orthologs, inparalogs, outparalogs, special inparalogs, and special outparalogs.

7. **Filtering of Gene Families**  
   Script: `result_filtering.py`  
   Input:  
   - All TSV files from `classify_genes.py`, `neighborhood.py`, and/or `process_tree.py`  
   Output: Filtered TSV files including only gene families with at least four unique conserved orthologs and copy gene relationships.

8. **Centroid Calculation for Expression Vectors**  
   Script: `get_centroids.R`  
   Input:  
   - Normalized TPM value file  
   - Conserved orthologs TSV file  
   Output: A TSV file containing orthogroup IDs and their expression centroids.

9. **Distance Calculation to Centroid**  
   Script: `calculate_distances.R`  
   Input:  
   - `centroid.tsv`  
   - TPM value file  
   - Conserved ortholog list  
   - Combined file of special inparalogs and outparalogs  
   Output:  
   - `distances.tsv`: Euclidean distance of each gene to its group's centroid  
   - Source-copy pairwise distance TSV  
   - Outlier classification (top 5%)

10. **Tandem Annotation for Expression Analysis**  
    Script: `distances_add_tandems.py`  
    Input:  
    - `Tandems.tsv` from `TandemsPaper.py`  
    - `distances.tsv` from `calculate_distances.R`  
    - `map.tsv` mapping Gene to Protein IDs  
    Output: Updated `distances.tsv` annotated with tandem gene information.

11. **Evolutionary Angle Calculations**  
    Script: `calculate_angles_optimized.py`  
    Input:  
    - Modified `distances.tsv` with tandem annotations  
    - `centroid.tsv`  
    - Pairwise distance TSV  
    Output: The input TSV files updated with angle calculations

