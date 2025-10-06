#include <Rcpp.h>

using namespace Rcpp;


void check_error_code(int ierr) {
    if (ierr == 0) return;
    
    std::string msg;
    switch(ierr) {
        // I/O errors
        case 101: msg = "Could not open file."; break;
        case 102: msg = "Could not read magic number."; break;
        case 103: msg = "Could not read type code."; break;
        case 104: msg = "Could not read number of dimensions."; break;
        case 105: msg = "Could not read array dimensions"; break;
        case 106: msg = "Could not read character length."; break;
        case 107: msg = "Could not read array data."; break;
        case 112: msg = "Could not write magic number"; break;
        case 113: msg = "Could not write type code"; break;
        case 114: msg = "Could not write number of dimensions"; break;
        case 115: msg = "Could not write dimensions"; break;
        case 116: msg = "Could not write character length"; break;
        case 117: msg = "Could not write array data"; break;
        
        // Format errors
        case 200: msg = "Invalid format detected."; break;
        case 201: msg = "Invalid input provided."; break;
        case 202: msg = "Empty input arrays provided."; break;
        case 203: msg = "Dimension mismatch detected."; break;
        case 204: msg = "NaN or Inf found in input data."; break;
        case 205: msg = "Unsupported data type encountered."; break;
        case 206: msg = "Array size mismatch detected"; break;
        
        // Memory errors
        case 301: msg = "Memory allocation failed."; break;
        case 302: msg = "Null pointer reference encountered."; break;
        
        // Fortran runtime errors
        case 5002: msg = "Fortran runtime error: unit not open / not connected."; break;
        
        // Internal errors
        case 9001: msg = "Internal error: unexpected state."; break;
        case 9999: msg = "Unknown error."; break;
        
        default: msg = "Unmapped error code: " + std::to_string(ierr);
    }
    stop(msg);
}


// ===================================================================
// FORTRAN FUNCTIONS
// ===================================================================

extern "C" {
    void euclidean_distance_c(double* vec1, double* vec2, int d, double* result);
    void distance_to_centroid_c(int n_genes, int n_families, double* genes, 
                                double* centroids, int* gene_to_fam, 
                                double* distances, int d);
    void compute_tissue_versatility_c(int n_axes, int n_vectors, 
                                      double* expression_vectors, 
                                      int* exp_vecs_selection_index,
                                      int n_selected_vectors, 
                                      int* axes_selection, 
                                      int n_selected_axes,
                                      double* tissue_versatilities, 
                                      double* tissue_angles_deg,
                                      int* ierr);
}

/**
 * Calculate Euclidean distance between two vectors
 * 
 * Computes the Euclidean distance between two vectors of the same dimension.
 * This function automatically checks for errors and throws informative exceptions.
 * 
 * @param vec1 First vector (numeric - accepts both integer and real types)
 * @param vec2 Second vector (numeric, same length as vec1)
 * 
 * @return Numeric value representing the Euclidean distance between the vectors
 * 
 * @throws stop() if vectors are not numeric, have different lengths, or are empty
 */
// [[Rcpp::export]]
double tox_euclidean_distance(RObject vec1, RObject vec2) {
    // Check if inputs are numeric (accept both REALSXP and INTSXP)
    if ((TYPEOF(vec1) != REALSXP && TYPEOF(vec1) != INTSXP) || 
        (TYPEOF(vec2) != REALSXP && TYPEOF(vec2) != INTSXP)) {
        stop("Both vectors must be numeric");
    }
    
    // Direct conversion to NumericVector
    NumericVector v1 = as<NumericVector>(vec1);
    NumericVector v2 = as<NumericVector>(vec2);
    
    if (v1.length() != v2.length()) {
        stop("Vectors must have the same length");
    }
    if (v1.length() == 0) {
        stop("Vectors cannot be empty");
    }
    
    int d = v1.length();
    double result = 0.0;
    
    // Call Fortran function: d by value, vectors and result by reference
    euclidean_distance_c(v1.begin(), v2.begin(), d, &result);
    return result;
}

/**
 * Calculate distances from genes to their family centroids
 * 
 * Computes the Euclidean distance from each gene to its corresponding family centroid.
 * This function automatically checks for errors and throws informative exceptions.
 * 
 * @param genes Matrix of gene expression data (genes as columns, dimensions as rows)
 * @param centroids Matrix of family centroids (families as columns, dimensions as rows) 
 * @param gene_to_fam Integer vector mapping each gene to its family index (1-based)
 * @param d Integer number of dimensions
 * 
 * @return Numeric vector of distances from each gene to its family centroid
 *         (returns -1 for genes without valid family assignment)
 * 
 * @throws stop() if inputs are not numeric, dimensions don't match, or indices are invalid
 */
// [[Rcpp::export]]
NumericVector tox_distance_to_centroid(RObject genes, RObject centroids, 
                                       RObject gene_to_fam, RObject d) {
    // Input type validation
    if ((TYPEOF(genes) != REALSXP && TYPEOF(genes) != INTSXP) || 
        (TYPEOF(centroids) != REALSXP && TYPEOF(centroids) != INTSXP)) {
        stop("genes and centroids must be numeric");
    }
    if (TYPEOF(gene_to_fam) != INTSXP && TYPEOF(gene_to_fam) != REALSXP) {
        stop("gene_to_fam must be numeric or integer");
    }
    if (TYPEOF(d) != INTSXP && TYPEOF(d) != REALSXP) {
        stop("d must be numeric or integer");
    }
    
    // Convert to appropriate types
    NumericVector genes_vec = as<NumericVector>(genes);
    NumericVector centroids_vec = as<NumericVector>(centroids);
    IntegerVector gene_to_fam_vec = as<IntegerVector>(gene_to_fam);
    int d_val = as<int>(d);
    
    // Calculate dimensions
    int n_genes = genes_vec.length() / d_val;
    int n_families = centroids_vec.length() / d_val;
    
    // Validate dimensions
    if (genes_vec.length() % d_val != 0) {
        stop("Length of genes must be divisible by d");
    }
    if (centroids_vec.length() % d_val != 0) {
        stop("Length of centroids must be divisible by d");
    }
    if (gene_to_fam_vec.length() != n_genes) {
        stop("Length of gene_to_fam must equal number of genes");
    }
    if (is_true(any(gene_to_fam_vec < 0))) {
        stop("gene_to_fam indices must be between 0 and n_families (0 = no family assignment)");
    }
    
    // Prepare output array
    NumericVector distances(n_genes);
    
    // Call Fortran function
    distance_to_centroid_c(n_genes, n_families, genes_vec.begin(), 
                          centroids_vec.begin(), gene_to_fam_vec.begin(), 
                          distances.begin(), d_val);
    
    // Fortran returns -1 for genes without valid family assignment
    return distances;
}

/**
 * Calculate Tissue Versatility
 * 
 * Computes normalized tissue versatility for selected expression vectors.
 * The metric is based on the angle between each gene expression vector and the space diagonal.
 * Versatility is normalized to [0, 1], where 0 means uniform expression and 1 means expression in only one axis.
 * This function automatically checks for errors and throws informative exceptions.
 * 
 * @param expression_vectors Matrix where each column is a gene expression vector (n_axes x n_vectors)
 * @param vector_selection Logical vector indicating which vectors to process (length n_vectors)
 * @param axis_selection Logical vector indicating which axes to include in calculation (length n_axes)
 * 
 * @return List containing:
 *   \item{tissue_versatilities}{Normalized tissue versatility values [0,1] for selected vectors}
 *   \item{tissue_angles_deg}{Angles in degrees [0,90] for selected vectors}
 *   \item{n_selected_vectors}{Number of vectors processed}
 *   \item{n_selected_axes}{Number of axes used in calculation}
 * 
 * @throws stop() if inputs are invalid, dimensions don't match, or Fortran encounters an error
 */
// [[Rcpp::export]]
List tox_calculate_tissue_versatility(RObject expression_vectors, 
                                      RObject vector_selection, 
                                      RObject axis_selection) {
    // Input validation
    if (!Rf_isMatrix(expression_vectors)) {
        stop("expression_vectors must be a matrix");
    }
    if (TYPEOF(vector_selection) != LGLSXP && TYPEOF(vector_selection) != INTSXP && TYPEOF(vector_selection) != REALSXP) {
        stop("vector_selection must be logical or numeric");
    }
    if (TYPEOF(axis_selection) != LGLSXP && TYPEOF(axis_selection) != INTSXP && TYPEOF(axis_selection) != REALSXP) {
        stop("axis_selection must be logical or numeric");
    }
    
    // Convert to appropriate types
    NumericMatrix expr_matrix = as<NumericMatrix>(expression_vectors);
    IntegerVector vector_sel = as<IntegerVector>(vector_selection);
    IntegerVector axis_sel = as<IntegerVector>(axis_selection);
    
    // Dimensions and counts
    int n_axes = expr_matrix.nrow();
    int n_vectors = expr_matrix.ncol();
    int n_selected_vectors = sum(vector_sel);
    int n_selected_axes = sum(axis_sel);

    // Validate dimensions
    if (vector_sel.length() != n_vectors) {
        stop("vector_selection length must match number of columns in expression_vectors");
    }
    if (axis_sel.length() != n_axes) {
        stop("axis_selection length must match number of rows in expression_vectors");
    }
    
    // Prepare output arrays
    NumericVector tissue_versatilities(n_selected_vectors);
    NumericVector tissue_angles_deg(n_selected_vectors);
    int ierr = 0;
    
    // Call Fortran function
    compute_tissue_versatility_c(n_axes, n_vectors, expr_matrix.begin(),
                                vector_sel.begin(), n_selected_vectors,
                                axis_sel.begin(), n_selected_axes,
                                tissue_versatilities.begin(), 
                                tissue_angles_deg.begin(), &ierr);
    
    // Check for errors and throw informative messages
    check_error_code(ierr);
    
    // Return structured result (no ierr since we checked for errors)
    return List::create(
        Named("tissue_versatilities") = tissue_versatilities,
        Named("tissue_angles_deg") = tissue_angles_deg,
        Named("n_selected_vectors") = n_selected_vectors,
        Named("n_selected_axes") = n_selected_axes
    );
}