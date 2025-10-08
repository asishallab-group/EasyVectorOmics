#include <Rcpp.h>

using namespace Rcpp;


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
 */
// [[Rcpp::export]]
double tox_euclidean_distance_rcpp(NumericVector vec1, NumericVector vec2) {
    int d = vec1.length();
    double result = 0.0;
    
    euclidean_distance_c(vec1.begin(), vec2.begin(), d, &result);
    return result;
}

/**
 * Calculate distances from genes to their family centroids
 */
// [[Rcpp::export]]
NumericVector tox_distance_to_centroid_rcpp(NumericVector genes, NumericVector centroids, 
                                       IntegerVector gene_to_fam, int d) {
    int n_genes = genes.length() / d;
    int n_families = centroids.length() / d;
    
    NumericVector distances(n_genes);
    
    distance_to_centroid_c(n_genes, n_families, genes.begin(), 
                          centroids.begin(), gene_to_fam.begin(), 
                          distances.begin(), d);
    
    return distances;
}

/**
 * Calculate Tissue Versatility
 */
// [[Rcpp::export]]
List tox_calculate_tissue_versatility_rcpp(NumericMatrix expression_vectors, 
                                      IntegerVector vector_selection, 
                                      IntegerVector axis_selection) {
    int n_axes = expression_vectors.nrow();
    int n_vectors = expression_vectors.ncol();
    int n_selected_vectors = sum(vector_selection);
    int n_selected_axes = sum(axis_selection);
    
    NumericVector tissue_versatilities(n_selected_vectors);
    NumericVector tissue_angles_deg(n_selected_vectors);
    int ierr = 0;
    
    compute_tissue_versatility_c(n_axes, n_vectors, expression_vectors.begin(),
                                vector_selection.begin(), n_selected_vectors,
                                axis_selection.begin(), n_selected_axes,
                                tissue_versatilities.begin(), 
                                tissue_angles_deg.begin(), &ierr);
    
    return List::create(
        Named("tissue_versatilities") = tissue_versatilities,
        Named("tissue_angles_deg") = tissue_angles_deg,
        Named("n_selected_vectors") = n_selected_vectors,
        Named("n_selected_axes") = n_selected_axes,
        Named("ierr") = ierr
    );
}