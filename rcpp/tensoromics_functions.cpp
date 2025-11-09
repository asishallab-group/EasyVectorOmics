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
    void compute_shift_vector_field_c(int d, int n_genes, int n_families,
                                      double* expression_vectors, double* family_centroids,
                                      int* gene_to_centroid, double* shift_vectors,
                                      int* ierr);

    void mean_vector_c(double* expression_vectors, int n_axes, int n_genes,
                       int* gene_indices, int n_selected_genes,
                       double* centroid_col, int* ierr);

    void group_centroid_c(double* expression_vectors, int n_axes, int n_genes,
                         int* gene_to_family, int n_families,
                         double* centroid_matrix, int use_all_mode,
                         int* ortholog_set, int* selected_indices, int selected_indices_len, int* ierr);
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

// [[Rcpp::export]]
List tox_compute_shift_vector_field_rcpp(NumericMatrix expression_vectors, NumericMatrix family_centroids, IntegerVector gene_to_centroid) {
    int n_axes_genes = expression_vectors.nrow();
    int n_vectors = expression_vectors.ncol();
    int n_axes_centroids = family_centroids.nrow();
    int n_families = family_centroids.ncol();

    NumericMatrix shift_vectors(2 * n_axes_genes, n_vectors);
    int ierr = 0;

    compute_shift_vector_field_c(n_axes_genes, n_vectors, n_families,
                                 expression_vectors.begin(), family_centroids.begin(),
                                 gene_to_centroid.begin(), shift_vectors.begin(), &ierr);

    NumericVector flat(shift_vectors.begin(), shift_vectors.end());

    return List::create(
        Named("shift_vectors") = flat,
        Named("ierr") = ierr
    );
}

// [[Rcpp::export]]
List tox_mean_vector_rcpp(NumericMatrix expression_vectors, IntegerVector gene_indices) {
    int n_axes = expression_vectors.nrow();
    int n_genes = expression_vectors.ncol();
    int n_selected_genes = gene_indices.length();

    NumericVector centroid_col(n_axes);
    int ierr = 0;

    mean_vector_c(expression_vectors.begin(), n_axes, n_genes,
                  gene_indices.begin(), n_selected_genes,
                  centroid_col.begin(), &ierr);

    return List::create(
        Named("centroid_col") = centroid_col,
        Named("ierr") = ierr
    );
}

// [[Rcpp::export]]
List tox_group_centroid_rcpp(NumericMatrix expression_vectors, IntegerVector gene_to_family, int n_families, IntegerVector ortholog_set, int use_all_mode) {
    int n_axes = expression_vectors.nrow();
    int n_genes = expression_vectors.ncol();

    NumericMatrix centroid_matrix(n_axes, n_families);
    IntegerVector selected_indices(n_genes);
    int selected_indices_len = n_genes;
    int ierr = 0;

        group_centroid_c(expression_vectors.begin(), n_axes, n_genes,
                         gene_to_family.begin(), n_families,
                         centroid_matrix.begin(), use_all_mode,
                         ortholog_set.begin(), selected_indices.begin(), selected_indices_len, &ierr);

    return List::create(
        Named("centroid_matrix") = centroid_matrix,
        Named("selected_indices") = selected_indices,
        Named("ierr") = ierr
    );
}
