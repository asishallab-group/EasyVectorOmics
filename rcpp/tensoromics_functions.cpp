#include <Rcpp.h>

using namespace Rcpp;


// ===================================================================
// FORTRAN FUNCTIONS
// ===================================================================


extern "C" {
void euclidean_distance_c(
    double* vec1,
    double* vec2,
    int d,
    double* result
);

    void distance_to_centroid_c(
        int n_genes, int n_families,
        double* genes, double* centroids,
        int* gene_to_fam,
        double* distances,
        int d
    );

    void compute_tissue_versatility_c(
        int n_axes, int n_vectors,
        double* expression_vectors,
        int* exp_vecs_selection_index,
        int n_selected_vectors,
        int* axes_selection,
        int n_selected_axes,
        double* tissue_versatilities,
        double* tissue_angles_deg,
        int* ierr
    );
                      
    void build_kd_index_C(
        double* points,
        int* num_dimensions,
        int* num_points,
        int* kd_indices,
        int* dimension_order,
        int* workspace,
        double* value_buffer,
        int* permutation,
        int* left_stack,
        int* right_stack,
        int* ierr
    );   
void which_c(
        int* mask,
        int n,
        int* idx_out,
        int m_max,
        int* m_out,
        int* ierr
);
void loess_smooth_2d_c(
        int n_total,
        int n_target,
        double* x_ref,
        double* y_ref,
        int* indices_used,
        int n_used,
        double* x_query,
        double kernel_sigma,
        double kernel_cutoff,
        double* y_out,
        int* ierr
);
     
void deserialize_int_C(
    int* arr_out, int arr_size, 
    int* filename_ascii, int fn_len, int* ierr);

void deserialize_real_C(
    double* arr_out, int arr_size, 
    int* filename_ascii, int fn_len, int* ierr);

void deserialize_char_flat_C(
    int* ascii_arr_out, int clen, 
    int total_array_size, int* filename_ascii, int fn_len, int* ierr);

void serialize_int_nd_C(
    void* arr_ptr, int* dims, 
    int ndim, int* filename_ascii, int fn_len, int* ierr);
void serialize_real_nd_C(
    void* arr_ptr, int* dims, 
    int ndim, int* filename_ascii, int fn_len, int* ierr);

void serialize_char_flat_C(
    int* ascii_arr, int* dims, 
    int ndim, int clen, int* filename_ascii, int fn_len, int* ierr);

void get_array_metadata_C(
    const int* filename_ascii, 
    int fn_len, int* dims_out, int dims_out_capacity, 
    int* ndims, int* ierr, int* clen);

void build_bst_index_C(
    const double* values, 
    int num_values, int* sorted_indices, int* left_stack, 
    int* right_stack, int* ierr);

void bst_range_query_C(
    const double* values, 
    const int* sorted_indices, int num_values, 
    double lower_bound, double upper_bound, int* output_indices, 
    int* num_matches, int* ierr);

void build_spherical_kd_C(
    const double* vectors, int num_dimensions, 
    int num_vectors, int* sphere_indices, int* dimension_order, 
    int* workspace, double* value_buffer, int* permutation, 
    int* left_stack, int* right_stack, int* ierr);
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


// [[Rcpp::export]]
List tox_calculate_tissue_versatility_rcpp(NumericMatrix expression_vectors, IntegerVector exp_vecs_selection_index, IntegerVector axes_selection) {
    int n_axes = expression_vectors.nrow();
    int n_vectors = expression_vectors.ncol();

    // Count selected vectors / axes (exp_vecs_selection_index & axes_selection are 0/1 ints)
    int n_selected_vectors = 0;
    for (int i = 0; i < exp_vecs_selection_index.size(); ++i) {
        if (exp_vecs_selection_index[i] != 0) ++n_selected_vectors;
    }
    int n_selected_axes = 0;
    for (int i = 0; i < axes_selection.size(); ++i) {
        if (axes_selection[i] != 0) ++n_selected_axes;
    }

    NumericVector tissue_versatilities(n_selected_vectors);
    NumericVector tissue_angles_deg(n_selected_vectors);
    int ierr = 0;

    compute_tissue_versatility_c(
        n_axes, n_vectors,
        expression_vectors.begin(),
        exp_vecs_selection_index.begin(),
        n_selected_vectors,
        axes_selection.begin(),
        n_selected_axes,
        tissue_versatilities.begin(),
        tissue_angles_deg.begin(),
        &ierr
    );

    return List::create(
        Named("tissue_versatilities") = tissue_versatilities,
        Named("tissue_angles_deg") = tissue_angles_deg,
        Named("n_selected_vectors") = n_selected_vectors,
        Named("n_selected_axes") = n_selected_axes,
        Named("ierr") = ierr
    );
}




// [[Rcpp::export]]
List tox_build_kd_index_rcpp(NumericMatrix points,
int num_dimensions,
int num_points,
IntegerVector dimension_order) {
        IntegerVector kd_indices(num_points);
        IntegerVector workspace(num_points);
        NumericVector value_buffer(num_points);
        IntegerVector permutation(num_points);
        IntegerVector left_stack(num_points);
        IntegerVector right_stack(num_points);
        int ierr = 0;
        build_kd_index_C(points.begin(), &num_dimensions, &num_points, kd_indices.begin(), dimension_order.begin(), workspace.begin(), value_buffer.begin(), permutation.begin(), left_stack.begin(), right_stack.begin(), &ierr);
        return List::create(
        Named("kd_indices") = kd_indices,
        Named("workspace") = workspace,
        Named("value_buffer") = value_buffer,
        Named("permutation") = permutation,
        Named("left_stack") = left_stack,
        Named("right_stack") = right_stack,
        Named("ierr") = ierr
        );
}

// [[Rcpp::export]]
List tox_which_rcpp(IntegerVector mask,int n,int m_max) {

        IntegerVector idx_out(m_max);
        int m_out = 0;
        int ierr = 0;
        which_c(mask.begin(), n, idx_out.begin(), m_max, &m_out, &ierr);
        return List::create(
        Named("idx_out") = idx_out,
        Named("m_out") = m_out,
        Named("ierr") = ierr
        );
}

// [[Rcpp::export]]
List tox_loess_smooth_2d_rcpp(int n_total,
int n_target,
NumericVector x_ref,
NumericVector y_ref,
IntegerVector indices_used,
int n_used,
NumericVector x_query,
double kernel_sigma,
double kernel_cutoff) {
        NumericVector y_out(n_target);
        int ierr = 0;
        loess_smooth_2d_c(n_total, n_target, x_ref.begin(), y_ref.begin(), indices_used.begin(), n_used, x_query.begin(), kernel_sigma, kernel_cutoff, y_out.begin(), &ierr);
        return List::create(
        Named("y_out") = y_out,
        Named("ierr") = ierr
        );
}



// [[Rcpp::export]]
int tox_deserialize_char_flat_rcpp(IntegerVector ascii_out, int clen, int total, IntegerVector filename_ascii, int fn_len) {
    int ierr = 0;
    deserialize_char_flat_C(ascii_out.begin(), clen, total, filename_ascii.begin(), fn_len, &ierr);
    return ierr;
}


// [[Rcpp::export]]
int tox_deserialize_int_rcpp(IntegerVector arr_out, IntegerVector filename_ascii, int fn_len) {
    int ierr = 0;
    int total = arr_out.size();
    deserialize_int_C(arr_out.begin(), total, filename_ascii.begin(), fn_len, &ierr);
    return ierr;
}


// [[Rcpp::export]]
int tox_deserialize_real_rcpp(NumericVector arr_out, IntegerVector filename_ascii, int fn_len) {
    int ierr = 0;
    int total = arr_out.size();
    deserialize_real_C(arr_out.begin(), total, filename_ascii.begin(), fn_len, &ierr);
    return ierr;
}


// [[Rcpp::export]]
int tox_serialize_int_nd_rcpp(IntegerVector arr, IntegerVector dims, int ndim, IntegerVector filename_ascii, int fn_len) {
    int ierr = 0;
    serialize_int_nd_C((void*)arr.begin(), dims.begin(), ndim, filename_ascii.begin(), fn_len, &ierr);
    return ierr;
}


// [[Rcpp::export]]
int tox_serialize_real_nd_rcpp(NumericVector arr, IntegerVector dims, int ndim, IntegerVector filename_ascii, int fn_len) {
    int ierr = 0;
    serialize_real_nd_C((void*)arr.begin(), dims.begin(), ndim, filename_ascii.begin(), fn_len, &ierr);
    return ierr;
}


// [[Rcpp::export]]
int tox_serialize_char_flat_rcpp(IntegerVector ascii_arr, IntegerVector dims, int ndim, int clen, IntegerVector filename_ascii, int fn_len) {
    int ierr = 0;
    serialize_char_flat_C(ascii_arr.begin(), dims.begin(), ndim, clen, filename_ascii.begin(), fn_len, &ierr);
    return ierr;
}
  
// [[Rcpp::export]]
List get_array_metadata_rcpp(IntegerVector filename_ascii, int fn_len, int dims_out_capacity = 5, bool with_clen = false) {
    IntegerVector dims_res(dims_out_capacity);
    int ndims = 0;
    int ierr = 0;
    int clen = 0;

    get_array_metadata_C(filename_ascii.begin(), fn_len, dims_res.begin(), dims_out_capacity, &ndims, &ierr, &clen);

    if (with_clen) {
        return List::create(Named("dims") = dims_res,
                            Named("ndim") = ndims,
                            Named("clen") = clen,
                            Named("ierr") = ierr);
    }

    return List::create(Named("dims") = dims_res,
                        Named("ndim") = ndims,
                        Named("ierr") = ierr);
}

// [[Rcpp::export]]
IntegerVector build_bst_index_rcpp(NumericVector values) {
    int num_values = static_cast<int>(values.size());
    IntegerVector sorted_indices(num_values);
    IntegerVector left_stack(num_values);
    IntegerVector right_stack(num_values);
    int ierr = 0;

    build_bst_index_C(values.begin(), num_values, sorted_indices.begin(), left_stack.begin(), right_stack.begin(), &ierr);

    // Return sorted indices (1-based Fortran indices preserved)
    return sorted_indices;
}


// [[Rcpp::export]]
List bst_range_query_rcpp(NumericVector values, IntegerVector sorted_indices, double lower_bound, double upper_bound) {
    int num_values = static_cast<int>(values.size());
    IntegerVector output_indices(num_values);
    int num_matches = 0;
    int ierr = 0;

    bst_range_query_C(values.begin(), sorted_indices.begin(), num_values, lower_bound, upper_bound, output_indices.begin(), &num_matches, &ierr);

    return List::create(Named("output_indices") = output_indices,
                                            Named("num_matches") = num_matches,
                                            Named("ierr") = ierr);
}

// [[Rcpp::export]]
List build_spherical_kd_rcpp(NumericMatrix vectors) {
    int num_dimensions = static_cast<int>(vectors.nrow());
    int num_vectors = static_cast<int>(vectors.ncol());

    IntegerVector sphere_indices(num_vectors);
    IntegerVector dimension_order(num_dimensions);
    IntegerVector workspace(num_vectors);
    NumericVector value_buffer(num_vectors);
    IntegerVector permutation(num_vectors);
    IntegerVector left_stack(num_vectors);
    IntegerVector right_stack(num_vectors);
    int ierr = 0;

    build_spherical_kd_C(vectors.begin(), num_dimensions, num_vectors, sphere_indices.begin(), dimension_order.begin(), workspace.begin(), value_buffer.begin(), permutation.begin(), left_stack.begin(), right_stack.begin(), &ierr);

    return List::create(Named("sphere_indices") = sphere_indices,
                                            Named("dimension_order") = dimension_order,
                                            Named("workspace") = workspace,
                                            Named("value_buffer") = value_buffer,
                                            Named("permutation") = permutation,
                                            Named("left_stack") = left_stack,
                                            Named("right_stack") = right_stack,
                                            Named("ierr") = ierr);
}