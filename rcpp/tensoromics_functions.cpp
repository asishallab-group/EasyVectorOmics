#include <Rcpp.h>
#include <vector>

using namespace Rcpp;


// ===================================================================
// FORTRAN FUNCTIONS
// ===================================================================


extern "C" {

    void normalize_by_std_dev_c(int n_genes, int n_tissues,
                                double *input_matrix, double *output_matrix, int *ierr);                                    
    void quantile_normalization_c(int n_genes, int n_tissues, double *input_matrix, double *output_matrix,
                                double *temp_col, double *rank_means,
                                  int *perm, int *stack_left, int *stack_right,
                                  int max_stack, int *ierr);
    void log2_transformation_c(int n_genes, int n_tissues,
                               double *input_matrix, double *output_matrix, int *ierr); 

    void calc_tiss_avg_c(int n_genes, int n_grps,int *group_s, int *group_c,
                         double *input_matrix, double *output_matrix, int *ierr);  
    void calc_fchange_c(int n_genes, int n_cols, int n_pairs,int *control_cols, int *cond_cols,
                        double *input_matrix, double *output_matrix, int *ierr);

    void normalization_pipeline_c(int n_genes, int n_tissues,
                                  double *input_matrix, double *buf_stddev, double *buf_quant,
                                  double *buf_avg, double *buf_log, double *temp_col,
                                  double *rank_means, int *perm, int *stack_left,
                                  int *stack_right, int max_stack,
                                  int *group_s, int *group_c, int n_grps, int *ierr);

      void compute_family_scaling_c(
        int n_genes, int n_families,
        double* distances, int* gene_to_fam,
        double* dscale,
        double* loess_x, double* loess_y, int* indices_used,
        int* ierr
      );
      
      void compute_family_scaling_expert_c(
        int n_genes, int n_families,
        double* distances, int* gene_to_fam,
        double* dscale,
        double* loess_x, double* loess_y, int* indices_used,
        int* perm_tmp, int* stack_left_tmp, int* stack_right_tmp,
        double* family_distances,
        int* ierr
      );
      
      void compute_rdi_c(
        int n_genes, int n_families,
        double* distances, int* gene_to_fam,
        double* dscale,
        double* rdi, double* sorted_rdi,
        int* perm, int* stack_left, int* stack_right
      );
      
      void identify_outliers_c(
        int n_genes,
        double* rdi, double* sorted_rdi,
        int* is_outlier_int,
        double* threshold,
        double percentile
      );
      
      void detect_outliers_c(
        int n_genes, int n_families,
        double* distances, int* gene_to_fam,
        double* work_array,
        int* perm, int* stack_left, int* stack_right,
        int* is_outlier_int,
        double* loess_x, double* loess_y, int* loess_n,
        int* ierr,
        double percentile
      );

      
      void euclidean_distance_c(double* vec1, double* vec2, int d, double* result);

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
      int num_dimensions,
      int num_points,
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
        int* arr, int arr_size, 
        int* filename_ascii, int fn_len, int* ierr);

    void deserialize_real_C(
        double* arr, int arr_size, 
        int* filename_ascii, int fn_len, int* ierr);

    void deserialize_char_flat_C(
        int* ascii_arr, int clen, 
        int total_array_size, int* filename_ascii, int fn_len, int* ierr);

    void serialize_int_nd_C(
        void* arr, int* dims, 
        int ndim, int* filename_ascii, int fn_len, int* ierr);
    void serialize_real_nd_C(
        void* arr, int* dims, 
        int ndim, int* filename_ascii, int fn_len, int* ierr);

    void serialize_char_flat_C(
        int* ascii_arr, int* dims, 
        int ndim, int clen, int* filename_ascii, int fn_len, int* ierr);
    void get_array_metadata_C(
      const int* filename_ascii,
      int fn_len, int* dims_out, int* dims_out_capacity,
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

/**
 * Calculate tissue versatility and angles
 */


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

/**
 * Calculate normalization by standard deviation
 */
// [[Rcpp::export]]
List tox_normalize_by_std_dev_rcpp(NumericMatrix input) {
    int n_genes = input.nrow();
    int n_tissues = input.ncol();
    NumericMatrix output(n_genes, n_tissues);
    int ierr = 0;

    normalize_by_std_dev_c(n_genes, n_tissues, input.begin(), output.begin(), &ierr);

    return List::create(
        Named("output_vector") = output,
        Named("ierr") = ierr
    );
}

/**
 * Perform quantile normalization
 */
// [[Rcpp::export]]
List tox_quantile_normalization_rcpp(NumericMatrix input) {
    int n_genes = input.nrow();
    int n_tissues = input.ncol();

    int max_stack = static_cast<int>(std::ceil(std::log2(n_genes)) + 10);

    NumericMatrix output(n_genes, n_tissues);
    NumericVector temp_col(n_genes);
    NumericVector rank_means(n_genes);
    IntegerVector perm(n_genes);
    IntegerVector stack_left(max_stack);
    IntegerVector stack_right(max_stack);
    int ierr = 0;

    quantile_normalization_c(n_genes, n_tissues, input.begin(), output.begin(),
                             temp_col.begin(), rank_means.begin(), perm.begin(),
                             stack_left.begin(), stack_right.begin(), max_stack, &ierr);

    return List::create(
        Named("output_vector")     = output,
        Named("rank_means") = rank_means,
        Named("perm")       = perm,
        Named("ierr")       = ierr
    );
}

/**
 * Perform log2 transformation
 */
// [[Rcpp::export]]
List tox_log2_transformation_rcpp(NumericMatrix input) {
    int n_genes = input.nrow();
    int n_tissues = input.ncol();
    NumericMatrix output(n_genes, n_tissues);
    int ierr = 0;

    log2_transformation_c(n_genes, n_tissues, input.begin(), output.begin(), &ierr);

    return List::create(
        Named("output_vector") = output,
        Named("ierr") = ierr
    );
}

/**
 * Calculate tissue averages
 */
// [[Rcpp::export]]
List tox_calc_tiss_avg_rcpp(NumericMatrix input, IntegerVector group_s, IntegerVector group_c) {
    int n_gene = input.nrow();
    int n_grps = group_s.size();
    NumericMatrix output(n_gene, n_grps);
    int ierr = 0;

    calc_tiss_avg_c(n_gene, n_grps, group_s.begin(), group_c.begin(), input.begin(), output.begin(), &ierr);
    return List::create(
        Named("output_vector") = output,
        Named("ierr") = ierr
    );
}

/**
 * Calculate fold change
 */
// [[Rcpp::export]]
List tox_calc_fchange_rcpp(NumericMatrix input, IntegerVector control_cols, IntegerVector cond_cols) {
    int n_genes = input.nrow();
    int n_cols = input.ncol();
    int n_pairs = control_cols.size();
    NumericMatrix output(n_genes, n_pairs);
    int ierr = 0;

    calc_fchange_c(n_genes, n_cols, n_pairs, control_cols.begin(), cond_cols.begin(), input.begin(), output.begin(), &ierr);
    return List::create(
        Named("output_vector") = output,
        Named("ierr") = ierr
    );
}

/**
 * Perform normalization pipeline
 */
// [[Rcpp::export]]
List tox_normalization_pipeline_rcpp(NumericMatrix input, IntegerVector group_s, IntegerVector group_c) {
    int n_genes = input.nrow();
    int n_tissues = input.ncol();
    int n_grps = group_s.size();

   int max_stack = static_cast<int>(std::ceil(std::log2(n_genes)) + 10);

    NumericMatrix buf_stddev(n_genes, n_tissues);
    NumericMatrix buf_quant(n_genes, n_tissues);
    NumericMatrix buf_avg(n_genes, n_grps);
    NumericMatrix buf_log(n_genes, n_grps);
    NumericVector temp_col(n_genes);
    NumericVector rank_means(n_genes);
    IntegerVector perm(n_genes);
    IntegerVector stack_left(max_stack);
    IntegerVector stack_right(max_stack);
    int ierr = 0;


     normalization_pipeline_c(n_genes, n_tissues, input.begin(),
                                     buf_stddev.begin(), buf_quant.begin(),
                                     buf_avg.begin(), buf_log.begin(),
                                     temp_col.begin(), rank_means.begin(), perm.begin(),
                                     stack_left.begin(), stack_right.begin(), max_stack,
                                     group_s.begin(), group_c.begin(), n_grps, &ierr);

    return List::create(
        Named("buf_stddev") = buf_stddev,
        Named("buf_quant") = buf_quant,
        Named("buf_avg") = buf_avg,
        Named("buf_log") = buf_log,
        Named("rank_means") = rank_means,
        Named("perm") = perm,
        Named("ierr") = ierr
    );
}




// ===================================================================
// OUTLIER DETECTION WRAPPERS
// ===================================================================

// [[Rcpp::export]]
List tox_compute_family_scaling_rcpp(NumericVector distances, IntegerVector gene_to_fam, int n_families) {
  int n_genes = distances.size();
  NumericVector dscale(n_families);
  NumericVector loess_x(n_families);
  NumericVector loess_y(n_families);
  IntegerVector indices_used(n_families);
  int ierr = 0;

  compute_family_scaling_c(
    n_genes, n_families,
    distances.begin(),
    gene_to_fam.begin(),
    dscale.begin(),
    loess_x.begin(),
    loess_y.begin(),
    indices_used.begin(),
    &ierr
  );

  return List::create(
    Named("dscale") = dscale,
    Named("loess_x") = loess_x,
    Named("loess_y") = loess_y,
    Named("indices_used") = indices_used,
    Named("ierr") = ierr
  );
}

// [[Rcpp::export]]
List tox_compute_family_scaling_expert_rcpp(int n_families, NumericVector distances, IntegerVector gene_to_fam,
                                            IntegerVector perm_tmp, IntegerVector stack_left_tmp, IntegerVector stack_right_tmp,
                                            NumericVector family_distances) {
  int n_genes = distances.size();
  NumericVector dscale(n_families);
  NumericVector loess_x(n_families);
  NumericVector loess_y(n_families);
  IntegerVector indices_used(n_families);
  int ierr = 0;

  compute_family_scaling_expert_c(
    n_genes, n_families,
    distances.begin(),
    gene_to_fam.begin(),
    dscale.begin(),
    loess_x.begin(), loess_y.begin(), indices_used.begin(),
    perm_tmp.begin(), stack_left_tmp.begin(), stack_right_tmp.begin(),
    family_distances.begin(),
    &ierr
  );

  return List::create(
    Named("dscale") = dscale,
    Named("loess_x") = loess_x,
    Named("loess_y") = loess_y,
    Named("indices_used") = indices_used,
    Named("perm_tmp") = perm_tmp,
    Named("stack_left_tmp") = stack_left_tmp,
    Named("stack_right_tmp") = stack_right_tmp,
    Named("family_distances") = family_distances,
    Named("ierr") = ierr
  );
}

// [[Rcpp::export]]
List tox_compute_rdi_rcpp(NumericVector distances, IntegerVector gene_to_fam, NumericVector dscale) {
  int n_genes = distances.size();
  int n_families = dscale.size();
  NumericVector rdi(n_genes);
  NumericVector sorted_rdi(n_genes);
  IntegerVector perm(n_genes);
  IntegerVector stack_left(n_genes);
  IntegerVector stack_right(n_genes);

  compute_rdi_c(
    n_genes, n_families,
    distances.begin(),
    gene_to_fam.begin(),
    dscale.begin(),
    rdi.begin(), sorted_rdi.begin(),
    perm.begin(), stack_left.begin(), stack_right.begin()
  );
  return List::create(
    Named("rdi") = rdi,
    Named("sorted_rdi") = sorted_rdi,
    Named("perm") = perm,
    Named("stack_left") = stack_left,
    Named("stack_right") = stack_right
  );
}

// [[Rcpp::export]]
List tox_identify_outliers_rcpp(NumericVector rdi, double percentile) {
  int n_genes = rdi.size();
  NumericVector sorted_rdi = clone(rdi);
  std::sort(sorted_rdi.begin(), sorted_rdi.end());
  IntegerVector is_outlier_int(n_genes);
  double threshold = 0.0;

  identify_outliers_c(
    n_genes,
    rdi.begin(), sorted_rdi.begin(),
    is_outlier_int.begin(),
    &threshold,
    percentile
  );

  // Convert integer 0/1 flags to logical vector for R
  LogicalVector is_outlier(n_genes);
  for (int i = 0; i < n_genes; ++i) {
    is_outlier[i] = (is_outlier_int[i] != 0);
  }

  return List::create(
    Named("is_outlier") = is_outlier,
    Named("threshold") = threshold
  );
}

// [[Rcpp::export]]
List tox_detect_outliers_rcpp(NumericVector distances, IntegerVector gene_to_fam, int n_families, double percentile) {
  int n_genes = distances.size();
  NumericVector work_array(n_genes);
  IntegerVector perm(n_genes);
  IntegerVector stack_left(n_genes);
  IntegerVector stack_right(n_genes);
  IntegerVector is_outlier_int(n_genes);
  NumericVector loess_x(n_families);
  NumericVector loess_y(n_families);
  IntegerVector loess_n(n_families);
  int ierr = 0;

  detect_outliers_c(
    n_genes, n_families,
    distances.begin(), gene_to_fam.begin(),
    work_array.begin(),
    perm.begin(), stack_left.begin(), stack_right.begin(),
    is_outlier_int.begin(),
    loess_x.begin(), loess_y.begin(), loess_n.begin(),
    &ierr,
    percentile
  );

  // Convert integer flags to logical vector for R
  LogicalVector is_outlier(n_genes);
  for (int i = 0; i < n_genes; ++i) {
    is_outlier[i] = (is_outlier_int[i] != 0);
  }

  return List::create(
    Named("is_outlier") = is_outlier,
    Named("loess_x") = loess_x,
    Named("loess_y") = loess_y,
    Named("loess_n") = loess_n,
    Named("ierr") = ierr
  );
}


// [[Rcpp::export]]
List tox_build_kd_index_rcpp(NumericMatrix X, IntegerVector dim_order) {
  int d = X.nrow();
    int n = X.ncol();
    
    IntegerVector kd_ix(n);
    IntegerVector work(n);
    NumericVector subarray(n);
    IntegerVector perm(n);
    IntegerVector stack_left(n);
    IntegerVector stack_right(n);
    int ierr = 0;
    
  build_kd_index_C(X.begin(), d, n, kd_ix.begin(), dim_order.begin(), work.begin(), subarray.begin(), perm.begin(), stack_left.begin(), stack_right.begin(), &ierr);
        return List::create(
        Named("kd_ix") = kd_ix,
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


static std::vector<int> filename_to_ascii(const std::string &filename) {
  std::vector<int> v;
  for (unsigned char c : filename) v.push_back((int)c);
  return v;
}

// --- Serialize wrappers (match R signatures: arr, filename) ---

// [[Rcpp::export]]
List tox_serialize_int_array_rcpp(IntegerVector arr, std::string filename) {
  IntegerVector dim = arr.hasAttribute("dim") ? as<IntegerVector>(arr.attr("dim")) : IntegerVector::create((int)arr.size());
  int ndim = dim.size();
  std::vector<int> dims(ndim);
  for (int i = 0; i < ndim; ++i) dims[i] = dim[i];

  auto fname = filename_to_ascii(filename);
  int fn_len = static_cast<int>(fname.size());
  int ierr = 0;
  serialize_int_nd_C((void*)arr.begin(), dims.data(), ndim, fname.data(), fn_len, &ierr);
  return List::create(Named("ierr") = ierr);
}

// [[Rcpp::export]]
List tox_serialize_real_array_rcpp(NumericVector arr, std::string filename) {
  IntegerVector dim = arr.hasAttribute("dim") ? as<IntegerVector>(arr.attr("dim")) : IntegerVector::create((int)arr.size());
  int ndim = dim.size();
  std::vector<int> dims(ndim);
  for (int i = 0; i < ndim; ++i) dims[i] = dim[i];

  auto fname = filename_to_ascii(filename);
  int fn_len = static_cast<int>(fname.size());
  int ierr = 0;
  serialize_real_nd_C((void*)arr.begin(), dims.data(), ndim, fname.data(), fn_len, &ierr);
  return List::create(Named("ierr") = ierr);
}

// [[Rcpp::export]]
List tox_serialize_char_array_rcpp(CharacterVector carr, std::string filename) {
  IntegerVector dim = carr.hasAttribute("dim") ? as<IntegerVector>(carr.attr("dim")) : IntegerVector::create((int)carr.size());
  int ndim = dim.size();
  std::vector<int> dims(ndim);
  for (int i = 0; i < ndim; ++i) dims[i] = dim[i];

  int total = 1;
  for (int i = 0; i < ndim; ++i) total *= dims[i];

  int clen = 0;
  for (int i = 0; i < total; ++i) {
    std::string s = as<std::string>(carr[i]);
    if ((int)s.size() > clen) clen = (int)s.size();
  }
  if (clen == 0) clen = 1;

  std::vector<int> ascii_flat(total * clen, 0);
  for (int i = 0; i < total; ++i) {
    std::string s = as<std::string>(carr[i]);
    for (int j = 0; j < (int)s.size() && j < clen; ++j) ascii_flat[i*clen + j] = (unsigned char)s[j];
  }

  auto fname = filename_to_ascii(filename);
  int fn_len = static_cast<int>(fname.size());
  int ierr = 0;
  serialize_char_flat_C(ascii_flat.data(), dims.data(), ndim, clen, fname.data(), fn_len, &ierr);
  return List::create(Named("ierr") = ierr);
}

// --- Deserialize wrappers (match R signatures: filename, max_dims=5) ---

// [[Rcpp::export]]
List tox_deserialize_int_array_rcpp(std::string filename, int max_dims = 5) {
  auto fname = filename_to_ascii(filename);
  int fn_len = static_cast<int>(fname.size());
  std::vector<int> dims_out(max_dims);
  int ndims = 0;
  int ierr = 0;
  int clen = 0;
  get_array_metadata_C(fname.data(), fn_len, dims_out.data(), &max_dims, &ndims, &ierr, &clen);
  if (ierr != 0) return List::create(Named("ierr") = ierr);

  int total = 1;
  for (int i = 0; i < ndims; ++i) total *= dims_out[i];

  IntegerVector out(total);
  deserialize_int_C(out.begin(), total, fname.data(), fn_len, &ierr);
  return List::create(Named("values") = out, Named("dims") = IntegerVector(dims_out.begin(), dims_out.begin()+ndims), Named("ndim") = ndims, Named("ierr") = ierr);
}

// [[Rcpp::export]]
List tox_deserialize_real_array_rcpp(std::string filename, int max_dims = 5) {
  auto fname = filename_to_ascii(filename);
  int fn_len = static_cast<int>(fname.size());
  std::vector<int> dims_out(max_dims);
  int ndims = 0;
  int ierr = 0;
  int clen = 0;
  get_array_metadata_C(fname.data(), fn_len, dims_out.data(), &max_dims, &ndims, &ierr, &clen);
  if (ierr != 0) return List::create(Named("ierr") = ierr);

  int total = 1;
  for (int i = 0; i < ndims; ++i) total *= dims_out[i];

  NumericVector out(total);
  deserialize_real_C(out.begin(), total, fname.data(), fn_len, &ierr);
  return List::create(Named("values") = out, Named("dims") = IntegerVector(dims_out.begin(), dims_out.begin()+ndims), Named("ndim") = ndims, Named("ierr") = ierr);
}

// [[Rcpp::export]]
List tox_deserialize_char_array_rcpp(std::string filename, int max_dims = 5) {
  auto fname = filename_to_ascii(filename);
  int fn_len = static_cast<int>(fname.size());
  std::vector<int> dims_out(max_dims);
  int ndims = 0;
  int ierr = 0;
  int clen = 0;
  get_array_metadata_C(fname.data(), fn_len, dims_out.data(), &max_dims, &ndims, &ierr, &clen);
  if (ierr != 0) return List::create(Named("ierr") = ierr);

  int total = 1;
  for (int i = 0; i < ndims; ++i) total *= dims_out[i];

  std::vector<int> ascii_out(total * clen);
  deserialize_char_flat_C(ascii_out.data(), clen, total, fname.data(), fn_len, &ierr);

  CharacterVector out(total);
  for (int i = 0; i < total; ++i) {
    std::string s;
    for (int j = 0; j < clen; ++j) {
      int ch = ascii_out[i*clen + j];
      if (ch == 0) break;
      s.push_back((char)ch);
    }
    out[i] = s;
  }
  return List::create(Named("values") = out, Named("dims") = IntegerVector(dims_out.begin(), dims_out.begin()+ndims), Named("ndim") = ndims, Named("ierr") = ierr);
}

// [[Rcpp::export]]
List tox_get_array_metadata_rcpp(std::string filename, int dims_out_capacity = 5, bool with_clen = false) {
  auto fname = filename_to_ascii(filename);
  int fn_len = static_cast<int>(fname.size());
  IntegerVector dims_res(dims_out_capacity);
  int ndims = 0;
  int ierr = 0;
  int clen = 0;
  get_array_metadata_C(fname.data(), fn_len, dims_res.begin(), &dims_out_capacity, &ndims, &ierr, &clen);
  if (with_clen) {
    return List::create(Named("dims") = dims_res, Named("ndim") = ndims, Named("clen") = clen, Named("ierr") = ierr);
  }
  return List::create(Named("dims") = dims_res, Named("ndim") = ndims, Named("ierr") = ierr);
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
List tox_build_spherical_kd_rcpp(NumericMatrix V,IntegerVector dim_order) {
    int d = V.nrow();
    int n = V.ncol();
    
    IntegerVector sphere_ix(n);
    IntegerVector work(n);
    NumericVector subarray(n);
    IntegerVector perm(n);
    IntegerVector stack_left(n);
    IntegerVector stack_right(n);
    int ierr = 0;
    

    build_spherical_kd_C(V.begin(), d, n, sphere_ix.begin(), dim_order.begin(), work.begin(), subarray.begin(), perm.begin(), stack_left.begin(), stack_right.begin(), &ierr);

    return List::create(Named("sphere_ix") = sphere_ix,
                        Named("ierr") = ierr
    );}