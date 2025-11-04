#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;

// ===================================================================
// FORTRAN FUNCTIONS
// ===================================================================

extern "C" {
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
List tox_euclidean_distance_rcpp(NumericVector vec1, NumericVector vec2) {
  int d = vec1.size();
  double result = 0.0;

  euclidean_distance_c(
    vec1.begin(), vec2.begin(),
    d,
    &result
  );

  return List::create(
    Named("distance") = result
  );
}

// [[Rcpp::export]]
List tox_distance_to_centroid_rcpp(NumericMatrix genes, NumericMatrix centroids, IntegerVector gene_to_fam, int d) {
  int n_genes = genes.ncol();
  int n_families = centroids.ncol();

  NumericVector distances(n_genes);

  distance_to_centroid_c(
    n_genes, n_families,
    genes.begin(), centroids.begin(),
    gene_to_fam.begin(),
    distances.begin(),
    d
  );

  return List::create(
    Named("distances") = distances
  );
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


