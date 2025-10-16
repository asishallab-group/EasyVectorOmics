#include <Rcpp.h>

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
    double percentile,
    int* ierr
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

 }

// ===================================================================
// OUTLIER DETECTION WRAPPERS
// ===================================================================

// [[Rcpp::export]]
List tox_compute_family_scaling_rcpp(
    int n_genes, int n_families,
    NumericVector distances,
    IntegerVector gene_to_fam,
    NumericVector dscale,
    NumericVector loess_x,
    NumericVector loess_y,
    IntegerVector indices_used) {

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
List tox_compute_family_scaling_expert_rcpp(
    int n_genes, int n_families,
    NumericVector distances,
    IntegerVector gene_to_fam,
    NumericVector dscale,
    NumericVector loess_x,
    NumericVector loess_y,
    IntegerVector indices_used,
    IntegerVector perm_tmp,
    IntegerVector stack_left_tmp,
    IntegerVector stack_right_tmp,
    NumericVector family_distances) {

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
    Named("ierr") = ierr
  );
}

// [[Rcpp::export]]
List tox_compute_rdi_rcpp(
    int n_genes, int n_families,
    NumericVector distances,
    IntegerVector gene_to_fam,
    NumericVector dscale,
    NumericVector rdi,
    NumericVector sorted_rdi,
    IntegerVector perm,
    IntegerVector stack_left,
    IntegerVector stack_right) {

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
    Named("perm") = perm
  );
}

// [[Rcpp::export]]
List tox_identify_outliers_rcpp(
    int n_genes,
    NumericVector rdi,
    NumericVector sorted_rdi,
    IntegerVector is_outlier_int,
    double percentile) {

  double threshold = 0.0;
  int ierr = 0;

  identify_outliers_c(
    n_genes,
    rdi.begin(), sorted_rdi.begin(),
    is_outlier_int.begin(),
    &threshold,
    percentile,
    &ierr
  );

  return List::create(
    Named("is_outlier_int") = is_outlier_int,
    Named("threshold") = threshold,
    Named("ierr") = ierr
  );
}

// [[Rcpp::export]]
List tox_detect_outliers_rcpp(
    int n_genes, int n_families,
    NumericVector distances,
    IntegerVector gene_to_fam,
    NumericVector work_array,
    IntegerVector perm,
    IntegerVector stack_left,
    IntegerVector stack_right,
    IntegerVector is_outlier_int,
    NumericVector loess_x,
    NumericVector loess_y,
    IntegerVector loess_n,
    double percentile) {

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

  return List::create(
    Named("is_outlier_int") = is_outlier_int,
    Named("ierr") = ierr
  );
}

