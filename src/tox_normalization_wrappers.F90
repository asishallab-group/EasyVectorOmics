! filepath: src/tox_normalization.F90
#include "precompiler_constants.F90"

!> Wrapper for R and Fortran usage (R .Fortran)
subroutine normalize_by_std_dev_r(n_genes, n_tissues, input_matrix, output_matrix)
  use tox_normalization
  integer, intent(in) :: n_genes, n_tissues
  real(8), intent(in)  :: input_matrix(n_genes, n_tissues)
  real(8), intent(out) :: output_matrix(n_genes, n_tissues)
  call normalize_by_std_dev_core(n_genes, n_tissues, input_matrix, output_matrix)
end subroutine normalize_by_std_dev_r

!> C/Python interface (bind(C)), expects flat arrays.
subroutine normalize_by_std_dev_c(n_genes, n_tissues, input_matrix, output_matrix) bind(C, name="normalize_by_std_dev_c")
  use iso_c_binding
  use tox_normalization
  integer(c_int), value :: n_genes, n_tissues
  real(c_double), intent(in), target :: input_matrix(*)
  real(c_double), intent(out), target :: output_matrix(*)
  real(c_double), pointer :: inmat(:,:), outmat(:,:)
  call c_f_pointer(c_loc(input_matrix(1)), inmat, [n_genes, n_tissues])
  call c_f_pointer(c_loc(output_matrix(1)), outmat, [n_genes, n_tissues])
  call normalize_by_std_dev_core(n_genes, n_tissues, inmat, outmat)
end subroutine normalize_by_std_dev_c

subroutine quantile_normalization_r(n_genes, n_tissues, input_matrix, output_matrix, &
                                        temp_col, rank_means, perm, stack_left, stack_right, max_stack)
  use tox_normalization
    integer, intent(in) :: n_genes, n_tissues, max_stack
  real(real64), intent(in)  :: input_matrix(n_genes, n_tissues)
  real(real64), intent(out) :: output_matrix(n_genes, n_tissues)
  real(real64), intent(inout) :: temp_col(n_genes)
  real(real64), intent(inout) :: rank_means(n_genes)
  integer, intent(inout) :: perm(n_genes)
  integer, intent(inout) :: stack_left(max_stack)
  integer, intent(inout) :: stack_right(max_stack)

  call quantile_normalization(n_genes, n_tissues, input_matrix, output_matrix, &
                            temp_col, rank_means, perm, stack_left, stack_right, max_stack)
end subroutine quantile_normalization_r

subroutine quantile_normalization_c(n_genes, n_tissues, input_matrix, output_matrix, &
                                    temp_col, rank_means, perm, stack_left, stack_right, max_stack) &
                                    bind(C, name="quantile_normalization_c")
  use iso_c_binding
  use tox_normalization
  integer(c_int), intent(in), value :: n_genes
  integer(c_int), intent(in), value :: n_tissues
  real(c_double), intent(in), target :: input_matrix(*)
  real(c_double), intent(out), target :: output_matrix(*)
  real(c_double), intent(inout), target :: temp_col(*)
  real(c_double), intent(inout), target :: rank_means(*)
  integer(c_int), intent(inout), target :: perm(*)
  integer(c_int), intent(inout), target :: stack_left(*)
  integer(c_int), intent(inout), target :: stack_right(*)
  integer(c_int), intent(in), value :: max_stack

  call quantile_normalization(n_genes, n_tissues, input_matrix, output_matrix, &
                            temp_col, rank_means, perm, stack_left, stack_right, max_stack)
end subroutine quantile_normalization_c

