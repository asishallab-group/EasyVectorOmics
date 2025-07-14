!> Module with normalization routines for tensor omics.
module tox_normalization
  use, intrinsic :: iso_fortran_env, only: real64
contains

  !> Normalizes each gene's expression vector using `sqrt(mean(x^2))`
  !> across tissues (not classical standard deviation).
  !>
  !> @param n_genes Number of genes (rows).<br>
  !> @param n_tissues Number of tissues (columns).<br>
  !> @param input_matrix Flattened matrix of gene expression values (column-major).<br>
  !> @param output_matrix Output normalized matrix (same shape as input).<br>
  pure subroutine normalize_by_std_dev(n_genes, n_tissues, input_matrix, output_matrix)
      implicit none

      ! Arguments
      integer, intent(in) :: n_genes, n_tissues
      real(real64), intent(in) :: input_matrix(n_genes * n_tissues)
      real(real64), intent(out) :: output_matrix(n_genes * n_tissues)

      ! Local variables
      integer :: i, j
      real(real64) :: std_dev, temp_sum

      ! Loop over each gene
      do i = 1, n_genes
          temp_sum = 0.0d0
          do j = 1, n_tissues
              temp_sum = temp_sum + input_matrix((j-1)*n_genes + i)**2
          end do

          std_dev = sqrt(temp_sum / dble(n_tissues))
          if (std_dev == 0.0d0) std_dev = 1.0d0

          do j = 1, n_tissues
              output_matrix((j-1)*n_genes + i) = input_matrix((j-1)*n_genes + i) / std_dev
          end do
      end do

      
  end subroutine normalize_by_std_dev

  !> Quantile normalization of a gene expression matrix (F42-compliant).
  !> Computes average expression per rank across tissues.
  !>
  !> @param n_genes       Number of genes (rows) <br>
  !> @param n_tissues     Number of tissues (columns)<br>
  !> @param input_matrix  Flattened input matrix (column-major)<br>
  !> @param output_matrix Flattened normalized output matrix<br>
  !> @param temp_col      Temporary vector for column sorting (size n_genes)<br>
  !> @param rank_means    Preallocated vector to store rank means (size n_genes)<br>
  !> @param perm          Permutation vector (size n_genes)<br>
  !> @param stack_left    Manual quicksort stack (≥ log2(n_genes) + 10)<br>
  !> @param stack_right   Manual quicksort stack (same size as stack_left)<br>
  !> @param max_stack     Stack size passed from R<br>
  pure subroutine quantile_normalization(n_genes, n_tissues, input_matrix, output_matrix, &
                                      temp_col, rank_means, perm, stack_left, stack_right, max_stack)
    use tox_sorting, only: sort_array

    implicit none

    integer, intent(in) :: n_genes         
    integer, intent(in) :: n_tissues       
    real(real64), intent(in) :: input_matrix(n_genes * n_tissues)  
    real(real64), intent(out) :: output_matrix(n_genes * n_tissues)  
    real(real64), intent(inout) :: temp_col(n_genes)      
    real(real64), intent(inout) :: rank_means(n_genes)    
    integer, intent(inout) :: perm(n_genes)         
    integer, intent(inout) :: stack_left(max_stack)     ! Add intent(inout)
    integer, intent(inout) :: stack_right(max_stack)    ! Add intent(inout)
    integer, intent(in) :: max_stack              

      ! Locals
      integer :: i, j

      ! Initialize rank means
      rank_means = 0.0d0

      ! === First pass: accumulate values by rank across tissues ===
      do j = 1, n_tissues
          ! Prepare current column and initialize permutation
          do i = 1, n_genes
              temp_col(i) = input_matrix((j - 1) * n_genes + i)
              perm(i) = i
          end do

          ! Sort current column with index tracking

          call sort_array(temp_col, perm, stack_left, stack_right)

          ! Accumulate values for each rank
          do i = 1, n_genes
              rank_means(i) = rank_means(i) + temp_col(perm(i))
          end do
      end do

      ! Average the rank values
      do i = 1, n_genes
          rank_means(i) = rank_means(i) / dble(n_tissues)
      end do

      ! === Second pass: assign averaged values by rank ===
      do j = 1, n_tissues
          ! Prepare column and reset permutation
          do i = 1, n_genes
              temp_col(i) = input_matrix((j - 1) * n_genes + i)
              perm(i) = i
          end do

          call sort_array(temp_col, perm, stack_left, stack_right)


          do i = 1, n_genes
              output_matrix((j - 1) * n_genes + perm(i)) = rank_means(i)
          end do
      end do
  end subroutine quantile_normalization

  !> Apply `log2(x + 1)` transformation to each element of the input matrix.
  !> This subroutine performs element-wise `log2(x + 1)` transformation on a
  !> matrix flattened in column-major order. The `log2` is computed via:
  !> `log(x + 1) / log(2)`, which is numerically equivalent and avoids the
  !> non-portable `log2` intrinsic for compatibility with WebAssembly (WASM).
  !>
  !> @param[in]  n_genes       Number of genes (rows)<br>
  !> @param[in]  n_tissues     Number of tissues (columns)<br>
  !> @param[in]  input_matrix  Flattened input matrix (size: n_genes * n_tissues)<br>
  !> @param[out] output_matrix Flattened output matrix (same size)<br>
  pure subroutine log2_transformation(n_genes, n_tissues, input_matrix, output_matrix)
      implicit none

      ! Arguments
      integer, intent(in) :: n_genes, n_tissues
      real(real64), intent(in) :: input_matrix(n_genes * n_tissues)
      real(real64), intent(out) :: output_matrix(n_genes * n_tissues)

      ! Locals
      integer :: i
      real(real64), parameter :: LOG2 = log(2.0d0)

      ! Loop through all elements in the flattened input matrix
      do i = 1, n_genes * n_tissues
          ! Apply the log2(x + 1) transformation
          output_matrix(i) = log(input_matrix(i) + 1.0d0) / LOG2
      end do
  end subroutine log2_transformation


  !> Calculate tissue averages by averaging replicates within each group.
  !> For each group of tissue replicates, this subroutine computes the average
  !> expression per gene. The input matrix is column-major, flattened as a 1D array.
  !>
  !> @param[in]  n_gene        Number of genes (rows)<br>
  !> @param[in]  n_grps        Number of tissue groups<br>
  !> @param[in]  group_s       Start column index for each group (length: n_grps)<br>
  !> @param[in]  group_c       Number of columns per group (length: n_grps)<br>
  !> @param[in]  input_matrix  Flattened input matrix (length: n_gene * n_col)<br>
  !> @param[out] output_matrix Flattened output matrix (length: n_gene * n_grps)<br>
  pure subroutine calc_tiss_avg(n_gene, n_grps, group_s, group_c, input_matrix, output_matrix)
      implicit none

      ! === Arguments ===
      integer, intent(in) :: n_gene, n_grps
      integer, intent(in) :: group_s(n_grps)
      integer, intent(in) :: group_c(n_grps)
      real(real64), intent(in) :: input_matrix(n_gene * sum(group_c))
      real(real64), intent(out) :: output_matrix(n_gene * n_grps)

      ! === Local variables ===
      integer :: i, j, g, col
      real(real64) :: sum_val
      integer :: start_idx, count_cols

      ! === Loop over each group ===
      do g = 1, n_grps
          start_idx = group_s(g)
          count_cols = group_c(g)

          do i = 1, n_gene
              sum_val = 0.0d0

              do j = 0, count_cols - 1
                  col = start_idx + j
                  sum_val = sum_val + input_matrix((col - 1) * n_gene + i)
              end do

              output_matrix((g - 1) * n_gene + i) = sum_val / dble(count_cols)
          end do
      end do
  end subroutine calc_tiss_avg


  !> Calculate `log2 fold changes` between condition and control columns.
  !> For each control-condition pair, this subroutine computes the `log2 fold change`
  !> by subtracting the expression value in the control column from the corresponding
  !> value in the condition column, for all genes.
  !>
  !> The input matrix must be column-major and flattened as a 1D array.
  !>
  !> @param n_genes        Number of genes (rows)<br>
  !> @param n_pairs        Number of condition-control column pairs<br>
  !> @param control_cols   Indices (1-based) of control columns (length n_pairs)<br>
  !> @param cond_cols      Indices (1-based) of condition columns (length n_pairs)<br>
  !> @param i_matrix       Input expression matrix, flattened (length: n_genes × N)<br>
  !> @param o_matrix       Output matrix for fold changes (length: n_genes × n_pairs)<br>
  pure subroutine calc_fchange(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix)

      implicit none

      ! === Arguments ===
      integer, intent(in) :: n_genes       !< Number of genes (rows)
      integer, intent(in) :: n_cols        !< Number of columns in the input matrix
      integer, intent(in) :: n_pairs       !< Number of control-condition pairs
      integer, intent(in) :: control_cols(n_pairs) !< Control column indices
      integer, intent(in) :: cond_cols(n_pairs)    !< Condition column indices
      real(real64), intent(in) :: i_matrix(n_genes * n_cols)     !< Input matrix (flattened)
      real(real64), intent(out) :: o_matrix(n_genes * n_pairs)   !< Output matrix (flattened)

      ! === Locals ===
      integer :: i, p
      integer :: control_col, cond_col

      ! === Loop over each pair ===
      do p = 1, n_pairs
          control_col = control_cols(p)
          cond_col    = cond_cols(p)

          do i = 1, n_genes
              o_matrix((p - 1) * n_genes + i) = &
                  i_matrix((cond_col - 1) * n_genes + i) - &
                  i_matrix((control_col - 1) * n_genes + i)
          end do
      end do

  end subroutine calc_fchange


end module tox_normalization


!> Wrapper for R and Fortran usage (R .Fortran)
subroutine normalize_by_std_dev_r(n_genes, n_tissues, input_matrix, output_matrix)
  use tox_normalization
  integer, intent(in) :: n_genes, n_tissues
  real(real64), intent(in)  :: input_matrix(n_genes, n_tissues)
  real(real64), intent(out) :: output_matrix(n_genes, n_tissues)
  call normalize_by_std_dev(n_genes, n_tissues, input_matrix, output_matrix)
end subroutine normalize_by_std_dev_r

!> C/Python interface (bind(C)), expects flat arrays.
subroutine normalize_by_std_dev_c(n_genes, n_tissues, input_matrix, output_matrix) bind(C, name="normalize_by_std_dev_c")
  use iso_c_binding
  use tox_normalization
  integer(c_int), value :: n_genes, n_tissues
  real(c_double), intent(in), target :: input_matrix(n_genes * n_tissues)
  real(c_double), intent(out), target :: output_matrix(n_genes * n_tissues)
  real(c_double), pointer :: inmat(:,:), outmat(:,:)
  call c_f_pointer(c_loc(input_matrix(1)), inmat, [n_genes, n_tissues])
  call c_f_pointer(c_loc(output_matrix(1)), outmat, [n_genes, n_tissues])
  call normalize_by_std_dev(n_genes, n_tissues, inmat, outmat)
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
  integer(c_int), intent(in), value :: max_stack
  real(c_double), intent(in), target :: input_matrix(n_genes, n_tissues)
  real(c_double), intent(out), target :: output_matrix(n_genes, n_tissues)
  real(c_double), intent(inout), target :: temp_col(n_genes)
  real(c_double), intent(inout), target :: rank_means(n_genes)
  integer(c_int), intent(inout), target :: perm(n_genes)
  integer(c_int), intent(inout), target :: stack_left(max_stack)
  integer(c_int), intent(inout), target :: stack_right(max_stack)

  call quantile_normalization(n_genes, n_tissues, input_matrix, output_matrix, &
                            temp_col, rank_means, perm, stack_left, stack_right, max_stack)
end subroutine quantile_normalization_c

subroutine log2_transformation_r(n_genes, n_tissues, input_matrix, output_matrix)
  use tox_normalization
  integer, intent(in) :: n_genes, n_tissues
  real(real64), intent(in) :: input_matrix(n_genes * n_tissues)
  real(real64), intent(out) :: output_matrix(n_genes * n_tissues)
  call log2_transformation(n_genes, n_tissues, input_matrix, output_matrix)
end subroutine log2_transformation_r

subroutine log2_transformation_c(n_genes, n_tissues, input_matrix, output_matrix) bind(C, name="log2_transformation_c")
  use iso_c_binding
  use tox_normalization
  integer(c_int), intent(in), value :: n_genes
  integer(c_int), intent(in), value :: n_tissues
  real(c_double), intent(in), target :: input_matrix(n_genes * n_tissues)
  real(c_double), intent(out), target :: output_matrix(n_genes * n_tissues)

  call log2_transformation(n_genes, n_tissues, input_matrix, output_matrix)
end subroutine log2_transformation_c

subroutine calc_tiss_avg_r(n_gene, n_grps, group_s, group_c, input_matrix, output_matrix)
  use tox_normalization
  integer, intent(in) :: n_gene, n_grps
  integer, intent(in) :: group_s(n_grps)
  integer, intent(in) :: group_c(n_grps)
  real(real64), intent(in) :: input_matrix(n_gene * sum(group_c))
  real(real64), intent(out) :: output_matrix(n_gene * n_grps)
  call calc_tiss_avg(n_gene, n_grps, group_s, group_c, input_matrix, output_matrix)
end subroutine calc_tiss_avg_r

subroutine calc_tiss_avg_c(n_gene, n_grps, group_s, group_c, input_matrix, output_matrix) bind(C, name="calc_tiss_avg_c")
  use iso_c_binding
  use tox_normalization
  integer(c_int), intent(in), value :: n_gene
  integer(c_int), intent(in), value :: n_grps
  integer(c_int), intent(in), target :: group_s(n_grps)
  integer(c_int), intent(in), target :: group_c(n_grps)
  real(c_double), intent(in), target :: input_matrix(n_gene * sum(group_c))
  real(c_double), intent(out) :: output_matrix(n_gene * n_grps)

  call calc_tiss_avg(n_gene, n_grps, group_s, group_c, input_matrix, output_matrix)
end subroutine calc_tiss_avg_c

subroutine calc_fchange_r(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix)
  use tox_normalization
  integer, intent(in) :: n_genes       !< Number of genes (rows)
  integer, intent(in) :: n_cols        !< Number of columns in the input matrix
  integer, intent(in) :: n_pairs       !< Number of control-condition pairs
  integer, intent(in) :: control_cols(n_pairs) !< Control column indices
  integer, intent(in) :: cond_cols(n_pairs)    !< Condition column indices
  real(real64), intent(in) :: i_matrix(n_genes * n_cols)     !< Input matrix (flattened)
  real(real64), intent(out) :: o_matrix(n_genes * n_pairs)   !< Output matrix (flattened)

  call calc_fchange(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix)
end subroutine calc_fchange_r

subroutine calc_fchange_c(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix) bind(C, name="calc_fchange_c")
  use iso_c_binding
  use tox_normalization
  integer(c_int), intent(in), value :: n_genes
  integer(c_int), intent(in), value :: n_cols
  integer(c_int), intent(in), value :: n_pairs
  integer(c_int), intent(in), target :: control_cols(n_pairs)
  integer(c_int), intent(in), target :: cond_cols(n_pairs)
  real(c_double), intent(in), target :: i_matrix(n_genes * n_cols)
  real(c_double), intent(out) :: o_matrix(n_genes * n_pairs) 

  call calc_fchange(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix)
end subroutine calc_fchange_c