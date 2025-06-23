!> Module with normalization routines for tensor omics.
module tox_normalization_mod
  use, intrinsic :: iso_fortran_env, only: real64
contains

  !> Normalizes each gene's expression vector by its root mean square (RMS) across tissues.
  !! This is not the classical standard deviation, but sqrt(mean(x^2)).
  !!
  !! @param n_genes      Number of genes (rows)
  !! @param n_tissues    Number of tissues (columns)
  !! @param input_matrix Input matrix (n_genes x n_tissues)
  !! @param output_matrix Output matrix (n_genes x n_tissues), normalized by RMS per gene
  subroutine normalize_by_std_dev_core(n_genes, n_tissues, input_matrix, output_matrix)
    integer, intent(in) :: n_genes, n_tissues
    real(real64), intent(in)  :: input_matrix(n_genes, n_tissues)
    real(real64), intent(out) :: output_matrix(n_genes, n_tissues)
    integer :: i, j
    real(real64) :: std_dev, temp_sum

    !$omp parallel do private(i, j, std_dev, temp_sum) schedule(static)
    do i = 1, n_genes
      temp_sum = 0.0_real64
      !$omp simd
      do j = 1, n_tissues
        temp_sum = temp_sum + input_matrix(i, j)**2
      end do

      std_dev = sqrt(temp_sum / real(n_tissues, real64))
      if (std_dev == 0.0_real64) std_dev = 1.0_real64

      !$omp simd
      do j = 1, n_tissues
        output_matrix(i, j) = input_matrix(i, j) / std_dev
      end do
    end do
    !$omp end parallel do
  end subroutine normalize_by_std_dev_core

  !> Performs quantile normalization across tissues for each gene.
  !! The algorithm sorts each gene's expression vector, computes the mean for each rank,
  !! and assigns the averaged value back according to the sorted ranks.
  !!
  !! Parallelization and vectorization follow project guidelines:
  !! - Only the outer loop (genes) is parallelized with OpenMP.
  !! - Inner loops (tissues) are vectorized with `!$omp simd`.
  !! - All temporary arrays are private to each thread.
  !!
  !! @param n_genes      Number of genes (rows)
  !! @param n_tissues    Number of tissues (columns)
  !! @param input_matrix Input matrix (n_genes x n_tissues)
  !! @param output_matrix Output matrix (n_genes x n_tissues), quantile normalized
  !! @param temp_col     Temporary column buffer (n_genes)
  !! @param rank_means   Buffer for mean per rank (n_genes)
  !! @param perm         Permutation buffer (n_genes)
  !! @param stack_left   Stack buffer for sorting (max_stack)
  !! @param stack_right  Stack buffer for sorting (max_stack)
  !! @param max_stack    Maximum stack size for sorting
  subroutine quantile_normalization(n_genes, n_tissues, input_matrix, output_matrix, &
                                    temp_col, rank_means, perm, stack_left, stack_right, max_stack)
    use, intrinsic :: iso_fortran_env, only: real64
    use tox_sorting, only: sort_array
    implicit none

    integer, intent(in) :: n_genes, n_tissues, max_stack
    real(real64), intent(in)  :: input_matrix(n_genes, n_tissues)
    real(real64), intent(out) :: output_matrix(n_genes, n_tissues)
    real(real64), intent(inout) :: temp_col(n_genes)
    real(real64), intent(inout) :: rank_means(n_genes)
    integer, intent(inout) :: perm(n_genes)
    integer, intent(inout) :: stack_left(max_stack)
    integer, intent(inout) :: stack_right(max_stack)

    integer :: i, j
    real(real64) :: local_temp_col(n_tissues)
    integer :: local_perm(n_tissues), local_stack_left(max_stack), local_stack_right(max_stack)

    !> Initialize rank means to zero
    !$omp simd
    do i = 1, n_genes
      rank_means(i) = 0.0_real64
    end do

    !> First pass: accumulate values by rank for each gene (parallel over genes)
    !$omp parallel do private(i, j, local_temp_col, local_perm, local_stack_left, local_stack_right) schedule(static)
    do i = 1, n_genes
      !> Prepare current gene's expression vector and permutation
      !$omp simd
      do j = 1, n_tissues
        local_temp_col(j) = input_matrix(i, j)
        local_perm(j) = j
      end do

      call sort_array(local_temp_col, local_perm, local_stack_left, local_stack_right)

      !> Accumulate values for each rank
      !$omp simd
      do j = 1, n_tissues
        rank_means(i) = rank_means(i) + local_temp_col(local_perm(j))
      end do
    end do
    !$omp end parallel do

    !> Average the rank values
    !$omp simd
    do i = 1, n_genes
      rank_means(i) = rank_means(i) / real(n_tissues, real64)
    end do

    !> Second pass: assign averaged values by rank (parallel over genes)
    !$omp parallel do private(i, j, local_temp_col, local_perm, local_stack_left, local_stack_right) schedule(static)
    do i = 1, n_genes
      !$omp simd
      do j = 1, n_tissues
        local_temp_col(j) = input_matrix(i, j)
        local_perm(j) = j
      end do

      call sort_array(local_temp_col, local_perm, local_stack_left, local_stack_right)

      !$omp simd
      do j = 1, n_tissues
        output_matrix(i, local_perm(j)) = rank_means(i)
      end do
    end do
    !$omp end parallel do

  end subroutine quantile_normalization

end module tox_normalization_mod