#include "macros.h"
!> # Global Jensen-Shannon-Divergence Compatibility Test (gJCT)
!!
!! This module implements the four main subroutines for the gJCT algorithm.
module tox_jenson_shannon_test
  use safeguard
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
  use f42_utils, only: heapsort_real, calc_percentile
  use tox_errors, only: validate_all_in_range_real, validate_in_range_int, is_err, set_ok
  implicit none
  
  private
  public :: compute_gene_means, compute_gene_means_helper, compute_residuals, compute_residuals_helper, pool_means_alloc, construct_neighborhoods, construct_neighborhoods_helper
  
contains

  pure subroutine compute_gene_means(n_genes, n_reps, expr, means, ierr)
    integer(int32), intent(in) :: n_genes
    !! Number of genes in the study
    integer(int32), intent(in) :: n_reps
    !! Number of biological replicates in the study
    real(real64), intent(in) :: expr(n_reps, n_genes)
    !! Expression matrix
    real(real64), intent(out) :: means(n_genes)
    !! Per-gene mean expression values
    integer(int32), intent(out) :: ierr
    !! Error code

    call set_ok(ierr)
    
    call validate_in_range_int(n_genes, ierr, min = 1_int32)
    call validate_in_range_int(n_reps, ierr, min = 1_int32)
    !! expression can contain NaN
    if(is_err(ierr)) return

    call compute_gene_means_helper(n_genes, n_reps, expr, means)
  end subroutine compute_gene_means
  
  !> Compute per-gene mean expression
  pure subroutine compute_gene_means_helper(n_genes, n_reps, expr, means)
    integer(int32), intent(in) :: n_genes
    !! Number of genes in the study
    integer(int32), intent(in) :: n_reps
    !! Number of biological replicates in the study
    real(real64), intent(in) :: expr(n_reps, n_genes)
    !! Expression matrix
    real(real64), intent(out) :: means(n_genes)
    !! Per-gene mean expression values
    
    integer(int32) :: g, i, valid_count
    real(real64) :: sum_val
    
    ! Use do concurrent for parallelization across genes
    do concurrent (g = 1:n_genes)
      sum_val = 0.0_real64
      valid_count = 0
      
      ! Count valid (non-NaN) replicates and compute sum
      do i = 1, n_reps
        if (.not. ieee_is_nan(expr(i, g))) then
          sum_val = sum_val + expr(i, g)
          valid_count = valid_count + 1
        end if
      end do
      
      ! Compute mean or set to NaN if no valid values
      if (valid_count > 0) then
        means(g) = sum_val / real(valid_count, real64)
      else
        means(g) = ieee_value(means(g), ieee_quiet_nan)
      end if
    end do
    
  end subroutine compute_gene_means_helper

  !> Compute signed residuals
  pure subroutine compute_residuals(n_genes, n_reps, expr, means, resid, ierr)
    integer(int32), intent(in) :: n_genes
    !! Number of genes in the study
    integer(int32), intent(in) :: n_reps
    !! Number of biological replicates in the study
    real(real64), intent(in) :: expr(n_reps, n_genes)
    !! Expression matrix containing
    real(real64), intent(in) :: means(n_genes)
    !! Per-gene mean expression values
    real(real64), intent(out) :: resid(n_reps, n_genes)
    !! Matrix of signed residuals
    integer(int32), intent(out) :: ierr
    !! Error code

    call set_ok(ierr)
    call validate_in_range_int(n_genes, ierr, min = 1_int32)
    call validate_in_range_int(n_reps, ierr, min = 1_int32)
    ! family means and expr containing NaN is expected behaviour
    if(is_err(ierr)) return

    call compute_residuals_helper(n_genes, n_reps, expr, means, resid)
  end subroutine compute_residuals

  !> Compute signed residuals helper
  pure subroutine compute_residuals_helper(n_genes, n_reps, expr, means, resid)
    integer(int32), intent(in) :: n_genes
    !! Number of genes in the study
    integer(int32), intent(in) :: n_reps
    !! Number of biological replicates in the study
    real(real64), intent(in) :: expr(n_reps, n_genes)
    !! Expression matrix containing
    real(real64), intent(in) :: means(n_genes)
    !! Per-gene mean expression values
    real(real64), intent(out) :: resid(n_reps, n_genes)
    !! Matrix of signed residuals
    
    integer(int32) :: g, i
    
    ! Use do concurrent for parallelization across genes and replicates
    do concurrent (g = 1:n_genes)
      do concurrent (i = 1:n_reps)
        if (.not. ieee_is_nan(expr(i, g))) then
          resid(i, g) = expr(i, g) - means(g)
        else
          resid(i, g) = ieee_value(resid(i, g), ieee_quiet_nan)
        end if
      end do
    end do
    
  end subroutine compute_residuals_helper

  !> Pool per-gene mean expression values across studies
  subroutine pool_means_alloc(n_genes_S1, mean_S1, n_genes_S2, mean_S2, n_points, N_pool, x_star, ierr)
    integer(int32), intent(in) :: n_genes_S1
    !! Number of genes in study S1
    integer(int32), intent(in) :: n_genes_S2
    !! Number of genes in study S2
    integer(int32), intent(in) :: n_points
    !! Number of reference points to define
    real(real64), intent(in) :: mean_S1(n_genes_S1)
    !! Per-gene mean expression values
    real(real64), intent(in) :: mean_S2(n_genes_S2)
    !! Per-gene mean expression values
    integer(int32), intent(out) :: N_pool
    !! Total number of valid (non-NaN) pooled mean-expression values
    real(real64), intent(out) :: x_star(n_points)
    !! Mean-expression reference points
    integer(int32), intent(out) :: ierr
    !! Error code
    
    real(real64), allocatable :: pooled_means(:)
    integer(int32), allocatable :: perm(:)
    integer(int32) :: i, valid_count

    call set_ok(ierr)

    call validate_in_range_int(n_genes_S1, ierr, min = 1_int32)
    call validate_in_range_int(n_genes_S2, ierr, min = 1_int32)
    call validate_in_range_int(n_points, ierr, min = 1_int32)
    ! mean values can contain NaN
    if(is_err(ierr)) return
    
    ! Count total non-NaN values
    valid_count = 0
    do i = 1, n_genes_S1
      if (.not. ieee_is_nan(mean_S1(i))) valid_count = valid_count + 1
    end do
    do i = 1, n_genes_S2
      if (.not. ieee_is_nan(mean_S2(i))) valid_count = valid_count + 1
    end do
    
    N_pool = valid_count
    
    ! Allocate arrays for pooled means
    allocate(pooled_means(N_pool))
    allocate(perm(N_pool))

    ! Initialize permutation array
    do i = 1, N_pool
      perm(i) = i
    end do

    call pool_means_helper(n_genes_S1, mean_S1, n_genes_S2, mean_S2, n_points, N_pool, pooled_means, perm, x_star, ierr)

  end subroutine pool_means_alloc
  
  !> Pool per-gene mean expression values across studies
  pure subroutine pool_means_helper(n_genes_S1, mean_S1, n_genes_S2, mean_S2, n_points, N_pool, pooled_means, perm, x_star, ierr)
    integer(int32), intent(in) :: n_genes_S1
    !! Number of genes in study S1
    integer(int32), intent(in) :: n_genes_S2
    !! Number of genes in study S2
    integer(int32), intent(in) :: n_points
    !! Number of reference points to define
    real(real64), intent(in) :: mean_S1(n_genes_S1)
    !! Per-gene mean expression values
    real(real64), intent(in) :: mean_S2(n_genes_S2)
    !! Per-gene mean expression values
    integer(int32), intent(in) :: N_pool
    !! Total number of valid (non-NaN) pooled mean-expression values
    real(real64), intent(inout) :: pooled_means(N_pool)
    !! Pooled means
    integer(int32), intent(inout) :: perm(N_pool)
    !! Sorting permutation
    real(real64), intent(out) :: x_star(n_points)
    !! Mean-expression reference points
    integer(int32), intent(out) :: ierr
    !! Error code
    
    integer(int32) :: i, j, idx, valid_count
    real(real64) :: pos, quantile_level
    
    ! Fill pooled_means with non-NaN values
    valid_count = 0
    do i = 1, n_genes_S1
      if (.not. ieee_is_nan(mean_S1(i))) then
        valid_count = valid_count + 1
        pooled_means(valid_count) = mean_S1(i)
      end if
    end do
    do i = 1, n_genes_S2
      if (.not. ieee_is_nan(mean_S2(i))) then
        valid_count = valid_count + 1
        pooled_means(valid_count) = mean_S2(i)
      end if
    end do
    
    call heapsort_real(pooled_means, perm)
        
    ! Compute reference points as empirical quantiles using the permutation
    do concurrent (j = 1:n_points)

      !!! CURRENT TOPIC OF DISCUSSION; PLEASE USE EITHER OR
      !!!==================
      !!! Interpolation method
      ! Convert j/(n_points+1) to a percentile (0-100)
      quantile_level = real(j, real64) / real(n_points + 1, real64) * 100.0_real64
      
      ! Use calc_percentile to compute the value
      call calc_percentile(pooled_means, perm, quantile_level, x_star(j), ierr)
      !!!==================

      !!!==================
      !!! Flooring approach
      !!!==================
      ! quantile_level = real(j, real64) / real(n_points + 1, real64)
      ! pos = quantile_level * real(N_pool, real64)
      ! idx = floor(pos)
      
      ! ! Clamp to valid indices (1-based indexing)
      ! if (idx < 1) idx = 1
      ! if (idx > N_pool) idx = N_pool
      
      ! ! Get value from original pooled_means using permutation index
      ! x_star(j) = pooled_means(perm(idx))
    end do    
  end subroutine pool_means_helper

!> Construct neighborhood-based residual sets (kNN)
pure subroutine construct_neighborhoods(n_points, x_star, n_genes_S, mean_S, n_reps_S, resid_S, distances, work_indices, &
                                   N_pool, k_x, neighborhood_residuals, neighborhood_indices, ierr, neighborhood_size)
  integer(int32), intent(in) :: n_points
  !! Number of reference points
  integer(int32), intent(in) :: n_genes_S
  !! Number of genes in the current study
  integer(int32), intent(in) :: n_reps_S
  !! Number of biological replicates in the study
  integer(int32), intent(in) :: N_pool
  !! Total number of pooled mean-expression values across both studies
  real(real64), intent(in) :: x_star(n_points)
  !! Mean-expression reference points
  real(real64), intent(in) :: mean_S(n_genes_S)
  !! Per-gene mean expression values
  real(real64), intent(in) :: resid_S(n_reps_S, n_genes_S)
  !! Matrix of signed residuals
  real(real64), intent(inout) :: distances(n_genes_S)
  !! Distances work array
  integer(int32), intent(inout) :: work_indices(n_genes_S)
  !! Work array for indices (not to be confused with permutation)
  integer(int32), intent(out) :: k_x
  !! Neighborhood size used (constant for all reference points)
  real(real64), intent(out) :: neighborhood_residuals(n_points, n_reps_S * 1000)
  !! Collection of residual vectors for each neighborhood
  integer(int32), intent(out) :: neighborhood_indices(n_points, 1000)
  !! Indices of selected neighborhood genes
  integer(int32), intent(out) :: ierr
  !! Error code
  integer(int32), intent(in), optional :: neighborhood_size
  !! Optional explicit neighborhood size, otherwise use default

  call set_ok(ierr)
  call validate_in_range_int(n_points, ierr, min = 1_int32)
  call validate_in_range_int(n_genes_S, ierr, min = 1_int32)
  call validate_in_range_int(n_reps_S, ierr, min = 1_int32)
  call validate_in_range_int(N_pool, ierr, min = 1_int32)
  if(present(neighborhood_size)) then
    call validate_in_range_int(neighborhood_size, ierr, min = 1_int32)
  end if
  if(is_err(ierr)) return

  call construct_neighborhoods_helper(n_points, x_star, n_genes_S, mean_S, n_reps_S, resid_S, distances, work_indices, N_pool, k_x, neighborhood_residuals, neighborhood_indices, neighborhood_size)

end subroutine construct_neighborhoods
  
!> Construct neighborhood-based residual sets (kNN)
pure subroutine construct_neighborhoods_helper(n_points, x_star, n_genes_S, mean_S, n_reps_S, resid_S, distances, work_indices, &
                                   N_pool, k_x, neighborhood_residuals, neighborhood_indices, neighborhood_size)
  integer(int32), intent(in) :: n_points
  !! Number of reference points
  integer(int32), intent(in) :: n_genes_S
  !! Number of genes in the current study
  integer(int32), intent(in) :: n_reps_S
  !! Number of biological replicates in the study
  integer(int32), intent(in) :: N_pool
  !! Total number of pooled mean-expression values across both studies
  real(real64), intent(in) :: x_star(n_points)
  !! Mean-expression reference points
  real(real64), intent(in) :: mean_S(n_genes_S)
  !! Per-gene mean expression values
  real(real64), intent(in) :: resid_S(n_reps_S, n_genes_S)
  !! Matrix of signed residuals
  real(real64), intent(inout) :: distances(n_genes_S)
  !! Distances work array
  integer(int32), intent(inout) :: work_indices(n_genes_S)
  !! Work array for indices 
  integer(int32), intent(out) :: k_x
  !! Neighborhood size used (constant for all reference points)
  real(real64), intent(out) :: neighborhood_residuals(n_points, n_reps_S * 1000)
  !! Collection of residual vectors for each neighborhood
  integer(int32), intent(out) :: neighborhood_indices(n_points, 1000)
  !! Indices of selected neighborhood genes
  integer(int32), intent(in), optional :: neighborhood_size
  !! Optional explicit neighborhood size
  
  integer(int32) :: j, g, i, m, residual_count
  integer(int32) :: perm(n_genes_S)  ! perm for heapsort
  
  ! Calculate neighborhood size
  if(present(neighborhood_size)) then 
    k_x = min(neighborhood_size, 1000)
  else
    k_x = min(max(100, N_pool / (2 * n_points)), 1000)
  end if
  ! Initialize output arrays
  neighborhood_residuals = ieee_value(0.0_real64, ieee_quiet_nan)
  neighborhood_indices = -1
  
  ! Process each reference point
  do concurrent (j = 1:n_points)
    ! Initialize permutation vector: [1, 2, 3, ..., n_genes_S]
    do g = 1, n_genes_S
      perm(g) = g
      ! Compute distances for ALL genes
      distances(g) = abs(mean_S(g) - x_star(j))  ! NaN if mean_S(g) is NaN
    end do
    
    ! Sort distances using quicksort
    ! heapsort_real will reorder perm so that distances(perm(1:n_genes_S)) is sorted
    call heapsort_real(distances, perm)
    
    ! Store the k_x nearest neighbor indices
    ! perm(1), perm(2), ..., perm(k_x) now contain indices of genes with smallest distances
    do m = 1, min(k_x, n_genes_S)
      g = perm(m)  ! Get the actual gene index from the permutation vector
      ! Only store if the distance is finite (not NaN)
      if (.not. ieee_is_nan(distances(perm(m)))) then
        neighborhood_indices(j, m) = g
      else
        ! Once we hit NaN, all remaining distances will be NaN
        exit
      end if
    end do
    
    ! Extract and pool residuals from k-nearest neighbors
    residual_count = 0
    do m = 1, min(k_x, n_genes_S)
      g = perm(m)
      if (g > 0 .and. .not. ieee_is_nan(mean_S(g))) then
        do i = 1, n_reps_S
          if (.not. ieee_is_nan(resid_S(i, g))) then
            residual_count = residual_count + 1
            if (residual_count <= n_reps_S * 1000) then
              neighborhood_residuals(j, residual_count) = resid_S(i, g)
            end if
          end if
        end do
      end if
    end do
  end do    
end subroutine construct_neighborhoods_helper
  
end module tox_jenson_shannon_test

!> C wrapper for compute_gene_means
subroutine compute_gene_means_c(n_genes, n_reps, expr, means, ierr) &
          bind(C, name="compute_gene_means_c")
  use, intrinsic :: iso_c_binding, only: c_int, c_double
  use tox_jenson_shannon_test, only: compute_gene_means
  M_USE_NULL_VALIDATION
  implicit none
  
  integer(c_int), intent(in), target :: n_genes
  !! Number of genes in the study
  integer(c_int), intent(in), target :: n_reps
  !! Number of biological replicates in the study
  real(c_double), intent(in), target :: expr(n_reps, n_genes)
  !! Expression matrix
  real(c_double), intent(out), target :: means(n_genes)
  !! Per-gene mean expression values
  integer(c_int), intent(out), target :: ierr
  !! Error code

  M_CHECK_IERR_NON_NULL
  M_CHECK_NON_NULL(n_genes)
  M_CHECK_NON_NULL(n_reps)
  M_CHECK_NON_NULL(expr)
  M_CHECK_NON_NULL(means)
  
  call compute_gene_means(n_genes, n_reps, expr, means, ierr)
end subroutine compute_gene_means_c

!> C wrapper for compute_residuals
subroutine compute_residuals_c(n_genes, n_reps, expr, means, resid, ierr) &
          bind(C, name="compute_residuals_c")
  use, intrinsic :: iso_c_binding, only: c_int, c_double
  use tox_jenson_shannon_test, only: compute_residuals
  M_USE_NULL_VALIDATION
  implicit none
  
  integer(c_int), intent(in), target :: n_genes
  !! Number of genes in the study
  integer(c_int), intent(in), target :: n_reps
  !! Number of biological replicates in the study
  real(c_double), intent(in), target :: expr(n_reps, n_genes)
  !! Expression matrix
  real(c_double), intent(in), target :: means(n_genes)
  !! Per-gene mean expression values
  real(c_double), intent(out), target :: resid(n_reps, n_genes)
  !! Matrix of signed residuals
  integer(c_int), intent(out), target :: ierr
  !! Error code
  
  M_CHECK_IERR_NON_NULL
  M_CHECK_NON_NULL(n_genes)
  M_CHECK_NON_NULL(n_reps)
  M_CHECK_NON_NULL(expr)
  M_CHECK_NON_NULL(means)
  M_CHECK_NON_NULL(resid)

  call compute_residuals(n_genes, n_reps, expr, means, resid, ierr)
end subroutine compute_residuals_c

!> C wrapper for pool_means_alloc
subroutine pool_means_c(n_genes_S1, mean_S1, n_genes_S2, mean_S2, &
                              n_points, N_pool, x_star, ierr) &
          bind(C, name="pool_means_c")
  use, intrinsic :: iso_c_binding, only: c_int, c_double
  use tox_jenson_shannon_test, only: pool_means_alloc
  M_USE_NULL_VALIDATION
  implicit none
  
  integer(c_int), intent(in), target :: n_genes_S1
  !! Number of genes in study S1
  integer(c_int), intent(in), target :: n_genes_S2
  !! Number of genes in study S2
  integer(c_int), intent(in), target :: n_points
  !! Number of reference points to define
  real(c_double), intent(in), target :: mean_S1(n_genes_S1)
  !! Per-gene mean expression values for S1
  real(c_double), intent(in), target :: mean_S2(n_genes_S2)
  !! Per-gene mean expression values for S2
  integer(c_int), intent(out), target :: N_pool
  !! Total number of valid (non-NaN) pooled mean-expression values
  real(c_double), intent(out), target :: x_star(n_points)
  !! Mean-expression reference points
  integer(c_int), intent(out), target :: ierr
  !! Error code

  M_CHECK_IERR_NON_NULL
  M_CHECK_NON_NULL(n_genes_S1)
  M_CHECK_NON_NULL(n_genes_S2)
  M_CHECK_NON_NULL(n_points)
  M_CHECK_NON_NULL(mean_S1)
  M_CHECK_NON_NULL(mean_S2)
  M_CHECK_NON_NULL(N_pool)
  M_CHECK_NON_NULL(x_star)
  
  call pool_means_alloc(n_genes_S1, mean_S1, &
                        n_genes_S2, mean_S2, &
                        n_points, N_pool, x_star, ierr)
end subroutine pool_means_c

!> C wrapper for construct_neighborhoods
subroutine construct_neighborhoods_c(n_points, x_star, n_genes_S, mean_S, &
                                    n_reps_S, resid_S, neighborhood_size, N_pool, k_x, &
                                    neighborhood_residuals, neighborhood_indices, ierr) &
          bind(C, name="construct_neighborhoods_c")
  use, intrinsic :: iso_c_binding, only: c_int, c_double
  use tox_jenson_shannon_test, only: construct_neighborhoods
  M_USE_NULL_VALIDATION
  implicit none
  
  integer(c_int), intent(in), target :: n_points
  !! Number of reference points
  integer(c_int), intent(in), target :: n_genes_S
  !! Number of genes in the current study
  integer(c_int), intent(in), target :: n_reps_S
  !! Number of biological replicates in the study
  integer(c_int), intent(in), target :: N_pool
  !! Total number of pooled mean-expression values across both studies
  real(c_double), intent(in), target :: x_star(n_points)
  !! Mean-expression reference points
  real(c_double), intent(in), target :: mean_S(n_genes_S)
  !! Per-gene mean expression values
  real(c_double), intent(in), target :: resid_S(n_reps_S, n_genes_S)
  !! Matrix of signed residuals
  integer(c_int), intent(in), target :: neighborhood_size
  !! Explicit size of the neighborhood
  integer(c_int), intent(out), target :: k_x
  !! Neighborhood size used (constant for all reference points)
  real(c_double), intent(out), target :: neighborhood_residuals(n_points, n_reps_S * 1000)
  !! Collection of residual vectors for each neighborhood
  integer(c_int), intent(out), target :: neighborhood_indices(n_points, 1000)
  !! Indices of selected neighborhood genes
  integer(c_int), intent(out), target :: ierr
  !! Error code
  
  real(c_double) :: distances(n_genes_S)
  integer(c_int) :: work_indices(n_genes_S)

  M_CHECK_IERR_NON_NULL
  M_CHECK_NON_NULL(n_points)
  M_CHECK_NON_NULL(n_genes_S)
  M_CHECK_NON_NULL(n_reps_S)
  M_CHECK_NON_NULL(N_pool)
  M_CHECK_NON_NULL(x_star)
  M_CHECK_NON_NULL(mean_S)
  M_CHECK_NON_NULL(resid_S)
  M_CHECK_NON_NULL(k_x)
  M_CHECK_NON_NULL(neighborhood_residuals)
  M_CHECK_NON_NULL(neighborhood_indices)
  M_CHECK_NON_NULL(neighborhood_size)

  if(neighborhood_size > 0) then
    call construct_neighborhoods(n_points, x_star, n_genes_S, mean_S, &
                              n_reps_S, resid_S, distances, work_indices, &
                              N_pool, k_x, neighborhood_residuals, neighborhood_indices, ierr, neighborhood_size)
  else
    call construct_neighborhoods(n_points, x_star, n_genes_S, mean_S, &
                              n_reps_S, resid_S, distances, work_indices, &
                              N_pool, k_x, neighborhood_residuals, neighborhood_indices, ierr)
  end if
end subroutine construct_neighborhoods_c