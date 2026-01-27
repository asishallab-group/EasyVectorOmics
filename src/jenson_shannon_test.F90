!> # Global Jensen-Shannon-Divergence Compatibility Test (gJCT)
!!
!! This module implements the four main subroutines for the gJCT algorithm.
module jenson_shannon_test
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
  use f42_utils, only: heapsort_real
  implicit none
  
  private
  public :: compute_gene_means, compute_residuals, pool_means, construct_neighborhoods
  
contains
  
  !> Compute per-gene mean expression
  subroutine compute_gene_means(n_genes, n_reps, expr, means)
    integer(int32), intent(in) :: n_genes
    !! Number of genes in the study
    integer(int32), intent(in) :: n_reps
    !! Number of biological replicates (samples) in the study
    real(real64), intent(in) :: expr(n_reps, n_genes)
    !! Expression matrix containing \( x_{i,g}^{(S)} \) values (samples × genes)
    real(real64), intent(out) :: means(n_genes)
    !! Per-gene mean expression values \( \bar{x}_g^{(S)} \)
    
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
    
  end subroutine compute_gene_means
  
  !> Compute signed residuals
  subroutine compute_residuals(n_genes, n_reps, expr, means, resid)
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
    do concurrent (g = 1:n_genes) (i = 1:n_reps)
      if (.not. ieee_is_nan(expr(i, g))) then
        resid(i, g) = expr(i, g) - means(g)
      else
        resid(i, g) = ieee_value(resid(i, g), ieee_quiet_nan)
      end if
    end do
    
  end subroutine compute_residuals
  
  !> Pool per-gene mean expression values across studies
  subroutine pool_means(n_genes_S1, mean_S1, n_genes_S2, mean_S2, n_points, N_pool, x_star)
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
    real(real64), intent(out) :: x_star(r)
    !! Mean-expression reference points
    
    real(real64), allocatable :: pooled_means(:)
    real(real64), allocatable :: sorted_means(:)
    integer(int32), allocatable :: perm(:)
    integer(int32) :: i, j, idx, valid_count
    real(real64) :: pos, quantile_level
    
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
    allocate(sorted_means(N_pool))
    allocate(perm(N_pool))
    
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
    
    ! Create permutation array and copy to sorted_means
    do i = 1, N_pool
      perm(i) = i
      sorted_means(i) = pooled_means(i)
    end do
    
    ! Sort pooled_means using quicksort
    call heapsort_real(pooled_means, perm)
    
    ! Reorder sorted_means according to permutation
    do i = 1, N_pool
      sorted_means(i) = pooled_means(perm(i))
    end do
    
    ! Compute reference points as empirical quantiles
    do concurrent (j = 1:n_points)
      quantile_level = real(j, real64) / real(n_points + 1, real64)
      pos = quantile_level * real(N_pool, real64)
      idx = floor(pos)
      
      ! Clamp to valid indices (1-based indexing)
      if (idx < 1) idx = 1
      if (idx > N_pool) idx = N_pool
      
      ! Get value from sorted means
      x_star(j) = sorted_means(idx)
    end do
    
    deallocate(pooled_means, sorted_means, perm)
    
  end subroutine pool_means
  
  !> Construct neighborhood-based residual sets (kNN)
  subroutine construct_neighborhoods(n_points, x_star, n_genes_S, mean_S, n_reps_S, resid_S, &
                                     N_pool, k_x, neighborhood_residuals, neighborhood_indices)
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
    integer(int32), intent(out) :: k_x
    !! Neighborhood size used (constant for all reference points)
    real(real64), intent(out) :: neighborhood_residuals(n_points, n_reps_S * 1000)
    !! Collection of residual vectors for each neighborhood
    integer(int32), intent(out) :: neighborhood_indices(n_points, 1000)
    !! Indices of selected neighborhood genes
    
    real(real64), allocatable :: distances(:)
    integer(int32), allocatable :: gene_indices(:)
    integer(int32) :: j, g, i, m, valid_count, neighbor_count, residual_count
    real(real64) :: dist
    integer(int32) :: temp_idx
    real(real64) :: temp_dist
    
    ! Calculate neighborhood size
    k_x = min(max(100, N_pool / (2 * n_points)), 1000)
    
    ! Pre-allocate arrays
    allocate(distances(n_genes_S))
    allocate(gene_indices(n_genes_S))
    
    ! Initialize output arrays
    neighborhood_residuals = ieee_value(0.0_real64, ieee_quiet_nan)
    neighborhood_indices = -1
    
    ! Process each reference point
    do j = 1, n_points
      valid_count = 0
      
      ! Compute distances for all valid genes (non-NaN means)
      do g = 1, n_genes_S
        if (.not. ieee_is_nan(mean_S(g))) then
          valid_count = valid_count + 1
          gene_indices(valid_count) = g
          distances(valid_count) = abs(mean_S(g) - x_star(j))
        end if
      end do
      
      ! Find k-nearest neighbors if there are valid genes
      if (valid_count > 0) then
        ! Simple selection sort for k_x smallest distances
        do m = 1, min(k_x, valid_count)
          ! Find minimum distance among remaining elements
          do i = m + 1, valid_count
            if (distances(i) < distances(m)) then
              ! Swap distances
              temp_dist = distances(m)
              distances(m) = distances(i)
              distances(i) = temp_dist
              
              ! Swap indices
              temp_idx = gene_indices(m)
              gene_indices(m) = gene_indices(i)
              gene_indices(i) = temp_idx
            end if
          end do
          
          ! Store the m-th nearest neighbor index
          neighborhood_indices(j, m) = gene_indices(m)
        end do
        
        ! Extract and pool residuals from k-nearest neighbors
        residual_count = 0
        do m = 1, min(k_x, valid_count)
          g = neighborhood_indices(j, m)
          do i = 1, n_reps_S
            if (.not. ieee_is_nan(resid_S(i, g))) then
              residual_count = residual_count + 1
              if (residual_count <= n_reps_S * 1000) then
                neighborhood_residuals(j, residual_count) = resid_S(i, g)
              end if
            end if
          end do
        end do
      end if
    end do
    
    deallocate(distances, gene_indices)
    
  end subroutine construct_neighborhoods
  
end module jenson_shannon_test