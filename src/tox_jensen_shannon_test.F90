#include "macros.h"

!> # Global Jensen-Shannon-Divergence Compatibility Test (gJCT)
!!
!! This module implements the four main subroutines for the gJCT algorithm.
module tox_jensen_shannon_test
    use safeguard
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
    use f42_utils, only: heapsort_real, calc_percentile_helper, clamp
    use tox_errors, only: validate_all_in_range_real, validate_in_range_int, is_err, set_ok, validate_dimension_size, ERR_ALLOC_FAIL, set_err
    implicit none

    private
    public :: compute_gene_means, compute_gene_means_helper, compute_residuals, compute_residuals_helper, pool_means_alloc, construct_neighborhoods, construct_neighborhoods_helper

contains

    !> Compute per-gene mean expression, ignoring NaN values
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

        call validate_dimension_size(n_genes, ierr)
        call validate_dimension_size(n_reps, ierr)
        ! expression can contain NaN
        if (is_err(ierr)) return

        call compute_gene_means_helper(n_genes, n_reps, expr, means)
    end subroutine compute_gene_means

    !> (no input validation) Compute per-gene mean expression, ignoring NaN values
    pure subroutine compute_gene_means_helper(n_genes, n_reps, expr, means)
        integer(int32), intent(in) :: n_genes
            !! Number of genes in the study
        integer(int32), intent(in) :: n_reps
            !! Number of biological replicates in the study
        real(real64), intent(in) :: expr(n_reps, n_genes)
            !! Expression matrix
        real(real64), intent(out) :: means(n_genes)
            !! Per-gene mean expression values

        integer(int32) :: i_gene, i_rep, n_included
        real(real64) :: sum_val

        ! Use do concurrent for parallelization across genes
        do concurrent(i_gene=1:n_genes) local(sum_val, n_included)
            sum_val = 0.0_real64
            n_included = 0

            ! Count valid (non-NaN) replicates and compute sum
            do concurrent(i_rep=1:n_reps) shared(expr) reduce(+:sum_val, n_included)
                if (.not. ieee_is_nan(expr(i_rep, i_gene))) then
                    sum_val = sum_val + expr(i_rep, i_gene)
                    n_included = n_included + 1
                end if
            end do

            means(i_gene) = sum_val/real(n_included, real64)
        end do
    end subroutine compute_gene_means_helper

    !> Compute signed residuals (centering by mean)
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
        call validate_dimension_size(n_genes, ierr)
        call validate_dimension_size(n_reps, ierr)
        ! family means and expr containing NaN is expected behaviour
        if (is_err(ierr)) return

        call compute_residuals_helper(n_genes, n_reps, expr, means, resid)
    end subroutine compute_residuals

    !> (no input validation) Compute signed residuals (centering by mean)
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

        integer(int32) :: i_gene, i_rep

        do concurrent(i_gene=1:n_genes)
            do concurrent(i_rep=1:n_reps) shared(expr, resid, means, i_gene)
                if (.not. ieee_is_nan(expr(i_rep, i_gene))) then
                    resid(i_rep, i_gene) = expr(i_rep, i_gene) - means(i_gene)
                else
                    resid(i_rep, i_gene) = ieee_value(resid(i_rep, i_gene), ieee_quiet_nan)
                end if
            end do
        end do

    end subroutine compute_residuals_helper

    !> Pool per-gene mean expression values across studies
    subroutine pool_means_alloc(n_genes_S1, mean_S1, n_genes_S2, mean_S2, n_points, n_pool, x_star, ierr)
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
        integer(int32), intent(out) :: n_pool
            !! Total number of included (non-NaN) pooled mean-expression values
        real(real64), intent(out) :: x_star(n_points)
            !! Mean-expression reference points
        integer(int32), intent(out) :: ierr
            !! Error code

        real(real64), allocatable :: pooled_means(:)
        integer(int32), allocatable :: perm(:)
        integer(int32) :: i_gene, pool_idx

        call set_ok(ierr)

        call validate_dimension_size(n_genes_S1, ierr)
        call validate_dimension_size(n_genes_S2, ierr)
        call validate_dimension_size(n_points, ierr)
        ! mean values can contain NaN
        if (is_err(ierr)) return

        ! Allocate arrays for pooled means
        M_ALLOCATE(pooled_means(n_genes_S1 + n_genes_S2))
        M_ALLOCATE(perm(n_genes_S1 + n_genes_S2))

        do concurrent(i_gene=1:n_genes_S1) shared(pooled_means, mean_S1)
            pooled_means(i_gene) = mean_S1(i_gene)
            perm(i_gene) = i_gene
        end do

        do concurrent(i_gene=1:n_genes_S2) local(pool_idx) shared(pooled_means, mean_S2)
            pool_idx = n_genes_S1 + i_gene
            pooled_means(pool_idx) = mean_S2(i_gene)
            perm(pool_idx) = pool_idx
        end do

        call heapsort_real(pooled_means, perm)

        call pool_means_helper(pooled_means, perm, n_genes_S1, n_genes_S2, n_points, n_pool, x_star)

    end subroutine pool_means_alloc

    !> Pool per-gene mean expression values across studies
    pure subroutine pool_means_helper(pooled_means, pooled_means_perm, n_genes_S1, n_genes_S2, n_points, n_pool, x_star)
        integer(int32), intent(in) :: n_genes_S1
            !! Number of genes in study S1
        integer(int32), intent(in) :: n_genes_S2
            !! Number of genes in study S2
        integer(int32), intent(in) :: n_points
            !! Number of reference points to define
        integer(int32), intent(out) :: n_pool
            !! Total number of included (non-NaN) pooled mean-expression values
        real(real64), intent(in) :: pooled_means(n_genes_S1 + n_genes_S2)
            !! Pooled means
        integer(int32), intent(in) :: pooled_means_perm(n_genes_S1 + n_genes_S2)
            !! Sorting permutation for `pooled_means`
        real(real64), intent(out) :: x_star(n_points)
            !! Mean-expression reference points

        integer(int32) :: i_gene, i_point
        real(real64) :: quantile_level

        n_pool = size(pooled_means, kind=int32)

        ! NaN is always last -> find last non-NaN index for percentile calculation
        do i_gene = n_pool, 1, -1
            if (ieee_is_nan(pooled_means(pooled_means_perm(i_gene)))) then
                n_pool = n_pool - 1
            else
                exit
            end if
        end do

        if (n_pool == 0) then
            x_star = ieee_value(x_star, ieee_quiet_nan)
        else
            ! Compute reference points as empirical quantiles using the permutation
            do concurrent(i_point=1:n_points) local(quantile_level) shared(n_points, pooled_means, pooled_means_perm, n_pool, x_star)
                quantile_level = real(i_point, real64)/real(n_points + 1, real64)*100.0_real64

                ! Use calc_percentile to compute the value
                call calc_percentile_helper(pooled_means(:n_pool), pooled_means_perm(:n_pool), quantile_level, x_star(i_point))
            end do
        end if
    end subroutine pool_means_helper

    pure integer(int32) function calc_neighborhood_size(n_pool, n_points, n_genes_S, mean_S, desired_size) result(n_neighbors)
        integer(int32), intent(in) :: n_pool
            !! Total number of pooled mean-expression values across both studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), intent(in) :: n_genes_S
            !! Number of genes in the current study
        real(real64), intent(in) :: mean_S(n_genes_S)
            !! Per-gene mean expression values
        integer(int32), intent(in), optional :: desired_size
            !! Optional desired neighborhood size, default=1000

        integer(int32) :: max_neighbors, i_gene, min_neighbors

        M_DEFAULT_VAL(desired_size, max_neighbors, 1000_int32)

        min_neighbors = 100_int32

        ! If less neighbors than default min are desired, take it
        if (max_neighbors < min_neighbors) then
            n_neighbors = max_neighbors

        ! Take at least `min_neighbors` neighbors. It could be lower if there are too many reference points -> low steps across x_star.
        else
            n_neighbors = int(clamp(n_pool/(2*n_points), min_val=min_neighbors, max_val=max_neighbors))
        end if

        ! don't take more neighbors than genes
        n_neighbors = min(n_neighbors, n_genes_S)

        ! Exclude genes with NaN mean, so those with all replicates NaN
        do concurrent(i_gene=1:n_genes_S) shared(mean_S) reduce(+:n_neighbors)
            if (ieee_is_nan(mean_S(i_gene))) then
                n_neighbors = n_neighbors - 1
            end if
        end do

        n_neighbors = max(0_int32, n_neighbors)
    end function calc_neighborhood_size

    !> Construct neighborhood-based residual sets (kNN)
    pure subroutine construct_neighborhoods_alloc(n_points, x_star, n_genes_S, mean_S, n_reps_S, resid_S, &
                                                  neighborhood_residuals, neighborhood_indices, n_neighbors, ierr)
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), intent(in) :: n_genes_S
            !! Number of genes in the current study
        integer(int32), intent(in) :: n_reps_S
            !! Number of biological replicates in the study
        real(real64), intent(in) :: x_star(n_points)
            !! Mean-expression reference points
        real(real64), intent(in) :: mean_S(n_genes_S)
            !! Per-gene mean expression values
        real(real64), intent(in) :: resid_S(n_reps_S, n_genes_S)
            !! Matrix of signed residuals
        real(real64), intent(out) :: neighborhood_residuals(n_reps_S, n_neighbors, n_points)
            !! Collection of residual vectors for each neighborhood
        integer(int32), intent(out) :: neighborhood_indices(n_neighbors, n_points)
            !! Indices of selected neighborhood genes
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors, **CALCULATE IT WITH [[tox_jensen_shannon_test(module):calc_neighborhood_size(function)]]**
        integer(int32), intent(out) :: ierr
            !! Error code

        real(real64), dimension(:), allocatable :: tmp_distances
        integer(int32), dimension(:), allocatable :: tmp_distances_perm

        call set_ok(ierr)
        call validate_dimension_size(n_points, ierr)
        call validate_dimension_size(n_genes_S, ierr)
        call validate_dimension_size(n_reps_S, ierr)
        call validate_dimension_size(n_neighbors, ierr)
        if (is_err(ierr)) return

        M_ALLOCATE(tmp_distances(n_neighbors))
        M_ALLOCATE(tmp_distances_perm(n_neighbors))

        call construct_neighborhoods_helper(n_points, x_star, n_genes_S, mean_S, n_reps_S, resid_S, tmp_distances, tmp_distances_perm, neighborhood_residuals, neighborhood_indices, n_neighbors)

    end subroutine construct_neighborhoods_alloc

    !> Construct neighborhood-based residual sets (kNN)
    pure subroutine construct_neighborhoods(n_points, x_star, n_genes_S, mean_S, n_reps_S, resid_S, tmp_distances, tmp_distances_perm, neighborhood_residuals, neighborhood_indices, n_neighbors, ierr)
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), intent(in) :: n_genes_S
            !! Number of genes in the current study
        integer(int32), intent(in) :: n_reps_S
            !! Number of biological replicates in the study
        real(real64), intent(in) :: x_star(n_points)
            !! Mean-expression reference points
        real(real64), intent(in) :: mean_S(n_genes_S)
            !! Per-gene mean expression values
        real(real64), intent(in) :: resid_S(n_reps_S, n_genes_S)
            !! Matrix of signed residuals
        real(real64), intent(out) :: tmp_distances(n_genes_S)
            !! Distances work array
        integer(int32), intent(out) :: tmp_distances_perm(n_genes_S)
            !! Work array for permutation vector to sort `tmp_distances_perm`
        real(real64), intent(out) :: neighborhood_residuals(n_reps_S, n_neighbors, n_points)
            !! Collection of residual vectors for each neighborhood
        integer(int32), intent(out) :: neighborhood_indices(n_neighbors, n_points)
            !! Indices of selected neighborhood genes
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors, **CALCULATE IT WITH [[tox_jensen_shannon_test(module):calc_neighborhood_size(function)]]**
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
        call validate_dimension_size(n_points, ierr)
        call validate_dimension_size(n_genes_S, ierr)
        call validate_dimension_size(n_reps_S, ierr)
        call validate_dimension_size(n_neighbors, ierr)
        call validate_in_range_int(n_neighbors, ierr, min=1_int32, max=n_genes_S)
        if (is_err(ierr)) return

        call construct_neighborhoods_helper(n_points, x_star, n_genes_S, mean_S, n_reps_S, resid_S, tmp_distances, tmp_distances_perm, neighborhood_residuals, neighborhood_indices, n_neighbors)

    end subroutine construct_neighborhoods

    !> (no input validation) Construct neighborhood-based residual sets (kNN)
    pure subroutine construct_neighborhoods_helper(n_points, x_star, n_genes_S, mean_S, n_reps_S, resid_S, tmp_distances, tmp_distances_perm, &
                                                   neighborhood_residuals, neighborhood_indices, n_neighbors)
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), intent(in) :: n_genes_S
            !! Number of genes in the current study
        integer(int32), intent(in) :: n_reps_S
            !! Number of biological replicates in the study
        real(real64), intent(in) :: x_star(n_points)
            !! Mean-expression reference points
        real(real64), intent(in) :: mean_S(n_genes_S)
            !! Per-gene mean expression values
        real(real64), intent(in) :: resid_S(n_reps_S, n_genes_S)
            !! Matrix of signed residuals
        real(real64), intent(out) :: tmp_distances(n_genes_S)
            !! Distances work array
        integer(int32), intent(out) :: tmp_distances_perm(n_genes_S)
            !! Work array for permutation vector to sort `tmp_distances_perm`
        real(real64), intent(out) :: neighborhood_residuals(n_reps_S, n_neighbors, n_points)
            !! Collection of residual vectors for each neighborhood
        integer(int32), intent(out) :: neighborhood_indices(n_neighbors, n_points)
            !! Indices of selected neighborhood genes
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors, **CALCULATE IT WITH [[tox_jensen_shannon_test(module):calc_neighborhood_size(function)]]**

        integer(int32) :: i_point, i_gene, gene_idx

        ! Process each reference point
        do i_point = 1, n_points

            ! Calculate distances.
            do concurrent(i_gene=1:n_genes_S) shared(tmp_distances, tmp_distances_perm, i_point, x_star, mean_S)
                tmp_distances(i_gene) = abs(mean_S(i_gene) - x_star(i_point))

                ! Initialize perm
                tmp_distances_perm(i_gene) = i_gene
            end do

            ! Sort distances using heapsort
            ! heapsort_real will reorder tmp_distances_perm so that tmp_distances(tmp_distances_perm(1:n_genes_S)) is sorted
            ! `calc_neighborhood_size` guarantees that the NaN `mean_S` indices are not included after sorting (they come after tmp_distances_perm(:n_neighbors))
            call heapsort_real(tmp_distances, tmp_distances_perm)

            ! Store the n_neighbors nearest neighbor indices
            ! tmp_distances_perm(1:n_neighbors) now contain indices of genes with smallest distances
            do concurrent(i_gene=1:n_neighbors) local(gene_idx) shared(tmp_distances_perm, neighborhood_indices, neighborhood_residuals, resid_S)
                gene_idx = tmp_distances_perm(i_gene)  ! Get the actual gene index from the permutation vector

                neighborhood_indices(i_gene, i_point) = gene_idx

                neighborhood_residuals(:, i_gene, i_point) = resid_S(:, gene_idx)
            end do
        end do
    end subroutine construct_neighborhoods_helper
end module tox_jensen_shannon_test
