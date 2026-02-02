#include "macros.h"

module tox_jensen_shannon_divergence
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use f42_utils, only: clamp, calc_percentile_helper, is_close, sort_array_heapsort, shuffle_vector, init_random
    use tox_errors, only: set_ok, set_err, is_err, ERR_ALLOC_FAIL, validate_dimension_size, validate_in_range_real, validate_all_in_range_real, validate_in_range_int, validate_all_in_range_int
    implicit none
contains

    !> Computes the shared residual range [-R, R] for the computed residuals from studies S1 and S2
    pure subroutine determine_shared_residual_range(abs_residual_pool, abs_residual_pool_perm, pool_size, shared_residual_range, ierr, residual_range_quantile)
        integer(int32), intent(in) :: pool_size
            !! Size of pool of residuals `abs_residual_pool`, usually `(n_reps_S1 + n_reps_2)*n_neighbors*n_points`
        real(real64), intent(in), optional :: residual_range_quantile
            !! Quantile for determining the residual range, default: 95.0
        real(real64), intent(out) :: shared_residual_range
            !! Computed residual range (R)
        real(real64), dimension(pool_size), intent(in) :: abs_residual_pool
            !! The absolute residual values of the concatenated S1,S2 residuals
        integer(int32), dimension(pool_size), intent(in) :: abs_residual_pool_perm
            !! The permutation vector that sorts `abs_residual_pool`
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(pool_size, ierr)
        call validate_in_range_real(residual_range_quantile, ierr, min=0.0_real64, max=100.0_real64)
        call validate_all_in_range_int(abs_residual_pool_perm, pool_size, ierr, min=1_int32, max=pool_size)

        if (is_err(ierr)) return

        call determine_shared_residual_range_helper(abs_residual_pool, abs_residual_pool_perm, pool_size, shared_residual_range, residual_range_quantile)
    end subroutine determine_shared_residual_range

    !> (no input validation) Computes the shared residual range [-R, R] for the computed residuals from studies S1 and S2
    pure subroutine determine_shared_residual_range_helper(abs_residual_pool, abs_residual_pool_perm, pool_size, shared_residual_range, residual_range_quantile)
        integer(int32), intent(in) :: pool_size
            !! Size of pool of residuals `abs_residual_pool`, usually `(n_reps_S1 + n_reps_2)*n_neighbors*n_points`
        real(real64), intent(in), optional :: residual_range_quantile
            !! Quantile for determining the residual range, default: 95.0
        real(real64), intent(out) :: shared_residual_range
            !! Computed residual range (R)
        real(real64), dimension(pool_size), intent(in) :: abs_residual_pool
            !! The absolute residual values of the concatenated S1,S2 residuals
        integer(int32), dimension(pool_size), intent(in) :: abs_residual_pool_perm
            !! The permutation vector that sorts `abs_residual_pool`

        integer(int32) :: i_pool, last_non_nan
        real(real64) :: actual_quantile

        M_DEFAULT_VAL(residual_range_quantile, actual_quantile, 95.0_real64)

        last_non_nan = pool_size
        ! NaN is always last -> find last non-NaN index for percentile calculation
        do i_pool = last_non_nan, 1, -1
            if (ieee_is_nan(abs_residual_pool(abs_residual_pool_perm(i_pool)))) then
                last_non_nan = last_non_nan - 1
            else
                exit
            end if
        end do

        if (last_non_nan == 0) then
            shared_residual_range = 0.0_real64
            return
        end if

        call calc_percentile_helper(abs_residual_pool(:last_non_nan), abs_residual_pool_perm(:last_non_nan), actual_quantile, shared_residual_range)
    end subroutine determine_shared_residual_range_helper

    !> Computes the shared residual range [-R, R] for the computed residuals from studies S1 and S2
    pure subroutine determine_shared_residual_range_alloc(neighborhood_residuals_S1, neighborhood_residuals_S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, shared_residual_range, ierr, residual_range_quantile)
        integer(int32), intent(in) :: n_reps_S1
            !! Number of replicates in study 1
        integer(int32), intent(in) :: n_reps_S2
            !! Number of replicates in study 2
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors in the studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1 ([[tox_jensen_shannon_test(module):construct_neighborhoods(subroutine)]]), NaN is explicitly allowed for missing values
        real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2 ([[tox_jensen_shannon_test(module):construct_neighborhoods(subroutine)]]), NaN is explicitly allowed for missing values
        real(real64), intent(in), optional :: residual_range_quantile
            !! Quantile for determining the residual range, default: 95.0
        real(real64), intent(out) :: shared_residual_range
            !! Computed residual range (R)
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i_rep, i_neighbor, pool_size, pool_idx, n_predecessors, i_point
        real(real64), dimension(:, :, :), allocatable, target :: abs_residual_pool
        real(real64), dimension(:), pointer :: abs_residual_pool_flat
        integer(int32), dimension(:, :, :), allocatable, target :: perm
        integer(int32), dimension(:), pointer :: perm_flat

        call set_ok(ierr)

        call validate_dimension_size(n_reps_S1, ierr)
        call validate_dimension_size(n_reps_S2, ierr)
        call validate_dimension_size(n_neighbors, ierr)
        call validate_dimension_size(n_points, ierr)

        if (is_err(ierr)) return

        M_ALLOCATE(abs_residual_pool(n_reps_S1 + n_reps_S2, n_neighbors, n_points))
        M_ALLOCATE(perm(n_reps_S1 + n_reps_S2, n_neighbors, n_points))

        pool_size = size(abs_residual_pool, kind=int32)

        ! Collect the absolute residual values and compute the percentile value
        do concurrent (i_point = 1:n_points)
            do concurrent (i_neighbor = 1:n_neighbors) local(n_predecessors) shared(i_point, n_reps_S1, n_reps_S2)
                n_predecessors = ((i_point - 1) * n_neighbors + (i_neighbor - 1)) * (n_reps_S1 + n_reps_S2)
                do concurrent (i_rep = 1:n_reps_S1) shared(n_predecessors, perm, neighborhood_residuals_S1, abs_residual_pool, i_point, i_neighbor)
                    abs_residual_pool(i_rep, i_neighbor, i_point) = abs(neighborhood_residuals_S1(i_rep, i_neighbor, i_point))

                    perm(i_rep, i_neighbor, i_point) = n_predecessors + i_rep
                end do

                do concurrent (i_rep = 1:n_reps_S2) local(pool_idx) shared(n_predecessors, perm, neighborhood_residuals_S2, abs_residual_pool, i_point, i_neighbor)
                    pool_idx = i_rep + n_reps_S1

                    abs_residual_pool(pool_idx, i_neighbor, i_point) = abs(neighborhood_residuals_S2(i_rep, i_neighbor, i_point))

                    perm(pool_idx, i_neighbor, i_point) = n_predecessors + pool_idx
                end do
            end do
        end do

        abs_residual_pool_flat(1:pool_size) => abs_residual_pool
        perm_flat(1:pool_size) => perm

        call sort_array_heapsort(abs_residual_pool_flat, perm_flat)

        call determine_shared_residual_range(abs_residual_pool, perm, pool_size, shared_residual_range, ierr, residual_range_quantile)
    end subroutine determine_shared_residual_range_alloc

    !> Summarizes the neighborhood residuals in absolute histogram counts and probability mass functions `pmf(residual, bin)` (actually a matrix)
    pure subroutine build_residual_histograms(neighborhood_residuals, n_reps, n_neighbors, n_points, shared_residual_range, n_bins, counts, pmf, included_n_reps, ierr, neighbor_mask)
        integer(int32), intent(in) :: n_reps
            !! Number of replicates of the study
        integer(int32), intent(in) :: n_neighbors
            !! Number of reference points (k)
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        real(real64), dimension(n_reps, n_neighbors, n_points), intent(in) :: neighborhood_residuals
            !! Computed neighborhood residuals for a study ([[tox_jensen_shannon_test(module):construct_neighborhoods(subroutine)]]), NaN is explicitly allowed for missing values
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range (R) from [[tox_jensen_shannon_divergence(module):determine_shared_residual_range_alloc(subroutine)]]
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins in range [-R,R]
        integer(int32), dimension(n_points, n_bins), intent(out) :: counts
            !! Absolute counts of a residual per bin
        real(real64), dimension(n_points, n_bins), intent(out) :: pmf
            !! `counts` normalized to `0 <= counts(:, i) <= 1` and `sum(counts(:, i)) == 1`
        integer(int32), dimension(n_points), intent(out) :: included_n_reps
            !! Stores the count of non-NaN replicates (included ones)
        integer(int32), intent(out) :: ierr
            !! Error code
        logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask
            !! Optional mask to exclude specific neighbors (e.g. for family-wise analysis)

        call set_ok(ierr)

        call validate_dimension_size(n_reps, ierr)
        call validate_dimension_size(n_neighbors, ierr)
        call validate_dimension_size(n_points, ierr)
        call validate_dimension_size(n_bins, ierr)
        call validate_in_range_real(shared_residual_range, ierr, min=0.0_real64)

        if (is_err(ierr)) return

        call build_residual_histograms_helper(neighborhood_residuals, n_reps, n_neighbors, n_points, shared_residual_range, n_bins, counts, pmf, included_n_reps, neighbor_mask)
    end subroutine build_residual_histograms

    !> (no input validation) Summarizes the neighborhood residuals in absolute histogram counts and probability mass functions `pmf(residual, bin)` (actually a matrix)
    pure subroutine build_residual_histograms_helper(neighborhood_residuals, n_reps, n_neighbors, n_points, shared_residual_range, n_bins, counts, pmf, included_n_reps, neighbor_mask)
        integer(int32), intent(in) :: n_reps
            !! Number of replicates of the study
        integer(int32), intent(in) :: n_neighbors
            !! Number of reference points (k)
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the study
        real(real64), dimension(n_reps, n_neighbors, n_points), intent(in) :: neighborhood_residuals
            !! Computed neighborhood residuals for a study ([[tox_jensen_shannon_test(module):construct_neighborhoods(subroutine)]]), NaN is explicitly allowed for missing values
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range (R) from [[tox_jensen_shannon_divergence(module):determine_shared_residual_range_alloc(subroutine)]]
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins in range [-R,R]
        integer(int32), dimension(n_points, n_bins), intent(out) :: counts
            !! Absolute counts of a residual per bin
        real(real64), dimension(n_points, n_bins), intent(out) :: pmf
            !! `counts` normalized to `0 <= counts(:, i) <= 1` and `sum(counts(:, i)) == 1`
        integer(int32), dimension(n_points), intent(out) :: included_n_reps
            !! Stores the count of non-NaN replicates (included ones)
        logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask
            !! Optional mask to exclude specific neighbors (e.g. for family-wise analysis)

        real(real64) :: bin_width, clamped_residual
        integer(int32) :: bin_idx, i_neighbor, i_rep, i_bin, included_reps, i_point
        logical :: filter_neighbors

        bin_width = 2.0_real64 * shared_residual_range / real(n_bins, real64)
        counts = 0_int32
        pmf = 0.0_real64

        filter_neighbors = present(neighbor_mask)

        ! 1. assign the bins to the residuals (increase the respective count)
        ! outer loop cannot be concurrent, as counts of same residuals and bins but different neighbors might be changed at the same time
        do concurrent (i_point = 1:n_points) local(included_reps) shared(included_n_reps)
            included_reps = 0_int32
            do concurrent (i_neighbor = 1:n_neighbors) &
                    local(i_rep, clamped_residual, bin_idx) &
                    shared(filter_neighbors, n_reps, counts, neighborhood_residuals, shared_residual_range, bin_width) &
                    reduce(+:included_reps)
                ! Exclude neighbor if desired
                if (filter_neighbors) then
                    if (.not. neighbor_mask(i_neighbor, i_point)) cycle
                end if

                ! Count non-NaNs and assign the to a bin
                do i_rep = 1, n_reps
                    if (.not. ieee_is_nan(neighborhood_residuals(i_rep, i_neighbor, i_point))) then
                        ! clamp residual to histogram range
                        clamped_residual = clamp(neighborhood_residuals(i_rep, i_neighbor, i_point), min_val=-shared_residual_range, max_val=shared_residual_range)

                        ! assign bin to residual
                        bin_idx = min(n_bins, int( (clamped_residual + shared_residual_range) / bin_width ) + 1)
                        counts(i_point, bin_idx) = counts(i_point, bin_idx) + 1

                        included_reps = included_reps + 1
                    end if
                end do
            end do
            included_n_reps(i_point) = included_reps
        end do

        ! 2. calculate pmf
        do concurrent (i_bin = 1:n_bins)
            do concurrent (i_point = 1:n_points) shared(pmf, i_bin, included_n_reps, counts)
                if (included_n_reps(i_point) == 0) then
                    pmf(i_point, i_bin) = 0.0_real64
                else
                    pmf(i_point, i_bin) = real(counts(i_point, i_bin), real64) / real(included_n_reps(i_point), real64)
                end if
            end do
        end do
    end subroutine build_residual_histograms_helper

    !> Having the probabilities `pmf` from [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]], this subroutine computes the Jensen-Shannon divergence per reference point/neighbor
    pure subroutine compute_divergence_per_reference_point(pmf_S1, pmf_S2, n_points, n_bins, js_divergences, ierr)
        integer(int32), intent(in) :: n_points
            !! Number of reference points (k)
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins in range [-R,R]
        real(real64), dimension(n_points, n_bins), intent(in) :: pmf_S1
            !! Computed normalized hostogram counts from [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]] for study 1
        real(real64), dimension(n_points, n_bins), intent(in) :: pmf_S2
            !! Computed normalized hostogram counts from [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]] for study 2
        real(real64), dimension(n_points), intent(out) :: js_divergences
            !! Jensen-Shannon divergence per reference point
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_points, ierr)
        call validate_dimension_size(n_bins, ierr)
        call validate_all_in_range_real(pmf_S1, size(pmf_S1, kind=int32), ierr, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_real(pmf_S2, size(pmf_S2, kind=int32), ierr, min=0.0_real64, max=1.0_real64)

        if (is_err(ierr)) return

        call compute_divergence_per_reference_point_helper(pmf_S1, pmf_S2, n_points, n_bins, js_divergences)
    end subroutine compute_divergence_per_reference_point

    !> (no input validation) Having the probabilities `pmf` from [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]], this subroutine computes the Jensen-Shannon divergence per reference point/neighbor
    pure subroutine compute_divergence_per_reference_point_helper(pmf_S1, pmf_S2, n_points, n_bins, js_divergences)
        integer(int32), intent(in) :: n_points
            !! Number of reference points (k)
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins in range [-R,R]
        real(real64), dimension(n_points, n_bins), intent(in) :: pmf_S1
            !! Computed normalized hostogram counts from [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]] for study 1
        real(real64), dimension(n_points, n_bins), intent(in) :: pmf_S2
            !! Computed normalized hostogram counts from [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]] for study 2
        real(real64), dimension(n_points), intent(out) :: js_divergences
            !! Jensen-Shannon divergence per reference point

        real(real64) :: S_mean, s1_val, s2_val
        integer(int32) :: i_bin, i_point

        js_divergences = 0.0_real64

        ! 1. compute the Knullback-Leibler (KL) divergences
        ! Note that the J-S divergence is defined as `0.5 * KL_S1 + 0.5 * KL_S2`, equivalent to `0.5 * (KL_S1 + KL_S2)`.
        ! Thus, instead of computing KL_S* independently, it accumulates directly in the `js_divergences` output.
        ! Another thing, switching the loops would enable both to run concurrently and the 0.5*js_divergences step could be done in one go,
        ! but cache locality still beats that, except for the case of thousands of neighbors, which might not be the common case.
        do i_bin = 1, n_bins
            do concurrent (i_point = 1:n_points) local(s1_val, s2_val, S_mean) shared(i_bin, pmf_S1, pmf_S2, js_divergences)
                s1_val = pmf_S1(i_point, i_bin)
                s2_val = pmf_S2(i_point, i_bin)
                S_mean = 0.5_real64 * (s1_val + s2_val)

                if (.not. is_close(S_mean, 0.0_real64)) then
                    if (s1_val > 0.0_real64) then
                        js_divergences(i_point) = js_divergences(i_point) + s1_val * log(s1_val / S_mean)
                    end if

                    if (s2_val > 0.0_real64) then
                        js_divergences(i_point) = js_divergences(i_point) + s2_val * log(s2_val / S_mean)
                    end if
                end if
            end do
        end do

        ! 2. Compute the js_divergences
        do concurrent (i_point = 1:n_points) shared(js_divergences)
            js_divergences(i_point) = 0.5_real64 * js_divergences(i_point)
        end do
    end subroutine compute_divergence_per_reference_point_helper

    !> Computes the global weighted Jensen-Shannon divergence from the per-neighbor divergences calculated by [[tox_jensen_shannon_divergence(module):compute_divergence_per_reference_point(subroutine)]]
    pure subroutine compute_weighted_global_divergence(js_divergences, n_points, included_n_reps_S1, included_n_reps_S2, global_js_divergence, weights, ierr)
        integer(int32), intent(in) :: n_points
            !! Number of reference points (k)
        real(real64), dimension(n_points), intent(in) :: js_divergences
            !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
        integer(int32), dimension(n_points), intent(in) :: included_n_reps_S1
            !! Count of non-NaN residuals (included ones) in study 1 (obtained from [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]])
        integer(int32), dimension(n_points), intent(in) :: included_n_reps_S2
            !! Count of non-NaN residuals (included ones) in study 2 (obtained from [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]])
        real(real64), intent(out) :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        real(real64), dimension(n_points), intent(out) :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_points, ierr)
        call validate_all_in_range_real(js_divergences, size(js_divergences, kind=int32), ierr, min=0.0_real64)
        call validate_all_in_range_int(included_n_reps_S1, size(included_n_reps_S1, kind=int32), ierr, min=0_int32)
        call validate_all_in_range_int(included_n_reps_S2, size(included_n_reps_S2, kind=int32), ierr, min=0_int32)

        if (is_err(ierr)) return

        call compute_weighted_global_divergence_helper(js_divergences, n_points, included_n_reps_S1, included_n_reps_S2, global_js_divergence, weights)
    end subroutine compute_weighted_global_divergence

    !> (no input validation) Computes the global weighted Jensen-Shannon divergence from the per-neighbor divergences calculated by [[tox_jensen_shannon_divergence(module):compute_divergence_per_reference_point(subroutine)]]
    pure subroutine compute_weighted_global_divergence_helper(js_divergences, n_points, included_n_reps_S1, included_n_reps_S2, global_js_divergence, weights)
        integer(int32), intent(in) :: n_points
            !! Number of reference points (k)
        real(real64), dimension(n_points), intent(in) :: js_divergences
            !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
        integer(int32), dimension(n_points), intent(in) :: included_n_reps_S1
            !! Count of non-NaN residuals (included ones) in study 1 (obtained from [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]])
        integer(int32), dimension(n_points), intent(in) :: included_n_reps_S2
            !! Count of non-NaN residuals (included ones) in study 2 (obtained from [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]])
        real(real64), intent(out) :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        real(real64), dimension(n_points), intent(out) :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`

        integer(int32) :: i_point, included_reps
        real(real64) :: total_sample_count

        global_js_divergence = 0.0_real64

        total_sample_count = real(sum(included_n_reps_S1) + sum(included_n_reps_S2), real64)

        if (is_close(total_sample_count, 0.0_real64)) then
            weights = 0.0_real64
        else
            ! Calculate the global Jensen-Shannon divergence
            do concurrent (i_point = 1:n_points) local(included_reps) shared(weights, included_n_reps_S1, included_n_reps_S2, total_sample_count) reduce(+:global_js_divergence)
                included_reps = included_n_reps_S1(i_point) + included_n_reps_S2(i_point)

                weights(i_point) = real(included_reps, real64) / total_sample_count

                global_js_divergence = global_js_divergence + weights(i_point) * js_divergences(i_point)
            end do
        end if
    end subroutine compute_weighted_global_divergence_helper

    !> Estimates how likely the observed divergence is to occur by chance under the null hypothesis that both studies are exchangeable
    subroutine gjct_permutation_test_alloc(neighborhood_residuals_S1, neighborhood_residuals_S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, global_jsd_observed, n_bins, shared_residual_range, n_permutations, jsd_null, p_value, ierr, random_seed, neighbor_mask_S1, neighbor_mask_S2)
        integer(int32), intent(in) :: n_reps_S1
            !! Number of replicates in study 1
        integer(int32), intent(in) :: n_reps_S2
            !! Number of replicates in study 2
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors in the studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1 ([[tox_jensen_shannon_test(module):construct_neighborhoods(subroutine)]]), NaN is explicitly allowed for missing values
        real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2 ([[tox_jensen_shannon_test(module):construct_neighborhoods(subroutine)]]), NaN is explicitly allowed for missing values
        real(real64), intent(in) :: global_jsd_observed
            !! Observed global JSD value for both studies (from [[tox_jensen_shannon_divergence(module):compute_weighted_global_divergence(subroutine)]])
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins used for the studies in [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]]
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range for both studies, from [[tox_jensen_shannon_divergence(module):determine_shared_residual_range(subroutine)]]
        integer(int32), intent(in) :: n_permutations
            !! Number of permutations to perform
        real(real64), dimension(n_permutations), intent(out) :: jsd_null
            !! Vector of global divergence values obtained under the null hypothesis
        real(real64), intent(out) :: p_value
            !! Empirical p-value of the permutation test: \( \frac{\text{count}(jsd\_null \ge global\_jsd\_observed) + 1}{n\_permutations} \)
        integer(int32), intent(out) :: ierr
            !! Error code
        integer(int32), intent(in), optional :: random_seed
            !! Seed to use for shuffling
        logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S1
            !! Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
        logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S2
            !! Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)

        real(real64), dimension(:, :, :), allocatable :: S1
        real(real64), dimension(:, :, :), allocatable :: S2
        real(real64), dimension(:, :), allocatable :: tmp_pool
        integer(int32), dimension(:, :), allocatable :: tmp_counts
        real(real64), dimension(:, :), allocatable :: tmp_pmf_S1
        real(real64), dimension(:, :), allocatable :: tmp_pmf_S2
        integer(int32), dimension(:), allocatable :: tmp_included_n_reps_S1
        integer(int32), dimension(:), allocatable :: tmp_included_n_reps_S2
        real(real64), dimension(:), allocatable :: tmp_js_divergences
        real(real64), dimension(:), allocatable :: tmp_weights

        call set_ok(ierr)

        call validate_dimension_size(n_reps_S1, ierr)
        call validate_dimension_size(n_reps_S2, ierr)
        call validate_dimension_size(n_neighbors, ierr)
        call validate_dimension_size(n_points, ierr)
        call validate_in_range_int(n_bins, ierr, min=1_int32)

        if (is_err(ierr)) return

        ! Allocate working arrays
        M_ALLOCATE(tmp_pool(n_reps_S1 + n_reps_S2, n_neighbors))
        M_ALLOCATE(tmp_counts(n_points, n_bins))
        M_ALLOCATE(tmp_pmf_S1(n_points, n_bins))
        M_ALLOCATE(tmp_pmf_S2(n_points, n_bins))
        M_ALLOCATE(tmp_included_n_reps_S1(n_points))
        M_ALLOCATE(tmp_included_n_reps_S2(n_points))
        M_ALLOCATE(tmp_js_divergences(n_points))
        M_ALLOCATE(tmp_weights(n_points))

        ! Allocate the residual copies
        M_ALLOCATE(S1(n_reps_S1, n_neighbors, n_points))
        M_ALLOCATE(S2(n_reps_S2, n_neighbors, n_points))

        S1 = neighborhood_residuals_S1
        S2 = neighborhood_residuals_S2

        call gjct_permutation_test(S1, S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, global_jsd_observed, n_bins, shared_residual_range, n_permutations, jsd_null, p_value, tmp_pool, tmp_pmf_S1, tmp_pmf_S2, tmp_counts, tmp_included_n_reps_S1, tmp_included_n_reps_S2, tmp_js_divergences, tmp_weights, ierr, random_seed, neighbor_mask_S1, neighbor_mask_S2)
    end subroutine gjct_permutation_test_alloc

    !> Estimates how likely the observed divergence is to occur by chance under the null hypothesis that both studies are exchangeable
    subroutine gjct_permutation_test( &
            neighborhood_residuals_S1_copy, neighborhood_residuals_S2_copy, n_reps_S1, n_reps_S2, n_neighbors, n_points, global_jsd_observed, n_bins, shared_residual_range, n_permutations, jsd_null, p_value, &
            tmp_pool, tmp_pmf_S1, tmp_pmf_S2, tmp_counts, tmp_included_n_reps_S1, tmp_included_n_reps_S2, tmp_js_divergences, tmp_weights, &
            ierr, random_seed, neighbor_mask_S1, neighbor_mask_S2 &
        )
        integer(int32), intent(in) :: n_reps_S1
            !! Number of replicates in study 1
        integer(int32), intent(in) :: n_reps_S2
            !! Number of replicates in study 2
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors in the studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(inout) :: neighborhood_residuals_S1_copy
            !! Copy (if wanted) of the computed neighborhood residuals for study 1, will be shuffled in-place
        real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(inout) :: neighborhood_residuals_S2_copy
            !! Copy (if wanted) of the computed neighborhood residuals for study 2, will be shuffled in-place
        real(real64), intent(in) :: global_jsd_observed
            !! Observed global JSD value for both studies (from [[tox_jensen_shannon_divergence(module):compute_weighted_global_divergence(subroutine)]])
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins used for the studies in [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]]
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range for both studies, from [[tox_jensen_shannon_divergence(module):determine_shared_residual_range(subroutine)]]
        integer(int32), intent(in) :: n_permutations
            !! Number of permutations to perform
        real(real64), dimension(n_permutations), intent(out) :: jsd_null
            !! Vector of global divergence values obtained under the null hypothesis
        real(real64), intent(out) :: p_value
            !! Empirical p-value of the permutation test: \( \frac{\text{count}(jsd\_null \ge global\_jsd\_observed) + 1}{n\_permutations} \)
        real(real64), dimension(n_reps_S1 + n_reps_S2, n_neighbors), intent(out) :: tmp_pool
            !! Working array for shuffling the concatenated residuals from both studies per reference point
        real(real64), dimension(n_points, n_bins), intent(out) :: tmp_pmf_S1
            !! Absolute counts of a residual per bin obtained from [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]]
        real(real64), dimension(n_points, n_bins), intent(out) :: tmp_pmf_S2
            !! Absolute counts of a residual per bin obtained from [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]]
        integer(int32), dimension(n_points, n_bins), intent(out) :: tmp_counts
            !! Working array for [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]]
        integer(int32), dimension(n_points), intent(out) :: tmp_included_n_reps_S1
            !! Working array for [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]]
        integer(int32), dimension(n_points), intent(out) :: tmp_included_n_reps_S2
            !! Working array for [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]]
        real(real64), dimension(n_points), intent(out) :: tmp_js_divergences
            !! Working array for [[tox_jensen_shannon_divergence(module):compute_divergence_per_reference_point(subroutine)]]
        real(real64), dimension(n_points), intent(out) :: tmp_weights
            !! Working array for [[tox_jensen_shannon_divergence(module):compute_weighted_global_divergence(subroutine)]]
        integer(int32), intent(out) :: ierr
            !! Error code
        integer(int32), intent(in), optional :: random_seed
            !! Seed to use for shuffling
        logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S1
            !! Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
        logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S2
            !! Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)

        call set_ok(ierr)

        call validate_dimension_size(n_reps_S1, ierr)
        call validate_dimension_size(n_reps_S2, ierr)
        call validate_dimension_size(n_neighbors, ierr)
        call validate_dimension_size(n_points, ierr)
        call validate_in_range_int(n_bins, ierr, min=1_int32)
        call validate_in_range_int(n_permutations, ierr, min=1_int32)
        call validate_in_range_real(shared_residual_range, ierr, min=0.0_real64)

        if (is_err(ierr)) return

        call gjct_permutation_test_helper(neighborhood_residuals_S1_copy, neighborhood_residuals_S2_copy, n_reps_S1, n_reps_S2, n_neighbors, n_points, global_jsd_observed, n_bins, shared_residual_range, n_permutations, jsd_null, p_value, tmp_pool, tmp_pmf_S1, tmp_pmf_S2, tmp_counts, tmp_included_n_reps_S1, tmp_included_n_reps_S2, tmp_js_divergences, tmp_weights, random_seed, neighbor_mask_S1, neighbor_mask_S2)
    end subroutine gjct_permutation_test

    !> (no input validation) Estimates how likely the observed divergence is to occur by chance under the null hypothesis that both studies are exchangeable
    subroutine gjct_permutation_test_helper( &
            neighborhood_residuals_S1_copy, neighborhood_residuals_S2_copy, n_reps_S1, n_reps_S2, n_neighbors, n_points, global_jsd_observed, n_bins, shared_residual_range, n_permutations, jsd_null, p_value, &
            tmp_pool, tmp_pmf_S1, tmp_pmf_S2, tmp_counts, tmp_included_n_reps_S1, tmp_included_n_reps_S2, tmp_js_divergences, tmp_weights, &
            random_seed, neighbor_mask_S1, neighbor_mask_S2 &
        )
        integer(int32), intent(in) :: n_reps_S1
            !! Number of replicates in study 1
        integer(int32), intent(in) :: n_reps_S2
            !! Number of replicates in study 2
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors in the studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(inout) :: neighborhood_residuals_S1_copy
            !! Copy (if wanted) of the computed neighborhood residuals for study 1, will be shuffled in-place
        real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(inout) :: neighborhood_residuals_S2_copy
            !! Copy (if wanted) of the computed neighborhood residuals for study 2, will be shuffled in-place
        real(real64), intent(in) :: global_jsd_observed
            !! Observed global JSD value for both studies (from [[tox_jensen_shannon_divergence(module):compute_weighted_global_divergence(subroutine)]])
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins used for the studies in [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]]
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range for both studies, from [[tox_jensen_shannon_divergence(module):determine_shared_residual_range(subroutine)]]
        integer(int32), intent(in) :: n_permutations
            !! Number of permutations to perform
        real(real64), dimension(n_permutations), intent(out) :: jsd_null
            !! Vector of global divergence values obtained under the null hypothesis
        real(real64), intent(out) :: p_value
            !! Empirical p-value of the permutation test: \( \frac{\text{count}(jsd\_null \ge global\_jsd\_observed) + 1}{n\_permutations} \)
        real(real64), dimension(n_reps_S1 + n_reps_S2, n_neighbors), intent(out) :: tmp_pool
            !! Working array for shuffling the concatenated residuals from both studies per reference point
        real(real64), dimension(n_points, n_bins), intent(out) :: tmp_pmf_S1
            !! Absolute counts of a residual per bin obtained from [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]]
        real(real64), dimension(n_points, n_bins), intent(out) :: tmp_pmf_S2
            !! Absolute counts of a residual per bin obtained from [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]]
        integer(int32), dimension(n_points, n_bins), intent(out) :: tmp_counts
            !! Working array for [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]]
        integer(int32), dimension(n_points), intent(out) :: tmp_included_n_reps_S1
            !! Working array for [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]]
        integer(int32), dimension(n_points), intent(out) :: tmp_included_n_reps_S2
            !! Working array for [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]]
        real(real64), dimension(n_points), intent(out) :: tmp_js_divergences
            !! Working array for [[tox_jensen_shannon_divergence(module):compute_divergence_per_reference_point(subroutine)]]
        real(real64), dimension(n_points), intent(out) :: tmp_weights
            !! Working array for [[tox_jensen_shannon_divergence(module):compute_weighted_global_divergence(subroutine)]]
        integer(int32), intent(in), optional :: random_seed
            !! Seed to use for shuffling
        logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S1
            !! Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
        logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S2
            !! Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)

        integer(int32) :: n_jsd_exceeding_observed, i_permutation, n_residuals_S1, n_residuals_S2, pool_size, i_point

        if (present(random_seed)) then
            call init_random(random_seed)
        end if

        n_jsd_exceeding_observed = 0_int32
        do i_permutation = 1, n_permutations
            ! 1. shuffle residuals
            do i_point = 1, n_points
                call shuffle_reference_point_helper( &
                    neighborhood_residuals_S1_copy(:, :, i_point), neighborhood_residuals_S2_copy(:, :, i_point), &
                    n_reps_S1, n_reps_S2, n_neighbors, tmp_pool &
                )
            end do

            ! 2. Pipeline to determine the global jsd for current permutation
            call jct_compute_jsd_pipeline_helper(neighborhood_residuals_S1_copy, neighborhood_residuals_S2_copy, n_reps_S1, n_reps_S2, n_neighbors, n_points, n_bins, shared_residual_range, tmp_js_divergences, tmp_included_n_reps_S1, tmp_included_n_reps_S2, jsd_null(i_permutation), tmp_weights, tmp_pmf_S1, tmp_pmf_S2, tmp_counts, neighbor_mask_S1, neighbor_mask_S2)

            if (jsd_null(i_permutation) >= global_jsd_observed) then
                n_jsd_exceeding_observed = n_jsd_exceeding_observed + 1
            end if
        end do

        p_value = real(n_jsd_exceeding_observed + 1, real64) / real(n_permutations + 1, real64)
    end subroutine gjct_permutation_test_helper

    !> Helper for [[tox_jensen_shannon_divergence(module):gjct_permutation_test_helper(subroutine)]] to shuffle reference points
    subroutine shuffle_reference_point_helper(reference_point_S1, reference_point_S2, n_reps_S1, n_reps_S2, n_neighbors, pool_flat)
        integer(int32), intent(in) :: n_reps_S1
            !! Number of replicates in study 1
        integer(int32), intent(in) :: n_reps_S2
            !! Number of replicates in study 2
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors in the studies
        real(real64), dimension(n_reps_S1 * n_neighbors), intent(inout) :: reference_point_S1
            !! Residuals for one reference point in study 1, will be shuffled in-place
        real(real64), dimension(n_reps_S2 * n_neighbors), intent(inout) :: reference_point_S2
            !! Residuals for one reference point in study 2, will be shuffled in-place
        real(real64), dimension((n_reps_S1 + n_reps_S2) * n_neighbors), intent(out) :: pool_flat
            !! Working array for shuffling the concatenated residuals from both studies per reference point

        integer(int32) :: pool_size, n_residuals_S1

        pool_size = size(pool_flat, kind=int32)
        n_residuals_S1 = size(reference_point_S1, kind=int32)

        pool_flat(1:n_residuals_S1) = reference_point_S1
        pool_flat(n_residuals_S1+1:pool_size) = reference_point_S2

        call shuffle_vector(pool_flat)

        reference_point_S1 = pool_flat(1:n_residuals_S1)
        reference_point_S2 = pool_flat(n_residuals_S1+1:pool_size)
    end subroutine shuffle_reference_point_helper

    !> TODO: #118
    pure subroutine fjct_compute_jsd_alloc(family_idx, gene_to_family_S1, gene_to_family_S2, n_genes_S1, n_genes_S2, neighborhood_residuals_S1, neighborhood_residuals_S2, &
            neighborhood_genes_S1, neighborhood_genes_S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, resid_S1, resid_S2, n_bins, shared_residual_range, js_divergences, &
            included_n_reps_S1, included_n_reps_S2, global_js_divergence, weights, ierr &
        )
        integer(int32), intent(in) :: n_genes_S1
            !! Number of genes in study 1
        integer(int32), intent(in) :: n_genes_S2
            !! Number of genes in study 2
        integer(int32), intent(in) :: n_reps_S1
            !! Number of replicates in study 1
        integer(int32), intent(in) :: n_reps_S2
            !! Number of replicates in study 2
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors in the studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        integer(int32), intent(in) :: family_idx
            !! Index of the family that should be analyzed
        integer(int32), dimension(n_genes_S1), intent(in) :: gene_to_family_S1
            !! Mapping for study 1: Each index (gene) holds the index of its family
        integer(int32), dimension(n_genes_S2), intent(in) :: gene_to_family_S2
            !! Mapping for study 2: Each index (gene) holds the index of its family
        real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1 ([[tox_jensen_shannon_test(module):construct_neighborhoods(subroutine)]]), NaN is explicitly allowed for missing values
        real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2 ([[tox_jensen_shannon_test(module):construct_neighborhoods(subroutine)]]), NaN is explicitly allowed for missing values
        integer(int32), dimension(n_neighbors, n_points), intent(in) :: neighborhood_genes_S1
            !! Indices of selected neighborhood genes, obtained from `neighborhood_indices` of [[tox_jensen_shannon_test(module):construct_neighborhoods(subroutine)]]
        integer(int32), dimension(n_neighbors, n_points), intent(in) :: neighborhood_genes_S2
            !! Indices of selected neighborhood genes, obtained from `neighborhood_indices` of [[tox_jensen_shannon_test(module):construct_neighborhoods(subroutine)]]
        real(real64), dimension(n_reps_S1, n_genes_S1), intent(in) :: resid_S1
            !! Matrix of signed residuals for study 1, from [[tox_jensen_shannon_test(module):compute_residuals(subroutine)]]
        real(real64), dimension(n_reps_S2, n_genes_S2), intent(in) :: resid_S2
            !! Matrix of signed residuals for study 1, from [[tox_jensen_shannon_test(module):compute_residuals(subroutine)]]
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins used for the studies in [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]]
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range for both studies, from [[tox_jensen_shannon_divergence(module):determine_shared_residual_range(subroutine)]]
        real(real64), dimension(n_points), intent(out) :: js_divergences
            !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
        integer(int32), dimension(n_points), intent(out) :: included_n_reps_S1
            !! Count of non-NaN residuals (included ones) in study 1 (obtained from [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]])
        integer(int32), dimension(n_points), intent(out) :: included_n_reps_S2
            !! Count of non-NaN residuals (included ones) in study 2 (obtained from [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]])
        real(real64), intent(out) :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        real(real64), dimension(n_points), intent(out) :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i_point, i_neighbor
        logical, dimension(:, :), allocatable :: neighbor_mask_S1, neighbor_mask_S2
        real(real64), dimension(:, :), allocatable :: pmf_S1, pmf_S2
        integer(int32), dimension(:, :), allocatable :: tmp_counts

        call set_ok(ierr)

        call validate_dimension_size(n_neighbors, ierr)
        call validate_in_range_int(family_idx, ierr, min=1_int32)
        call validate_all_in_range_int(gene_to_family_S1, n_genes_S1, ierr, min=0_int32)
        call validate_all_in_range_int(gene_to_family_S2, n_genes_S2, ierr, min=0_int32)
        call validate_all_in_range_int(neighborhood_genes_S1, size(neighborhood_genes_S1, kind=int32), ierr, min=1_int32, max=n_genes_S1)
        call validate_all_in_range_int(neighborhood_genes_S2, size(neighborhood_genes_S2, kind=int32), ierr, min=1_int32, max=n_genes_S2)

        if (is_err(ierr)) return

        M_ALLOCATE(neighbor_mask_S1(n_neighbors, n_points))
        M_ALLOCATE(neighbor_mask_S2(n_neighbors, n_points))
        M_ALLOCATE(pmf_S1(n_points, n_bins))
        M_ALLOCATE(pmf_S2(n_points, n_bins))
        M_ALLOCATE(tmp_counts(n_points, n_bins))

        ! Set up mask for filtered analysis -> only include neighbors being part of the family
        do concurrent (i_point = 1:n_points)
            do concurrent (i_neighbor = 1:n_neighbors) shared(family_idx, neighbor_mask_S1, neighbor_mask_S2, gene_to_family_S1, neighborhood_genes_S1, gene_to_family_S2, neighborhood_genes_S2)
                neighbor_mask_S1(i_neighbor, i_point) = gene_to_family_S1(neighborhood_genes_S1(i_neighbor, i_point)) == family_idx
                neighbor_mask_S2(i_neighbor, i_point) = gene_to_family_S2(neighborhood_genes_S2(i_neighbor, i_point)) == family_idx
            end do
        end do

        call fjct_compute_jsd(neighborhood_residuals_S1, neighborhood_residuals_S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, neighbor_mask_S1, neighbor_mask_S2, n_bins, shared_residual_range, js_divergences, included_n_reps_S1, included_n_reps_S2, global_js_divergence, weights, pmf_S1, pmf_S2, tmp_counts, ierr)
    end subroutine fjct_compute_jsd_alloc

    !> TODO: #118
    pure subroutine fjct_compute_jsd(neighborhood_residuals_S1, neighborhood_residuals_S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, neighbor_mask_S1, neighbor_mask_S2, n_bins, shared_residual_range, js_divergences, included_n_reps_S1, included_n_reps_S2, global_js_divergence, weights, pmf_S1, pmf_S2, tmp_counts, ierr)
        integer(int32), intent(in) :: n_reps_S1
            !! Number of replicates in study 1
        integer(int32), intent(in) :: n_reps_S2
            !! Number of replicates in study 2
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors in the studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1 ([[tox_jensen_shannon_test(module):construct_neighborhoods(subroutine)]]), NaN is explicitly allowed for missing values
        real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2 ([[tox_jensen_shannon_test(module):construct_neighborhoods(subroutine)]]), NaN is explicitly allowed for missing values
        logical, dimension(n_neighbors, n_points), intent(in) :: neighbor_mask_S1
            !! Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
        logical, dimension(n_neighbors, n_points), intent(in) :: neighbor_mask_S2
            !! Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins used for the studies in [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]]
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range for both studies, from [[tox_jensen_shannon_divergence(module):determine_shared_residual_range(subroutine)]]
        real(real64), dimension(n_points), intent(out) :: js_divergences
            !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
        integer(int32), dimension(n_points), intent(out) :: included_n_reps_S1
            !! Count of non-NaN residuals (included ones) in study 1 (obtained from [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]])
        integer(int32), dimension(n_points), intent(out) :: included_n_reps_S2
            !! Count of non-NaN residuals (included ones) in study 2 (obtained from [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]])
        real(real64), intent(out) :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        real(real64), dimension(n_points), intent(out) :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
        real(real64), dimension(n_points, n_bins), intent(out) :: pmf_S1
            !! Absolute counts of a residual per bin obtained from [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]]
        real(real64), dimension(n_points, n_bins), intent(out) :: pmf_S2
            !! Absolute counts of a residual per bin obtained from [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]]
        integer(int32), dimension(n_points, n_bins), intent(out) :: tmp_counts
            !! Working array for [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]]
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_reps_S1, ierr)
        call validate_dimension_size(n_reps_S2, ierr)
        call validate_dimension_size(n_neighbors, ierr)
        call validate_dimension_size(n_points, ierr)
        call validate_in_range_int(n_bins, ierr, min=1_int32)
        call validate_in_range_real(shared_residual_range, ierr, min=0.0_real64)

        if (is_err(ierr)) return

        call jct_compute_jsd_pipeline_helper(neighborhood_residuals_S1, neighborhood_residuals_S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, n_bins, shared_residual_range, js_divergences, included_n_reps_S1, included_n_reps_S2, global_js_divergence, weights, pmf_S1, pmf_S2, tmp_counts, neighbor_mask_S1, neighbor_mask_S2)
    end subroutine fjct_compute_jsd

    !> TODO: #118
    pure subroutine jct_compute_jsd_pipeline_helper(neighborhood_residuals_S1, neighborhood_residuals_S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, n_bins, shared_residual_range, js_divergences, included_n_reps_S1, included_n_reps_S2, global_js_divergence, weights, pmf_S1, pmf_S2, tmp_counts, neighbor_mask_S1, neighbor_mask_S2)
        integer(int32), intent(in) :: n_reps_S1
            !! Number of replicates in study 1
        integer(int32), intent(in) :: n_reps_S2
            !! Number of replicates in study 2
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors in the studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1 ([[tox_jensen_shannon_test(module):construct_neighborhoods(subroutine)]]), NaN is explicitly allowed for missing values
        real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2 ([[tox_jensen_shannon_test(module):construct_neighborhoods(subroutine)]]), NaN is explicitly allowed for missing values
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins used for the studies in [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]]
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range for both studies, from [[tox_jensen_shannon_divergence(module):determine_shared_residual_range(subroutine)]]
        real(real64), dimension(n_points), intent(out) :: js_divergences
            !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
        integer(int32), dimension(n_points), intent(out) :: included_n_reps_S1
            !! Count of non-NaN residuals (included ones) in study 1 (obtained from [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]])
        integer(int32), dimension(n_points), intent(out) :: included_n_reps_S2
            !! Count of non-NaN residuals (included ones) in study 2 (obtained from [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]])
        real(real64), intent(out) :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        real(real64), dimension(n_points), intent(out) :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
        real(real64), dimension(n_points, n_bins), intent(out) :: pmf_S1
            !! Absolute counts of a residual per bin obtained from [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]]
        real(real64), dimension(n_points, n_bins), intent(out) :: pmf_S2
            !! Absolute counts of a residual per bin obtained from [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]]
        integer(int32), dimension(n_points, n_bins), intent(out) :: tmp_counts
            !! Working array for [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]]
        logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S1
            !! Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
        logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S2
            !! Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)

        call build_residual_histograms_helper(neighborhood_residuals_S1, n_reps_S1, n_neighbors, n_points, shared_residual_range, n_bins, tmp_counts, pmf_S1, included_n_reps_S1, neighbor_mask_S1)
        call build_residual_histograms_helper(neighborhood_residuals_S2, n_reps_S2, n_neighbors, n_points, shared_residual_range, n_bins, tmp_counts, pmf_S2, included_n_reps_S2, neighbor_mask_S2)
        call compute_divergence_per_reference_point_helper(pmf_S1, pmf_S2, n_points, n_bins, js_divergences)
        call compute_weighted_global_divergence_helper(js_divergences, n_points, included_n_reps_S1, included_n_reps_S2, global_js_divergence, weights)
    end subroutine jct_compute_jsd_pipeline_helper

end module tox_jensen_shannon_divergence

!> C-compatible wrapper for [[tox_jensen_shannon_divergence(module):determine_shared_residual_range(subroutine)]]
pure subroutine determine_shared_residual_range_expert_c( &
    abs_residual_pool, abs_residual_pool_perm, pool_size, &
    residual_range_quantile, shared_residual_range, ierr &
    ) bind(C, name="determine_shared_residual_range_expert_c")

    use tox_jensen_shannon_divergence, only: determine_shared_residual_range
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: pool_size
        !! Size of pool of residuals `abs_residual_pool`, usually `(n_reps_S1 + n_reps_2)*n_neighbors*n_points`
    real(c_double), intent(in), target :: residual_range_quantile
        !! Quantile for determining the residual range, default: 95.0
    real(c_double), intent(out), target :: shared_residual_range
        !! Computed residual range (R)
    real(c_double), dimension(pool_size), intent(in), target :: abs_residual_pool
        !! The absolute residual values of the concatenated S1,S2 residuals
    integer(c_int), dimension(pool_size), intent(in), target :: abs_residual_pool_perm
        !! The permutation vector that sorts `abs_residual_pool`
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(pool_size)
    M_CHECK_NON_NULL(residual_range_quantile)
    M_CHECK_NON_NULL(shared_residual_range)
    M_CHECK_NON_NULL(abs_residual_pool)
    M_CHECK_NON_NULL(abs_residual_pool_perm)

    call determine_shared_residual_range( &
        abs_residual_pool, abs_residual_pool_perm, pool_size, &
        shared_residual_range, ierr, residual_range_quantile )

end subroutine determine_shared_residual_range_expert_c

!> C-compatible wrapper for [[tox_jensen_shannon_divergence(module):determine_shared_residual_range_alloc(subroutine)]]
pure subroutine determine_shared_residual_range_c( &
    neighborhood_residuals_S1, neighborhood_residuals_S2, &
    n_reps_S1, n_reps_S2, n_neighbors, n_points, &
    residual_range_quantile, &
    shared_residual_range, ierr ) &
    bind(C, name="determine_shared_residual_range_c")

    use tox_jensen_shannon_divergence, only: determine_shared_residual_range_alloc
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_reps_S1
        !! Number of replicates in study 1
    integer(c_int), intent(in), target :: n_reps_S2
        !! Number of replicates in study 2
    integer(c_int), intent(in), target :: n_neighbors
        !! Number of reference points (k)
    integer(c_int), intent(in), target :: n_points
        !! Number of reference points in the studies
    real(c_double), dimension(n_reps_S1, n_neighbors, n_points), intent(in), target :: neighborhood_residuals_S1
        !! Computed neighborhood residuals for study 1 ([[tox_jensen_shannon_test(module):construct_neighborhoods(subroutine)]]), NaN is explicitly allowed for missing values
    real(c_double), dimension(n_reps_S2, n_neighbors, n_points), intent(in), target :: neighborhood_residuals_S2
        !! Computed neighborhood residuals for study 2 ([[tox_jensen_shannon_test(module):construct_neighborhoods(subroutine)]]), NaN is explicitly allowed for missing values
    real(c_double), intent(in), target :: residual_range_quantile
        !! Quantile for determining the residual range, default: 95.0
    real(c_double), intent(out), target :: shared_residual_range
        !! Computed residual range (R)
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_reps_S1)
    M_CHECK_NON_NULL(n_reps_S2)
    M_CHECK_NON_NULL(n_neighbors)
    M_CHECK_NON_NULL(n_points)
    M_CHECK_NON_NULL(neighborhood_residuals_S1)
    M_CHECK_NON_NULL(neighborhood_residuals_S2)
    M_CHECK_NON_NULL(residual_range_quantile)
    M_CHECK_NON_NULL(shared_residual_range)

    call determine_shared_residual_range_alloc( &
        neighborhood_residuals_S1, neighborhood_residuals_S2, &
        n_reps_S1, n_reps_S2, n_neighbors, n_points, &
        shared_residual_range, ierr, residual_range_quantile )

end subroutine determine_shared_residual_range_c

!> C-compatible wrapper for [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]]
pure subroutine build_residual_histograms_c( &
    neighborhood_residuals, &
    n_reps, n_neighbors, n_points, &
    shared_residual_range, &
    n_bins, &
    counts, pmf, included_n_reps, &
    ierr ) &
    bind(C, name="build_residual_histograms_c")

    use tox_jensen_shannon_divergence, only: build_residual_histograms
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_reps
        !! Number of replicates in the study
    integer(c_int), intent(in), target :: n_neighbors
        !! Number of reference points (k)
    integer(c_int), intent(in), target :: n_points
        !! Number of reference points in the studies
    real(c_double), dimension(n_reps, n_neighbors, n_points), intent(in), target :: neighborhood_residuals
        !! Computed neighborhood residuals for a study ([[tox_jensen_shannon_test(module):construct_neighborhoods(subroutine)]]), NaN is explicitly allowed for missing values
    real(c_double), intent(in), target :: shared_residual_range
        !! Computed residual range (R) from [[tox_jensen_shannon_divergence(module):determine_shared_residual_range_alloc(subroutine)]]
    integer(c_int), intent(in), target :: n_bins
        !! Number of equally sized histogram bins in range [-R,R]
    integer(c_int), dimension(n_points, n_bins), intent(out), target :: counts
        !! Absolute counts of a residual per bin
    real(c_double), dimension(n_points, n_bins), intent(out), target :: pmf
        !! `counts` normalized to `0 <= counts(:, i) <= 1` and `sum(counts(:, i)) == 1`
    integer(c_int), dimension(n_points), intent(out), target :: included_n_reps
        !! Stores the count of non-NaN replicates (included ones)
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(neighborhood_residuals)
    M_CHECK_NON_NULL(n_reps)
    M_CHECK_NON_NULL(n_neighbors)
    M_CHECK_NON_NULL(n_points)
    M_CHECK_NON_NULL(shared_residual_range)
    M_CHECK_NON_NULL(n_bins)
    M_CHECK_NON_NULL(counts)
    M_CHECK_NON_NULL(pmf)
    M_CHECK_NON_NULL(included_n_reps)

    call build_residual_histograms( &
        neighborhood_residuals, &
        n_reps, n_neighbors, n_points, &
        shared_residual_range, &
        n_bins, &
        counts, pmf, included_n_reps, &
        ierr )

end subroutine build_residual_histograms_c

!> C-compatible wrapper for [[tox_jensen_shannon_divergence(module):compute_divergence_per_reference_point(subroutine)]]
pure subroutine compute_divergence_per_reference_point_c( &
    pmf_S1, pmf_S2, &
    n_points, n_bins, &
    js_divergences, ierr ) &
    bind(C, name="compute_divergence_per_reference_point_c")

    use tox_jensen_shannon_divergence, only: compute_divergence_per_reference_point
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_points
        !! Number of reference points (k)
    integer(c_int), intent(in), target :: n_bins
        !! Number of equally sized histogram bins in range [-R,R]
    real(c_double), dimension(n_points, n_bins), intent(in), target :: pmf_S1
        !! Computed normalized hostogram counts from [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]] for study 1
    real(c_double), dimension(n_points, n_bins), intent(in), target :: pmf_S2
        !! Computed normalized hostogram counts from [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]] for study 2
    real(c_double), dimension(n_points), intent(out), target :: js_divergences
        !! Jensen-Shannon divergence per reference point
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_points)
    M_CHECK_NON_NULL(n_bins)
    M_CHECK_NON_NULL(pmf_S1)
    M_CHECK_NON_NULL(pmf_S2)
    M_CHECK_NON_NULL(js_divergences)

    call compute_divergence_per_reference_point( &
        pmf_S1, pmf_S2, &
        n_points, n_bins, &
        js_divergences, ierr )

end subroutine compute_divergence_per_reference_point_c

!> C-compatible wrapper for [[tox_jensen_shannon_divergence(module):compute_weighted_global_divergence(subroutine)]]
pure subroutine compute_weighted_global_divergence_c( &
    js_divergences, &
    n_points, &
    included_n_reps_S1, included_n_reps_S2, &
    global_js_divergence, weights, &
    ierr ) &
    bind(C, name="compute_weighted_global_divergence_c")

    use tox_jensen_shannon_divergence, only: compute_weighted_global_divergence
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_points
        !! Number of reference points (k)
    real(c_double), dimension(n_points), intent(in), target :: js_divergences
        !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
    integer(c_int), dimension(n_points), intent(in), target :: included_n_reps_S1
        !! Count of non-NaN residuals (included ones) in study 1
    integer(c_int), dimension(n_points), intent(in), target :: included_n_reps_S2
        !! Count of non-NaN residuals (included ones) in study 2
    real(c_double), intent(out), target :: global_js_divergence
        !! Weighted global Jensen-Shannon divergence
    real(c_double), dimension(n_points), intent(out), target :: weights
        !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_points)
    M_CHECK_NON_NULL(js_divergences)
    M_CHECK_NON_NULL(included_n_reps_S1)
    M_CHECK_NON_NULL(included_n_reps_S2)
    M_CHECK_NON_NULL(global_js_divergence)
    M_CHECK_NON_NULL(weights)

    call compute_weighted_global_divergence( &
        js_divergences, n_points, &
        included_n_reps_S1, included_n_reps_S2, &
        global_js_divergence, weights, ierr )

end subroutine compute_weighted_global_divergence_c

!> C-compatible wrapper for [[tox_jensen_shannon_divergence(module):gjct_permutation_test_alloc(subroutine)]]
subroutine gjct_permutation_test_c( &
    neighborhood_residuals_S1, neighborhood_residuals_S2, &
    n_reps_S1, n_reps_S2, n_neighbors, n_points, &
    global_jsd_observed, n_bins, shared_residual_range, n_permutations, &
    jsd_null, p_value, ierr, random_seed) &
    bind(C, name="gjct_permutation_test_c")

    use tox_jensen_shannon_divergence, only: gjct_permutation_test_alloc
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_reps_S1
        !! Number of replicates in study 1
    integer(c_int), intent(in), target :: n_reps_S2
        !! Number of replicates in study 2
    integer(c_int), intent(in), target :: n_neighbors
        !! Number of neighbors in the studies
    integer(c_int), intent(in), target :: n_points
        !! Number of reference points in the studies
    real(c_double), dimension(n_reps_S1, n_neighbors, n_points), intent(in), target :: neighborhood_residuals_S1
        !! Computed neighborhood residuals for study 1 ([[tox_jensen_shannon_test(module):construct_neighborhoods(subroutine)]]), NaN is explicitly allowed for missing values
    real(c_double), dimension(n_reps_S2, n_neighbors, n_points), intent(in), target :: neighborhood_residuals_S2
        !! Computed neighborhood residuals for study 2 ([[tox_jensen_shannon_test(module):construct_neighborhoods(subroutine)]]), NaN is explicitly allowed for missing values
    real(c_double), intent(in), target :: global_jsd_observed
        !! Observed global JSD value for both studies (from [[tox_jensen_shannon_divergence(module):compute_weighted_global_divergence(subroutine)]])
    integer(c_int), intent(in), target :: n_bins
        !! Number of equally sized histogram bins used for the studies in [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]]
    real(c_double), intent(in), target :: shared_residual_range
        !! Computed residual range for both studies, from [[tox_jensen_shannon_divergence(module):determine_shared_residual_range(subroutine)]]
    integer(c_int), intent(in), target :: n_permutations
        !! Number of permutations to perform
    real(c_double), dimension(n_permutations), intent(out), target :: jsd_null
        !! Vector of global divergence values obtained under the null hypothesis
    real(c_double), intent(out), target :: p_value
        !! Empirical p-value of the permutation test: \( \frac{\text{count}(jsd\_null \ge global\_jsd\_observed) + 1}{n\_permutations} \)
    integer(c_int), intent(out), target :: ierr
        !! Error code
    integer(c_int), intent(in), target :: random_seed
        !! Seed to use for shuffling

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_reps_S1)
    M_CHECK_NON_NULL(n_reps_S2)
    M_CHECK_NON_NULL(n_neighbors)
    M_CHECK_NON_NULL(n_points)
    M_CHECK_NON_NULL(neighborhood_residuals_S1)
    M_CHECK_NON_NULL(neighborhood_residuals_S2)
    M_CHECK_NON_NULL(global_jsd_observed)
    M_CHECK_NON_NULL(n_bins)
    M_CHECK_NON_NULL(shared_residual_range)
    M_CHECK_NON_NULL(n_permutations)
    M_CHECK_NON_NULL(jsd_null)
    M_CHECK_NON_NULL(p_value)
    M_CHECK_NON_NULL(random_seed)

    call gjct_permutation_test_alloc( &
        neighborhood_residuals_S1, neighborhood_residuals_S2, &
        n_reps_S1, n_reps_S2, n_neighbors, n_points, &
        global_jsd_observed, n_bins, shared_residual_range, n_permutations, &
        jsd_null, p_value, ierr, random_seed)

end subroutine gjct_permutation_test_c
