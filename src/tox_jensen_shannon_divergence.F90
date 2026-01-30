#include "macros.h"

module tox_jensen_shannon_divergence
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use f42_utils, only: clamp, calc_percentile_helper, is_close, sort_array_heapsort
    use tox_errors, only: set_ok, set_err, is_err, ERR_ALLOC_FAIL, validate_dimension_size, validate_in_range_real, validate_all_in_range_real, validate_in_range_int, validate_all_in_range_int
    implicit none
contains

    !> Computes the shared residual range [-R, R] for the computed residuals from studies S1 and S2
    pure subroutine determine_shared_residual_range(neighborhood_residuals_S1, neighborhood_residuals_S2, n_residuals, n_neighbors, shared_residual_range, tmp_abs_residual_pool, tmp_abs_residual_pool_perm, ierr, residual_range_quantile)
        integer(int32), intent(in) :: n_residuals
            !! Number of residuals
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors (k)
        real(real64), dimension(n_residuals, n_neighbors), intent(in) :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1 (kNN), NaN is explicitly allowed for missing values
        real(real64), dimension(n_residuals, n_neighbors), intent(in) :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2 (kNN), NaN is explicitly allowed for missing values
        real(real64), intent(in), optional :: residual_range_quantile
            !! Quantile for determining the residual range, default: 95.0
        real(real64), intent(out) :: shared_residual_range
            !! Computed residual range (R)
        real(real64), dimension(2 * size(neighborhood_residuals_S1, kind=int32)), intent(out) :: tmp_abs_residual_pool
            !! Work array for the absolute residual values of the concatenated S1,S2 residuals
        integer(int32), dimension(size(tmp_abs_residual_pool, kind=int32)), intent(out) :: tmp_abs_residual_pool_perm
            !! Work array for the sorting permutation of `tmp_abs_residual_pool`
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_residuals, ierr)
        call validate_in_range_real(residual_range_quantile, ierr, min=0.0_real64, max=100.0_real64)

        if (is_err(ierr)) return

        call determine_shared_residual_range_helper(neighborhood_residuals_S1, neighborhood_residuals_S2, n_residuals, n_neighbors, shared_residual_range, tmp_abs_residual_pool, tmp_abs_residual_pool_perm, residual_range_quantile)
    end subroutine determine_shared_residual_range

    !> (no input validation) Computes the shared residual range [-R, R] for the computed residuals from studies S1 and S2
    pure subroutine determine_shared_residual_range_helper(neighborhood_residuals_S1, neighborhood_residuals_S2, n_residuals, n_neighbors, shared_residual_range, tmp_abs_residual_pool, tmp_abs_residual_pool_perm, residual_range_quantile)
        integer(int32), intent(in) :: n_residuals
            !! Number of residuals
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors (k)
        real(real64), dimension(n_residuals, n_neighbors), intent(in) :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1 (kNN), NaN is explicitly allowed for missing values
        real(real64), dimension(n_residuals, n_neighbors), intent(in) :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2 (kNN), NaN is explicitly allowed for missing values
        real(real64), intent(in), optional :: residual_range_quantile
            !! Quantile for determining the residual range, default: 95.0
        real(real64), intent(out) :: shared_residual_range
            !! Computed residual range (R)
        real(real64), dimension(2 * size(neighborhood_residuals_S1, kind=int32)), intent(out) :: tmp_abs_residual_pool
            !! Work array for the absolute residual values of the concatenated S1,S2 residuals
        integer(int32), dimension(size(tmp_abs_residual_pool, kind=int32)), intent(out) :: tmp_abs_residual_pool_perm
            !! Work array for the sorting permutation of `tmp_abs_residual_pool`

        integer(int32) :: i_residual, i_neighbor, n_pool, n_predecessors, pool_idx
        real(real64) :: actual_quantile

        M_DEFAULT_VAL(residual_range_quantile, actual_quantile, 95.0_real64)

        ! Collect the absolute residual values and compute the percentile value
        do concurrent (i_neighbor = 1:n_neighbors, i_residual = 1:n_residuals) local(n_predecessors, pool_idx) shared(tmp_abs_residual_pool_perm, neighborhood_residuals_S1, neighborhood_residuals_S2, tmp_abs_residual_pool)
            n_predecessors = (i_neighbor - 1) * n_residuals * 2
            pool_idx = n_predecessors + i_residual * 2
            tmp_abs_residual_pool(pool_idx - 1) = abs(neighborhood_residuals_S1(i_residual, i_neighbor))
            tmp_abs_residual_pool(pool_idx) = abs(neighborhood_residuals_S2(i_residual, i_neighbor))
            
            tmp_abs_residual_pool_perm(pool_idx) = pool_idx
            tmp_abs_residual_pool_perm(pool_idx - 1) = pool_idx - 1
        end do

        call sort_array_heapsort(tmp_abs_residual_pool, tmp_abs_residual_pool_perm)

        n_pool = size(tmp_abs_residual_pool, kind=int32)
        ! NaN is always last -> find last non-NaN index for percentile calculation
        do i_residual = n_pool, 1, -1
            if (ieee_is_nan(tmp_abs_residual_pool(tmp_abs_residual_pool_perm(i_residual)))) then
                n_pool = n_pool - 1
            else
                exit
            end if
        end do

        call calc_percentile_helper(tmp_abs_residual_pool(:n_pool), tmp_abs_residual_pool_perm(:n_pool), actual_quantile, shared_residual_range)
    end subroutine determine_shared_residual_range_helper

    !> Computes the shared residual range [-R, R] for the computed residuals from studies S1 and S2
    pure subroutine determine_shared_residual_range_alloc(neighborhood_residuals_S1, neighborhood_residuals_S2, n_residuals, n_neighbors, shared_residual_range, ierr, residual_range_quantile)
        integer(int32), intent(in) :: n_residuals
            !! Number of residuals
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors (k)
        real(real64), dimension(n_residuals, n_neighbors), intent(in) :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1 (kNN), NaN is explicitly allowed for missing values
        real(real64), dimension(n_residuals, n_neighbors), intent(in) :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2 (kNN), NaN is explicitly allowed for missing values
        real(real64), intent(in), optional :: residual_range_quantile
            !! Quantile for determining the residual range, default: 95.0
        real(real64), intent(out) :: shared_residual_range
            !! Computed residual range (R)
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: n_pool
        real(real64), dimension(:), allocatable :: abs_residual_pool
        integer(int32), dimension(:), allocatable :: perm

        call set_ok(ierr)

        call validate_dimension_size(n_residuals, ierr)

        if (is_err(ierr)) return

        n_pool = 2_int32 * size(neighborhood_residuals_S1, kind=int32)

        M_ALLOCATE(abs_residual_pool(n_pool))
        M_ALLOCATE(perm(n_pool))

        call determine_shared_residual_range(neighborhood_residuals_S1, neighborhood_residuals_S2, n_residuals, n_neighbors, shared_residual_range, abs_residual_pool, perm, ierr, residual_range_quantile)
    end subroutine determine_shared_residual_range_alloc

    !> Summarizes the neighborhood residuals in absolute histogram counts and probability mass functions `pmf(residual, bin)` (actually a matrix)
    pure subroutine build_residual_histograms(neighborhood_residuals, n_residuals, n_neighbors, shared_residual_range, n_bins, counts, pmf, included_n_residuals, ierr)
        integer(int32), intent(in) :: n_residuals
            !! Number of residuals
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors (k)
        real(real64), dimension(n_residuals, n_neighbors), intent(in) :: neighborhood_residuals
            !! Computed neighborhood residuals for a study (kNN), NaN is explicitly allowed for missing values
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range (R) from [[tox_jensen_shannon_divergence(module):determine_shared_residual_range_alloc(subroutine)]]
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins in range [-R,R]
        integer(int32), dimension(n_neighbors, n_bins), intent(out) :: counts
            !! Absolute counts of a residual per bin
        real(real64), dimension(n_neighbors, n_bins), intent(out) :: pmf
            !! `counts` normalized to `0 <= counts(:, i) <= 1` and `sum(counts(:, i)) == 1`
        integer(int32), dimension(n_neighbors), intent(out) :: included_n_residuals
            !! Stores the count of non-NaN residuals (included ones)
        integer(int32), intent(out) :: ierr
            !! Error code

        real(real64) :: bin_width, clamped_residual
        integer(int32) :: bin_idx, i_neighbor, i_residual, i_bin, included_residuals

        call set_ok(ierr)

        call validate_dimension_size(n_residuals, ierr)
        call validate_dimension_size(n_neighbors, ierr)
        call validate_in_range_int(n_bins, ierr, min=1_int32)
        call validate_in_range_real(shared_residual_range, ierr, min=0.0_real64)

        if (is_err(ierr)) return

        call build_residual_histograms_helper(neighborhood_residuals, n_residuals, n_neighbors, shared_residual_range, n_bins, counts, pmf, included_n_residuals)
    end subroutine build_residual_histograms

    !> (no input validation) Summarizes the neighborhood residuals in absolute histogram counts and probability mass functions `pmf(residual, bin)` (actually a matrix)
    pure subroutine build_residual_histograms_helper(neighborhood_residuals, n_residuals, n_neighbors, shared_residual_range, n_bins, counts, pmf, included_n_residuals)
        integer(int32), intent(in) :: n_residuals
            !! Number of residuals
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors (k)
        real(real64), dimension(n_residuals, n_neighbors), intent(in) :: neighborhood_residuals
            !! Computed neighborhood residuals for a study (kNN), NaN is explicitly allowed for missing values
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range (R) from [[tox_jensen_shannon_divergence(module):determine_shared_residual_range_alloc(subroutine)]]
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins in range [-R,R]
        integer(int32), dimension(n_neighbors, n_bins), intent(out) :: counts
            !! Absolute counts of a residual per bin
        real(real64), dimension(n_neighbors, n_bins), intent(out) :: pmf
            !! `counts` normalized to `0 <= counts(:, i) <= 1` and `sum(counts(:, i)) == 1`
        integer(int32), dimension(n_neighbors), intent(out) :: included_n_residuals
            !! Stores the count of non-NaN residuals (included ones)

        real(real64) :: bin_width, clamped_residual
        integer(int32) :: bin_idx, i_neighbor, i_residual, i_bin, included_residuals

        bin_width = 2.0_real64 * shared_residual_range / real(n_bins, real64)
        counts = 0_int32
        pmf = 0.0_real64

        ! 1. assign the bins to the residuals (increase the respective count)
        ! outer loop cannot be concurrent, as counts of same residuals and bins but different neighbors might be changed at the same time
        do concurrent (i_neighbor = 1:n_neighbors) local(included_residuals, i_residual, clamped_residual, bin_idx) shared(n_residuals, counts, neighborhood_residuals, shared_residual_range, bin_width)
            included_residuals = 0_int32

            do i_residual = 1, n_residuals
                if (.not. ieee_is_nan(neighborhood_residuals(i_residual, i_neighbor))) then
                    ! clamp residual to histogram range
                    clamped_residual = clamp(neighborhood_residuals(i_residual, i_neighbor), min_val=-shared_residual_range, max_val=shared_residual_range)

                    ! assign bin to residual
                    bin_idx = min(n_bins, int( (clamped_residual + shared_residual_range) / bin_width ) + 1)
                    counts(i_neighbor, bin_idx) = counts(i_neighbor, bin_idx) + 1

                    included_residuals = included_residuals + 1
                end if
            end do

            included_n_residuals(i_neighbor) = included_residuals
        end do

        ! 2. calculate pmf
        do concurrent (i_bin = 1:n_bins) local(i_neighbor) shared(counts, pmf, included_n_residuals)
            do concurrent (i_neighbor = 1:n_neighbors) shared(pmf, i_bin, included_n_residuals, counts)
                if (included_n_residuals(i_neighbor) == 0) then
                    pmf(i_neighbor, i_bin) = 0.0_real64
                else
                    pmf(i_neighbor, i_bin) = real(counts(i_neighbor, i_bin), real64) / real(included_n_residuals(i_neighbor), real64)
                end if
            end do
        end do
    end subroutine build_residual_histograms_helper

    !> Having the probabilities `pmf` from [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]], this subroutine computes the Jensen-Shannon divergence per reference point/neighbor
    pure subroutine compute_divergence_per_reference_point(pmf_S1, pmf_S2, n_neighbors, n_bins, js_divergences, ierr)
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors (k)
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins in range [-R,R]
        real(real64), dimension(n_neighbors, n_bins), intent(in) :: pmf_S1
            !! Computed normalized hostogram counts from [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]] for study 1
        real(real64), dimension(n_neighbors, n_bins), intent(in) :: pmf_S2
            !! Computed normalized hostogram counts from [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]] for study 2
        real(real64), dimension(n_neighbors), intent(out) :: js_divergences
            !! Jensen-Shannon divergence per neighbor
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_neighbors, ierr)
        call validate_dimension_size(n_bins, ierr)
        call validate_all_in_range_real(pmf_S1, size(pmf_S1, kind=int32), ierr, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_real(pmf_S2, size(pmf_S2, kind=int32), ierr, min=0.0_real64, max=1.0_real64)

        if (is_err(ierr)) return

        call compute_divergence_per_reference_point_helper(pmf_S1, pmf_S2, n_neighbors, n_bins, js_divergences)
    end subroutine compute_divergence_per_reference_point

    !> (no input validation) Having the probabilities `pmf` from [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]], this subroutine computes the Jensen-Shannon divergence per reference point/neighbor
    pure subroutine compute_divergence_per_reference_point_helper(pmf_S1, pmf_S2, n_neighbors, n_bins, js_divergences)
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors (k)
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins in range [-R,R]
        real(real64), dimension(n_neighbors, n_bins), intent(in) :: pmf_S1
            !! Computed normalized hostogram counts from [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]] for study 1
        real(real64), dimension(n_neighbors, n_bins), intent(in) :: pmf_S2
            !! Computed normalized hostogram counts from [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]] for study 2
        real(real64), dimension(n_neighbors), intent(out) :: js_divergences
            !! Jensen-Shannon divergence per neighbor

        real(real64) :: S_mean, s1_val, s2_val
        integer(int32) :: i_bin, i_neighbor

        js_divergences = 0.0_real64

        ! 1. compute the Knullback-Leibler (KL) divergences
        ! Note that the J-S divergence is defined as `0.5 * KL_S1 + 0.5 * KL_S2`, equivalent to `0.5 * (KL_S1 + KL_S2)`.
        ! Thus, instead of computing KL_S* independently, it accumulates directly in the `js_divergences` output.
        ! Another thing, switching the loops would enable both to run concurrently and the 0.5*js_divergences step could be done in one go,
        ! but cache locality still beats that, except for the case of thousands of neighbors, which might not be the common case.
        do i_bin = 1, n_bins
            do concurrent (i_neighbor = 1:n_neighbors) local(s1_val, s2_val, S_mean) shared(i_bin, pmf_S1, pmf_S2, js_divergences)
                s1_val = pmf_S1(i_neighbor, i_bin)
                s2_val = pmf_S2(i_neighbor, i_bin)
                S_mean = 0.5_real64 * (s1_val + s2_val)

                if (.not. is_close(S_mean, 0.0_real64)) then
                    if (s1_val > 0.0_real64) then
                        js_divergences(i_neighbor) = js_divergences(i_neighbor) + s1_val * log(s1_val / S_mean)
                    end if

                    if (s2_val > 0.0_real64) then
                        js_divergences(i_neighbor) = js_divergences(i_neighbor) + s2_val * log(s2_val / S_mean)
                    end if
                end if
            end do
        end do

        ! 2. Compute the js_divergences
        do concurrent (i_neighbor = 1:n_neighbors) shared(js_divergences)
            js_divergences(i_neighbor) = 0.5_real64 * js_divergences(i_neighbor)
        end do
    end subroutine compute_divergence_per_reference_point_helper

    !> Computes the global weighted Jensen-Shannon divergence from the per-neighbor divergences calculated by [[tox_jensen_shannon_divergence(module):compute_divergence_per_reference_point(subroutine)]]
    pure subroutine compute_weighted_global_divergence(js_divergences, n_neighbors, included_n_residuals_S1, included_n_residuals_S2, global_js_divergence, weights, ierr)
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors (k)
        real(real64), dimension(n_neighbors), intent(in) :: js_divergences
            !! Jensen-Shannon divergence per neighbor, computed for studies S1 and S2
        integer(int32), dimension(n_neighbors), intent(in) :: included_n_residuals_S1
            !! Count of non-NaN residuals (included ones) in study 1
        integer(int32), dimension(n_neighbors), intent(in) :: included_n_residuals_S2
            !! Count of non-NaN residuals (included ones) in study 2
        real(real64), intent(out) :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        real(real64), dimension(n_neighbors), intent(out) :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_neighbors, ierr)
        call validate_all_in_range_real(js_divergences, size(js_divergences, kind=int32), ierr, min=0.0_real64)
        call validate_all_in_range_int(included_n_residuals_S1, size(included_n_residuals_S1, kind=int32), ierr, min=0_int32)
        call validate_all_in_range_int(included_n_residuals_S2, size(included_n_residuals_S2, kind=int32), ierr, min=0_int32)

        if (is_err(ierr)) return

        call compute_weighted_global_divergence_helper(js_divergences, n_neighbors, included_n_residuals_S1, included_n_residuals_S2, global_js_divergence, weights)
    end subroutine compute_weighted_global_divergence

    !> (no input validation) Computes the global weighted Jensen-Shannon divergence from the per-neighbor divergences calculated by [[tox_jensen_shannon_divergence(module):compute_divergence_per_reference_point(subroutine)]]
    pure subroutine compute_weighted_global_divergence_helper(js_divergences, n_neighbors, included_n_residuals_S1, included_n_residuals_S2, global_js_divergence, weights)
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors (k)
        real(real64), dimension(n_neighbors), intent(in) :: js_divergences
            !! Jensen-Shannon divergence per neighbor, computed for studies S1 and S2
        integer(int32), dimension(n_neighbors), intent(in) :: included_n_residuals_S1
            !! Count of non-NaN residuals (included ones) in study 1 (obtained from [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]])
        integer(int32), dimension(n_neighbors), intent(in) :: included_n_residuals_S2
            !! Count of non-NaN residuals (included ones) in study 2 (obtained from [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]])
        real(real64), intent(out) :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        real(real64), dimension(n_neighbors), intent(out) :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`

        integer(int32) :: i_neighbor, included_residuals
        real(real64) :: total_sample_count

        global_js_divergence = 0.0_real64

        total_sample_count = real(sum(included_n_residuals_S1) + sum(included_n_residuals_S2), real64)

        ! Calculate the global Jensen-Shannon divergence
        do concurrent (i_neighbor = 1:n_neighbors) local(included_residuals) shared(weights, included_n_residuals_S1, included_n_residuals_S2, total_sample_count) reduce(+:global_js_divergence)
            included_residuals = included_n_residuals_S1(i_neighbor) + included_n_residuals_S2(i_neighbor)

            weights(i_neighbor) = real(included_residuals, real64) / total_sample_count

            global_js_divergence = global_js_divergence + weights(i_neighbor) * js_divergences(i_neighbor)
        end do
    end subroutine compute_weighted_global_divergence_helper

end module tox_jensen_shannon_divergence

!> C-compatible wrapper for [[tox_jensen_shannon_divergence(module):determine_shared_residual_range(subroutine)]]
pure subroutine determine_shared_residual_range_expert_c( &
    neighborhood_residuals_S1, neighborhood_residuals_S2, &
    n_residuals, n_neighbors, &
    residual_range_quantile, &
    shared_residual_range, &
    tmp_abs_residual_pool, tmp_abs_residual_pool_perm, &
    ierr ) &
    bind(C, name="determine_shared_residual_range_expert_c")

    use tox_jensen_shannon_divergence, only: determine_shared_residual_range
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_residuals
        !! Number of residuals
    integer(c_int), intent(in), target :: n_neighbors
        !! Number of neighbors (k)
    real(c_double), dimension(n_residuals, n_neighbors), intent(in), target :: neighborhood_residuals_S1
        !! Residuals for study 1
    real(c_double), dimension(n_residuals, n_neighbors), intent(in), target :: neighborhood_residuals_S2
        !! Residuals for study 2
    real(c_double), intent(in), target :: residual_range_quantile
        !! Quantile for determining the residual range
    real(c_double), intent(out), target :: shared_residual_range
        !! Computed residual range (R)
    real(c_double), dimension(2*n_residuals*n_neighbors), intent(out), target :: tmp_abs_residual_pool
        !! Work array for absolute residuals
    integer(c_int), dimension(2*n_residuals*n_neighbors), intent(out), target :: tmp_abs_residual_pool_perm
        !! Work array for sorting permutation
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_residuals)
    M_CHECK_NON_NULL(n_neighbors)
    M_CHECK_NON_NULL(neighborhood_residuals_S1)
    M_CHECK_NON_NULL(neighborhood_residuals_S2)
    M_CHECK_NON_NULL(residual_range_quantile)
    M_CHECK_NON_NULL(shared_residual_range)
    M_CHECK_NON_NULL(tmp_abs_residual_pool)
    M_CHECK_NON_NULL(tmp_abs_residual_pool_perm)

    call determine_shared_residual_range( &
        neighborhood_residuals_S1, neighborhood_residuals_S2, &
        n_residuals, n_neighbors, &
        shared_residual_range, &
        tmp_abs_residual_pool, tmp_abs_residual_pool_perm, &
        ierr, residual_range_quantile )

end subroutine determine_shared_residual_range_expert_c

!> C-compatible wrapper for [[tox_jensen_shannon_divergence(module):determine_shared_residual_range_alloc(subroutine)]]
pure subroutine determine_shared_residual_range_c( &
    neighborhood_residuals_S1, neighborhood_residuals_S2, &
    n_residuals, n_neighbors, &
    residual_range_quantile, &
    shared_residual_range, ierr ) &
    bind(C, name="determine_shared_residual_range_c")

    use tox_jensen_shannon_divergence, only: determine_shared_residual_range_alloc
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_residuals
        !! Number of residuals
    integer(c_int), intent(in), target :: n_neighbors
        !! Number of neighbors (k)
    real(c_double), dimension(n_residuals, n_neighbors), intent(in), target :: neighborhood_residuals_S1
        !! Residuals for study 1
    real(c_double), dimension(n_residuals, n_neighbors), intent(in), target :: neighborhood_residuals_S2
        !! Residuals for study 2
    real(c_double), intent(in), target :: residual_range_quantile
        !! Quantile for determining the residual range
    real(c_double), intent(out), target :: shared_residual_range
        !! Computed residual range (R)
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_residuals)
    M_CHECK_NON_NULL(n_neighbors)
    M_CHECK_NON_NULL(neighborhood_residuals_S1)
    M_CHECK_NON_NULL(neighborhood_residuals_S2)
    M_CHECK_NON_NULL(residual_range_quantile)
    M_CHECK_NON_NULL(shared_residual_range)

    call determine_shared_residual_range_alloc( &
        neighborhood_residuals_S1, neighborhood_residuals_S2, &
        n_residuals, n_neighbors, &
        shared_residual_range, ierr, residual_range_quantile )

end subroutine determine_shared_residual_range_c

!> C-compatible wrapper for [[tox_jensen_shannon_divergence(module):build_residual_histograms(subroutine)]]
pure subroutine build_residual_histograms_c( &
    neighborhood_residuals, &
    n_residuals, n_neighbors, &
    shared_residual_range, &
    n_bins, &
    counts, pmf, included_n_residuals, &
    ierr ) &
    bind(C, name="build_residual_histograms_c")

    use tox_jensen_shannon_divergence, only: build_residual_histograms
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_residuals
        !! Number of residuals
    integer(c_int), intent(in), target :: n_neighbors
        !! Number of neighbors (k)
    real(c_double), dimension(n_residuals, n_neighbors), intent(in), target :: neighborhood_residuals
        !! Residuals for a study (kNN), NaN allowed
    real(c_double), intent(in), target :: shared_residual_range
        !! Shared residual range R
    integer(c_int), intent(in), target :: n_bins
        !! Number of histogram bins
    integer(c_int), dimension(n_neighbors, n_bins), intent(out), target :: counts
        !! Histogram counts
    real(c_double), dimension(n_neighbors, n_bins), intent(out), target :: pmf
        !! Normalized PMF
    integer(c_int), dimension(n_neighbors), intent(out), target :: included_n_residuals
        !! Count of non-NaN residuals per neighbor
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(neighborhood_residuals)
    M_CHECK_NON_NULL(n_residuals)
    M_CHECK_NON_NULL(n_neighbors)
    M_CHECK_NON_NULL(shared_residual_range)
    M_CHECK_NON_NULL(n_bins)
    M_CHECK_NON_NULL(counts)
    M_CHECK_NON_NULL(pmf)
    M_CHECK_NON_NULL(included_n_residuals)

    call build_residual_histograms( &
        neighborhood_residuals, &
        n_residuals, n_neighbors, &
        shared_residual_range, &
        n_bins, &
        counts, pmf, included_n_residuals, &
        ierr )

end subroutine build_residual_histograms_c

!> C-compatible wrapper for [[tox_jensen_shannon_divergence(module):compute_divergence_per_reference_point(subroutine)]]
pure subroutine compute_divergence_per_reference_point_c( &
    pmf_S1, pmf_S2, &
    n_neighbors, n_bins, &
    js_divergences, ierr ) &
    bind(C, name="compute_divergence_per_reference_point_c")

    use tox_jensen_shannon_divergence, only: compute_divergence_per_reference_point
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_neighbors
        !! Number of neighbors (k)
    integer(c_int), intent(in), target :: n_bins
        !! Number of histogram bins
    real(c_double), dimension(n_neighbors, n_bins), intent(in), target :: pmf_S1
        !! PMF for study S1
    real(c_double), dimension(n_neighbors, n_bins), intent(in), target :: pmf_S2
        !! PMF for study S2
    real(c_double), dimension(n_neighbors), intent(out), target :: js_divergences
        !! Jensen–Shannon divergence per neighbor
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_neighbors)
    M_CHECK_NON_NULL(n_bins)
    M_CHECK_NON_NULL(pmf_S1)
    M_CHECK_NON_NULL(pmf_S2)
    M_CHECK_NON_NULL(js_divergences)

    call compute_divergence_per_reference_point( &
        pmf_S1, pmf_S2, &
        n_neighbors, n_bins, &
        js_divergences, ierr )

end subroutine compute_divergence_per_reference_point_c

!> C-compatible wrapper for [[tox_jensen_shannon_divergence(module):compute_weighted_global_divergence(subroutine)]]
pure subroutine compute_weighted_global_divergence_c( &
    js_divergences, &
    n_neighbors, &
    included_n_residuals_S1, included_n_residuals_S2, &
    global_js_divergence, weights, &
    ierr ) &
    bind(C, name="compute_weighted_global_divergence_c")

    use tox_jensen_shannon_divergence, only: compute_weighted_global_divergence
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_neighbors
        !! Number of neighbors (k)
    real(c_double), dimension(n_neighbors), intent(in), target :: js_divergences
        !! Per-neighbor JSD values
    integer(c_int), dimension(n_neighbors), intent(in), target :: included_n_residuals_S1
        !! Residual counts for study S1
    integer(c_int), dimension(n_neighbors), intent(in), target :: included_n_residuals_S2
        !! Residual counts for study S2
    real(c_double), intent(out), target :: global_js_divergence
        !! Weighted global JSD
    real(c_double), dimension(n_neighbors), intent(out), target :: weights
        !! Normalized weights
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_neighbors)
    M_CHECK_NON_NULL(js_divergences)
    M_CHECK_NON_NULL(included_n_residuals_S1)
    M_CHECK_NON_NULL(included_n_residuals_S2)
    M_CHECK_NON_NULL(global_js_divergence)
    M_CHECK_NON_NULL(weights)

    call compute_weighted_global_divergence( &
        js_divergences, n_neighbors, &
        included_n_residuals_S1, included_n_residuals_S2, &
        global_js_divergence, weights, ierr )

end subroutine compute_weighted_global_divergence_c
