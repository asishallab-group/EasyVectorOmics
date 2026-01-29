#include "macros.h"

module tox_jensen_shannon_divergence
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use f42_utils, only: clamp, calc_percentile, is_close, sort_array_heapsort
    use tox_errors, only: set_ok, set_err, is_err, ERR_ALLOC_FAIL, validate_dimension_size, validate_in_range_real, validate_all_in_range_real, validate_in_range_int
    implicit none
contains

    !> Computes the shared residual range [-R, R] for the computed residuals from studies S1 and S2
    pure subroutine determine_shared_residual_range_alloc(neighborhood_residuals_S1, neighborhood_residuals_S2, n_residuals, n_neighbors, shared_residual_range, residual_range_quantile, ierr)
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

        integer(int32) :: i_residual, i_neighbor, n_pool, n_predecessors, pool_idx
        real(real64) :: actual_quantile
        real(real64), dimension(:), allocatable :: abs_residual_pool
        integer(int32), dimension(:), allocatable :: perm

        M_DEFAULT_VAL(residual_range_quantile, actual_quantile, 95.0_real64)

        call set_ok(ierr)

        call validate_dimension_size(n_residuals, ierr)
        call validate_in_range_real(actual_quantile, ierr, min=0.0_real64, max=100.0_real64)

        if (is_err(ierr)) return

        n_pool = 2_int32 * size(neighborhood_residuals_S1, kind=int32)

        M_ALLOCATE(abs_residual_pool(n_pool))
        M_ALLOCATE(perm(n_pool))

        ! Collect the absolute residual values and compute the percentile value
        do concurrent (i_neighbor = 1:n_neighbors, i_residual = 1:n_pool) local(n_predecessors, pool_idx) shared(perm, neighborhood_residuals_S1, neighborhood_residuals_S2, abs_residual_pool)
            n_predecessors = (i_neighbor - 1) * n_residuals * 2
            pool_idx = n_predecessors + i_residual * 2
            abs_residual_pool(pool_idx - 1) = abs(neighborhood_residuals_S1(i_residual, i_neighbor))
            abs_residual_pool(pool_idx) = abs(neighborhood_residuals_S2(i_residual, i_neighbor))
            
            perm(pool_idx) = pool_idx
            perm(pool_idx - 1) = pool_idx - 1
        end do

        call sort_array_heapsort(abs_residual_pool, perm)

        ! NaN is always last -> find last non-NaN index for percentile calculation
        do i_residual = n_pool, 1, -1
            if (ieee_is_nan(abs_residual_pool(perm(i_residual)))) then
                n_pool = n_pool - 1
            else
                exit
            end if
        end do

        call calc_percentile(abs_residual_pool(:n_pool), perm(:n_pool), actual_quantile, shared_residual_range, ierr)
    end subroutine determine_shared_residual_range_alloc

    !> Summarizes the neighborhood residuals in absolute histogram counts and probability mass functions `pmf(residual, bin)` (actually a matrix)
    pure subroutine build_residual_histograms(neighborhood_residuals, n_residuals, n_neighbors, shared_residual_range, n_bins, counts, pmf, ierr)
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
        integer(int32), intent(out) :: ierr
            !! Error code

        real(real64) :: bin_width, clamped_residual, bin_residuals
        integer(int32) :: bin_idx, i_neighbor, i_residual, i_bin

        call set_ok(ierr)

        call validate_dimension_size(n_residuals, ierr)
        call validate_dimension_size(n_neighbors, ierr)
        call validate_in_range_int(n_bins, ierr, min=1_int32)
        call validate_in_range_real(shared_residual_range, ierr, min=0.0_real64)

        if (is_err(ierr)) return

        bin_width = 2.0_real64 * shared_residual_range / real(n_bins, real64)
        counts = 0_int32

        ! 1. assign the bins to the residuals (increase the respective count)
        ! outer loop cannot be concurrent, as counts of same residuals and bins but different neighbors might be changed at the same time
        do concurrent (i_neighbor = 1:n_neighbors) local(i_residual, clamped_residual, bin_idx) shared(n_residuals, counts, neighborhood_residuals, shared_residual_range, bin_width)
            do i_residual = 1, n_residuals
                if (.not. ieee_is_nan(neighborhood_residuals(i_residual, i_neighbor))) then
                    ! clamp residual to histogram range
                    clamped_residual = clamp(neighborhood_residuals(i_residual, i_neighbor), min_val=-shared_residual_range, max_val=shared_residual_range)

                    ! assign bin to residual
                    bin_idx = int( (clamped_residual + shared_residual_range) / bin_width )
                    counts(i_neighbor, bin_idx) = counts(i_neighbor, bin_idx) + 1
                end if
            end do
        end do

        ! 2. calculate pmf
        do concurrent (i_bin = 1:n_bins) local(bin_residuals, i_residual) shared(counts, pmf)
            bin_residuals = real(sum(counts(:, i_bin)), real64)
            if (is_close(bin_residuals, 0.0_real64)) then
                pmf(:, i_bin) = 0.0_real64
            else
                do concurrent (i_neighbor = 1:n_neighbors) shared(pmf, bin_residuals, i_bin)
                    pmf(i_neighbor, i_bin) = real(counts(i_neighbor, i_bin), real64) / bin_residuals
                end do
            end if
        end do
    end subroutine build_residual_histograms

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
            !! Jensen-Shannon divergence per residual
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
            !! Jensen-Shannon divergence per residual

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

                if (s1_val > 0.0_real64) then
                    js_divergences(i_neighbor) = js_divergences(i_neighbor) + s1_val * log(s1_val / S_mean)
                end if

                if (s2_val > 0.0_real64) then
                    js_divergences(i_neighbor) = js_divergences(i_neighbor) + s2_val * log(s2_val / S_mean)
                end if
            end do
        end do

        ! 2. Compute the js_divergences
        do concurrent (i_neighbor = 1:n_neighbors) shared(js_divergences)
            js_divergences(i_neighbor) = 0.5_real64 * js_divergences(i_neighbor)
        end do
    end subroutine compute_divergence_per_reference_point_helper

    !> Computes the global weighted Jensen-Shannon divergence from the per-neighbor divergences calculated by [[tox_jensen_shannon_divergence(module):compute_divergence_per_reference_point(subroutine)]]
    pure subroutine compute_weighted_global_divergence(js_divergences, neighborhood_residuals_S1, neighborhood_residuals_S2, n_residuals, n_neighbors, global_js_divergence, included_n_residuals, weights, ierr)
        integer(int32), intent(in) :: n_residuals
            !! Number of residuals
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors (k)
        real(real64), dimension(n_neighbors), intent(in) :: js_divergences
            !! Jensen-Shannon divergence per residual, computed for `neighborhood_residuals_S1` and `neighborhood_residuals_S2`
        real(real64), dimension(n_residuals, n_neighbors), intent(in) :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1 (kNN), NaN is explicitly allowed for missing values
        real(real64), dimension(n_residuals, n_neighbors), intent(in) :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2 (kNN), NaN is explicitly allowed for missing values
        real(real64), intent(out) :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        integer(int32), dimension(n_neighbors), intent(out) :: included_n_residuals
            !! Stores the count of non-NaN residuals (included ones)
        real(real64), dimension(n_neighbors), intent(out) :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i_neighbor, i_residual
        real(real64) :: total_sample_count

        call set_ok(ierr)

        call validate_dimension_size(n_neighbors, ierr)
        call validate_dimension_size(n_residuals, ierr)
        call validate_all_in_range_real(js_divergences, size(js_divergences, kind=int32), ierr, min=0.0_real64)

        if (is_err(ierr)) return

        call compute_weighted_global_divergence_helper(js_divergences, neighborhood_residuals_S1, neighborhood_residuals_S2, n_residuals, n_neighbors, global_js_divergence, included_n_residuals, weights)
    end subroutine compute_weighted_global_divergence

    !> (no input validation) Computes the global weighted Jensen-Shannon divergence from the per-neighbor divergences calculated by [[tox_jensen_shannon_divergence(module):compute_divergence_per_reference_point(subroutine)]]
    pure subroutine compute_weighted_global_divergence_helper(js_divergences, neighborhood_residuals_S1, neighborhood_residuals_S2, n_residuals, n_neighbors, global_js_divergence, included_n_residuals, weights)
        integer(int32), intent(in) :: n_residuals
            !! Number of residuals
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors (k)
        real(real64), dimension(n_neighbors), intent(in) :: js_divergences
            !! Jensen-Shannon divergence per residual, computed for `neighborhood_residuals_S1` and `neighborhood_residuals_S2`
        real(real64), dimension(n_residuals, n_neighbors), intent(in) :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1 (kNN), NaN is explicitly allowed for missing values
        real(real64), dimension(n_residuals, n_neighbors), intent(in) :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2 (kNN), NaN is explicitly allowed for missing values
        real(real64), intent(out) :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        integer(int32), dimension(n_neighbors), intent(out) :: included_n_residuals
            !! Stores the count of non-NaN residuals (included ones)
        real(real64), dimension(n_neighbors), intent(out) :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`

        integer(int32) :: i_neighbor, i_residual
        real(real64) :: total_sample_count

        global_js_divergence = 0.0_real64
        included_n_residuals = 0_int32

        ! Count included residuals
        do concurrent (i_neighbor = 1:n_neighbors) local(i_residual) shared(n_residuals, neighborhood_residuals_S1, neighborhood_residuals_S2, included_n_residuals)
            do i_residual = 1, n_residuals
                if (ieee_is_nan(neighborhood_residuals_S1(i_residual, i_neighbor)) .or. ieee_is_nan(neighborhood_residuals_S2(i_residual, i_neighbor))) then
                    included_n_residuals(i_neighbor) = included_n_residuals(i_neighbor) + 1
                end if
            end do
        end do

        total_sample_count = real(sum(included_n_residuals), real64)

        ! Calculate the global Jensen-Shannon divergence
        do concurrent (i_neighbor = 1:n_neighbors) shared(weights, included_n_residuals, total_sample_count) reduce(+:global_js_divergence)
            weights(i_neighbor) = real(included_n_residuals(i_neighbor), real64) / total_sample_count

            global_js_divergence = global_js_divergence + weights(i_neighbor)
        end do
    end subroutine compute_weighted_global_divergence_helper

end module tox_jensen_shannon_divergence