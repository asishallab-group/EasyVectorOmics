#include "macros.h"

module tox_trajectory_contribution_analysis
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use f42_utils, only: is_close
    use tox_errors, only: set_ok, is_err, set_err, ERR_IDX_OUT_OF_BOUNDS, ERR_DIVISION_BY_ZERO, ERR_INVALID_INPUT, ERR_NAN_INF, ERR_ALLOC_FAIL, validate_dimension_size, validate_all_in_range_real, validate_all_in_range_int, validate_in_range_real
    implicit none

    integer(int32), parameter :: MODE_NORMAL = 1
    integer(int32), parameter :: MODE_RAP    = 2

    ! Baseline computation modes
    integer(int32), parameter :: BASELINE_RAW  = 1
    integer(int32), parameter :: BASELINE_MIN  = 2
    integer(int32), parameter :: BASELINE_MEAN = 3
contains

    !> Calculates trajectory contribution using cosine similarity.
    pure subroutine trajectory_contribution(factor, dependent, n_timepoints, mode, contribution, ierr)
        integer(int32), intent(in) :: n_timepoints
        real(real64), dimension(n_timepoints), intent(in) :: factor
        real(real64), dimension(n_timepoints), intent(in) :: dependent
        integer(int32), intent(in) :: mode
        real(real64), intent(out) :: contribution
        integer(int32), intent(out) :: ierr

        integer(int32) :: i_timepoint
        real(real64) :: norm_factor, norm_dependent, magnitude, dot_prod

        call set_ok(ierr)

        dot_prod = 0.0_real64
        norm_factor = 0.0_real64
        norm_dependent = 0.0_real64
        do i_timepoint = 1, n_timepoints
            norm_factor = norm_factor + factor(i_timepoint) ** 2
            norm_dependent = norm_dependent + dependent(i_timepoint) ** 2
            dot_prod = dot_prod + factor(i_timepoint) * dependent(i_timepoint)
        end do

        if (ieee_is_nan(dot_prod)) then
            call set_err(ierr, ERR_NAN_INF)
            return
        end if

        magnitude = sqrt(norm_factor) * sqrt(norm_dependent)
        if (is_close(magnitude, 0.0_real64)) then
            call set_err(ierr, ERR_DIVISION_BY_ZERO)
            return
        end if

        select case (mode)
        case (MODE_NORMAL)
            contribution = dot_prod / magnitude
        case (MODE_RAP)
            contribution = acos(dot_prod / magnitude)
        case default
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end select
    end subroutine trajectory_contribution

    !> Calculates spike contributions using cosine similarity for specific indices.
    pure subroutine spike_contribution(factor, dependent, n_timepoints, mode, contribution, ierr)
        integer(int32), intent(in) :: n_timepoints
        real(real64), dimension(n_timepoints), intent(in) :: factor
        real(real64), dimension(n_timepoints), intent(in) :: dependent
        integer(int32), intent(in) :: mode
        real(real64), dimension(n_timepoints), intent(out) :: contribution
        integer(int32), intent(out) :: ierr

        integer(int32) :: i_timepoint
        real(real64) :: norm_factor, norm_dependent, magnitude

        call set_ok(ierr)

        norm_factor = 0.0_real64
        norm_dependent = 0.0_real64
        do i_timepoint = 1, n_timepoints
            norm_factor = norm_factor + factor(i_timepoint) ** 2
            norm_dependent = norm_dependent + dependent(i_timepoint) ** 2
        end do

        magnitude = sqrt(norm_factor) * sqrt(norm_dependent)
        if (ieee_is_nan(magnitude)) then
            call set_err(ierr, ERR_NAN_INF)
            return
        end if

        if (is_close(magnitude, 0.0_real64)) then
            call set_err(ierr, ERR_DIVISION_BY_ZERO)
            return
        end if

        select case (mode)
        case (MODE_NORMAL)
            do i_timepoint = 1, n_timepoints
                contribution(i_timepoint) = (factor(i_timepoint) * dependent(i_timepoint)) / magnitude
            end do
        case (MODE_RAP)
            do i_timepoint = 1, n_timepoints
                contribution(i_timepoint) = acos((factor(i_timepoint) * dependent(i_timepoint)) / magnitude)
            end do
        case default
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end select
    end subroutine spike_contribution

    !> Calculates both spike and integrated contributions for a given factor against a dependent variable.
    pure subroutine calc_contributions(trajectories, n_factors, n_samples, n_timepoints, i_factor, dependent_idx, mode, spike_contribs, integrated_contribs, temp_factor_vector, temp_dependent_vector, ierr)
        integer(int32), intent(in) :: n_factors
        integer(int32), intent(in) :: n_samples
        integer(int32), intent(in) :: n_timepoints
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
        integer(int32), intent(in) :: i_factor
        integer(int32), intent(in) :: dependent_idx
        integer(int32), intent(in) :: mode
        real(real64), dimension(n_samples), intent(out) :: integrated_contribs
        real(real64), dimension(n_timepoints, n_samples), intent(out) :: spike_contribs
        real(real64), dimension(n_timepoints), intent(out) :: temp_factor_vector
        real(real64), dimension(n_timepoints), intent(out) :: temp_dependent_vector
        integer(int32), intent(out) :: ierr

        integer(int32) :: i_sample

        call set_ok(ierr)

        do i_sample = 1, n_samples
            call get_vec_across_timepoints(trajectories, n_factors, n_samples, n_timepoints, dependent_idx, i_sample, temp_dependent_vector, ierr)
            if (is_err(ierr)) exit
            call get_vec_across_timepoints(trajectories, n_factors, n_samples, n_timepoints, i_factor, i_sample, temp_factor_vector, ierr)
            if (is_err(ierr)) exit

            call trajectory_contribution(temp_factor_vector, temp_dependent_vector, n_timepoints, mode, integrated_contribs(i_sample), ierr)
            if (is_err(ierr)) exit
            call spike_contribution(temp_factor_vector, temp_dependent_vector, n_timepoints, mode, spike_contribs(:, i_sample), ierr)
            if (is_err(ierr)) exit
        end do
    end subroutine calc_contributions

    !> Allocation-friendly wrapper around calc_contributions.
    subroutine calc_contributions_alloc(trajectories, n_factors, n_samples, n_timepoints, i_factor, dependent_idx, mode, spike_contribs, integrated_contribs, ierr)
        integer(int32), intent(in) :: n_factors
        integer(int32), intent(in) :: n_samples
        integer(int32), intent(in) :: n_timepoints
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
        integer(int32), intent(in) :: i_factor
        integer(int32), intent(in) :: dependent_idx
        integer(int32), intent(in) :: mode
        real(real64), dimension(n_timepoints, n_samples), intent(out) :: spike_contribs
        real(real64), dimension(n_samples), intent(out) :: integrated_contribs
        integer(int32), intent(out) :: ierr

        real(real64), allocatable :: temp_factor_vector(:)
        real(real64), allocatable :: temp_dependent_vector(:)
        integer :: stat_alloc

        call set_ok(ierr)

        call validate_dimension_size(n_timepoints, ierr)
        if (is_err(ierr)) return

        allocate(temp_factor_vector(n_timepoints), temp_dependent_vector(n_timepoints), stat=stat_alloc)
        if (stat_alloc /= 0) then
            call set_err(ierr, ERR_ALLOC_FAIL)
            return
        end if

        call calc_contributions(trajectories, n_factors, n_samples, n_timepoints, i_factor, dependent_idx, mode, spike_contribs, integrated_contribs, temp_factor_vector, temp_dependent_vector, ierr)

        if (allocated(temp_dependent_vector)) deallocate(temp_dependent_vector)
        if (allocated(temp_factor_vector))    deallocate(temp_factor_vector)
    end subroutine calc_contributions_alloc

    !> Extract a sample vector across all samples for a given factor/timepoint.
    pure subroutine get_vec_across_samples(trajectories, n_factors, n_samples, n_timepoints, i_factor, i_timepoint, result_vector, ierr)
        integer(int32), intent(in) :: n_factors
        integer(int32), intent(in) :: n_samples
        integer(int32), intent(in) :: n_timepoints
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
        integer(int32), intent(in) :: i_factor
        integer(int32), intent(in) :: i_timepoint
        real(real64), dimension(n_samples), intent(out) :: result_vector
        integer(int32), intent(out) :: ierr

        integer(int32) :: i_sample

        call set_ok(ierr)

        call validate_dimension_size(n_samples, ierr)
        call validate_dimension_size(n_factors, ierr)
        call validate_dimension_size(n_timepoints, ierr)
        if (is_err(ierr)) return

        if (i_factor < 1 .or. i_factor > n_factors .or. i_timepoint < 1 .or. i_timepoint > n_timepoints) then
            call set_err(ierr, ERR_IDX_OUT_OF_BOUNDS)
            return
        end if

        do i_sample = 1, n_samples
            result_vector(i_sample) = trajectories(i_factor, i_sample, i_timepoint)
        end do
    end subroutine get_vec_across_samples

    !> Extract a trajectory vector across all timepoints for a given factor/sample.
    pure subroutine get_vec_across_timepoints(trajectories, n_factors, n_samples, n_timepoints, i_factor, i_sample, result_vector, ierr)
        integer(int32), intent(in) :: n_factors
        integer(int32), intent(in) :: n_samples
        integer(int32), intent(in) :: n_timepoints
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
        integer(int32), intent(in) :: i_factor
        integer(int32), intent(in) :: i_sample
        real(real64), dimension(n_timepoints), intent(out) :: result_vector
        integer(int32), intent(out) :: ierr

        integer(int32) :: i_timepoint

        call set_ok(ierr)

        call validate_dimension_size(n_samples, ierr)
        call validate_dimension_size(n_factors, ierr)
        call validate_dimension_size(n_timepoints, ierr)
        if (is_err(ierr)) return

        if (i_factor < 1 .or. i_factor > n_factors .or. i_sample < 1 .or. i_sample > n_samples) then
            call set_err(ierr, ERR_IDX_OUT_OF_BOUNDS)
            return
        end if

        do i_timepoint = 1, n_timepoints
            result_vector(i_timepoint) = trajectories(i_factor, i_sample, i_timepoint)
        end do
    end subroutine get_vec_across_timepoints

    !> Process trajectories computing thresholds and outliers using per-timepoint percentiles.
    subroutine process_trajectories_alloc(trajectories, n_factors, n_samples, n_timepoints, factor_mask, factor_mask_count, dependent_idx, mode, percentile, &
                                integrated_contribs, spike_contribs, thresholds_integrated_contrib, &
                                outliers_integrated_contrib, thresholds_spike_contrib, &
                                outliers_spike_contrib, ierr)
        use tox_trajectory_contribution_analysis_outlier_detection
        integer(int32), intent(in) :: n_factors
        integer(int32), intent(in) :: n_samples
        integer(int32), intent(in) :: n_timepoints
        integer(int32), intent(in) :: factor_mask_count
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
        logical, dimension(n_factors), intent(in) :: factor_mask
        integer(int32), intent(in) :: dependent_idx
        integer(int32), intent(in) :: mode
        real(real64), intent(in) :: percentile
        real(real64), dimension(n_samples, factor_mask_count), intent(out) :: integrated_contribs
        real(real64), dimension(n_timepoints, n_samples, factor_mask_count), intent(out) :: spike_contribs
        real(real64), dimension(factor_mask_count), intent(out) :: thresholds_integrated_contrib
        real(real64), dimension(n_timepoints, factor_mask_count), intent(out) :: thresholds_spike_contrib
        logical, dimension(n_samples, factor_mask_count), intent(out) :: outliers_integrated_contrib
        logical, dimension(n_timepoints, n_samples, factor_mask_count), intent(out) :: outliers_spike_contrib
        integer(int32), intent(out) :: ierr

        integer(int32) :: i_factor, i_included_factor

        call set_ok(ierr)

        call validate_dimension_size(n_samples, ierr)
        call validate_dimension_size(n_factors, ierr)
        call validate_dimension_size(n_timepoints, ierr)
        if (is_err(ierr)) return

        if (dependent_idx < 1 .or. dependent_idx > n_factors) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if

        if (factor_mask(dependent_idx)) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if

        i_included_factor = 0
        do i_factor = 1, n_factors
            if (factor_mask(i_factor)) then
                i_included_factor = i_included_factor + 1

                call calc_contributions_alloc(trajectories, n_factors, n_samples, n_timepoints, &
                    i_factor, dependent_idx, mode, spike_contribs(:, :, i_included_factor), integrated_contribs(:, i_included_factor), ierr)
                if (is_err(ierr)) exit

                call calc_integrated_threshold_alloc(integrated_contribs(:, i_included_factor), n_samples, percentile, thresholds_integrated_contrib(i_included_factor), ierr)
                if (is_err(ierr)) exit

                call detect_outliers_integrated(integrated_contribs(:, i_included_factor), n_samples, thresholds_integrated_contrib(i_included_factor), &
                    outliers_integrated_contrib(:, i_included_factor), ierr)
                if (is_err(ierr)) exit

                call calc_spike_thresholds_alloc(spike_contribs(:, :, i_included_factor), n_timepoints, n_samples, percentile, thresholds_spike_contrib(:, i_included_factor), ierr)
                if (is_err(ierr)) exit

                call detect_outliers_spike(spike_contribs(:, :, i_included_factor), n_timepoints, n_samples, thresholds_spike_contrib(:, i_included_factor), &
                    outliers_spike_contrib(:, :, i_included_factor), ierr)
                if (is_err(ierr)) exit
            end if
        end do
    end subroutine process_trajectories_alloc

    !> Process trajectories using a global percentile threshold for spike contributions.
    subroutine process_trajectories_flat_alloc(trajectories, n_factors, n_samples, n_timepoints, factor_mask, factor_mask_count, dependent_idx, mode, percentile, &
                                    integrated_contribs, spike_contribs, thresholds_integrated_contrib, &
                                    outliers_integrated_contrib, thresholds_spike_contrib, &
                                    outliers_spike_contrib, ierr)
        use tox_trajectory_contribution_analysis_outlier_detection
        integer(int32), intent(in) :: n_factors
        integer(int32), intent(in) :: n_samples
        integer(int32), intent(in) :: n_timepoints
        integer(int32), intent(in) :: factor_mask_count
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
        logical, dimension(n_factors), intent(in) :: factor_mask
        integer(int32), intent(in) :: dependent_idx
        integer(int32), intent(in) :: mode
        real(real64), intent(in) :: percentile
        real(real64), dimension(n_samples, factor_mask_count), intent(out) :: integrated_contribs
        real(real64), dimension(n_timepoints, n_samples, factor_mask_count), intent(out) :: spike_contribs
        real(real64), dimension(factor_mask_count), intent(out) :: thresholds_integrated_contrib
        real(real64), dimension(factor_mask_count), intent(out) :: thresholds_spike_contrib
        logical, dimension(n_samples, factor_mask_count), intent(out) :: outliers_integrated_contrib
        logical, dimension(n_timepoints, n_samples, factor_mask_count), intent(out) :: outliers_spike_contrib
        integer(int32), intent(out) :: ierr

        integer(int32) :: i_factor, i_included_factor
        real(real64), allocatable :: flattened_spikes(:)
        logical, allocatable :: flattened_mask(:)
        integer :: stat_alloc

        call set_ok(ierr)

        call validate_dimension_size(n_samples, ierr)
        call validate_dimension_size(n_factors, ierr)
        call validate_dimension_size(n_timepoints, ierr)
        if (is_err(ierr)) return

        if (dependent_idx < 1 .or. dependent_idx > n_factors) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if

        if (factor_mask(dependent_idx)) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if

        allocate(flattened_spikes(n_timepoints * n_samples), stat=stat_alloc)
        if (stat_alloc /= 0) then
            call set_err(ierr, ERR_ALLOC_FAIL)
            return
        end if

        allocate(flattened_mask(n_timepoints * n_samples), stat=stat_alloc)
        if (stat_alloc /= 0) then
            call set_err(ierr, ERR_ALLOC_FAIL)
            deallocate(flattened_spikes)
            return
        end if

        i_included_factor = 0
        do i_factor = 1, n_factors
            if (factor_mask(i_factor)) then
                i_included_factor = i_included_factor + 1

                call calc_contributions_alloc(trajectories, n_factors, n_samples, n_timepoints, &
                    i_factor, dependent_idx, mode, spike_contribs(:, :, i_included_factor), integrated_contribs(:, i_included_factor), ierr)
                if (is_err(ierr)) exit

                call calc_integrated_threshold_alloc(integrated_contribs(:, i_included_factor), n_samples, percentile, thresholds_integrated_contrib(i_included_factor), ierr)
                if (is_err(ierr)) exit

                call detect_outliers_integrated(integrated_contribs(:, i_included_factor), n_samples, thresholds_integrated_contrib(i_included_factor), &
                    outliers_integrated_contrib(:, i_included_factor), ierr)
                if (is_err(ierr)) exit

                flattened_spikes = reshape(spike_contribs(:, :, i_included_factor), [n_timepoints * n_samples])

                call calc_integrated_threshold_alloc(flattened_spikes, n_timepoints * n_samples, percentile, thresholds_spike_contrib(i_included_factor), ierr)
                if (is_err(ierr)) exit

                call detect_outliers_integrated(flattened_spikes, n_timepoints * n_samples, thresholds_spike_contrib(i_included_factor), flattened_mask, ierr)
                if (is_err(ierr)) exit

                outliers_spike_contrib(:, :, i_included_factor) = reshape(flattened_mask, [n_timepoints, n_samples])
            end if
        end do

        if (allocated(flattened_mask)) deallocate(flattened_mask)
        if (allocated(flattened_spikes)) deallocate(flattened_spikes)
    end subroutine process_trajectories_flat_alloc

    !> This routine performs contribution analysis for a specific factor–dependent pair, no input validation
    pure subroutine compute_contributions_helper(factor, dependent, n_dims, mode, local_contributions, total_contribution, ierr)
        integer(int32), intent(in) :: n_dims
            !! Number of elements in `factor` and `dependent`
        real(real64), dimension(n_dims), intent(in) :: factor
            !! Factor time series, length n_timepoints
        real(real64), dimension(n_dims), intent(in) :: dependent
            !! Dependent variable time series, length n_timepoints
        integer(int32), intent(in) :: mode
            !! Baseline mode: 1=RAW, 2=MIN, 3=MEAN
        real(real64), dimension(n_dims), intent(out) :: local_contributions
            !! Per-element contributions
        real(real64), intent(out) :: total_contribution
            !! Total contribution (`sum(local_contributions)`)
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i_dim
        real(real64) :: factor_baseline, dependent_baseline

        call set_ok(ierr)

        call compute_baselines_factor_dependent(n_dims, factor, dependent, mode, factor_baseline, dependent_baseline, ierr)
        if (is_err(ierr)) return

        total_contribution = 0.0_real64
        do i_dim = 1, n_dims
            local_contributions(i_dim) = (factor(i_dim) - factor_baseline) * (dependent(i_dim) - dependent_baseline)
            total_contribution = total_contribution + local_contributions(i_dim)
        end do
    end subroutine compute_contributions_helper

    !> This routine performs contribution analysis for a specific factor–dependent pair, including input validation
    pure subroutine compute_contributions(factor, dependent, n_dims, mode, local_contributions, total_contribution, ierr)
        integer(int32), intent(in) :: n_dims
            !! Number of elements in `factor` and `dependent`
        real(real64), dimension(n_dims), intent(in) :: factor
            !! Factor time series, length n_timepoints
        real(real64), dimension(n_dims), intent(in) :: dependent
            !! Dependent variable time series, length n_timepoints
        integer(int32), intent(in) :: mode
            !! Baseline mode: 1=RAW, 2=MIN, 3=MEAN
        real(real64), dimension(n_dims), intent(out) :: local_contributions
            !! Per-element contributions
        real(real64), intent(out) :: total_contribution
            !! Total contribution (`sum(local_contributions)`)
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_dims, ierr)
        call validate_all_in_range_real(factor, n_dims, ierr)
        call validate_all_in_range_real(dependent, n_dims, ierr)

        if (is_err(ierr)) return

        call compute_contributions_helper(factor, dependent, n_dims, mode, local_contributions, total_contribution, ierr)
    end subroutine compute_contributions

    !> This routine performs contribution analysis for every selected factor–dependent pair
    pure subroutine compute_all_contributions(trajectories, n_factors, n_samples, n_timepoints, factor_indices, n_selected_factors, dependent_indices, n_selected_dependents, mode, local_contributions, total_contributions, temp_factors, temp_dependent, ierr)
        integer(int32), intent(in) :: n_factors
            !! number of factors
        integer(int32), intent(in) :: n_samples
            !! number of samples
        integer(int32), intent(in) :: n_timepoints
            !! number of timepoints
        integer(int32), intent(in) :: n_selected_factors
            !! number of selected factors in `factor_indices`
        integer(int32), intent(in) :: n_selected_dependents
            !! number of selected dependents in `dependent_indices`
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
        integer(int32), dimension(n_selected_factors), intent(in) :: factor_indices
            !! indices of factors to compute the contributions for
        integer(int32), dimension(n_selected_dependents), intent(in) :: dependent_indices
            !! indices of dependents to compute the contributions for
        integer(int32), intent(in) :: mode
            !! Baseline mode: 1=RAW, 2=MIN, 3=MEAN
        real(real64), dimension(n_timepoints, n_selected_factors, n_selected_dependents, n_samples), intent(out) :: local_contributions
            !! Per-timepoint contributions per sample-dependent-factor combination
        real(real64), dimension(n_selected_factors, n_selected_dependents, n_samples), intent(out) :: total_contributions
            !! Total contribution (`sum(local_contributions)`) per sample-dependent-factor combination
        real(real64), dimension(n_timepoints, n_selected_factors), intent(out) :: temp_factors
            !! Working array to hold the currently handled sample's factors in contiguous memory
        real(real64), dimension(n_timepoints), intent(out) :: temp_dependent
            !! Working array to hold the currently handled dependent in contiguous memory
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i_timepoint, i_dependent, i_factor, i_sel_factor, i_sel_dependent, i_sample

        call set_ok(ierr)

        call validate_dimension_size(n_factors, ierr)
        call validate_dimension_size(n_samples, ierr)
        call validate_dimension_size(n_timepoints, ierr)
        call validate_dimension_size(n_selected_factors, ierr)
        call validate_dimension_size(n_selected_dependents, ierr)
        call validate_all_in_range_real(trajectories, size(trajectories, kind=int32), ierr)
        call validate_all_in_range_int(factor_indices, n_selected_factors, ierr, min=1, max=n_factors)
        call validate_all_in_range_int(dependent_indices, n_selected_dependents, ierr, min=1, max=n_factors)

        if (is_err(ierr)) return

        do i_sample = 1, n_samples
            ! create factor vectors for current sample
            do i_sel_factor = 1, n_selected_factors
                i_factor = factor_indices(i_sel_factor)
                do i_timepoint = 1, n_timepoints
                    temp_factors(i_timepoint, i_sel_factor) = trajectories(i_factor, i_sample, i_timepoint)
                end do
            end do

            ! calculate contributions for each factor-dependent combination
            do i_sel_dependent = 1, n_selected_dependents
                ! create dependent vector for current sample
                i_dependent = dependent_indices(i_sel_dependent)
                do i_timepoint = 1, n_timepoints
                    temp_dependent(i_timepoint) = trajectories(i_dependent, i_sample, i_timepoint)
                end do

                do i_sel_factor = 1, n_selected_factors
                    call compute_contributions_helper(temp_factors(:, i_sel_factor), temp_dependent, n_timepoints, mode, local_contributions(:, i_sel_factor, i_sel_dependent, i_sample), total_contributions(i_sel_factor, i_sel_dependent, i_sample), ierr)
                    if (is_err(ierr)) return
                end do
            end do
        end do
    end subroutine compute_all_contributions

   !> Compute scalar baselines for a factor and dependent variable time series.
    pure subroutine compute_baselines_factor_dependent(n_timepoints, factor, dependent, mode, &
                                                       factor_baseline, dependent_baseline, ierr)
        
        integer(int32), intent(in) :: n_timepoints
            !! Number of timepoints in both factor and dependent arrays
        real(real64), intent(in)  :: factor(n_timepoints)
            !! Factor time series, length n_timepoints
        real(real64), intent(in)  :: dependent(n_timepoints)
            !! Dependent variable time series, length n_timepoints
        integer(int32), intent(in) :: mode
            !! Baseline mode: 1=RAW, 2=MIN, 3=MEAN
        real(real64), intent(out) :: factor_baseline
            !! Computed baseline for factor
        real(real64), intent(out) :: dependent_baseline
            !! Computed baseline for dependent variable
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
        factor_baseline = 0.0_real64
        dependent_baseline = 0.0_real64

        ! Validate that n_timepoints > 0
        call validate_dimension_size(n_timepoints, ierr)
        if (is_err(ierr)) return

        select case (mode)

        case (BASELINE_RAW)
            ! Raw contributions: no centering
            factor_baseline = 0.0_real64
            dependent_baseline = 0.0_real64

        case (BASELINE_MIN)
            ! Minimum-centered contributions
            factor_baseline = minval(factor)
            dependent_baseline = minval(dependent)

        case (BASELINE_MEAN)
            ! Mean-centered contributions
            factor_baseline = sum(factor) / real(n_timepoints, kind=real64)
            dependent_baseline = sum(dependent) / real(n_timepoints, kind=real64)

        case default
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end select

        ! Validate that baselines are finite (non-NaN, non-Inf)
        call validate_in_range_real(factor_baseline, ierr)
        call validate_in_range_real(dependent_baseline, ierr)

    end subroutine compute_baselines_factor_dependent

    !> Helper to map baseline mode string ("min", "mean", "raw") to integer constant
    pure subroutine get_baseline_mode(str, mode, ierr)
        character(len=*), intent(in) :: str
            !! mode string ("min", "mean", "raw")
        integer(int32), intent(out) :: mode
            !! integer representation for the mode passed by `str`
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        select case (trim(str))
            case ("raw")
                mode = BASELINE_RAW
            case ("min")
                mode = BASELINE_MIN
            case ("mean")
                mode = BASELINE_MEAN
            case default
                call set_err(ierr, ERR_INVALID_INPUT)
        end select
    end subroutine get_baseline_mode

!> Compute velocity trajectories from position trajectories

subroutine compute_velocity_trajectories(trajectories, velocity, &
                                         n_samples, n_timepoints, n_variables, ierr)

    integer(int32), intent(in)  :: n_samples, n_timepoints, n_variables
    integer(int32), intent(out) :: ierr
    real(real64), intent(in)  :: trajectories(n_samples, n_timepoints, n_variables)
    real(real64), intent(out) :: velocity(n_samples, n_timepoints, n_variables)

    integer(int32) :: sample, var, t

    call set_ok(ierr)
    call validate_dimension_size(n_samples, ierr)
    call validate_dimension_size(n_timepoints, ierr)
    call validate_dimension_size(n_variables, ierr)
    if (is_err(ierr)) return

    velocity = 0.0_real64
    do sample = 1, n_samples
        do var = 1, n_variables
            do t = 2, n_timepoints
                velocity(sample, t, var) = trajectories(sample, t, var) - &
                                           trajectories(sample, t - 1, var)
            end do
        end do
    end do
end subroutine compute_velocity_trajectories

!> Compute acceleration trajectories from velocity trajectories
subroutine compute_acceleration_from_velocity(velocity, acceleration, &
                                              n_samples, n_timepoints, n_variables, ierr)

    integer(int32), intent(in)  :: n_samples, n_timepoints, n_variables
    integer(int32), intent(out) :: ierr
    real(real64), intent(in)  :: velocity(n_samples, n_timepoints, n_variables)
    real(real64), intent(out) :: acceleration(n_samples, n_timepoints, n_variables)

    integer(int32) :: sample, var, t

    call set_ok(ierr)

    call validate_dimension_size(n_samples, ierr)
    call validate_dimension_size(n_timepoints, ierr)
    call validate_dimension_size(n_variables, ierr)
    if (is_err(ierr)) return

    acceleration = 0.0_real64
    if (n_timepoints <= 2) return

    do sample = 1, n_samples
        do var = 1, n_variables
            do t = 3, n_timepoints
                acceleration(sample, t, var) = velocity(sample, t, var) - velocity(sample, t - 1, var)
            end do
        end do
    end do
end subroutine compute_acceleration_from_velocity

!> Compute velocity and acceleration contributions for all variable pairs in the trajectories
subroutine compute_velocity_acceleration_contributions(trajectories, n_samples, n_timepoints, n_variables, mode, velocity, acceleration, &
    factor_velocity, dependent_velocity, velocity_contributions, &
    factor_acceleration, dependent_acceleration, acceleration_contributions, &
    C_velocity, velocity_contribution_series, &
    C_acceleration, acceleration_contribution_series, ierr)

    integer(int32), intent(in) :: n_samples, n_timepoints, n_variables
    integer(int32), intent(in) :: mode
    real(real64),   intent(in) :: trajectories(n_samples, n_timepoints, n_variables)

    ! Workspace (preallocated by caller)
    real(real64), intent(out) :: velocity(n_samples, n_timepoints, n_variables)
    real(real64), intent(out) :: acceleration(n_samples, n_timepoints, n_variables)

    real(real64), intent(inout) :: factor_velocity(n_timepoints-1)
    real(real64), intent(inout) :: dependent_velocity(n_timepoints-1)
    real(real64), intent(inout) :: velocity_contributions(n_timepoints-1)

    real(real64), intent(inout) :: factor_acceleration(n_timepoints-2)
    real(real64), intent(inout) :: dependent_acceleration(n_timepoints-2)
    real(real64), intent(inout) :: acceleration_contributions(n_timepoints-2)

    ! Outputs
    real(real64), intent(out) :: C_velocity(n_samples, n_variables, n_variables)
    real(real64), intent(out) :: velocity_contribution_series(n_samples, n_variables, n_variables, n_timepoints)
    real(real64), intent(out) :: C_acceleration(n_samples, n_variables, n_variables)
    real(real64), intent(out) :: acceleration_contribution_series(n_samples, n_variables, n_variables, n_timepoints)

    integer(int32), intent(out) :: ierr

    integer(int32) :: sample, factor_index, dependent_index, time_index
    real(real64)   :: total_velocity_contribution, total_acceleration_contribution
    integer(int32) :: n_vel, n_acc

    call set_ok(ierr)

    call validate_dimension_size(n_samples, ierr)
    call validate_dimension_size(n_timepoints, ierr)
    call validate_dimension_size(n_variables, ierr)
    if (is_err(ierr)) return

    C_velocity                     = 0.0_real64
    velocity_contribution_series   = 0.0_real64
    C_acceleration                 = 0.0_real64
    acceleration_contribution_series = 0.0_real64

    ! ---- Step 1: velocity ----
    call compute_velocity_trajectories(trajectories, velocity, n_samples, n_timepoints, n_variables, ierr)
    if (is_err(ierr)) return

    ! ---- Step 2: acceleration from velocity ----
    call compute_acceleration_from_velocity(velocity, acceleration, n_samples, n_timepoints, n_variables, ierr)
    if (is_err(ierr)) return

    n_vel = n_timepoints - 1_int32
    n_acc = n_timepoints - 2_int32

    do sample = 1, n_samples
        do factor_index = 1, n_variables
            do dependent_index = 1, n_variables

                ! -------- Velocity contributions (t = 2..T) --------
                if (n_vel > 0) then
                    do time_index = 2, n_timepoints
                        factor_velocity(time_index-1)    = velocity(sample, time_index, factor_index)
                        dependent_velocity(time_index-1) = velocity(sample, time_index, dependent_index)
                    end do

                    call compute_contributions( &
                        factor_velocity, dependent_velocity, n_vel, mode, &
                        velocity_contributions, total_velocity_contribution, ierr)
                    if (is_err(ierr)) return

                    C_velocity(sample, factor_index, dependent_index) = total_velocity_contribution

                    do time_index = 2, n_timepoints
                        velocity_contribution_series(sample, factor_index, dependent_index, time_index) = &
                            velocity_contributions(time_index-1)
                    end do
                end if

                ! -------- Acceleration contributions (t = 3..T) --------
                if (n_acc > 0) then
                    do time_index = 3, n_timepoints
                        factor_acceleration(time_index-2)    = acceleration(sample, time_index, factor_index)
                        dependent_acceleration(time_index-2) = acceleration(sample, time_index, dependent_index)
                    end do

                    call compute_contributions( &
                        factor_acceleration, dependent_acceleration, n_acc, mode, &
                        acceleration_contributions, total_acceleration_contribution, ierr)
                    if (is_err(ierr)) return

                    C_acceleration(sample, factor_index, dependent_index) = total_acceleration_contribution

                    do time_index = 3, n_timepoints
                        acceleration_contribution_series(sample, factor_index, dependent_index, time_index) = &
                            acceleration_contributions(time_index-2)
                    end do
                end if

            end do
        end do
    end do
end subroutine compute_velocity_acceleration_contributions

subroutine compute_velocity_acceleration_contributions_alloc(trajectories, n_samples, n_timepoints, n_variables, mode, &
    C_velocity, velocity_contribution_series, &
    C_acceleration, acceleration_contribution_series, ierr)

    integer(int32), intent(in) :: n_samples, n_timepoints, n_variables
    integer(int32), intent(in) :: mode
    real(real64),   intent(in) :: trajectories(n_samples, n_timepoints, n_variables)

    real(real64), intent(out) :: C_velocity(n_samples, n_variables, n_variables)
    real(real64), intent(out) :: velocity_contribution_series(n_samples, n_variables, n_variables, n_timepoints)
    real(real64), intent(out) :: C_acceleration(n_samples, n_variables, n_variables)
    real(real64), intent(out) :: acceleration_contribution_series(n_samples, n_variables, n_variables, n_timepoints)
    integer(int32), intent(out) :: ierr

    ! Workspace (allocated once here)
    real(real64), allocatable :: velocity(:,:,:)
    real(real64), allocatable :: acceleration(:,:,:)

    real(real64), allocatable :: factor_velocity(:), dependent_velocity(:), velocity_contributions(:)
    real(real64), allocatable :: factor_acceleration(:), dependent_acceleration(:), acceleration_contributions(:)

    integer :: stat_alloc

    call set_ok(ierr)

    call validate_dimension_size(n_samples, ierr)
    call validate_dimension_size(n_timepoints, ierr)
    call validate_dimension_size(n_variables, ierr)
    if (is_err(ierr)) return

    ! Allocate big work arrays once
    allocate(velocity(n_samples, n_timepoints, n_variables), stat=stat_alloc)
    if (stat_alloc /= 0) then
        call set_err(ierr, ERR_ALLOC_FAIL)
        return
    end if

    allocate(acceleration(n_samples, n_timepoints, n_variables), stat=stat_alloc)
    if (stat_alloc /= 0) then
        call set_err(ierr, ERR_ALLOC_FAIL)
        deallocate(velocity)
        return
    end if

    ! Allocate 1D reusable work vectors once
    if (n_timepoints > 1) then
        allocate(factor_velocity(n_timepoints-1), dependent_velocity(n_timepoints-1), velocity_contributions(n_timepoints-1), stat=stat_alloc)
        if (stat_alloc /= 0) then
            call set_err(ierr, ERR_ALLOC_FAIL)
            deallocate(acceleration, velocity)
            return
        end if
    end if

    if (n_timepoints > 2) then
        allocate(factor_acceleration(n_timepoints-2), dependent_acceleration(n_timepoints-2), acceleration_contributions(n_timepoints-2), stat=stat_alloc)
        if (stat_alloc /= 0) then
            call set_err(ierr, ERR_ALLOC_FAIL)
            if (allocated(factor_velocity)) deallocate(factor_velocity, dependent_velocity, velocity_contributions)
            deallocate(acceleration, velocity)
            return
        end if
    end if

    ! Call the SK routine (no allocation inside)
    call compute_velocity_acceleration_contributions(trajectories, n_samples, n_timepoints, n_variables, mode, &
        velocity, acceleration, &
        factor_velocity, dependent_velocity, velocity_contributions, &
        factor_acceleration, dependent_acceleration, acceleration_contributions, &
        C_velocity, velocity_contribution_series, &
        C_acceleration, acceleration_contribution_series, ierr)

    ! Cleanup once (even on error)
    if (allocated(factor_acceleration)) deallocate(factor_acceleration, dependent_acceleration, acceleration_contributions)
    if (allocated(factor_velocity))     deallocate(factor_velocity, dependent_velocity, velocity_contributions)
    if (allocated(acceleration))        deallocate(acceleration)
    if (allocated(velocity))            deallocate(velocity)

end subroutine compute_velocity_acceleration_contributions_alloc


end module tox_trajectory_contribution_analysis

!> C-compatible wrapper for [[tox_trajectory_contribution_analysis(module):compute_all_contributions(subroutine)]]
pure subroutine compute_all_contributions_c(trajectories, n_factors, n_samples, n_timepoints, &
    factor_indices, n_selected_factors, dependent_indices, n_selected_dependents, mode, &
    local_contributions, total_contributions, temp_factors, temp_dependent, ierr) &
    bind(C, name="compute_all_contributions_c")

    use, intrinsic :: iso_fortran_env, only: int32
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char
    use tox_trajectory_contribution_analysis, only: compute_all_contributions, get_baseline_mode
    use tox_conversions, only: c_char_1d_as_string
    use tox_errors, only: is_err
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_factors
        !! number of factors
    integer(c_int), intent(in), target :: n_samples
        !! number of samples
    integer(c_int), intent(in), target :: n_timepoints
        !! number of timepoints
    integer(c_int), intent(in), target :: n_selected_factors
        !! number of selected factors in `factor_indices`
    integer(c_int), intent(in), target :: n_selected_dependents
        !! number of selected dependents in `dependent_indices`
    real(c_double), dimension(n_factors, n_samples, n_timepoints), intent(in), target :: trajectories
        !! trajectories array: (n_factors, n_samples, n_timepoints)
    integer(c_int), dimension(n_selected_factors), intent(in), target :: factor_indices
        !! indices of factors to compute the contributions for
    integer(c_int), dimension(n_selected_dependents), intent(in), target :: dependent_indices
        !! indices of dependents to compute the contributions for
    character(len=1, kind=c_char), dimension(4), intent(in), target :: mode
        !! Baseline mode: "raw", "min", "mean"
    real(c_double), dimension(n_timepoints, n_selected_factors, n_selected_dependents, n_samples), intent(out), target :: local_contributions
        !! Per-timepoint contributions per sample-dependent-factor combination
    real(c_double), dimension(n_selected_factors, n_selected_dependents, n_samples), intent(out), target :: total_contributions
        !! Total contribution (`sum(local_contributions)`) per sample-dependent-factor combination
    real(c_double), dimension(n_timepoints, n_selected_factors), intent(out), target :: temp_factors
        !! Working array to hold the currently handled sample's factors in contiguous memory
    real(c_double), dimension(n_timepoints), intent(out), target :: temp_dependent
        !! Working array to hold the currently handled dependent in contiguous memory
    integer(c_int), intent(out), target :: ierr
        !! Error code

    character(len=:), allocatable :: mode_f
    integer(int32) :: mode_int

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_factors)
    M_CHECK_NON_NULL(n_samples)
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(n_selected_factors)
    M_CHECK_NON_NULL(n_selected_dependents)
    M_CHECK_NON_NULL(trajectories)
    M_CHECK_NON_NULL(factor_indices)
    M_CHECK_NON_NULL(dependent_indices)
    M_CHECK_NON_NULL(mode)
    M_CHECK_NON_NULL(local_contributions)
    M_CHECK_NON_NULL(total_contributions)
    M_CHECK_NON_NULL(temp_factors)
    M_CHECK_NON_NULL(temp_dependent)

    call c_char_1d_as_string(mode, mode_f, ierr)
    if (is_err(ierr)) return

    call get_baseline_mode(mode_f, mode_int, ierr)
    if (is_err(ierr)) return
    
    call compute_all_contributions(trajectories, n_factors, n_samples, n_timepoints, &
        factor_indices, n_selected_factors, dependent_indices, n_selected_dependents, mode_int, &
        local_contributions, total_contributions, temp_factors, temp_dependent, ierr)
end subroutine compute_all_contributions_c

!> C-compatible wrapper for [[tox_trajectory_contribution_analysis(module):trajectory_contribution(subroutine)]]
pure subroutine trajectory_contribution_c(factor, dependent, n_timepoints, mode, contribution, ierr) &
    bind(C, name="trajectory_contribution_c")

    use tox_trajectory_contribution_analysis, only: trajectory_contribution
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_errors, only: is_err
    M_USE_NULL_VALIDATION
    implicit none

    real(c_double), intent(in), target :: factor(n_timepoints)
    real(c_double), intent(in), target :: dependent(n_timepoints)
    integer(c_int), intent(in), target :: n_timepoints
    integer(c_int), intent(in), target :: mode
    real(c_double), intent(out), target :: contribution
    integer(c_int), intent(out), target :: ierr

    integer(int32) :: n_tp, mode_int
    integer(int32) :: ierr_local
    real(real64) :: contrib_local

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(factor)
    M_CHECK_NON_NULL(dependent)
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(mode)
    M_CHECK_NON_NULL(contribution)

    n_tp = int(n_timepoints, int32)
    mode_int = int(mode, int32)

    contrib_local = 0.0_real64
    ierr_local = 0_int32

    call trajectory_contribution(factor, dependent, n_tp, mode_int, contrib_local, ierr_local)

    contribution = real(contrib_local, kind=c_double)
    ierr = int(ierr_local, kind=c_int)
end subroutine trajectory_contribution_c

!> C-compatible wrapper for [[tox_trajectory_contribution_analysis(module):spike_contribution(subroutine)]]
pure subroutine spike_contribution_c(factor, dependent, n_timepoints, mode, contribution, ierr) &
    bind(C, name="spike_contribution_c")

    use tox_trajectory_contribution_analysis, only: spike_contribution
    use, intrinsic :: iso_fortran_env, only: int32
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_errors, only: is_err
    M_USE_NULL_VALIDATION
    implicit none

    real(c_double), intent(in), target :: factor(n_timepoints)
    real(c_double), intent(in), target :: dependent(n_timepoints)
    integer(c_int), intent(in), target :: n_timepoints
    integer(c_int), intent(in), target :: mode
    real(c_double), intent(out), target :: contribution(n_timepoints)
    integer(c_int), intent(out), target :: ierr

    integer(int32) :: n_tp, mode_int
    integer(int32) :: ierr_local

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(factor)
    M_CHECK_NON_NULL(dependent)
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(mode)
    M_CHECK_NON_NULL(contribution)

    n_tp = int(n_timepoints, int32)
    mode_int = int(mode, int32)

    ierr_local = 0_int32
    contribution = 0.0_c_double

    call spike_contribution(factor, dependent, n_tp, mode_int, contribution, ierr_local)

    ierr = int(ierr_local, kind=c_int)
end subroutine spike_contribution_c

!> C-compatible wrapper for [[tox_trajectory_contribution_analysis(module):calc_contributions_alloc(subroutine)]]
subroutine calc_contributions_c(trajectories, n_factors, n_samples, n_timepoints, i_factor, dependent_idx, mode, &
                               spike_contribs, integrated_contribs, ierr) &
    bind(C, name="calc_contributions_c")

    use tox_trajectory_contribution_analysis, only: calc_contributions_alloc
    use, intrinsic :: iso_fortran_env, only: int32
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_errors, only: set_ok, set_err, is_err, ERR_INVALID_INPUT
    M_USE_NULL_VALIDATION
    implicit none

    ! Dimensions and indices
    integer(c_int), intent(in), target :: n_factors
    integer(c_int), intent(in), target :: n_samples
    integer(c_int), intent(in), target :: n_timepoints
    integer(c_int), intent(in), target :: i_factor
    integer(c_int), intent(in), target :: dependent_idx
    integer(c_int), intent(in), target :: mode

    ! Arrays
    real(c_double), intent(in), target :: trajectories(n_factors, n_samples, n_timepoints)
    real(c_double), intent(out), target :: spike_contribs(n_timepoints, n_samples)
    real(c_double), intent(out), target :: integrated_contribs(n_samples)

    ! Error code
    integer(c_int), intent(out), target :: ierr

    integer(int32) :: nf, ns, nt
    integer(int32) :: i_factor_f, dependent_f, mode_f
    integer(int32) :: ierr_local

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_factors)
    M_CHECK_NON_NULL(n_samples)
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(i_factor)
    M_CHECK_NON_NULL(dependent_idx)
    M_CHECK_NON_NULL(mode)
    M_CHECK_NON_NULL(trajectories)
    M_CHECK_NON_NULL(spike_contribs)
    M_CHECK_NON_NULL(integrated_contribs)

    nf = int(n_factors, int32)
    ns = int(n_samples, int32)
    nt = int(n_timepoints, int32)
    i_factor_f = int(i_factor, int32) + 1_int32
    dependent_f = int(dependent_idx, int32) + 1_int32
    mode_f = int(mode, int32)

    call set_ok(ierr_local)

    if (nf < 1 .or. ns < 1 .or. nt < 1) then
        call set_err(ierr_local, ERR_INVALID_INPUT)
    else if (i_factor_f < 1 .or. i_factor_f > nf .or. dependent_f < 1 .or. dependent_f > nf) then
        call set_err(ierr_local, ERR_INVALID_INPUT)
    end if

    if (.not. is_err(ierr_local)) then
        call calc_contributions_alloc(trajectories, nf, ns, nt, i_factor_f, dependent_f, mode_f, spike_contribs, integrated_contribs, ierr_local)
    end if

    ierr = int(ierr_local, c_int)
end subroutine calc_contributions_c

!> C-compatible wrapper for [[tox_trajectory_contribution_analysis(module):calc_contributions(subroutine)]] with caller work buffers
subroutine calc_contributions_expert_c(trajectories, n_factors, n_samples, n_timepoints, i_factor, dependent_idx, mode, &
                                       spike_contribs, integrated_contribs, temp_factor_vector, temp_dependent_vector, ierr) &
    bind(C, name="calc_contributions_expert_c")

    use tox_trajectory_contribution_analysis, only: calc_contributions
    use, intrinsic :: iso_fortran_env, only: int32
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_errors, only: set_ok, set_err, is_err, ERR_INVALID_INPUT
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_factors
    integer(c_int), intent(in), target :: n_samples
    integer(c_int), intent(in), target :: n_timepoints
    integer(c_int), intent(in), target :: i_factor
    integer(c_int), intent(in), target :: dependent_idx
    integer(c_int), intent(in), target :: mode
    real(c_double), intent(in), target :: trajectories(n_factors, n_samples, n_timepoints)
    real(c_double), intent(out), target :: spike_contribs(n_timepoints, n_samples)
    real(c_double), intent(out), target :: integrated_contribs(n_samples)
    real(c_double), intent(inout), target :: temp_factor_vector(n_timepoints)
    real(c_double), intent(inout), target :: temp_dependent_vector(n_timepoints)
    integer(c_int), intent(out), target :: ierr

    integer(int32) :: nf, ns, nt
    integer(int32) :: i_factor_f, dependent_f, mode_f
    integer(int32) :: ierr_local

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_factors)
    M_CHECK_NON_NULL(n_samples)
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(i_factor)
    M_CHECK_NON_NULL(dependent_idx)
    M_CHECK_NON_NULL(mode)
    M_CHECK_NON_NULL(trajectories)
    M_CHECK_NON_NULL(spike_contribs)
    M_CHECK_NON_NULL(integrated_contribs)
    M_CHECK_NON_NULL(temp_factor_vector)
    M_CHECK_NON_NULL(temp_dependent_vector)

    nf = int(n_factors, int32)
    ns = int(n_samples, int32)
    nt = int(n_timepoints, int32)
    i_factor_f = int(i_factor, int32) + 1_int32
    dependent_f = int(dependent_idx, int32) + 1_int32
    mode_f = int(mode, int32)

    call set_ok(ierr_local)

    if (nf < 1 .or. ns < 1 .or. nt < 1) then
        call set_err(ierr_local, ERR_INVALID_INPUT)
    else if (i_factor_f < 1 .or. i_factor_f > nf .or. dependent_f < 1 .or. dependent_f > nf) then
        call set_err(ierr_local, ERR_INVALID_INPUT)
    end if

    if (.not. is_err(ierr_local)) then
        call calc_contributions(trajectories, nf, ns, nt, i_factor_f, dependent_f, mode_f, spike_contribs, integrated_contribs, temp_factor_vector, temp_dependent_vector, ierr_local)
    end if

    ierr = int(ierr_local, c_int)
end subroutine calc_contributions_expert_c

!> C-compatible wrapper for [[tox_trajectory_contribution_analysis(module):compute_contributions(subroutine)]]
pure subroutine compute_contributions_c(factor, dependent, n_dims, mode, local_contributions, total_contribution, ierr) &
    bind(C, name="compute_contributions_c")
    use, intrinsic :: iso_fortran_env, only: int32
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char
    use tox_trajectory_contribution_analysis, only: compute_contributions, get_baseline_mode
    use tox_conversions, only: c_char_1d_as_string
    use tox_errors, only: is_err
    M_USE_NULL_VALIDATION
    implicit none

    ! Arguments mapped to C types
    integer(c_int), intent(in), target :: n_dims
        !! Number of elements in `factor` and `dependent`
    real(c_double), dimension(n_dims), intent(in), target :: factor
        !! Factor time series
    real(c_double), dimension(n_dims), intent(in), target :: dependent
        !! Dependent variable time series
    character(len=1, kind=c_char), dimension(4), intent(in), target :: mode
        !! Baseline mode: "raw", "min", "mean"
    real(c_double), dimension(n_dims), intent(out), target :: local_contributions
        !! Per-element contributions
    real(c_double), intent(out), target :: total_contribution
        !! Total contribution
    integer(c_int), intent(out), target :: ierr
        !! Error code

    character(len=:), allocatable :: mode_f
    integer(int32) :: mode_int

    ! Null checks
    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_dims)
    M_CHECK_NON_NULL(factor)
    M_CHECK_NON_NULL(dependent)
    M_CHECK_NON_NULL(mode)
    M_CHECK_NON_NULL(local_contributions)
    M_CHECK_NON_NULL(total_contribution)

    call c_char_1d_as_string(mode, mode_f, ierr)
    if (is_err(ierr)) return

    call get_baseline_mode(mode_f, mode_int, ierr)
    if (is_err(ierr)) return

    call compute_contributions(factor, dependent, n_dims, mode_int, local_contributions, total_contribution, ierr)
end subroutine compute_contributions_c

!> C-compatible wrapper for [[tox_trajectory_contribution_analysis(module):compute_baselines_factor_dependent(subroutine)]]
subroutine compute_baselines_factor_dependent_c(factor, dependent, n_timepoints, mode, &
                                               factor_baseline, dependent_baseline, ierr) &
    bind(C, name="tox_compute_baselines_factor_dependent")

    use tox_trajectory_contribution_analysis, only: compute_baselines_factor_dependent, BASELINE_RAW, BASELINE_MIN, BASELINE_MEAN
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use tox_errors, only: set_ok, is_err, set_err, ERR_INVALID_INPUT
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), value :: n_timepoints
        !! Number of timepoints in both factor and dependent arrays
    integer(c_int), value :: mode
        !! Baseline mode: 1=RAW, 2=MIN, 3=MEAN
    real(c_double), intent(in),  target :: factor(n_timepoints)
        !! Factor time series, length n_timepoints
    real(c_double), intent(in),  target :: dependent(n_timepoints)
        !! Dependent variable time series, length n_timepoints
    real(c_double), intent(out), target :: factor_baseline
        !! Computed baseline for factor
    real(c_double), intent(out), target :: dependent_baseline
        !! Computed baseline for dependent variable
    integer(c_int), intent(out), target :: ierr
        !! Error code

    integer(int32) :: n_tp, mode_int
    integer(int32) :: ierr_local
    real(real64) :: factor_baseline_local
    real(real64) :: dependent_baseline_local

    !! Null-pointer validation 
    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(factor)
    M_CHECK_NON_NULL(dependent)
    M_CHECK_NON_NULL(factor_baseline)
    M_CHECK_NON_NULL(dependent_baseline)

    n_tp = int(n_timepoints, int32)
    mode_int = int(mode, int32)

    call set_ok(ierr_local)

    if (n_tp < 1) then
        call set_err(ierr_local, ERR_INVALID_INPUT)
    else if (mode_int < BASELINE_RAW .or. mode_int > BASELINE_MEAN) then
        call set_err(ierr_local, ERR_INVALID_INPUT)
    end if

    if (.not. is_err(ierr_local)) then
        factor_baseline_local = 0.0_real64
        dependent_baseline_local = 0.0_real64
        call compute_baselines_factor_dependent(n_tp, factor, dependent, mode_int, factor_baseline_local, dependent_baseline_local, ierr_local)
        factor_baseline = real(factor_baseline_local, kind=c_double)
        dependent_baseline = real(dependent_baseline_local, kind=c_double)
    else
        factor_baseline = 0.0_c_double
        dependent_baseline = 0.0_c_double
    end if

    ierr = int(ierr_local, kind=c_int)
end subroutine compute_baselines_factor_dependent_c

subroutine process_trajectories_c(trajectories, n_factors, n_samples, n_timepoints, n_processed, factor_mask, &
                                 dependent_idx, mode, percentile, integrated_contribs, spike_contribs, &
                                 thresholds_integrated, outliers_integrated, thresholds_spike, outliers_spike, ierr) &
    bind(C, name="process_trajectories_C")

    use tox_trajectory_contribution_analysis, only: process_trajectories_alloc
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_errors, only: set_ok, set_err, is_err, ERR_INVALID_INPUT, ERR_ALLOC_FAIL
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_factors
    integer(c_int), intent(in), target :: n_samples
    integer(c_int), intent(in), target :: n_timepoints
    integer(c_int), intent(in), target :: n_processed
    real(c_double), intent(in), target :: trajectories(n_factors, n_samples, n_timepoints)
    integer(c_int), intent(in), target :: factor_mask(n_factors)
    integer(c_int), intent(in), target :: dependent_idx
    integer(c_int), intent(in), target :: mode
    real(c_double), intent(in), target :: percentile
    real(c_double), intent(out), target :: integrated_contribs(n_samples, n_processed)
    real(c_double), intent(out), target :: spike_contribs(n_timepoints, n_samples, n_processed)
    real(c_double), intent(out), target :: thresholds_integrated(n_processed)
    integer(c_int), intent(out), target :: outliers_integrated(n_samples, n_processed)
    real(c_double), intent(out), target :: thresholds_spike(n_timepoints, n_processed)
    integer(c_int), intent(out), target :: outliers_spike(n_timepoints, n_samples, n_processed)
    integer(c_int), intent(out), target :: ierr

    integer(int32) :: nf, ns, nt, n_proc
    integer(int32) :: dependent_f, mode_f
    real(real64) :: percentile_f
    integer(int32) :: ierr_local
    integer(int32) :: mask_true_count
    integer :: stat_alloc
    logical, allocatable :: factor_mask_f(:)
    logical, allocatable :: outliers_integrated_f(:,:)
    logical, allocatable :: outliers_spike_f(:,:,:)

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_factors)
    M_CHECK_NON_NULL(n_samples)
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(n_processed)
    M_CHECK_NON_NULL(trajectories)
    M_CHECK_NON_NULL(factor_mask)
    M_CHECK_NON_NULL(dependent_idx)
    M_CHECK_NON_NULL(mode)
    M_CHECK_NON_NULL(percentile)
    M_CHECK_NON_NULL(integrated_contribs)
    M_CHECK_NON_NULL(spike_contribs)
    M_CHECK_NON_NULL(thresholds_integrated)
    M_CHECK_NON_NULL(outliers_integrated)
    M_CHECK_NON_NULL(thresholds_spike)
    M_CHECK_NON_NULL(outliers_spike)

    integrated_contribs = 0.0_c_double
    spike_contribs = 0.0_c_double
    thresholds_integrated = 0.0_c_double
    thresholds_spike = 0.0_c_double
    outliers_integrated = 0_c_int
    outliers_spike = 0_c_int

    nf = int(n_factors, int32)
    ns = int(n_samples, int32)
    nt = int(n_timepoints, int32)
    n_proc = int(n_processed, int32)
    dependent_f = int(dependent_idx, int32)
    mode_f = int(mode, int32)
    percentile_f = real(percentile, real64)

    call set_ok(ierr_local)

    if (nf < 1 .or. ns < 1 .or. nt < 1) then
        call set_err(ierr_local, ERR_INVALID_INPUT)
    else if (n_proc < 1 .or. n_proc > nf) then
        call set_err(ierr_local, ERR_INVALID_INPUT)
    else if (dependent_f < 1 .or. dependent_f > nf) then
        call set_err(ierr_local, ERR_INVALID_INPUT)
    end if

    if (.not. is_err(ierr_local)) then
        allocate(factor_mask_f(nf), stat=stat_alloc)
        if (stat_alloc /= 0) then
            call set_err(ierr_local, ERR_ALLOC_FAIL)
        else
            factor_mask_f = .false.
        end if
    end if

    if (.not. is_err(ierr_local)) then
        factor_mask_f = (factor_mask /= 0_c_int)
        mask_true_count = int(count(factor_mask_f), int32)
        if (mask_true_count /= n_proc) then
            call set_err(ierr_local, ERR_INVALID_INPUT)
        end if
    end if

    if (.not. is_err(ierr_local)) then
        allocate(outliers_integrated_f(ns, n_proc), stat=stat_alloc)
        if (stat_alloc /= 0) then
            call set_err(ierr_local, ERR_ALLOC_FAIL)
        else
            outliers_integrated_f = .false.
        end if
    end if

    if (.not. is_err(ierr_local)) then
        allocate(outliers_spike_f(nt, ns, n_proc), stat=stat_alloc)
        if (stat_alloc /= 0) then
            call set_err(ierr_local, ERR_ALLOC_FAIL)
        else
            outliers_spike_f = .false.
        end if
    end if

    if (.not. is_err(ierr_local)) then
        call process_trajectories_alloc(trajectories, nf, ns, nt, factor_mask_f, n_proc, dependent_f, mode_f, percentile_f, &
            integrated_contribs, spike_contribs, thresholds_integrated, outliers_integrated_f, thresholds_spike, outliers_spike_f, ierr_local)
    end if

    if (.not. is_err(ierr_local)) then
        outliers_integrated = merge(1_c_int, 0_c_int, outliers_integrated_f)
        outliers_spike = merge(1_c_int, 0_c_int, outliers_spike_f)
    else
        outliers_integrated = 0_c_int
        outliers_spike = 0_c_int
    end if

    if (allocated(outliers_spike_f)) deallocate(outliers_spike_f)
    if (allocated(outliers_integrated_f)) deallocate(outliers_integrated_f)
    if (allocated(factor_mask_f)) deallocate(factor_mask_f)

    ierr = int(ierr_local, c_int)
end subroutine process_trajectories_c

!> C-compatible wrapper for [[tox_trajectory_contribution_analysis(module):process_trajectories_flat_alloc(subroutine)]]
subroutine process_trajectories_flat_c(trajectories, n_factors, n_samples, n_timepoints, n_processed, factor_mask, &
                                       dependent_idx, mode, percentile, integrated_contribs, spike_contribs, &
                                       thresholds_integrated, outliers_integrated, thresholds_spike, outliers_spike, ierr) &
    bind(C, name="process_trajectories_flat_C")

    use tox_trajectory_contribution_analysis, only: process_trajectories_flat_alloc
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_errors, only: set_ok, set_err, is_err, ERR_INVALID_INPUT, ERR_ALLOC_FAIL
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_factors
    integer(c_int), intent(in), target :: n_samples
    integer(c_int), intent(in), target :: n_timepoints
    integer(c_int), intent(in), target :: n_processed
    real(c_double), intent(in), target :: trajectories(n_factors, n_samples, n_timepoints)
    integer(c_int), intent(in), target :: factor_mask(n_factors)
    integer(c_int), intent(in), target :: dependent_idx
    integer(c_int), intent(in), target :: mode
    real(c_double), intent(in), target :: percentile
    real(c_double), intent(out), target :: integrated_contribs(n_samples, n_processed)
    real(c_double), intent(out), target :: spike_contribs(n_timepoints, n_samples, n_processed)
    real(c_double), intent(out), target :: thresholds_integrated(n_processed)
    integer(c_int), intent(out), target :: outliers_integrated(n_samples, n_processed)
    real(c_double), intent(out), target :: thresholds_spike(n_processed)
    integer(c_int), intent(out), target :: outliers_spike(n_timepoints, n_samples, n_processed)
    integer(c_int), intent(out), target :: ierr

    integer(int32) :: nf, ns, nt, n_proc
    integer(int32) :: dependent_f, mode_f
    real(real64) :: percentile_f
    integer(int32) :: ierr_local
    integer(int32) :: mask_true_count
    integer :: stat_alloc
    logical, allocatable :: factor_mask_f(:)
    logical, allocatable :: outliers_integrated_f(:,:)
    logical, allocatable :: outliers_spike_f(:,:,:)

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_factors)
    M_CHECK_NON_NULL(n_samples)
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(n_processed)
    M_CHECK_NON_NULL(trajectories)
    M_CHECK_NON_NULL(factor_mask)
    M_CHECK_NON_NULL(dependent_idx)
    M_CHECK_NON_NULL(mode)
    M_CHECK_NON_NULL(percentile)
    M_CHECK_NON_NULL(integrated_contribs)
    M_CHECK_NON_NULL(spike_contribs)
    M_CHECK_NON_NULL(thresholds_integrated)
    M_CHECK_NON_NULL(outliers_integrated)
    M_CHECK_NON_NULL(thresholds_spike)
    M_CHECK_NON_NULL(outliers_spike)

    integrated_contribs = 0.0_c_double
    spike_contribs = 0.0_c_double
    thresholds_integrated = 0.0_c_double
    thresholds_spike = 0.0_c_double
    outliers_integrated = 0_c_int
    outliers_spike = 0_c_int

    nf = int(n_factors, int32)
    ns = int(n_samples, int32)
    nt = int(n_timepoints, int32)
    n_proc = int(n_processed, int32)
    dependent_f = int(dependent_idx, int32)
    mode_f = int(mode, int32)
    percentile_f = real(percentile, real64)

    call set_ok(ierr_local)

    if (nf < 1 .or. ns < 1 .or. nt < 1) then
        call set_err(ierr_local, ERR_INVALID_INPUT)
    else if (n_proc < 1 .or. n_proc > nf) then
        call set_err(ierr_local, ERR_INVALID_INPUT)
    else if (dependent_f < 1 .or. dependent_f > nf) then
        call set_err(ierr_local, ERR_INVALID_INPUT)
    end if

    if (.not. is_err(ierr_local)) then
        allocate(factor_mask_f(nf), stat=stat_alloc)
        if (stat_alloc /= 0) then
            call set_err(ierr_local, ERR_ALLOC_FAIL)
        else
            factor_mask_f = .false.
        end if
    end if

    if (.not. is_err(ierr_local)) then
        factor_mask_f = (factor_mask /= 0_c_int)
        mask_true_count = int(count(factor_mask_f), int32)
        if (mask_true_count /= n_proc) then
            call set_err(ierr_local, ERR_INVALID_INPUT)
        end if
    end if

    if (.not. is_err(ierr_local)) then
        allocate(outliers_integrated_f(ns, n_proc), stat=stat_alloc)
        if (stat_alloc /= 0) then
            call set_err(ierr_local, ERR_ALLOC_FAIL)
        else
            outliers_integrated_f = .false.
        end if
    end if

    if (.not. is_err(ierr_local)) then
        allocate(outliers_spike_f(nt, ns, n_proc), stat=stat_alloc)
        if (stat_alloc /= 0) then
            call set_err(ierr_local, ERR_ALLOC_FAIL)
        else
            outliers_spike_f = .false.
        end if
    end if

    if (.not. is_err(ierr_local)) then
        call process_trajectories_flat_alloc(trajectories, nf, ns, nt, factor_mask_f, n_proc, dependent_f, mode_f, percentile_f, &
            integrated_contribs, spike_contribs, thresholds_integrated, outliers_integrated_f, thresholds_spike, outliers_spike_f, ierr_local)
    end if

    if (.not. is_err(ierr_local)) then
        outliers_integrated = merge(1_c_int, 0_c_int, outliers_integrated_f)
        outliers_spike = merge(1_c_int, 0_c_int, outliers_spike_f)
    else
        outliers_integrated = 0_c_int
        outliers_spike = 0_c_int
    end if

    if (allocated(outliers_spike_f)) deallocate(outliers_spike_f)
    if (allocated(outliers_integrated_f)) deallocate(outliers_integrated_f)
    if (allocated(factor_mask_f)) deallocate(factor_mask_f)

    ierr = int(ierr_local, c_int)
end subroutine process_trajectories_flat_c

!> C wrapper for compute_velocity_trajectories
subroutine tox_compute_velocity_trajectories_c(trajectories, n_samples, n_timepoints, n_variables, &
                                               velocity, ierr) &
    bind(C, name="tox_compute_velocity_trajectories")
    use tox_trajectory_contribution_analysis, only : compute_velocity_trajectories
    use, intrinsic :: iso_fortran_env, only: int32
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use tox_errors, only: is_err
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in),  target :: n_samples
    integer(c_int), intent(in),  target :: n_timepoints
    integer(c_int), intent(in),  target :: n_variables
    real(c_double), intent(in),  target :: trajectories(n_samples, n_timepoints, n_variables)
    real(c_double), intent(out), target :: velocity(n_samples, n_timepoints, n_variables)
    integer(c_int), intent(out), target :: ierr

    integer(int32) :: ns, nt, nv

    !! Null-pointer validation 
    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_samples)
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(n_variables)
    M_CHECK_NON_NULL(trajectories)
    M_CHECK_NON_NULL(velocity)

    ns = int(n_samples, int32)
    nt = int(n_timepoints, int32)
    nv = int(n_variables, int32)

    call compute_velocity_trajectories(trajectories, velocity, ns, nt, nv, ierr)
end subroutine tox_compute_velocity_trajectories_c

!> C wrapper for compute_acceleration_from_velocity
subroutine tox_compute_acceleration_from_velocity_c(velocity, n_samples, n_timepoints, n_variables, &
                                                    acceleration, ierr) &
    bind(C, name="tox_compute_acceleration_from_velocity")

    use tox_trajectory_contribution_analysis, only : compute_acceleration_from_velocity
    use, intrinsic :: iso_fortran_env, only: int32
    use, intrinsic :: iso_c_binding, only : c_int, c_double
    use tox_errors, only: is_err

    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in),  target :: n_samples
    integer(c_int), intent(in),  target :: n_timepoints
    integer(c_int), intent(in),  target :: n_variables
    real(c_double), intent(in),  target :: velocity(n_samples, n_timepoints, n_variables)
    real(c_double), intent(out), target :: acceleration(n_samples, n_timepoints, n_variables)
    integer(c_int), intent(out), target :: ierr

    integer(int32) :: ns, nt, nv

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_samples)
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(n_variables)
    M_CHECK_NON_NULL(velocity)
    M_CHECK_NON_NULL(acceleration)

    ns = int(n_samples, int32)
    nt = int(n_timepoints, int32)
    nv = int(n_variables, int32)

    call compute_acceleration_from_velocity(velocity, acceleration, ns, nt, nv, ierr)
end subroutine tox_compute_acceleration_from_velocity_c


!> C wrapper for compute_velocity_and_acceleration_contributions
subroutine tox_compute_velocity_acceleration_contributions_c(trajectories, n_samples, n_timepoints, n_variables, mode, &
    velocity, acceleration, &
    factor_velocity, dependent_velocity, velocity_contributions, &
    factor_acceleration, dependent_acceleration, acceleration_contributions, &
    C_velocity, velocity_contribution_series, &
    C_acceleration, acceleration_contribution_series, ierr) &
    bind(C, name="tox_compute_velocity_acceleration_contributions")
    use tox_trajectory_contribution_analysis, only: compute_velocity_acceleration_contributions

    use, intrinsic :: iso_c_binding, only : c_int, c_double
    use, intrinsic :: iso_fortran_env, only: int32
    use tox_errors, only: is_err
    M_USE_NULL_VALIDATION
    implicit none
    integer(c_int), intent(in),  target :: n_samples
    integer(c_int), intent(in),  target :: n_timepoints
    integer(c_int), intent(in),  target :: n_variables
    integer(c_int), intent(in),  target :: mode
    integer(c_int), intent(out), target :: ierr
    real(c_double), intent(in),  target :: trajectories(n_samples, n_timepoints, n_variables)

    ! ---- Workspace (passed in from C) ----
    real(c_double), intent(out), target :: velocity(n_samples, n_timepoints, n_variables)
    real(c_double), intent(out), target :: acceleration(n_samples, n_timepoints, n_variables)

    ! 1D work vectors (length depends on n_timepoints)
    ! Caller must allocate:
    !   factor_velocity, dependent_velocity, velocity_contributions : length (n_timepoints-1) if n_timepoints>1
    !   factor_acceleration, dependent_acceleration, acceleration_contributions : length (n_timepoints-2) if n_timepoints>2
    real(c_double), intent(inout), target :: factor_velocity(*)
    real(c_double), intent(inout), target :: dependent_velocity(*)
    real(c_double), intent(inout), target :: velocity_contributions(*)

    real(c_double), intent(inout), target :: factor_acceleration(*)
    real(c_double), intent(inout), target :: dependent_acceleration(*)
    real(c_double), intent(inout), target :: acceleration_contributions(*)

    real(c_double), intent(out), target :: C_velocity(n_samples, n_variables, n_variables)
    real(c_double), intent(out), target :: velocity_contribution_series(n_samples, n_variables, n_variables, n_timepoints)
    real(c_double), intent(out), target :: C_acceleration(n_samples, n_variables, n_variables)
    real(c_double), intent(out), target :: acceleration_contribution_series(n_samples, n_variables, n_variables, n_timepoints)

   
    integer(int32) :: ns, nt, nv, mode_int 

    ! ---- Null checks ----
    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_samples)
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(n_variables)
    M_CHECK_NON_NULL(mode)
    M_CHECK_NON_NULL(trajectories)
    M_CHECK_NON_NULL(velocity)
    M_CHECK_NON_NULL(acceleration)
    M_CHECK_NON_NULL(factor_velocity)
    M_CHECK_NON_NULL(dependent_velocity)
    M_CHECK_NON_NULL(velocity_contributions)
    M_CHECK_NON_NULL(factor_acceleration)
    M_CHECK_NON_NULL(dependent_acceleration)
    M_CHECK_NON_NULL(acceleration_contributions)
    M_CHECK_NON_NULL(C_velocity)
    M_CHECK_NON_NULL(velocity_contribution_series)
    M_CHECK_NON_NULL(C_acceleration)
    M_CHECK_NON_NULL(acceleration_contribution_series)

    ! Cast C ints to int32 used internally
    ns     = int(n_samples,  int32)
    nt     = int(n_timepoints,int32)
    nv     = int(n_variables,int32)
    mode_int = int(mode,      int32)

    call compute_velocity_acceleration_contributions(trajectories, ns, nt, nv, mode_int, &
        velocity, acceleration, &
        factor_velocity, dependent_velocity, velocity_contributions, &
        factor_acceleration, dependent_acceleration, acceleration_contributions, &
        C_velocity, velocity_contribution_series, &
        C_acceleration, acceleration_contribution_series, ierr)

end subroutine tox_compute_velocity_acceleration_contributions_c

!> C wrapper for compute_velocity_acceleration_contributions_alloc 
subroutine tox_compute_velocity_acceleration_contributions_alloc_c(trajectories, n_samples, n_timepoints, n_variables, mode, &
    C_velocity, velocity_contribution_series, &
    C_acceleration, acceleration_contribution_series, ierr) &
    bind(C, name="tox_compute_velocity_acceleration_contributions_alloc")

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding,  only: c_int, c_double
    use tox_errors, only: is_err
    use tox_trajectory_contribution_analysis, only: compute_velocity_acceleration_contributions_alloc
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in),  target :: n_samples
    integer(c_int), intent(in),  target :: n_timepoints
    integer(c_int), intent(in),  target :: n_variables
    integer(c_int), intent(in),  target :: mode

    real(c_double), intent(in),  target :: trajectories(n_samples, n_timepoints, n_variables)
    real(c_double), intent(out), target :: C_velocity(n_samples, n_variables, n_variables)
    real(c_double), intent(out), target :: velocity_contribution_series(n_samples, n_variables, n_variables, n_timepoints)
    real(c_double), intent(out), target :: C_acceleration(n_samples, n_variables, n_variables)
    real(c_double), intent(out), target :: acceleration_contribution_series(n_samples, n_variables, n_variables, n_timepoints)

    integer(c_int), intent(out), target :: ierr

    integer(int32) :: ns, nt, nv, mode_int

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_samples)
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(n_variables)
    M_CHECK_NON_NULL(mode)
    M_CHECK_NON_NULL(trajectories)
    M_CHECK_NON_NULL(C_velocity)
    M_CHECK_NON_NULL(velocity_contribution_series)
    M_CHECK_NON_NULL(C_acceleration)
    M_CHECK_NON_NULL(acceleration_contribution_series)

    ns     = int(n_samples,  int32)
    nt     = int(n_timepoints,int32)
    nv     = int(n_variables,int32)
    mode_int = int(mode,      int32)

    call compute_velocity_acceleration_contributions_alloc( &
        trajectories, ns, nt, nv, mode_int, &
        C_velocity, velocity_contribution_series, &
        C_acceleration, acceleration_contribution_series, ierr)

end subroutine tox_compute_velocity_acceleration_contributions_alloc_c

