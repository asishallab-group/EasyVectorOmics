#include "macros.h"

module tox_trajectory_contribution_analysis
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use f42_utils, only: is_close
    use tox_errors, only: set_ok, set_err, is_err, ERR_IDX_OUT_OF_BOUNDS, ERR_DIVISION_BY_ZERO, &
                          ERR_INVALID_INPUT, ERR_NAN_INF, ERR_ALLOC_FAIL, validate_dimension_size, &
                          validate_in_range_real
    implicit none

    integer(int32), parameter :: MODE_NORMAL = 1
    integer(int32), parameter :: MODE_RAP    = 2

    ! Baseline computation modes
    integer(int32), parameter :: BASELINE_RAW  = 1
    integer(int32), parameter :: BASELINE_MIN  = 2
    integer(int32), parameter :: BASELINE_MEAN = 3
contains

    !> Calculates trajectory contribution using cosine similarity
    !! $$ contribution = \frac{factor \cdot dependent}{\left \Vert factor \right \Vert \cdot \left \Vert dependent \right \Vert} $$
    !! There are two modes:
    !!
    !! 1. **Normal Mode**: meant for normal vectors, it returns the cosine similarity
    !! 2. **RAP Mode**: meant for RAP projected vectors, it focusses only on angular similarity and returns the `arccos` of the cosine similarity
    pure subroutine trajectory_contribution(factor, dependent, n_timepoints, mode, contribution, ierr)
        integer(int32), intent(in) :: n_timepoints
            !! Number of timepoints included in `factor` and `dependent`
        real(real64), dimension(n_timepoints), intent(in) :: factor
            !! Vector of independent variable that contributes to `dependent`
        real(real64), dimension(n_timepoints), intent(in) :: dependent
            !! Vector of dependent variable
        integer(int32), intent(in) :: mode
            !! | Mode   | Value |
            !! |--------|-------|
            !! | Normal | 1     |
            !! | RAP    | 2     |
        real(real64), intent(out) :: contribution
            !! Contribution value
        integer(int32), intent(out) :: ierr
            !! Error code

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
            ! -1 <= x <= 1 ; the closer to one, the more directionally similar
            contribution = dot_prod / magnitude
        case (MODE_RAP)
            ! the smaller the angle, the more aligned
            contribution = acos(dot_prod / magnitude)
        case default
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end select
    end subroutine trajectory_contribution

    !> Calculates spike contributions using cosine similarity for specific indices
    !! $$ contributions(i) =
    !!         \frac{factor(i) \cdot dependent(i)}
    !!              {\left \Vert factor \right \Vert \cdot \left \Vert dependent \right \Vert}
    !!         \; | \; 1 \le i \le n\_timepoints $$
    !! There are two modes:
    !!
    !! 1. **Normal Mode**: meant for normal vectors, it returns the cosine similarity
    !! 2. **RAP Mode**: meant for RAP projected vectors, it focusses only on angular similarity and returns the `arccos` of the cosine similarity
    pure subroutine spike_contribution(factor, dependent, n_timepoints, mode, contribution, ierr)
        integer(int32), intent(in) :: n_timepoints
            !! Number of timepoints included in `factor` and `dependent`
        real(real64), dimension(n_timepoints), intent(in) :: factor
            !! Vector of independent variable that contributes to `dependent`
        real(real64), dimension(n_timepoints), intent(in) :: dependent
            !! Vector of independent variable that contributes  `dependent`
        integer(int32), intent(in) :: mode
            !! | Mode   | Value |
            !! |--------|-------|
            !! | Normal | 1     |
            !! | RAP    | 2     |
        real(real64), dimension(n_timepoints), intent(out) :: contribution
            !! Contribution value
        integer(int32), intent(out) :: ierr
            !! Error code

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
                ! -1 <= x <= 1 ; the closer to one, the more directionally similar for specific timepoint
                contribution(i_timepoint) = (factor(i_timepoint) * dependent(i_timepoint)) / magnitude
            end do
        case (MODE_RAP)
            do i_timepoint = 1, n_timepoints
                ! the smaller the angle, the more aligned for specific timepoint
                contribution(i_timepoint) = acos((factor(i_timepoint) * dependent(i_timepoint)) / magnitude)
            end do
        case default
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end select
    end subroutine spike_contribution

    !> Calculates both spike and integrated contributions for a given factor against a dependent variable
    pure subroutine calc_contributions(trajectories, n_factors, n_samples, n_timepoints, i_factor, dependent_idx, mode, spike_contribs, integrated_contribs, temp_factor_vector, temp_dependent_vector, ierr)
        integer(int32), intent(in) :: n_factors
            !! Number of factors in `trajectories`
        integer(int32), intent(in) :: n_samples
            !! Number of samples in `trajectories`
        integer(int32), intent(in) :: n_timepoints
            !! Number of timepoints in `trajectories`
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
            !! Trajectories of expression data
        integer(int32), intent(in) :: i_factor
            !! Factor (independent variable) index the contributions should be calculated for
        integer(int32), intent(in) :: dependent_idx
            !! Index of the dependent variable the contributions should be calculated for
        integer(int32), intent(in) :: mode
            !! | Mode   | Value |
            !! |--------|-------|
            !! | Normal | 1     |
            !! | RAP    | 2     |
        real(real64), dimension(n_samples), intent(out) :: integrated_contribs
            !! Output array to hold the contribution per sample
        real(real64), dimension(n_timepoints, n_samples), intent(out) :: spike_contribs
            !! Output array to hold the timpoint-wise contributions per sample
        real(real64), dimension(n_timepoints), intent(out) :: temp_factor_vector
            !! Work array to hold the vector of `i_factor` per sample
        real(real64), dimension(n_timepoints), intent(out) :: temp_dependent_vector
            !! Work array to hold the vector of `dependent_idx` per sample
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i_sample

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

    !> Calculates both spike and integrated contributions for a given factor against a dependent variable
    pure subroutine calc_contributions_alloc(trajectories, n_factors, n_samples, n_timepoints, i_factor, dependent_idx, mode, spike_contribs, integrated_contribs, ierr)
        integer(int32), intent(in) :: n_factors
            !! Number of factors in `trajectories`
        integer(int32), intent(in) :: n_samples
            !! Number of samples in `trajectories`
        integer(int32), intent(in) :: n_timepoints
            !! Number of timepoints in `trajectories`
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
            !! Trajectories of expression data
        integer(int32), intent(in) :: i_factor
            !! Factor (independent variable) index the contributions should be calculated for
        integer(int32), intent(in) :: dependent_idx
            !! Index of the dependent variable the contributions should be calculated for
        integer(int32), intent(in) :: mode
            !! | Mode   | Value |
            !! |--------|-------|
            !! | Normal | 1     |
            !! | RAP    | 2     |
        real(real64), dimension(n_samples), intent(out) :: integrated_contribs
            !! Output array to hold the contribution per sample
        real(real64), dimension(n_timepoints, n_samples), intent(out) :: spike_contribs
            !! Output array to hold the timpoint-wise contributions per sample
        integer(int32), intent(out) :: ierr
            !! Error code

        real(real64), dimension(:), allocatable :: temp_factor_vector, temp_dependent_vector

        call set_ok(ierr)
        allocate(temp_factor_vector(n_timepoints), temp_dependent_vector(n_timepoints), stat=ierr)
        if (is_err(ierr)) then
            call set_err(ierr, ERR_ALLOC_FAIL)
            return
        end if

        call calc_contributions(trajectories, n_factors, n_samples, n_timepoints, i_factor, dependent_idx, mode, spike_contribs, integrated_contribs, temp_factor_vector, temp_dependent_vector, ierr)
    end subroutine calc_contributions_alloc

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



    !> Accessor routine to extract sample data for a specific factor and timepoint
    pure subroutine get_vec_across_samples(trajectories, n_factors, n_samples, n_timepoints, i_factor, i_timepoint, result_vector, ierr)
        integer(int32), intent(in) :: n_factors
            !! Number of factors in `trajectories`
        integer(int32), intent(in) :: n_samples
            !! Number of samples in `trajectories`
        integer(int32), intent(in) :: n_timepoints
            !! Number of timepoints in `trajectories`
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
            !! Trajectories of expression data
        integer(int32), intent(in) :: i_factor
            !! Factor index the retrieved samples should be related to
        integer(int32), intent(in) :: i_timepoint
            !! Timepoint index the retrieved samples should be related to
        real(real64), dimension(n_samples), intent(out) :: result_vector
            !! Result vector, equivalent to `trajectories(i_factor, :, i_timepoint)`
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i_sample

        call set_ok(ierr)

        call validate_dimension_size(n_samples, ierr)
        call validate_dimension_size(n_factors, ierr)
        call validate_dimension_size(n_timepoints, ierr)
        if(is_err(ierr)) return

        if (i_factor < 1 .or. i_timepoint < 1 .or. i_factor > n_factors .or. i_timepoint > n_timepoints) then
            call set_err(ierr, ERR_IDX_OUT_OF_BOUNDS)
            return
        end if

        do i_sample = 1, n_samples
            result_vector(i_sample) = trajectories(i_factor, i_sample, i_timepoint)
        end do
    end subroutine get_vec_across_samples

    !> Accessor routine to extract timepoint data for a specific factor and sample
    pure subroutine get_vec_across_timepoints(trajectories, n_factors, n_samples, n_timepoints, i_factor, i_sample, result_vector, ierr)
        integer(int32), intent(in) :: n_factors
            !! Number of factors in `trajectories`
        integer(int32), intent(in) :: n_samples
            !! Number of samples in `trajectories`
        integer(int32), intent(in) :: n_timepoints
            !! Number of timepoints in `trajectories`
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
            !! Trajectories of expression data
        integer(int32), intent(in) :: i_factor
            !! Factor index the retrieved timepoints should be related to
        integer(int32), intent(in) :: i_sample
            !! Timepoint index the retrieved timepoints should be related to
        real(real64), dimension(n_timepoints), intent(out) :: result_vector
            !! Result vector, equivalent to `trajectories(i_factor, i_sample, :)`
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i_timepoint

        call set_ok(ierr)

        call validate_dimension_size(n_samples, ierr)
        call validate_dimension_size(n_factors, ierr)
        call validate_dimension_size(n_timepoints, ierr)
        if(is_err(ierr)) return

        if (i_factor < 1 .or. i_sample < 1 .or. i_factor > n_factors .or. i_sample > n_samples) then
            call set_err(ierr, ERR_IDX_OUT_OF_BOUNDS)
            return
        end if

        do i_timepoint = 1, n_timepoints
            result_vector(i_timepoint) = trajectories(i_factor, i_sample, i_timepoint)
        end do
    end subroutine get_vec_across_timepoints

    !> Process trajectories with one percentile per timepoint for spike contributions.
    !! This subroutine uses multiple allocations. If you want to avoid this, you need to do the full pipeline manually step-by-step.
    subroutine process_trajectories_alloc(trajectories, n_factors, n_samples, n_timepoints, factor_mask, factor_mask_count, dependent_idx, mode, percentile, &
                                integrated_contribs, spike_contribs, thresholds_integrated_contrib, &
                                outliers_integrated_contrib, thresholds_spike_contrib, &
                                outliers_spike_contrib, ierr)
        use tox_trajectory_contribution_analysis_outlier_detection
        integer(int32), intent(in) :: n_factors
        !! Number of factors
        integer(int32), intent(in) :: n_samples
        !! Number of samples
        integer(int32), intent(in) :: n_timepoints
        !! Number of timepoints
        integer(int32), intent(in) :: factor_mask_count
        !! Number of true values in the factor mask
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
        !! Trajectories array
        logical, dimension(n_factors), intent(in) :: factor_mask
        !! Mask for what factors should be included
        integer(int32), intent(in) :: dependent_idx
        !! Index of the dependend factor
        integer(int32), intent(in) :: mode
        !! Mode for contribution calculation
        real(real64), intent(in) :: percentile
        !! Percentile (0-100)
        
        real(real64), dimension(n_samples, factor_mask_count), intent(out) :: integrated_contribs
        !! Integrated contributions
        real(real64), dimension(n_timepoints, n_samples, factor_mask_count), intent(out) :: spike_contribs
        !! Spike contributions
        real(real64), dimension(factor_mask_count), intent(out) :: thresholds_integrated_contrib
        !! Thresholds of integrated contributions based on the provided percentile
        real(real64), dimension(n_timepoints, factor_mask_count), intent(out) :: thresholds_spike_contrib
        !! Thresholds of spike contributions based on the provided percentile
        logical, dimension(n_samples, factor_mask_count), intent(out) :: outliers_integrated_contrib
        !! Detected outliers of the integrated contributions
        logical, dimension(n_timepoints, n_samples, factor_mask_count), intent(out) :: outliers_spike_contrib
        !! Detected outliers of the spike contributions

        integer(int32), intent(out) :: ierr
        integer(int32) :: i_factor, i_included_factor
        
        call set_ok(ierr)
        
        ! Input validation
        call validate_dimension_size(n_samples, ierr)
        call validate_dimension_size(n_factors, ierr)
        call validate_dimension_size(n_timepoints, ierr)
        if(is_err(ierr)) return
        
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
                
                ! Calculate contributions
                call calc_contributions_alloc(trajectories, n_factors, n_samples, n_timepoints, &
                                    i_factor, dependent_idx, mode, &
                                    spike_contribs(:, :, i_included_factor), &
                                    integrated_contribs(:, i_included_factor), &
                                    ierr)
                if (is_err(ierr)) exit
                
                ! Calculate integrated threshold and outliers
                call calc_integrated_threshold_alloc(integrated_contribs(:, i_included_factor), &
                                            n_samples, percentile, thresholds_integrated_contrib(i_included_factor), &
                                            ierr)
                if (is_err(ierr)) exit

                call detect_outliers_integrated(integrated_contribs(:, i_included_factor), n_samples, &
                                            thresholds_integrated_contrib(i_included_factor), &
                                            outliers_integrated_contrib(:, i_included_factor), ierr)
                if (is_err(ierr)) exit

                ! Calculate spike thresholds and outliers (one per timepoint)
                call calc_spike_thresholds_alloc(spike_contribs(:, :, i_included_factor), &
                                        n_timepoints, n_samples, percentile, & 
                                        thresholds_spike_contrib(:, i_included_factor), ierr)
                if (is_err(ierr)) exit
                
                call detect_outliers_spike(spike_contribs(:, :, i_included_factor), n_timepoints, &
                                        n_samples, thresholds_spike_contrib(:, i_included_factor), &
                                        outliers_spike_contrib(:, :, i_included_factor), ierr)
                if (is_err(ierr)) exit
            end if
        end do
    end subroutine process_trajectories_alloc

    !> Process trajectories with global percentile for spike contributions (flattened)
    !! This subroutine uses multiple allocations. If you want to avoid this, you need to do the full pipeline manually step-by-step.
    subroutine process_trajectories_flat_alloc(trajectories, n_factors, n_samples, n_timepoints, factor_mask, factor_mask_count, dependent_idx, mode, percentile, &
                                    integrated_contribs, spike_contribs, thresholds_integrated_contrib, &
                                    outliers_integrated_contrib, thresholds_spike_contrib, &
                                    outliers_spike_contrib, ierr)
        use tox_trajectory_contribution_analysis_outlier_detection
        integer(int32), intent(in) :: n_factors
        !! Number of factors
        integer(int32), intent(in) :: n_samples
        !! Number of samples
        integer(int32), intent(in) :: n_timepoints
        !! Number of timepoints
        integer(int32), intent(in) :: factor_mask_count
        !! Number of true values in the factor mask
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
        !! Trajectories array
        logical, dimension(n_factors), intent(in) :: factor_mask
        !! Logical mask for factors
        integer(int32), intent(in) :: dependent_idx
        !! Index of the dependent factor
        integer(int32), intent(in) :: mode
        !! Mode to use for contribution calculation
        real(real64), intent(in) :: percentile
        !! Percentile value (0.0 - 100.0)
        
        real(real64), dimension(n_samples, factor_mask_count), intent(out) :: integrated_contribs
        !! Integrated contributions
        real(real64), dimension(n_timepoints, n_samples, factor_mask_count), intent(out) :: spike_contribs
        !! Spike contributions
        real(real64), dimension(factor_mask_count), intent(out) :: thresholds_integrated_contrib
        !! Thresholds of integrated contributions based on provided percentile
        real(real64), dimension(factor_mask_count), intent(out) :: thresholds_spike_contrib
        !! Thresholds of spike contributions based on provided percentile
        logical, dimension(n_samples, factor_mask_count), intent(out) :: outliers_integrated_contrib
        !! Outliers of integrated contributions based on provided percentile
        logical, dimension(n_timepoints, n_samples, factor_mask_count), intent(out) :: outliers_spike_contrib
        !! Outliers of spike contributions based on provided percentile
        integer(int32), intent(out) :: ierr
        !! Error code
        
        integer(int32) :: i_factor, i_included_factor
        
        call set_ok(ierr)
        
        ! Input validation
        call validate_dimension_size(n_samples, ierr)
        call validate_dimension_size(n_factors, ierr)
        call validate_dimension_size(n_timepoints, ierr)
        if(is_err(ierr)) return
        
        if (dependent_idx < 1 .or. dependent_idx > n_factors) then
            call set_err(ierr, ERR_IDX_OUT_OF_BOUNDS)
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
                
                ! Calculate contributions
                call calc_contributions_alloc(trajectories, n_factors, n_samples, n_timepoints, &
                                    i_factor, dependent_idx, mode, &
                                    spike_contribs(:, :, i_included_factor), &
                                    integrated_contribs(:, i_included_factor), ierr)
                if (is_err(ierr)) exit
                
                ! Calculate integrated threshold and outliers
                call calc_integrated_threshold_alloc(integrated_contribs(:, i_included_factor), n_samples, &
                                            percentile, thresholds_integrated_contrib(i_included_factor),ierr)
                if (is_err(ierr)) exit
                
                call detect_outliers_integrated(integrated_contribs(:, i_included_factor), n_samples, &
                                            thresholds_integrated_contrib(i_included_factor), &
                                            outliers_integrated_contrib(:, i_included_factor), ierr)
                if (is_err(ierr)) exit
                
                ! Use flattened spike contributions with integrated routines
                call calc_integrated_threshold_alloc(spike_contribs(:, :, i_included_factor), n_timepoints * n_samples, &
                                            percentile, thresholds_spike_contrib(i_included_factor), ierr)
                if (is_err(ierr)) exit
                
                call detect_outliers_integrated(spike_contribs(:, :, i_included_factor), n_timepoints * n_samples, &
                                            thresholds_spike_contrib(i_included_factor), &
                                            outliers_spike_contrib(:, :, i_included_factor), ierr)
                if (is_err(ierr)) exit
            end if
        end do
        
    end subroutine process_trajectories_flat_alloc
    !> Compute time-resolved contributions and total contribution from factor and dependent variable time series.
    subroutine compute_contributions( factor_vector, dependent_vector, n_timepoints, mode, &
        contributions, total_contribution, ierr )

    use tox_errors, only : set_ok, set_err, is_err, ERR_ALLOC_FAIL, validate_dimension_size
    implicit none

    integer(int32), intent(in) :: n_timepoints
    integer(int32), intent(in) :: mode
    real(real64), intent(in)  :: factor_vector(n_timepoints)
    real(real64), intent(in)  :: dependent_vector(n_timepoints)
    real(real64), intent(out) :: contributions(n_timepoints)
    real(real64), intent(out) :: total_contribution
    integer(int32), intent(out) :: ierr

    integer(int32) :: time_index

    !> Baseline value for the factor variable.
    real(real64) :: baseline_factor
    !> Baseline value for the dependent variable.
    real(real64) :: baseline_dependent


    call set_ok(ierr)

  
    ! Ensure the dependent vector has the same length.
    if (size(dependent_vector) /= n_timepoints .or. size(factor_vector) /= n_timepoints) then
        ! There is no dedicated "dimension mismatch" error code in this context,
        ! so we reuse ERR_ALLOC_FAIL as a generic failure.
        call set_err(ierr, ERR_ALLOC_FAIL)
        return
    end if

    ! Use the common dimension validation utility.
    call validate_dimension_size(n_timepoints, ierr)
    if (is_err(ierr)) return


    baseline_factor    = factor_vector(1)
    baseline_dependent = dependent_vector(1)

    contributions      = 0.0_real64
    total_contribution = 0.0_real64

    do time_index = 1, n_timepoints

        contributions(time_index) = ( factor_vector(time_index)    - baseline_factor ) &
                                  * ( dependent_vector(time_index) - baseline_dependent )

        total_contribution = total_contribution + contributions(time_index)

    end do

end subroutine compute_contributions


!> Compute velocity trajectories from position trajectories

subroutine compute_velocity_trajectories(trajectories, velocity, &
                                         n_samples, n_timepoints, n_variables, ierr)
    use tox_errors, only : set_ok, validate_dimension_size
    implicit none

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
    use tox_errors, only : set_ok, validate_dimension_size
    implicit none

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
                acceleration(sample, t, var) = velocity(sample, t, var) - &
                                               2.0_real64 * velocity(sample, t - 1, var) + &
                                               velocity(sample, t - 2, var)
            end do
        end do
    end do
end subroutine compute_acceleration_from_velocity

!> Compute velocity and acceleration contributions for all variable pairs in the trajectories
subroutine compute_velocity_acceleration_contributions( trajectories, n_samples, n_timepoints, n_variables, mode, &
        C_velocity, velocity_contribution_series, &
        C_acceleration, acceleration_contribution_series, ierr)

    use tox_errors, only : set_ok, set_err, is_err, ERR_ALLOC_FAIL, validate_dimension_size
    implicit none

    integer(int32), intent(in)  :: n_samples, n_timepoints, n_variables
    integer(int32), intent(in)  :: mode
    real(real64),   intent(in)  :: trajectories(n_samples, n_timepoints, n_variables)
    real(real64),   intent(out) :: C_velocity(n_samples, n_variables, n_variables)
    real(real64),   intent(out) :: velocity_contribution_series(n_samples, n_variables, n_variables, n_timepoints)
    real(real64),   intent(out) :: C_acceleration(n_samples, n_variables, n_variables)
    real(real64),   intent(out) :: acceleration_contribution_series(n_samples, n_variables, n_variables, n_timepoints)
    integer(int32), intent(out) :: ierr

    real(real64), allocatable :: velocity(:,:,:)
    real(real64), allocatable :: acceleration(:,:,:)
    real(real64), allocatable :: factor_velocity(:)
    real(real64), allocatable :: dependent_velocity(:)
    real(real64), allocatable :: velocity_contributions(:)
    real(real64), allocatable :: factor_acceleration(:)
    real(real64), allocatable :: dependent_acceleration(:)
    real(real64), allocatable :: acceleration_contributions(:)
    real(real64) :: total_velocity_contribution
    real(real64) :: total_acceleration_contribution
    integer(int32) :: sample, factor_index, dependent_index, time_index
    integer :: stat_alloc

    call set_ok(ierr)

    call validate_dimension_size(n_samples, ierr)
    call validate_dimension_size(n_timepoints, ierr)
    call validate_dimension_size(n_variables, ierr)
    if (is_err(ierr)) return

    C_velocity                     = 0.0_real64
    velocity_contribution_series   = 0.0_real64
    C_acceleration                 = 0.0_real64
    acceleration_contribution_series = 0.0_real64

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

    call compute_velocity_trajectories(trajectories, velocity, n_samples, n_timepoints, n_variables, ierr)
    if (is_err(ierr)) then
        deallocate(acceleration, velocity)
        return
    end if

    call compute_acceleration_from_velocity(velocity, acceleration, n_samples, n_timepoints, n_variables, ierr)
    if (is_err(ierr)) then
        deallocate(acceleration, velocity)
        return
    end if

    do sample = 1, n_samples
        do factor_index = 1, n_variables
            do dependent_index = 1, n_variables

                if (n_timepoints > 1) then
                    allocate(factor_velocity(n_timepoints-1), dependent_velocity(n_timepoints-1), &
                             velocity_contributions(n_timepoints-1), stat=stat_alloc)
                    if (stat_alloc /= 0) then
                        call set_err(ierr, ERR_ALLOC_FAIL)
                        deallocate(acceleration, velocity)
                        return
                    end if

                    do time_index = 2, n_timepoints
                        factor_velocity(time_index-1)    = velocity(sample, time_index, factor_index)
                        dependent_velocity(time_index-1) = velocity(sample, time_index, dependent_index)
                    end do

                    call compute_contributions(factor_velocity, dependent_velocity, size(factor_velocity), mode, &
                                               velocity_contributions, total_velocity_contribution, ierr)
                    if (is_err(ierr)) then
                        deallocate(factor_velocity, dependent_velocity, velocity_contributions)
                        deallocate(acceleration, velocity)
                        return
                    end if

                    C_velocity(sample, factor_index, dependent_index) = total_velocity_contribution
                    do time_index = 2, n_timepoints
                        velocity_contribution_series(sample, factor_index, dependent_index, time_index) = &
                            velocity_contributions(time_index-1)
                    end do

                    deallocate(factor_velocity, dependent_velocity, velocity_contributions)
                end if

                if (n_timepoints > 2) then
                    allocate(factor_acceleration(n_timepoints-2), dependent_acceleration(n_timepoints-2), &
                             acceleration_contributions(n_timepoints-2), stat=stat_alloc)
                    if (stat_alloc /= 0) then
                        call set_err(ierr, ERR_ALLOC_FAIL)
                        deallocate(acceleration, velocity)
                        return
                    end if

                    do time_index = 3, n_timepoints
                        factor_acceleration(time_index-2)    = acceleration(sample, time_index, factor_index)
                        dependent_acceleration(time_index-2) = acceleration(sample, time_index, dependent_index)
                    end do

                    call compute_contributions(factor_acceleration, dependent_acceleration, size(factor_acceleration), mode, &
                                               acceleration_contributions, total_acceleration_contribution, ierr)
                    if (is_err(ierr)) then
                        deallocate(factor_acceleration, dependent_acceleration, acceleration_contributions)
                        deallocate(acceleration, velocity)
                        return
                    end if

                    C_acceleration(sample, factor_index, dependent_index) = total_acceleration_contribution
                    do time_index = 3, n_timepoints
                        acceleration_contribution_series(sample, factor_index, dependent_index, time_index) = &
                            acceleration_contributions(time_index-2)
                    end do

                    deallocate(factor_acceleration, dependent_acceleration, acceleration_contributions)
                end if

            end do
        end do
    end do

    deallocate(acceleration, velocity)
end subroutine compute_velocity_acceleration_contributions


end module tox_trajectory_contribution_analysis

pure subroutine trajectory_contribution_c(factor, dependent, n_timepoints, mode, contribution, ierr) bind(C, name="trajectory_contribution_c")
    use tox_trajectory_contribution_analysis, only: trajectory_contribution
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_timepoints
        !! Number of timepoints included in `factor` and `dependent`
    real(c_double), dimension(n_timepoints), intent(in), target :: factor
        !! Vector of independent variable that contributes to `dependent`
    real(c_double), dimension(n_timepoints), intent(in), target :: dependent
        !! Vector of dependent variable
    integer(c_int), intent(in), target :: mode
        !! | Mode   | Value |
        !! |--------|-------|
        !! | Normal | 1     |
        !! | RAP    | 2     |
    real(c_double), intent(out), target :: contribution
        !! Contribution value
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(factor)
    M_CHECK_NON_NULL(dependent)
    M_CHECK_NON_NULL(mode)
    M_CHECK_NON_NULL(contribution)

    call trajectory_contribution(factor, dependent, n_timepoints, mode, contribution, ierr)
end subroutine trajectory_contribution_c

!> Calculates spike contributions using cosine similarity for specific indices
!! $$ contributions(i) =
!!         \frac{factor(i) \cdot dependent(i)}
!!              {\left \Vert factor \right \Vert \cdot \left \Vert dependent \right \Vert}
!!         \; | \; 1 \le i \le n\_timepoints $$
!! There are two modes:
!!
!! 1. **Normal Mode**: meant for normal vectors, it returns the cosine similarity
!! 2. **RAP Mode**: meant for RAP projected vectors, it focusses only on angular similarity and returns the `arccos` of the cosine similarity
pure subroutine spike_contribution_c(factor, dependent, n_timepoints, mode, contribution, ierr) bind(C, name="spike_contribution_c")
    use tox_trajectory_contribution_analysis, only: spike_contribution
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_timepoints
        !! Number of timepoints included in `factor` and `dependent`
    real(c_double), dimension(n_timepoints), intent(in), target :: factor
        !! Vector of independent variable that contributes to `dependent`
    real(c_double), dimension(n_timepoints), intent(in), target :: dependent
        !! Vector of independent variable that contributes  `dependent`
    integer(c_int), intent(in), target :: mode
        !! | Mode   | Value |
        !! |--------|-------|
        !! | Normal | 1     |
        !! | RAP    | 2     |
    real(c_double), dimension(n_timepoints), intent(out), target :: contribution
        !! Contribution value
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(factor)
    M_CHECK_NON_NULL(dependent)
    M_CHECK_NON_NULL(mode)
    M_CHECK_NON_NULL(contribution)

    call spike_contribution(factor, dependent, n_timepoints, mode, contribution, ierr)
end subroutine spike_contribution_c

!> Calculates both spike and integrated contributions for a given factor against a dependent variable
pure subroutine calc_contributions_expert_c(trajectories, n_factors, n_samples, n_timepoints, i_factor, dependent_idx, mode, spike_contribs, integrated_contribs, temp_factor_vector, temp_dependent_vector, ierr) bind(C, name="calc_contributions_expert_c")
    use tox_trajectory_contribution_analysis, only: calc_contributions
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_factors
        !! Number of factors in `trajectories`
    integer(c_int), intent(in), target :: n_samples
        !! Number of samples in `trajectories`
    integer(c_int), intent(in), target :: n_timepoints
        !! Number of timepoints in `trajectories`
    real(c_double), dimension(n_factors, n_samples, n_timepoints), intent(in), target :: trajectories
        !! Trajectories of expression data
    integer(c_int), intent(in), target :: i_factor
        !! Factor (independent variable) index the contributions should be calculated for, starting from 0
    integer(c_int), intent(in), target :: dependent_idx
        !! Index of the dependent variable the contributions should be calculated for, starting from 0
    integer(c_int), intent(in), target :: mode
        !! | Mode   | Value |
        !! |--------|-------|
        !! | Normal | 1     |
        !! | RAP    | 2     |
    real(c_double), dimension(n_samples), intent(out), target :: integrated_contribs
        !! Output array to hold the contribution per sample
    real(c_double), dimension(n_timepoints, n_samples), intent(out), target :: spike_contribs
        !! Output array to hold the timpoint-wise contributions per sample
    real(c_double), dimension(n_timepoints), intent(out), target :: temp_factor_vector
        !! Work array to hold the vector of `i_factor` per sample
    real(c_double), dimension(n_timepoints), intent(out), target :: temp_dependent_vector
        !! Work array to hold the vector of `dependent_idx` per sample
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_factors)
    M_CHECK_NON_NULL(n_samples)
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(trajectories)
    M_CHECK_NON_NULL(i_factor)
    M_CHECK_NON_NULL(dependent_idx)
    M_CHECK_NON_NULL(mode)
    M_CHECK_NON_NULL(integrated_contribs)
    M_CHECK_NON_NULL(spike_contribs)
    M_CHECK_NON_NULL(temp_factor_vector)
    M_CHECK_NON_NULL(temp_dependent_vector)

    call calc_contributions(trajectories, n_factors, n_samples, n_timepoints, i_factor + 1, dependent_idx + 1, mode, spike_contribs, integrated_contribs, temp_factor_vector, temp_dependent_vector, ierr)
end subroutine calc_contributions_expert_c

!> Calculates both spike and integrated contributions for a given factor against a dependent variable
pure subroutine calc_contributions_c(trajectories, n_factors, n_samples, n_timepoints, i_factor, dependent_idx, mode, spike_contribs, integrated_contribs, ierr) bind(C, name="calc_contributions_c")
    use tox_trajectory_contribution_analysis, only: calc_contributions_alloc
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_factors
        !! Number of factors in `trajectories`
    integer(c_int), intent(in), target :: n_samples
        !! Number of samples in `trajectories`
    integer(c_int), intent(in), target :: n_timepoints
        !! Number of timepoints in `trajectories`
    real(c_double), dimension(n_factors, n_samples, n_timepoints), intent(in), target :: trajectories
        !! Trajectories of expression data
    integer(c_int), intent(in), target :: i_factor
        !! Factor (independent variable) index the contributions should be calculated for, starting from 0
    integer(c_int), intent(in), target :: dependent_idx
        !! Index of the dependent variable the contributions should be calculated for, starting from 0
    integer(c_int), intent(in), target :: mode
        !! | Mode   | Value |
        !! |--------|-------|
        !! | Normal | 1     |
        !! | RAP    | 2     |
    real(c_double), dimension(n_samples), intent(out), target :: integrated_contribs
        !! Output array to hold the contribution per sample
    real(c_double), dimension(n_timepoints, n_samples), intent(out), target :: spike_contribs
        !! Output array to hold the timpoint-wise contributions per sample
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_factors)
    M_CHECK_NON_NULL(n_samples)
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(trajectories)
    M_CHECK_NON_NULL(i_factor)
    M_CHECK_NON_NULL(dependent_idx)
    M_CHECK_NON_NULL(mode)
    M_CHECK_NON_NULL(integrated_contribs)
    M_CHECK_NON_NULL(spike_contribs)

    call calc_contributions_alloc(trajectories, n_factors, n_samples, n_timepoints, i_factor + 1, dependent_idx + 1, mode, spike_contribs, integrated_contribs, ierr)
end subroutine calc_contributions_c

!> C wrapper for process_trajectories_alloc
!> @Note Indices are 1 based, make sure all variables align with this
subroutine process_trajectories_C(trajectories, n_factors, n_samples, n_timepoints, n_processed, &
                                      factor_mask_int, dependent_idx, mode, percentile, &
                                      integrated_contribs, spike_contribs, &
                                      thresholds_integrated_contrib, outliers_integrated_contrib_int, &
                                      thresholds_spike_contrib, outliers_spike_contrib_int, ierr) &
                                      bind(C, name="process_trajectories_C")   

    use tox_conversions, only: logical_as_c_int, c_int_as_logical
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use tox_errors, only: set_ok, set_err, is_err, ERR_ALLOC_FAIL, validate_dimension_size
    use tox_trajectory_contribution_analysis, only: process_trajectories_alloc
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_factors
        !! Number of factors
    integer(c_int), intent(in), target :: n_samples
        !! Number of samples
    integer(c_int), intent(in), target :: n_timepoints
        !! Number of timepoints
    integer(c_int), intent(in), target :: n_processed
        !! Number of processed factors
    real(c_double), intent(in), target :: trajectories(n_factors, n_samples, n_timepoints)
        !! Trajectories array
    integer(c_int), intent(in), target :: factor_mask_int(n_factors)
        !! Mask for factors as C integers (0=false, 1=true)
    integer(c_int), intent(in), target :: dependent_idx
        !! Index of the dependent factor (1-based from C)
    integer(c_int), intent(in), target :: mode
        !! Mode for contribution calculation
    real(c_double), intent(in), target :: percentile
        !! Percentile (0-100)
    
    real(c_double), intent(out), target :: integrated_contribs(n_samples, n_processed)
        !! Integrated contributions
    real(c_double), intent(out), target :: spike_contribs(n_timepoints, n_samples, n_processed)
        !! Spike contributions
    real(c_double), intent(out), target :: thresholds_integrated_contrib(n_processed)
        !! Thresholds for integrated contributions
    integer(c_int), intent(out), target :: outliers_integrated_contrib_int(n_samples, n_processed)
        !! Outliers for integrated contributions as C integers
    real(c_double), intent(out), target :: thresholds_spike_contrib(n_timepoints, n_processed)
        !! Thresholds for spike contributions
    integer(c_int), intent(out), target :: outliers_spike_contrib_int(n_timepoints, n_samples, n_processed)
        !! Outliers for spike contributions as C integers
    integer(c_int), intent(out), target :: ierr
        !! Error code
    
    logical, allocatable :: factor_mask(:)
    logical, allocatable :: outliers_integrated_contrib(:,:)
    logical, allocatable :: outliers_spike_contrib(:,:,:)


    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_factors)
    M_CHECK_NON_NULL(n_samples)
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(n_processed)
    M_CHECK_NON_NULL(trajectories)
    M_CHECK_NON_NULL(factor_mask_int)
    M_CHECK_NON_NULL(dependent_idx)
    M_CHECK_NON_NULL(mode)
    M_CHECK_NON_NULL(percentile)
    M_CHECK_NON_NULL(integrated_contribs)
    M_CHECK_NON_NULL(spike_contribs)
    M_CHECK_NON_NULL(thresholds_integrated_contrib)
    M_CHECK_NON_NULL(outliers_integrated_contrib_int)
    M_CHECK_NON_NULL(thresholds_spike_contrib)
    M_CHECK_NON_NULL(outliers_spike_contrib_int)    

    ! Initialize error
    call set_ok(ierr)
    
    ! Check for valid dimensions
    call validate_dimension_size(n_samples, ierr)
    call validate_dimension_size(n_factors, ierr)
    call validate_dimension_size(n_timepoints, ierr)
    if(is_err(ierr)) return
    
    ! Convert C integer mask to Fortran logical
    allocate(factor_mask(n_factors), stat=ierr)
    if (is_err(ierr)) then
        call set_err(ierr, ERR_ALLOC_FAIL)
        return
    end if

    call c_int_as_logical(factor_mask_int, factor_mask)
    
    ! Allocate temporary arrays for logical outputs
    allocate(outliers_integrated_contrib(n_samples, n_processed), stat=ierr)
    if (is_err(ierr)) then
        call set_err(ierr, ERR_ALLOC_FAIL)
        deallocate(factor_mask)
        return
    end if
    
    allocate(outliers_spike_contrib(n_timepoints, n_samples, n_processed), stat=ierr)
    if (is_err(ierr)) then
        call set_err(ierr, ERR_ALLOC_FAIL)
        deallocate(factor_mask, outliers_integrated_contrib)
        return
    end if
    
    ! Call Fortran subroutine
    call process_trajectories_alloc(trajectories, n_factors, n_samples, n_timepoints, &
                                  factor_mask, n_processed, dependent_idx, mode, percentile, &
                                  integrated_contribs, spike_contribs, &
                                  thresholds_integrated_contrib, outliers_integrated_contrib, &
                                  thresholds_spike_contrib, outliers_spike_contrib, ierr)
    
    if (.not. is_err(ierr)) then
        ! Convert Fortran logical outputs to C integers
        call logical_as_c_int(outliers_integrated_contrib, outliers_integrated_contrib_int)
        call logical_as_c_int(outliers_spike_contrib, outliers_spike_contrib_int)
    end if
    
end subroutine process_trajectories_C

!> C wrapper for process_trajectories_flat_alloc
!> @Note Indices are 1 based, make sure all variables align with this
subroutine process_trajectories_flat_C(trajectories, n_factors, n_samples, n_timepoints, n_processed, &
                                           factor_mask_int, dependent_idx, mode, percentile, &
                                           integrated_contribs, spike_contribs, &
                                           thresholds_integrated_contrib, outliers_integrated_contrib_int, &
                                           thresholds_spike_contrib, outliers_spike_contrib_int, ierr) &
                                           bind(C, name="process_trajectories_flat_C")
    use tox_conversions, only: logical_as_c_int, c_int_as_logical
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use tox_errors, only: set_ok, set_err, is_err, ERR_ALLOC_FAIL, validate_dimension_size
    use tox_trajectory_contribution_analysis, only: process_trajectories_flat_alloc
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_factors
        !! Number of factors
    integer(c_int), intent(in), target :: n_samples
        !! Number of samples
    integer(c_int), intent(in), target :: n_timepoints
        !! Number of timepoints
    integer(c_int), intent(in), target :: n_processed
        !! Number of processed factors
    real(c_double), intent(in), target :: trajectories(n_factors, n_samples, n_timepoints)
        !! Trajectories array
    integer(c_int), intent(in), target :: factor_mask_int(n_factors)
        !! Mask for factors as C integers (0=false, 1=true)
    integer(c_int), intent(in), target :: dependent_idx
        !! Index of the dependent factor (0-based from C)
    integer(c_int), intent(in), target :: mode
        !! Mode for contribution calculation
    real(c_double), intent(in), target :: percentile
        !! Percentile (0-100)
    
    real(c_double), intent(out), target :: integrated_contribs(n_samples, n_processed)
        !! Integrated contributions
    real(c_double), intent(out), target :: spike_contribs(n_timepoints, n_samples, n_processed)
        !! Spike contributions
    real(c_double), intent(out), target :: thresholds_integrated_contrib(n_processed)
        !! Thresholds for integrated contributions
    integer(c_int), intent(out), target :: outliers_integrated_contrib_int(n_samples, n_processed)
        !! Outliers for integrated contributions as C integers
    real(c_double), intent(out), target :: thresholds_spike_contrib(n_processed)
        !! Thresholds for spike contributions (1D for flat version)
    integer(c_int), intent(out), target :: outliers_spike_contrib_int(n_timepoints, n_samples, n_processed)
        !! Outliers for spike contributions as C integers
    integer(c_int), intent(out), target :: ierr
        !! Error code
    
    logical, allocatable :: factor_mask(:)
    logical, allocatable :: outliers_integrated_contrib(:,:)
    logical, allocatable :: outliers_spike_contrib(:,:,:)

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_factors)
    M_CHECK_NON_NULL(n_samples)
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(n_processed)
    M_CHECK_NON_NULL(trajectories)
    M_CHECK_NON_NULL(factor_mask_int)
    M_CHECK_NON_NULL(dependent_idx)
    M_CHECK_NON_NULL(mode)
    M_CHECK_NON_NULL(percentile)
    M_CHECK_NON_NULL(integrated_contribs)
    M_CHECK_NON_NULL(spike_contribs)
    M_CHECK_NON_NULL(thresholds_integrated_contrib)
    M_CHECK_NON_NULL(outliers_integrated_contrib_int)
    M_CHECK_NON_NULL(thresholds_spike_contrib)
    M_CHECK_NON_NULL(outliers_spike_contrib_int)

    ! Initialize error
    call set_ok(ierr)
    
    ! Check for valid dimensions
    call validate_dimension_size(n_samples, ierr)
    call validate_dimension_size(n_factors, ierr)
    call validate_dimension_size(n_timepoints, ierr)
    if(is_err(ierr)) return
    
    ! Convert C integer mask to Fortran logical
    allocate(factor_mask(n_factors), stat=ierr)
    if (is_err(ierr)) then
        call set_err(ierr, ERR_ALLOC_FAIL)
        return
    end if
    
    call c_int_as_logical(factor_mask_int, factor_mask)
    
    ! Allocate temporary arrays for logical outputs
    allocate(outliers_integrated_contrib(n_samples, n_processed), stat=ierr)
    if (is_err(ierr)) then
        call set_err(ierr, ERR_ALLOC_FAIL)
        deallocate(factor_mask)
        return
    end if
    
    allocate(outliers_spike_contrib(n_timepoints, n_samples, n_processed), stat=ierr)
    if (is_err(ierr)) then
        call set_err(ierr, ERR_ALLOC_FAIL)
        deallocate(factor_mask, outliers_integrated_contrib)
        return
    end if
    
    ! Call Fortran subroutine
    call process_trajectories_flat_alloc(trajectories, n_factors, n_samples, n_timepoints, &
                                       factor_mask, n_processed, dependent_idx, mode, percentile, &
                                       integrated_contribs, spike_contribs, &
                                       thresholds_integrated_contrib, outliers_integrated_contrib, &
                                       thresholds_spike_contrib, outliers_spike_contrib, ierr)
    
    if (.not. is_err(ierr)) then
        ! Convert Fortran logical outputs to C integers
        call logical_as_c_int(outliers_integrated_contrib, outliers_integrated_contrib_int)
        call logical_as_c_int(outliers_spike_contrib, outliers_spike_contrib_int)
    end if    
end subroutine process_trajectories_flat_C
!> C wrapper for compute_baselines_factor_dependent
subroutine compute_baselines_factor_dependent_c(factor, dependent, n_timepoints, mode, &
                                               factor_baseline, dependent_baseline, ierr) &
    bind(C, name="tox_compute_baselines_factor_dependent")

    use tox_trajectory_contribution_analysis, only: compute_baselines_factor_dependent
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    M_USE_NULL_VALIDATION
    implicit none


    integer(c_int), intent(in),  target :: n_timepoints
    real(c_double), intent(in),  target :: factor(n_timepoints)
    real(c_double), intent(in),  target :: dependent(n_timepoints)
    integer(c_int), intent(in),  target :: mode
    real(c_double), intent(out), target :: factor_baseline
    real(c_double), intent(out), target :: dependent_baseline
    integer(c_int), intent(out), target :: ierr

    !! Null-pointer validation 
    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(factor)
    M_CHECK_NON_NULL(dependent)
    M_CHECK_NON_NULL(mode)
    M_CHECK_NON_NULL(factor_baseline)
    M_CHECK_NON_NULL(dependent_baseline)

    !! Call Fortran subroutine directly (no type conversion needed)
    call compute_baselines_factor_dependent(n_timepoints, factor, dependent, mode, &
                                           factor_baseline, dependent_baseline, ierr)

end subroutine compute_baselines_factor_dependent_c
!> C wrapper for compute_velocity_trajectories
subroutine tox_compute_velocity_trajectories_c(trajectories, n_samples, n_timepoints, n_variables, &
                                               velocity, ierr) &
    bind(C, name="tox_compute_velocity_trajectories")
    use tox_trajectory_contribution_analysis, only : compute_velocity_trajectories
    use, intrinsic :: iso_c_binding, only : c_int, c_double
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in),  target :: n_samples
    integer(c_int), intent(in),  target :: n_timepoints
    integer(c_int), intent(in),  target :: n_variables
    real(c_double), intent(in),  target :: trajectories(n_samples, n_timepoints, n_variables)
    real(c_double), intent(out), target :: velocity(n_samples, n_timepoints, n_variables)
    integer(c_int), intent(out), target :: ierr

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_samples)
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(n_variables)
    M_CHECK_NON_NULL(trajectories)
    M_CHECK_NON_NULL(velocity)

    call compute_velocity_trajectories(trajectories, velocity, &
                                       n_samples, n_timepoints, n_variables, ierr)
end subroutine tox_compute_velocity_trajectories_c

!> C wrapper for compute_acceleration_from_velocity
subroutine tox_compute_acceleration_from_velocity_c(velocity, n_samples, n_timepoints, n_variables, &
                                                    acceleration, ierr) &
    bind(C, name="tox_compute_acceleration_from_velocity")
    use tox_trajectory_contribution_analysis, only : compute_acceleration_from_velocity
    use, intrinsic :: iso_c_binding, only : c_int, c_double
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in),  target :: n_samples
    integer(c_int), intent(in),  target :: n_timepoints
    integer(c_int), intent(in),  target :: n_variables
    real(c_double), intent(in),  target :: velocity(n_samples, n_timepoints, n_variables)
    real(c_double), intent(out), target :: acceleration(n_samples, n_timepoints, n_variables)
    integer(c_int), intent(out), target :: ierr

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_samples)
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(n_variables)
    M_CHECK_NON_NULL(velocity)
    M_CHECK_NON_NULL(acceleration)

    call compute_acceleration_from_velocity(velocity, acceleration, &
                                            n_samples, n_timepoints, n_variables, ierr)
end subroutine tox_compute_acceleration_from_velocity_c
!> C wrapper for compute_contributions
subroutine tox_compute_contributions_c(factor_vector, dependent_vector, n_timepoints, mode, &
                                       contributions, total_contribution, ierr) &
    bind(C, name="tox_compute_contributions")
    use tox_trajectory_contribution_analysis, only : compute_contributions
    use, intrinsic :: iso_c_binding, only : c_int, c_double
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in),  target :: n_timepoints
    real(c_double), intent(in),  target :: factor_vector(n_timepoints)
    real(c_double), intent(in),  target :: dependent_vector(n_timepoints)
    integer(c_int), intent(in),  target :: mode
    real(c_double), intent(out), target :: contributions(n_timepoints)
    real(c_double), intent(out), target :: total_contribution
    integer(c_int), intent(out), target :: ierr

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(factor_vector)
    M_CHECK_NON_NULL(dependent_vector)
    M_CHECK_NON_NULL(mode)
    M_CHECK_NON_NULL(contributions)
    M_CHECK_NON_NULL(total_contribution)

    call compute_contributions(factor_vector, dependent_vector, n_timepoints, mode, &
                               contributions, total_contribution, ierr)
end subroutine tox_compute_contributions_c

!> C wrapper for compute_velocity_and_acceleration_contributions
subroutine tox_compute_velocity_acceleration_contributions_c( &
    trajectories, n_samples, n_timepoints, n_variables, mode, &
    C_velocity, velocity_contribution_series, &
    C_acceleration, acceleration_contribution_series, ierr) &
    bind(C, name="tox_compute_velocity_acceleration_contributions")

    use tox_trajectory_contribution_analysis, only: compute_velocity_acceleration_contributions
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in),  target :: n_samples
    integer(c_int), intent(in),  target :: n_timepoints
    integer(c_int), intent(in),  target :: n_variables
    real(c_double), intent(in),  target :: trajectories(n_samples, n_timepoints, n_variables)
    integer(c_int), intent(in),  target :: mode
    real(c_double), intent(out), target :: C_velocity(n_samples, n_variables, n_variables)
    real(c_double), intent(out), target :: velocity_contribution_series(n_samples, n_variables, n_variables, n_timepoints)
    real(c_double), intent(out), target :: C_acceleration(n_samples, n_variables, n_variables)
    real(c_double), intent(out), target :: acceleration_contribution_series(n_samples, n_variables, n_variables, n_timepoints)
    integer(c_int), intent(out), target :: ierr

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_samples)
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(n_variables)
    M_CHECK_NON_NULL(trajectories)
    M_CHECK_NON_NULL(mode)
    M_CHECK_NON_NULL(C_velocity)
    M_CHECK_NON_NULL(velocity_contribution_series)
    M_CHECK_NON_NULL(C_acceleration)
    M_CHECK_NON_NULL(acceleration_contribution_series)

    call compute_velocity_acceleration_contributions(trajectories, n_samples, n_timepoints, n_variables, mode, &
        C_velocity, velocity_contribution_series, C_acceleration, acceleration_contribution_series, ierr)
end subroutine tox_compute_velocity_acceleration_contributions_c