#include "macros.h"

module tox_trajectory_contribution_analysis
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, is_err, set_err, ERR_ALLOC_FAIL, ERR_INVALID_INPUT, validate_dimension_size, validate_all_in_range_real, validate_all_in_range_int, validate_in_range_real, validate_in_range_int
    use f42_utils, only: init_random, rand_range
    implicit none

    ! Baseline computation modes
    integer(int32), parameter :: BASELINE_RAW  = 1
    integer(int32), parameter :: BASELINE_MIN  = 2
    integer(int32), parameter :: BASELINE_MEAN = 3
contains

    !> Selects a random sample different to `current_sample`. For reproducibility call [[f42_utils(module):init_random(subroutine)]] beforehand.
    subroutine select_random_sample_helper(n_samples, current_sample, random_sample, ierr)
        integer(int32), intent(in) :: n_samples
            !! number of samples
        integer(int32), intent(in) :: current_sample
            !! sample that won't be selected
        integer(int32), intent(out) :: random_sample
            !! selected random sample
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_in_range_int(n_samples, ierr, min=2_int32)
        call validate_in_range_int(current_sample, ierr, min=1_int32, max=n_samples)

        if (is_err(ierr)) return

        ! pick random sample in range [1,n_samples-1], so a set without `current_sample`
        random_sample = int(rand_range(1.0_real64, real(n_samples, real64)), kind=int32)

        ! adjust picked sample for excluded `current_sample`
        if (random_sample >= current_sample) then
            random_sample = random_sample + 1
        end if
    end subroutine select_random_sample_helper

    !> For a given factor-dependent pair, this subroutine calculates the contributions by taking the same dependent but from a random different sample.
    subroutine perform_permutation_test(trajectories, n_factors, n_samples, n_timepoints, factor_idx, dependent_idx, sample_idx, mode, n_permutations, local_contributions, total_contributions, temp_factor, temp_dependent, ierr, random_seed)
        integer(int32), intent(in) :: n_factors
            !! number of factors
        integer(int32), intent(in) :: n_samples
            !! number of samples
        integer(int32), intent(in) :: n_timepoints
            !! number of timepoints
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
            !! expression vectors across different samples over time
        integer(int32), intent(in) :: factor_idx
            !! index of factor to compute the permutation contributions for
        integer(int32), intent(in) :: dependent_idx
            !! index of dependent to compute the permutation contributions for
        integer(int32), intent(in) :: sample_idx
            !! index of sample to compute the permutation contributions for
        integer(int32), intent(in) :: mode
            !! Baseline mode: 1=RAW, 2=MIN, 3=MEAN
        integer(int32), intent(in) :: n_permutations
            !! number of permutations to perform
        real(real64), dimension(n_timepoints, n_permutations), intent(out) :: local_contributions
            !! Per-timepoint contributions per permutation
        real(real64), dimension(n_permutations), intent(out) :: total_contributions
            !! Total contribution (`sum(local_contributions)`) per permutation
        real(real64), dimension(n_timepoints), intent(out) :: temp_factor
            !! Working array to hold the factor in contiguous memory
        real(real64), dimension(n_timepoints), intent(out) :: temp_dependent
            !! Working array to hold the random dependent in contiguous memory
        integer(int32), intent(out) :: ierr
            !! Error code
        integer(int32), intent(in), optional :: random_seed
            !! Seed to use for random number generation.

        integer(int32) :: random_sample, i_perm, i_timepoint

        call set_ok(ierr)

        call validate_dimension_size(n_factors, ierr)
        call validate_dimension_size(n_samples, ierr)
        call validate_dimension_size(n_timepoints, ierr)
        call validate_dimension_size(n_permutations, ierr)
        call validate_all_in_range_real(trajectories, size(trajectories, kind=int32), ierr)
        call validate_in_range_int(factor_idx, ierr, min=1, max=n_factors)
        call validate_in_range_int(dependent_idx, ierr, min=1, max=n_factors)
        call validate_in_range_int(sample_idx, ierr, min=1, max=n_samples)

        if (is_err(ierr)) return

        if (present(random_seed)) then
            call init_random(random_seed)
        end if

        do i_timepoint = 1, n_timepoints
            temp_factor(i_timepoint) = trajectories(factor_idx, sample_idx, i_timepoint)
        end do

        do i_perm = 1, n_permutations
            call select_random_sample_helper(n_samples, sample_idx, random_sample, ierr)
            if (is_err(ierr)) return

            do i_timepoint = 1, n_timepoints
                temp_dependent(i_timepoint) = trajectories(dependent_idx, random_sample, i_timepoint)
            end do

            call compute_contributions_helper(temp_factor, temp_dependent, n_timepoints, mode, local_contributions(:, i_perm), total_contributions(i_perm), ierr)
            if (is_err(ierr)) return
        end do
    end subroutine perform_permutation_test

    !> Once the permutation tests are calculated ([[tox_trajectory_contribution_analysis(module):perform_permutation_test(subroutine)]]),
    !| this subroutine calculates the p values for the contributions, i.e. how many of the permutation contributions were at least as high as the real contributions.
    pure subroutine compute_p_values(local_contributions_observed, total_contribution_observed, local_contributions_perm, total_contributions_perm, n_timepoints, n_permutations, local_p_values, total_p_value, ierr)
        integer(int32), intent(in) :: n_timepoints
            !! number of timepoints
        integer(int32), intent(in) :: n_permutations
            !! number of permutations to perform
        real(real64), dimension(n_timepoints), intent(in) :: local_contributions_observed
            !! Per-timepoint contributions for the observed factor-dependent-sample combination
        real(real64), intent(in) :: total_contribution_observed
            !! Total contribution (`sum(local_contributions)`) for the observed factor-dependent-sample combination
        real(real64), dimension(n_timepoints, n_permutations), intent(in) :: local_contributions_perm
            !! Per-timepoint contributions for the factor-dependent-random_sample combinations from [[tox_trajectory_contribution_analysis(module):perform_permutation_test(subroutine)]]
        real(real64), dimension(n_permutations), intent(in) :: total_contributions_perm
            !! Total contribution (`sum(local_contributions)`) for the factor-dependent-random_sample combinations from [[tox_trajectory_contribution_analysis(module):perform_permutation_test(subroutine)]]
        real(real64), dimension(n_timepoints), intent(out) :: local_p_values
            !! calculated p values for local contributions, like: `(local_contributions_perm >= local_contributions_observed)/n_permutations`
        real(real64), intent(out) :: total_p_value
            !! calculated p values for total contributions, like: `(total_contributions_perm >= total_contribution_observed)/n_permutations`
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i_perm, i_timepoint

        call set_ok(ierr)

        call validate_dimension_size(n_timepoints, ierr)
        call validate_dimension_size(n_permutations, ierr)
        call validate_all_in_range_real(local_contributions_observed, size(local_contributions_observed, kind=int32), ierr)
        call validate_in_range_real(total_contribution_observed, ierr)
        call validate_all_in_range_real(local_contributions_perm, size(local_contributions_perm, kind=int32), ierr)
        call validate_all_in_range_real(total_contributions_perm, size(total_contributions_perm, kind=int32), ierr)

        if (is_err(ierr)) return

        total_p_value = 0.0_real64
        local_p_values = 0.0_real64

        do i_perm = 1, n_permutations
            if (total_contributions_perm(i_perm) >= total_contribution_observed) then
                total_p_value = total_p_value + 1.0_real64
            end if

            do i_timepoint = 1, n_timepoints
                if (local_contributions_perm(i_timepoint, i_perm) >= local_contributions_observed(i_timepoint)) then
                    local_p_values(i_timepoint) = local_p_values(i_timepoint) + 1.0_real64
                end if
            end do
        end do

        do i_timepoint = 1, n_timepoints
            local_p_values(i_timepoint) = real(nint(local_p_values(i_timepoint), kind=int32), kind=real64) / real(n_permutations, kind=real64)
        end do

        total_p_value = real(nint(total_p_value, kind=int32), kind=real64) / real(n_permutations, kind=real64)
    end subroutine compute_p_values

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
            !! expression vectors across different samples over time
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
    pure subroutine get_baseline_mode(c_mode_str, mode, ierr)
        use, intrinsic :: iso_c_binding, only: c_char
        use tox_conversions, only: c_char_1d_as_string
        character(len=1, kind=c_char), dimension(4), intent(in) :: c_mode_str
            !! mode string ("min", "mean", "raw")
        integer(int32), intent(out) :: mode
            !! integer representation for the mode passed by `c_mode_str`
        integer(int32), intent(out) :: ierr
            !! Error code

        character(len=:), allocatable :: mode_str_f

        call set_ok(ierr)

        call c_char_1d_as_string(c_mode_str, mode_str_f, ierr)
        if (is_err(ierr)) return

        select case (trim(mode_str_f))
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
    pure subroutine compute_velocity_trajectories(trajectories, velocity, &
                                             n_samples, n_timepoints, n_variables, ierr)

        integer(int32), intent(in)  :: n_samples, n_timepoints, n_variables
        !! number of samples, timepoints, and variables
        integer(int32), intent(out) :: ierr
        !! Error code
        real(real64), intent(in)  :: trajectories(n_samples, n_timepoints, n_variables)
        !! input position trajectories
        real(real64), intent(out) :: velocity(n_samples, n_timepoints, n_variables)
        !! output velocity trajectories

        integer(int32) :: sample, var

        call set_ok(ierr)
        call validate_dimension_size(n_samples, ierr)
        call validate_dimension_size(n_timepoints, ierr)
        call validate_dimension_size(n_variables, ierr)
        if (is_err(ierr)) return

        velocity = 0.0_real64
        if (n_timepoints <= 1) return

        do sample = 1, n_samples
            do var = 1, n_variables
                call compute_velocity_trajectory(trajectories(sample, :, var), &
                                                velocity(sample, :, var), &
                                                n_timepoints, ierr)
                if (is_err(ierr)) return
            end do
        end do
    end subroutine compute_velocity_trajectories

    !> Compute acceleration trajectories from velocity trajectories
    pure subroutine compute_acceleration_from_velocity(velocity, acceleration, &
                                                  n_samples, n_timepoints, n_variables, ierr)

        integer(int32), intent(in)  :: n_samples, n_timepoints, n_variables
        !! number of samples, timepoints, and variables
        integer(int32), intent(out) :: ierr
        !! Error code
        real(real64), intent(in)  :: velocity(n_samples, n_timepoints, n_variables)
        !! input velocity trajectories
        real(real64), intent(out) :: acceleration(n_samples, n_timepoints, n_variables)
        !! output acceleration trajectories

        integer(int32) :: sample, var

        call set_ok(ierr)

        call validate_dimension_size(n_samples, ierr)
        call validate_dimension_size(n_timepoints, ierr)
        call validate_dimension_size(n_variables, ierr)
        if (is_err(ierr)) return

        acceleration = 0.0_real64
        if (n_timepoints <= 2) return

        do sample = 1, n_samples
            do var = 1, n_variables
                call compute_acceleration_from_velocity_trajectory(velocity(sample, :, var), &
                                                                  acceleration(sample, :, var), &
                                                                  n_timepoints, ierr)
                if (is_err(ierr)) return
            end do
        end do
    end subroutine compute_acceleration_from_velocity

    !> Compute velocity and acceleration contributions for all variable pairs in the trajectories
    pure subroutine compute_velocity_acceleration_contributions(trajectories, n_samples, n_timepoints, n_variables, mode, velocity, acceleration, &
        factor_velocity, dependent_velocity, velocity_contributions, &
        factor_acceleration, dependent_acceleration, acceleration_contributions, &
        C_velocity, velocity_contribution_series, &
        C_acceleration, acceleration_contribution_series, ierr)

        integer(int32), intent(in) :: n_samples, n_timepoints, n_variables
        !! number of samples, timepoints, and variables
        integer(int32), intent(in) :: mode
        !! Baseline mode: 1=RAW, 2=MIN, 3=MEAN
        real(real64),   intent(in) :: trajectories(n_samples, n_timepoints, n_variables)
        !! input position trajectories

        ! Workspace (preallocated by caller)
        real(real64), intent(out) :: velocity(n_samples, n_timepoints, n_variables)
        !! output velocity trajectories
        real(real64), intent(out) :: acceleration(n_samples, n_timepoints, n_variables)
        !! output acceleration trajectories

        real(real64), intent(inout) :: factor_velocity(n_timepoints-1, n_variables)
        !! velocity factor workspace
        real(real64), intent(inout) :: dependent_velocity(n_timepoints-1)
        !! velocity dependent workspace
        real(real64), intent(inout) :: velocity_contributions(n_timepoints-1)
        !! velocity contributions workspace

        real(real64), intent(inout) :: factor_acceleration(n_timepoints-2, n_variables)
        !! acceleration factor workspace
        real(real64), intent(inout) :: dependent_acceleration(n_timepoints-2)
        !! acceleration dependent workspace
        real(real64), intent(inout) :: acceleration_contributions(n_timepoints-2)
        !! acceleration contributions workspace

        ! Outputs
        real(real64), intent(out) :: C_velocity(n_samples, n_variables, n_variables)
        !! output velocity contributions
        real(real64), intent(out) :: velocity_contribution_series(n_samples, n_variables, n_variables, n_timepoints)
        !! output velocity contribution series
        real(real64), intent(out) :: C_acceleration(n_samples, n_variables, n_variables)
        !! output acceleration contributions
        real(real64), intent(out) :: acceleration_contribution_series(n_samples, n_variables, n_variables, n_timepoints)
        !! output acceleration contribution series

        integer(int32), intent(out) :: ierr
        !! Error code

        integer(int32) :: sample, factor_index, dependent_index, time_index
        !! Loop indices for samples, factors, dependents, and time points
        real(real64)   :: total_velocity_contribution, total_acceleration_contribution
        !! Total contributions for velocity and acceleration
        integer(int32) :: n_vel, n_acc
        !! Number of valid velocity and acceleration time points

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

        if (n_vel > 0) then
            do sample = 1, n_samples
                ! Extract ALL factor velocities once per sample into 2D array
                do factor_index = 1, n_variables
                    do time_index = 2, n_timepoints
                        factor_velocity(time_index-1, factor_index) = velocity(sample, time_index, factor_index)
                    end do
                end do
            
                do dependent_index = 1, n_variables
                    ! Extract dependent velocity data once per dependent
                    do time_index = 2, n_timepoints
                        dependent_velocity(time_index-1) = velocity(sample, time_index, dependent_index)
                    end do
                    
                    ! Process all factors for this dependent
                    do factor_index = 1, n_variables
                        ! Set first velocity to zero (no change from predecessor)
                        velocity_contribution_series(sample, factor_index, dependent_index, 1) = 0.0_real64

                        call compute_contributions( &
                            factor_velocity(:, factor_index), dependent_velocity, n_vel, mode, &
                            velocity_contributions, total_velocity_contribution, ierr)
                        if (is_err(ierr)) return

                        C_velocity(sample, factor_index, dependent_index) = total_velocity_contribution

                        do time_index = 2, n_timepoints
                            velocity_contribution_series(sample, factor_index, dependent_index, time_index) = &
                                velocity_contributions(time_index-1)
                        end do
                    end do
                end do
            end do
        end if

        if (n_acc > 0) then
            do sample = 1, n_samples
                ! Extract ALL factor accelerations once per sample into 2D array
                do factor_index = 1, n_variables
                    do time_index = 3, n_timepoints
                        factor_acceleration(time_index-2, factor_index) = acceleration(sample, time_index, factor_index)
                    end do
                end do
                
                ! Now iterate through dependents
                do dependent_index = 1, n_variables
                    ! Extract dependent acceleration data once per dependent
                    do time_index = 3, n_timepoints
                        dependent_acceleration(time_index-2) = acceleration(sample, time_index, dependent_index)
                    end do
                    
                    ! Process all factors for this dependent
                    do factor_index = 1, n_variables
                        ! Set first two accelerations to zero (no change for missing predecessors)
                        acceleration_contribution_series(sample, factor_index, dependent_index, 1) = 0.0_real64
                        acceleration_contribution_series(sample, factor_index, dependent_index, 2) = 0.0_real64

                        call compute_contributions( &
                            factor_acceleration(:, factor_index), dependent_acceleration, n_acc, mode, &
                            acceleration_contributions, total_acceleration_contribution, ierr)
                        if (is_err(ierr)) return

                        C_acceleration(sample, factor_index, dependent_index) = total_acceleration_contribution

                        do time_index = 3, n_timepoints
                            acceleration_contribution_series(sample, factor_index, dependent_index, time_index) = &
                                acceleration_contributions(time_index-2)
                        end do
                    end do
                end do
            end do
        end if
    end subroutine compute_velocity_acceleration_contributions

    subroutine compute_velocity_acceleration_contributions_alloc(trajectories, n_samples, n_timepoints, n_variables, mode, &
        C_velocity, velocity_contribution_series, &
        C_acceleration, acceleration_contribution_series, ierr)

        integer(int32), intent(in) :: n_samples, n_timepoints, n_variables
        !! number of samples
        !! number of timepoints
        !! number of variables
        integer(int32), intent(in) :: mode
        !! Baseline mode: 1=RAW, 2=MIN, 3=MEAN
        real(real64),   intent(in) :: trajectories(n_samples, n_timepoints, n_variables)
         !! input position trajectories

        real(real64), intent(out) :: C_velocity(n_samples, n_variables, n_variables)
        !! output velocity contributions
        real(real64), intent(out) :: velocity_contribution_series(n_samples, n_variables, n_variables, n_timepoints)
        !! output acceleration contributions
        real(real64), intent(out) :: C_acceleration(n_samples, n_variables, n_variables)
        !! output acceleration contributions
        real(real64), intent(out) :: acceleration_contribution_series(n_samples, n_variables, n_variables, n_timepoints)
        !! output acceleration contributions
        integer(int32), intent(out) :: ierr
        !! Error code

        ! Workspace (allocated once here)
        real(real64), allocatable :: velocity(:,:,:)
        !! velocity trajectory
        real(real64), allocatable :: acceleration(:,:,:)
        !! acceleration trajectory

        real(real64), allocatable :: factor_velocity(:, :), dependent_velocity(:), velocity_contributions(:)
        !! velocity contributions workspace
        real(real64), allocatable :: factor_acceleration(:, :), dependent_acceleration(:), acceleration_contributions(:)
        !! acceleration contributions workspace

        call set_ok(ierr)

        call validate_dimension_size(n_samples, ierr)
        call validate_dimension_size(n_timepoints, ierr)
        call validate_dimension_size(n_variables, ierr)
        if (is_err(ierr)) return

        ! Allocate big work arrays once
        M_ALLOCATE(velocity(n_samples, n_timepoints, n_variables))
        M_ALLOCATE(acceleration(n_samples, n_timepoints, n_variables))

        ! Allocate 1D reusable work vectors once
        if (n_timepoints > 1) then
            M_ALLOCATE(factor_velocity(n_timepoints-1, n_variables))
            M_ALLOCATE(dependent_velocity(n_timepoints-1))
            M_ALLOCATE(velocity_contributions(n_timepoints-1))
        end if

        if (n_timepoints > 2) then
            M_ALLOCATE(factor_acceleration(n_timepoints-2, n_variables))
            M_ALLOCATE(dependent_acceleration(n_timepoints-2))
            M_ALLOCATE(acceleration_contributions(n_timepoints-2))
        end if

        ! Call the SK routine (no allocation inside)
        call compute_velocity_acceleration_contributions(trajectories, n_samples, n_timepoints, n_variables, mode, &
            velocity, acceleration, &
            factor_velocity, dependent_velocity, velocity_contributions, &
            factor_acceleration, dependent_acceleration, acceleration_contributions, &
            C_velocity, velocity_contribution_series, &
            C_acceleration, acceleration_contribution_series, ierr)

        ! Explicitly free work arrays (Fortran would normally clean these on exit, but being explicit helps
        ! when interoping with ctypes and mixed runtimes).
        if (allocated(velocity))                deallocate(velocity)
        if (allocated(acceleration))            deallocate(acceleration)
        if (allocated(factor_velocity))         deallocate(factor_velocity)
        if (allocated(dependent_velocity))      deallocate(dependent_velocity)
        if (allocated(velocity_contributions))  deallocate(velocity_contributions)
        if (allocated(factor_acceleration))     deallocate(factor_acceleration)
        if (allocated(dependent_acceleration))  deallocate(dependent_acceleration)
        if (allocated(acceleration_contributions)) deallocate(acceleration_contributions)
    end subroutine compute_velocity_acceleration_contributions_alloc

    !> Compute velocity trajectory from a single position trajectory
    !!
    !! Velocity at time t is computed as: v(t) = x(t) - x(t-1)
    !! The first timepoint (t=1) has zero velocity since there is no predecessor.
    pure subroutine compute_velocity_trajectory(trajectory, velocity, n_timepoints, ierr)

        integer(int32), intent(in)  :: n_timepoints
        !! number of timepoints
        integer(int32), intent(out) :: ierr
        !! Error code
        real(real64), intent(in)  :: trajectory(n_timepoints)
        !! input position trajectory
        real(real64), intent(out) :: velocity(n_timepoints)
        !! output velocity trajectory

        integer(int32) :: t

        call set_ok(ierr)
        call validate_dimension_size(n_timepoints, ierr)
        if (is_err(ierr)) return

        velocity = 0.0_real64
        if (n_timepoints <= 1) return

        do t = 2, n_timepoints
            velocity(t) = trajectory(t) - trajectory(t - 1)
        end do
    end subroutine compute_velocity_trajectory

    !> Compute acceleration trajectory from a single velocity trajectory
    !!
    !! Acceleration at time t is computed as: a(t) = v(t) - v(t-1)
    !! The first two timepoints (t=1,2) have zero acceleration since there are insufficient predecessors.
    pure subroutine compute_acceleration_from_velocity_trajectory(velocity_traj, acceleration, &
                                                                 n_timepoints, ierr)

        integer(int32), intent(in)  :: n_timepoints
        !! number of timepoints
        integer(int32), intent(out) :: ierr
        !! Error code
        real(real64), intent(in)  :: velocity_traj(n_timepoints)
        !! velocity trajectory
        real(real64), intent(out) :: acceleration(n_timepoints)
        !! acceleration trajectory

        integer(int32) :: t

        call set_ok(ierr)
        call validate_dimension_size(n_timepoints, ierr)
        if (is_err(ierr)) return

        acceleration = 0.0_real64
        if (n_timepoints <= 2) return

        do t = 3, n_timepoints
            acceleration(t) = velocity_traj(t) - velocity_traj(t - 1)
        end do
    end subroutine compute_acceleration_from_velocity_trajectory
end module tox_trajectory_contribution_analysis

!> C-compatible wrapper for [[tox_trajectory_contribution_analysis(module):compute_all_contributions(subroutine)]]
pure subroutine compute_all_contributions_c(trajectories, n_factors, n_samples, n_timepoints, &
    factor_indices, n_selected_factors, dependent_indices, n_selected_dependents, mode, &
    local_contributions, total_contributions, temp_factors, temp_dependent, ierr) &
    bind(C, name="compute_all_contributions_c")

    use, intrinsic :: iso_fortran_env, only: int32
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char
    use tox_trajectory_contribution_analysis, only: compute_all_contributions, get_baseline_mode
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
    character(len=1, kind=c_char), dimension(*), intent(in), target :: mode
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

    call get_baseline_mode(mode, mode_int, ierr)
    if (is_err(ierr)) return
    
    call compute_all_contributions(trajectories, n_factors, n_samples, n_timepoints, &
        factor_indices, n_selected_factors, dependent_indices, n_selected_dependents, mode_int, &
        local_contributions, total_contributions, temp_factors, temp_dependent, ierr)
end subroutine compute_all_contributions_c

!> C-compatible wrapper for [[tox_trajectory_contribution_analysis(module):compute_contributions(subroutine)]]
pure subroutine compute_contributions_c(factor, dependent, n_dims, mode, local_contributions, total_contribution, ierr) &
    bind(C, name="compute_contributions_c")
    use, intrinsic :: iso_fortran_env, only: int32
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char
    use tox_trajectory_contribution_analysis, only: compute_contributions, get_baseline_mode
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
    character(len=1, kind=c_char), dimension(*), intent(in), target :: mode
        !! Baseline mode: "raw", "min", "mean"
    real(c_double), dimension(n_dims), intent(out), target :: local_contributions
        !! Per-element contributions
    real(c_double), intent(out), target :: total_contribution
        !! Total contribution
    integer(c_int), intent(out), target :: ierr
        !! Error code

    integer(int32) :: mode_int

    ! Null checks
    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_dims)
    M_CHECK_NON_NULL(factor)
    M_CHECK_NON_NULL(dependent)
    M_CHECK_NON_NULL(mode)
    M_CHECK_NON_NULL(local_contributions)
    M_CHECK_NON_NULL(total_contribution)

    call get_baseline_mode(mode, mode_int, ierr)
    if (is_err(ierr)) return

    call compute_contributions(factor, dependent, n_dims, mode_int, local_contributions, total_contribution, ierr)
end subroutine compute_contributions_c

!> C-compatible wrapper for [[tox_trajectory_contribution_analysis(module):compute_baselines_factor_dependent(subroutine)]]
pure subroutine compute_baselines_factor_dependent_c(factor, dependent, n_timepoints, mode, &
                                               factor_baseline, dependent_baseline, ierr) &
    bind(C, name="tox_compute_baselines_factor_dependent_c")

    use tox_trajectory_contribution_analysis, only: compute_baselines_factor_dependent, get_baseline_mode
    use, intrinsic :: iso_fortran_env, only: int32
    use, intrinsic :: iso_c_binding, only: c_double, c_int, c_char
    use tox_errors, only: is_err
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in),  target :: n_timepoints
        !! Number of timepoints in both factor and dependent arrays
    real(c_double), intent(in),  target :: factor(n_timepoints)
        !! Factor time series, length n_timepoints
    real(c_double), intent(in),  target :: dependent(n_timepoints)
        !! Dependent variable time series, length n_timepoints
    character(len=1, kind=c_char), dimension(*), intent(in), target :: mode
        !! Baseline mode: "raw", "min", "mean"
    real(c_double), intent(out), target :: factor_baseline
        !! Computed baseline for factor
    real(c_double), intent(out), target :: dependent_baseline
        !! Computed baseline for dependent variable
    integer(c_int), intent(out), target :: ierr
        !! Error code

    integer(int32) :: mode_int

    !! Null-pointer validation 
    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(factor)
    M_CHECK_NON_NULL(dependent)
    M_CHECK_NON_NULL(mode)
    M_CHECK_NON_NULL(factor_baseline)
    M_CHECK_NON_NULL(dependent_baseline)

    call get_baseline_mode(mode, mode_int, ierr)
    if (is_err(ierr)) return

    call compute_baselines_factor_dependent(n_timepoints, factor, dependent, mode_int, factor_baseline, dependent_baseline, ierr)
end subroutine compute_baselines_factor_dependent_c

!> C-compatible wrapper for [[tox_trajectory_contribution_analysis(module):perform_permutation_test(subroutine)]]
subroutine perform_permutation_test_c(trajectories, n_factors, n_samples, n_timepoints, &
    factor_idx, dependent_idx, sample_idx, mode, n_permutations, &
    local_contributions, total_contributions, temp_factor, temp_dependent, ierr, random_seed) &
    bind(C, name="perform_permutation_test_c")

    use, intrinsic :: iso_fortran_env, only: int32
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char
    use tox_trajectory_contribution_analysis, only: perform_permutation_test, get_baseline_mode
    M_USE_NULL_VALIDATION

    integer(c_int), intent(in), target :: n_factors
        !! number of factors
    integer(c_int), intent(in), target :: n_samples
        !! number of samples
    integer(c_int), intent(in), target :: n_timepoints
        !! number of timepoints
    real(c_double), dimension(n_factors, n_samples, n_timepoints), intent(in), target :: trajectories
        !! expression vectors across different samples over time
    integer(c_int), intent(in), target :: factor_idx
        !! index of factor to compute the permutation contributions for
    integer(c_int), intent(in), target :: dependent_idx
        !! index of dependent to compute the permutation contributions for
    integer(c_int), intent(in), target :: sample_idx
        !! index of sample to compute the permutation contributions for
    character(len=1, kind=c_char), dimension(*), intent(in), target :: mode
        !! Baseline mode: 1=RAW, 2=MIN, 3=MEAN
    integer(c_int), intent(in), target :: n_permutations
        !! number of permutations to perform
    real(c_double), dimension(n_timepoints, n_permutations), intent(out), target :: local_contributions
        !! Per-timepoint contributions per permutation
    real(c_double), dimension(n_permutations), intent(out), target :: total_contributions
        !! Total contribution (`sum(local_contributions)`) per permutation
    real(c_double), dimension(n_timepoints), intent(out), target :: temp_factor
        !! Working array to hold the factor in contiguous memory
    real(c_double), dimension(n_timepoints), intent(out), target :: temp_dependent
        !! Working array to hold the random dependent in contiguous memory
    integer(c_int), intent(out), target :: ierr
        !! Error code
    integer(c_int), intent(in), target :: random_seed
        !! Seed to use for random number generation.

    integer(int32) :: mode_int

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_samples)
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(trajectories)
    M_CHECK_NON_NULL(factor_idx)
    M_CHECK_NON_NULL(dependent_idx)
    M_CHECK_NON_NULL(sample_idx)
    M_CHECK_NON_NULL(mode)
    M_CHECK_NON_NULL(n_permutations)
    M_CHECK_NON_NULL(local_contributions)
    M_CHECK_NON_NULL(total_contributions)
    M_CHECK_NON_NULL(temp_factor)
    M_CHECK_NON_NULL(temp_dependent)
    M_CHECK_NON_NULL(random_seed)

    call get_baseline_mode(mode, mode_int, ierr)

    call perform_permutation_test(trajectories, n_factors, n_samples, n_timepoints, &
        factor_idx, dependent_idx, sample_idx, mode_int, n_permutations, &
        local_contributions, total_contributions, temp_factor, temp_dependent, ierr, random_seed)
end subroutine perform_permutation_test_c

!> C-compatible wrapper for [[tox_trajectory_contribution_analysis(module):compute_p_values(subroutine)]]
pure subroutine compute_p_values_c(local_contributions_observed, total_contribution_observed, &
    local_contributions_perm, total_contributions_perm, n_timepoints, n_permutations, &
    local_p_values, total_p_value, ierr) bind(C, name="compute_p_values_c")

    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_trajectory_contribution_analysis, only: compute_p_values
    M_USE_NULL_VALIDATION

    integer(c_int), intent(in), target :: n_timepoints
        !! number of timepoints
    integer(c_int), intent(in), target :: n_permutations
        !! number of permutations
    real(c_double), dimension(n_timepoints), intent(in), target :: local_contributions_observed
        !! observed local contributions
    real(c_double), intent(in), target :: total_contribution_observed
        !! observed total contribution
    real(c_double), dimension(n_timepoints, n_permutations), intent(in), target :: local_contributions_perm
        !! permutation local contributions
    real(c_double), dimension(n_permutations), intent(in), target :: total_contributions_perm
        !! permutation total contributions
    real(c_double), dimension(n_timepoints), intent(out), target :: local_p_values
        !! output local p-values
    real(c_double), intent(out), target :: total_p_value
        !! output total p-value
    integer(c_int), intent(out), target :: ierr
        !! error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(n_permutations)
    M_CHECK_NON_NULL(local_contributions_observed)
    M_CHECK_NON_NULL(total_contribution_observed)
    M_CHECK_NON_NULL(local_contributions_perm)
    M_CHECK_NON_NULL(total_contributions_perm)
    M_CHECK_NON_NULL(local_p_values)
    M_CHECK_NON_NULL(total_p_value)

    call compute_p_values(local_contributions_observed, total_contribution_observed, &
        local_contributions_perm, total_contributions_perm, n_timepoints, n_permutations, &
        local_p_values, total_p_value, ierr)
end subroutine compute_p_values_c


!> C wrapper for compute_velocity_trajectories
subroutine tox_compute_velocity_trajectories_c(trajectories, n_samples, n_timepoints, n_variables, &
                                               velocity, ierr) &
    bind(C, name="tox_compute_velocity_trajectories_c")
    use tox_trajectory_contribution_analysis, only : compute_velocity_trajectories
    use, intrinsic :: iso_fortran_env, only: int32
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use tox_errors, only: is_err
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in),  target :: n_samples
    !! number of samples
    integer(c_int), intent(in),  target :: n_timepoints
    !! number of timepoints
    integer(c_int), intent(in),  target :: n_variables
    !! number of variables
    real(c_double), intent(in),  target :: trajectories(n_samples, n_timepoints, n_variables)
    !! input trajectories
    real(c_double), intent(out), target :: velocity(n_samples, n_timepoints, n_variables)
    !! output velocity trajectories
    integer(c_int), intent(out), target :: ierr
    !! error code

    !! Null-pointer validation 
    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_samples)
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(n_variables)
    M_CHECK_NON_NULL(trajectories)
    M_CHECK_NON_NULL(velocity)

    call compute_velocity_trajectories(trajectories, velocity, n_samples, n_timepoints, n_variables, ierr)
end subroutine tox_compute_velocity_trajectories_c

!> C wrapper for compute_acceleration_from_velocity
subroutine tox_compute_acceleration_from_velocity_c(velocity, n_samples, n_timepoints, n_variables, &
                                                    acceleration, ierr) &
    bind(C, name="tox_compute_acceleration_from_velocity_c")

    use tox_trajectory_contribution_analysis, only : compute_acceleration_from_velocity
    use, intrinsic :: iso_fortran_env, only: int32
    use, intrinsic :: iso_c_binding, only : c_int, c_double
    use tox_errors, only: is_err

    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in),  target :: n_samples
    !! number of samples
    integer(c_int), intent(in),  target :: n_timepoints
    !! number of timepoints
    integer(c_int), intent(in),  target :: n_variables
    !! number of variables
    real(c_double), intent(in),  target :: velocity(n_samples, n_timepoints, n_variables)
    !! input velocity trajectories
    real(c_double), intent(out), target :: acceleration(n_samples, n_timepoints, n_variables)
    !! output acceleration trajectories
    integer(c_int), intent(out), target :: ierr
    !! error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_samples)
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(n_variables)
    M_CHECK_NON_NULL(velocity)
    M_CHECK_NON_NULL(acceleration)

    call compute_acceleration_from_velocity(velocity, acceleration, n_samples, n_timepoints, n_variables, ierr)
end subroutine tox_compute_acceleration_from_velocity_c

!> C wrapper for compute_velocity_and_acceleration_contributions
subroutine tox_compute_velocity_acceleration_contributions_c(trajectories, n_samples, n_timepoints, n_variables, mode, &
    velocity, acceleration, &
    factor_velocity, dependent_velocity, velocity_contributions, &
    factor_acceleration, dependent_acceleration, acceleration_contributions, &
    C_velocity, velocity_contribution_series, &
    C_acceleration, acceleration_contribution_series, ierr) &
    bind(C, name="tox_compute_velocity_acceleration_contributions_c")

    use tox_trajectory_contribution_analysis, only: compute_velocity_acceleration_contributions, get_baseline_mode
    use, intrinsic :: iso_c_binding, only : c_int, c_double, c_char
    use, intrinsic :: iso_fortran_env, only: int32
    use tox_errors, only: is_err
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in),  target :: n_samples
    !! number of samples
    integer(c_int), intent(in),  target :: n_timepoints
    !! number of timepoints
    integer(c_int), intent(in),  target :: n_variables
    !! number of variables
    character(len=1, kind=c_char), dimension(*), intent(in), target :: mode
    !! Baseline mode: "raw", "min", "mean"
    integer(c_int), intent(out), target :: ierr
    !! error code
    real(c_double), intent(in),  target :: trajectories(n_samples, n_timepoints, n_variables)
    !! input trajectories

    ! ---- Workspace (passed in from C) ----
    real(c_double), intent(out), target :: velocity(n_samples, n_timepoints, n_variables)
    !! output velocity trajectories
    real(c_double), intent(out), target :: acceleration(n_samples, n_timepoints, n_variables)
    !! output acceleration trajectories

    ! 1D work vectors (length depends on n_timepoints)
    ! Caller must allocate:
    !   factor_velocity, dependent_velocity, velocity_contributions : length (n_timepoints-1) if n_timepoints>1
    !   factor_acceleration, dependent_acceleration, acceleration_contributions : length (n_timepoints-2) if n_timepoints>2
    real(c_double), intent(inout), target :: factor_velocity(n_timepoints - 1, n_variables)
    !! factor for velocity contributions (2D workspace)
    real(c_double), intent(inout), target :: dependent_velocity(n_timepoints - 1)
    !! dependent variable for velocity contributions
    real(c_double), intent(inout), target :: velocity_contributions(n_timepoints - 1)
    !! velocity contributions

    real(c_double), intent(inout), target :: factor_acceleration(n_timepoints - 2, n_variables)
    !! factor for acceleration contributions (2D workspace)
    real(c_double), intent(inout), target :: dependent_acceleration(n_timepoints - 2)
    !! dependent variable for acceleration contributions
    real(c_double), intent(inout), target :: acceleration_contributions(n_timepoints - 2)
    !! acceleration contributions

    real(c_double), intent(out), target :: C_velocity(n_samples, n_variables, n_variables)
    !! velocity covariance matrix
    real(c_double), intent(out), target :: velocity_contribution_series(n_samples, n_variables, n_variables, n_timepoints)
    !! velocity contribution series
    real(c_double), intent(out), target :: C_acceleration(n_samples, n_variables, n_variables)
    !! acceleration covariance matrix
    real(c_double), intent(out), target :: acceleration_contribution_series(n_samples, n_variables, n_variables, n_timepoints)
    !! acceleration contribution series
    integer(int32) :: mode_int
    !! integer representation of baseline mode

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

    call get_baseline_mode(mode, mode_int, ierr)
    if (is_err(ierr)) return

    call compute_velocity_acceleration_contributions(trajectories, n_samples, n_timepoints, n_variables, mode_int, &
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
    bind(C, name="tox_compute_velocity_acceleration_contributions_alloc_c")

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding,  only: c_int, c_double, c_char
    use tox_errors, only: is_err
    use tox_trajectory_contribution_analysis, only: compute_velocity_acceleration_contributions_alloc, get_baseline_mode
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in),  target :: n_samples
    !! number of samples
    integer(c_int), intent(in),  target :: n_timepoints
    !! number of timepoints
    integer(c_int), intent(in),  target :: n_variables
    !! number of variables
    character(len=1, kind=c_char), dimension(*), intent(in), target :: mode
    !! Baseline mode: "raw", "min", "mean"

    real(c_double), intent(in),  target :: trajectories(n_samples, n_timepoints, n_variables)
    !! input trajectories
    real(c_double), intent(out), target :: C_velocity(n_samples, n_variables, n_variables)
    !! velocity covariance matrix
    real(c_double), intent(out), target :: velocity_contribution_series(n_samples, n_variables, n_variables, n_timepoints)
    !! velocity contribution series
    real(c_double), intent(out), target :: C_acceleration(n_samples, n_variables, n_variables)
    !! acceleration covariance matrix
    real(c_double), intent(out), target :: acceleration_contribution_series(n_samples, n_variables, n_variables, n_timepoints)
    !! acceleration contribution series
    integer(c_int), intent(out), target :: ierr
    !! error code

    integer(int32) :: mode_int
    ! ---- Null checks ----

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

    call get_baseline_mode(mode, mode_int, ierr)
    if (is_err(ierr)) return

    call compute_velocity_acceleration_contributions_alloc( &
        trajectories, n_samples, n_timepoints, n_variables, mode_int, &
        C_velocity, velocity_contribution_series, &
        C_acceleration, acceleration_contribution_series, ierr)

end subroutine tox_compute_velocity_acceleration_contributions_alloc_c

!> C wrapper for compute_velocity_trajectory
subroutine tox_compute_velocity_trajectory_c(trajectory, n_timepoints, velocity, ierr) &
    bind(C, name="tox_compute_velocity_trajectory_c")
    use tox_trajectory_contribution_analysis, only : compute_velocity_trajectory
    use, intrinsic :: iso_fortran_env, only: int32
    use, intrinsic :: iso_c_binding,  only: c_int, c_double
    use tox_errors, only: is_err
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in),  target :: n_timepoints
    !! number of timepoints
    real(c_double), intent(in),  target :: trajectory(n_timepoints)
    !! input position trajectory
    real(c_double), intent(out), target :: velocity(n_timepoints)
    !! output velocity trajectory
    integer(c_int), intent(out), target :: ierr
    !! error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(trajectory)
    M_CHECK_NON_NULL(velocity)

    call compute_velocity_trajectory(trajectory, velocity, n_timepoints, ierr)
end subroutine tox_compute_velocity_trajectory_c

!> C wrapper for compute_acceleration_from_velocity_trajectory
subroutine tox_compute_acceleration_from_velocity_trajectory_c(velocity, n_timepoints, acceleration, ierr) &
    bind(C, name="tox_compute_acceleration_from_velocity_trajectory_c")
    use tox_trajectory_contribution_analysis, only : compute_acceleration_from_velocity_trajectory
    use, intrinsic :: iso_fortran_env, only: int32
    use, intrinsic :: iso_c_binding,  only: c_int, c_double
    use tox_errors, only: is_err
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in),  target :: n_timepoints
    !! number of timepoints
    real(c_double), intent(in),  target :: velocity(n_timepoints)
    !! input velocity trajectory
    real(c_double), intent(out), target :: acceleration(n_timepoints)
    !! output acceleration trajectory
    integer(c_int), intent(out), target :: ierr
    !! error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(velocity)
    M_CHECK_NON_NULL(acceleration)

    call compute_acceleration_from_velocity_trajectory(velocity, acceleration, n_timepoints, ierr)
end subroutine tox_compute_acceleration_from_velocity_trajectory_c