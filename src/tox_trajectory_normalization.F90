#include "macros.h"
module tox_trajectory_normalization
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_errors, only: ERR_NAN_INF, set_ok, set_err, is_err, validate_dimension_size, validate_all_in_range_real
    use safeguard
    implicit none
    
    private
    public :: normalize_variable_timeseries, &
              normalize_single_trajectory, &
              normalize_all_trajectories
    
contains
    
    !> Normalize a single variable across time using min-max scaling
    pure subroutine normalize_variable_timeseries(v, v_norm, n_points, ierr)
        integer(int32), intent(in) :: n_points
        !! Vector length (number of time points)
        real(real64), intent(in) :: v(n_points)
        !! Original time series
        real(real64), intent(out) :: v_norm(n_points)
        !! Normalized time series        
        integer(int32), intent(out) :: ierr
        !! Error code
        
        real(real64) :: min_val, max_val, denominator, epsilon_val
        integer(int32) :: i
        
        ! Initialize
        call set_ok(ierr)

        call validate_all_in_range_real(v, n_points, ierr)
        if(is_err(ierr)) return
        
        ! Check for empty array
        call validate_dimension_size(n_points, ierr)
        if (is_err(ierr)) return

        epsilon_val = 1.0e-12_real64
        
        ! Find min and max values
        min_val = minval(v)
        max_val = maxval(v)
        
        ! Calculate denominator
        denominator = max_val - min_val
        
        ! Check for division by zero
        if (abs(denominator) < epsilon_val) then
            denominator = denominator + epsilon_val
        end if
        
        ! Apply min-max normalization
        do i = 1, n_points
            v_norm(i) = (v(i) - min_val) / denominator
        end do
        
    end subroutine normalize_variable_timeseries
    
    !> Normalize all factors in a single trajectory independently across time
    !! Input: trajectory(n_factors, n_timepoints) for ONE sample/entity
    pure subroutine normalize_single_trajectory(trajectory, trajectory_norm, n_factors, n_timepoints, ierr)
        integer(int32), intent(in) :: n_factors
        !! Number of factors/variables
        integer(int32), intent(in) :: n_timepoints
        !! Number of time points
        real(real64), intent(in) :: trajectory(n_factors, n_timepoints)
        !! Original trajectory for one sample
        real(real64), intent(out) :: trajectory_norm(n_factors, n_timepoints)
        !! Normalized trajectory for one sample
        integer(int32), intent(out) :: ierr
        !! Error code
        
        integer(int32) :: i_factor

        real(real64), dimension(n_timepoints) :: temp_series, temp_series_norm
        
        ! Initialize
        call set_ok(ierr)
        
        call validate_dimension_size(n_factors, ierr)
        if (is_err(ierr)) return
        
        ! Normalize each factor independently across time
        do i_factor = 1, n_factors
            temp_series = trajectory(i_factor, :)
            call normalize_variable_timeseries( &
                temp_series, &           ! Time series for this factor
                temp_series_norm, &      ! Normalized time series
                n_timepoints, ierr)
            trajectory_norm(i_factor, :) = temp_series_norm
            
            if (is_err(ierr)) return
        end do
        
    end subroutine normalize_single_trajectory
    
    !> Normalize all trajectories across multiple entities
    !! Input: trajectories(n_factors, n_samples, n_timepoints)
    !! Normalizes each factor independently across time for each sample
    pure subroutine normalize_all_trajectories(trajectories, trajectories_norm, &
                                          n_factors, n_samples, n_timepoints, ierr)
        integer(int32), intent(in) :: n_factors
        !! Number of factors
        integer(int32), intent(in) :: n_samples
        !! Number of samples/entities
        integer(int32), intent(in) :: n_timepoints
        !! Number of time points
        real(real64), intent(in) :: trajectories(n_factors, n_samples, n_timepoints)
        !! Original trajectories
        real(real64), intent(out) :: trajectories_norm(n_factors, n_samples, n_timepoints)
        !! Normalized trajectories
        integer(int32), intent(out) :: ierr
        !! Error code
        
        integer(int32) :: i_sample, i_factor

        real(real64), dimension(n_timepoints) :: temp_series, temp_series_norm
        
        call set_ok(ierr)

        call validate_dimension_size(n_samples, ierr)
        if (is_err(ierr)) return
        
        ! Normalize each sample/entity independently
        do i_sample = 1, n_samples
            do i_factor = 1, n_factors
                temp_series = trajectories(i_factor, i_sample, :)

                call normalize_variable_timeseries(temp_series, temp_series_norm, n_timepoints, ierr)
                
                if (is_err(ierr)) return
                trajectories_norm(i_factor, i_sample, :) = temp_series_norm
            end do
        end do
        
    end subroutine normalize_all_trajectories

end module tox_trajectory_normalization

pure subroutine normalize_variable_timeseries_C(v, v_norm, n_points, ierr) bind(C, name="normalize_variable_timeseries_C")
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use tox_trajectory_normalization, only: normalize_variable_timeseries
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_points
    !! Vector length
    real(c_double), dimension(n_points), intent(in), target :: v
    !! Original time series
    real(c_double), dimension(n_points), intent(out), target :: v_norm
    !! Normalized time series        
    integer(c_int), intent(out), target :: ierr
    !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_points)
    M_CHECK_NON_NULL(v)
    M_CHECK_NON_NULL(v_norm)

    call normalize_variable_timeseries(v, v_norm, n_points, ierr)

end subroutine normalize_variable_timeseries_C

pure subroutine normalize_single_trajectory_C(trajectory, trajectory_norm, n_factors, n_timepoints, ierr) bind(C, name="normalize_single_trajectory_C")
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use tox_trajectory_normalization, only: normalize_single_trajectory
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_factors
    !! Number of factors
    integer(c_int), intent(in), target :: n_timepoints
    !! Number of time points
    real(c_double), dimension(n_factors, n_timepoints), intent(in), target :: trajectory
    !! Original trajectory for one sample
    real(c_double), dimension(n_factors, n_timepoints), intent(out), target :: trajectory_norm
    !! Normalized trajectory for one sample
    integer(c_int), intent(out), target :: ierr
    !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_factors)
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(trajectory)
    M_CHECK_NON_NULL(trajectory_norm)

    call normalize_single_trajectory(trajectory, trajectory_norm, n_factors, n_timepoints, ierr)

end subroutine normalize_single_trajectory_C

pure subroutine normalize_all_trajectories_C(trajectories, trajectories_norm, n_factors, n_samples, n_timepoints, ierr) bind(C, name="normalize_all_trajectories_C")
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use tox_trajectory_normalization, only: normalize_all_trajectories
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_factors
    !! Number of factors
    integer(c_int), intent(in), target :: n_samples
    !! Number of samples
    integer(c_int), intent(in), target :: n_timepoints
    !! Number of time points
    real(c_double), dimension(n_factors, n_samples, n_timepoints), intent(in), target :: trajectories
    !! Original trajectories
    real(c_double), dimension(n_factors, n_samples, n_timepoints), intent(out), target :: trajectories_norm
    !! Normalized trajectories
    integer(c_int), intent(out), target :: ierr
    !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_factors)
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(n_samples)
    M_CHECK_NON_NULL(trajectories)
    M_CHECK_NON_NULL(trajectories_norm)

    call normalize_all_trajectories(trajectories, trajectories_norm, n_factors, n_samples, n_timepoints, ierr)

end subroutine normalize_all_trajectories_C