module tox_trajectory_normalization
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_errors
    implicit none
    
    private
    public :: normalize_variable_timeseries, &
              normalize_single_trajectory, &
              normalize_all_trajectories
    
contains
    
    !> Normalize a single variable across time using min-max scaling
    subroutine normalize_variable_timeseries(v, v_norm, n_points, ierr)
        integer(int32), intent(in) :: n_points
        !! Vector length
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
        
        ! Check for empty array
        call validate_dimension_size(n_points, ierr)
        if (is_err(ierr)) return

        epsilon_val = 1.0e-12_real64
        
        ! Find min and max values
        min_val = minval(v)
        max_val = maxval(v)
        
        ! Calculate denominator with epsilon for numerical stability
        denominator = max_val - min_val + epsilon_val
        
        ! Check for division by zero (should be prevented by epsilon)
        if (abs(denominator) < tiny(1.0_real64)) then
            call set_err(ierr, ERR_DIVISION_BY_ZERO)
            v_norm = 0.0_real64
            return
        end if
        
        ! Apply min-max normalization
        do i = 1, n_points
            v_norm(i) = (v(i) - min_val) / denominator
        end do
        
    end subroutine normalize_variable_timeseries
    
    !> Normalize all variables in a single trajectory independently
    subroutine normalize_single_trajectory(trajectory, trajectory_norm, n_factors, n_samples, ierr)
        integer(int32), intent(in) :: n_factors
        !! Number of factors
        integer(int32), intent(in) :: n_samples
        !! Number of samples
        real(real64), intent(in) :: trajectory(n_factors, n_samples)
        !! Original trajectory
        real(real64), intent(out) :: trajectory_norm(n_factors, n_samples)
        !! Normalized trajectory
        integer(int32), intent(out) :: ierr
        !! Error code
        
        integer(int32) :: j
        
        ! Initialize
        call set_ok(ierr)
        
        call validate_dimension_size(n_factors, ierr)
        call validate_dimension_size(n_samples, ierr)
        if (is_err(ierr)) return
        
        ! Normalize each variable independently
        do j = 1, n_samples
            call normalize_variable_timeseries(trajectory(:, j), trajectory_norm(:, j), n_factors, ierr)
            
            if(is_err(ierr)) return
        end do
        
    end subroutine normalize_single_trajectory
    
    !===================================================================
    ! Subroutine: normalize_all_trajectories
    ! Purpose: Normalize all trajectories across multiple entities
    ! Inputs:
    !   trajectories(:,:,:) - Original 3D array (time × variables × entities)
    !   eps - Small constant to prevent division by zero
    ! Outputs:
    !   trajectories_norm(:,:,:) - Normalized trajectories
    !   err - Error code (0 = success)
    !===================================================================
    subroutine normalize_all_trajectories(trajectories, trajectories_norm, n_factors, n_samples, n_timepoints, ierr)
        integer(int32), intent(in) :: n_factors
        !! Number of factors
        integer(int32), intent(in) :: n_samples
        !! Number of samples
        integer(int32), intent(in) :: n_timepoints
        !! Number of time points
        real(real64), intent(in) :: trajectories(n_factors, n_samples, n_timepoints)
        !! Original trajectories
        real(real64), intent(out) :: trajectories_norm(n_factors, n_samples, n_timepoints)
        !! Normalized trajectories
        integer(int32), intent(out) :: ierr
        !! Error code
        
        integer(int32) :: entity        
        ! Initialize
        call set_ok(ierr)

        call validate_dimension_size(n_factors, ierr)
        call validate_dimension_size(n_samples, ierr)
        call validate_dimension_size(n_timepoints, ierr)
        if (is_err(ierr)) return
        
        ! Normalize each entity independently
        do entity = 1, n_timepoints
            call normalize_single_trajectory(trajectories(:,:,entity), trajectories_norm(:,:,entity), n_factors, n_samples, ierr)
            
            if(is_err(ierr)) return
        end do
        
    end subroutine normalize_all_trajectories

end module tox_trajectory_normalization

subroutine normalize_variable_timeseries_C(v, v_norm, n_points, ierr) bind(C, name="normalize_variable_timeseries_C")
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use safeguard
    use tox_trajectory_normalization, only: normalize_variable_timeseries
    implicit none

    integer(c_int) :: n_points
    !! Vector length
    real(c_double), dimension(n_points), intent(in) :: v
    !! Original time series
    real(c_double), dimension(n_points), intent(out) :: v_norm
    !! Normalized time series        
    integer(c_int), intent(out) :: ierr
    !! Error code

    call normalize_variable_timeseries(v, v_norm, n_points, ierr)

end subroutine normalize_variable_timeseries_C

subroutine normalize_single_trajectory_C(trajectory, trajectory_norm, n_factors, n_samples, ierr) bind(C, name="normalize_single_trajectory_C")
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use safeguard
    use tox_trajectory_normalization, only: normalize_single_trajectory
    implicit none

    integer(c_int) :: n_factors
    !! Number of factors
    integer(c_int) :: n_samples
    !! Number of samples
    real(c_double), dimension(n_factors, n_samples), intent(in) :: trajectory
    !! Original trajectory
    real(c_double), dimension(n_factors, n_samples), intent(out) :: trajectory_norm
    !! Normalized trajectory
    integer(c_int), intent(out) :: ierr
    !! Error code

    call normalize_single_trajectory(trajectory, trajectory_norm, n_factors, n_samples, ierr)

end subroutine normalize_single_trajectory_C

subroutine normalize_all_trajectories_C(trajectories, trajectories_norm, n_factors, n_samples, n_timepoints, ierr) bind(C, name="normalize_all_trajectories_C")
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use safeguard
    use tox_trajectory_normalization, only: normalize_all_trajectories
    implicit none

    integer(c_int) :: n_factors
    !! Number of factors
    integer(c_int) :: n_samples
    !! Number of samples
    integer(c_int) :: n_timepoints
    !! Number of time points
    real(c_double), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
    !! Original trajectories
    real(c_double), dimension(n_factors, n_samples, n_timepoints), intent(out) :: trajectories_norm
    !! Normalized trajectories
    integer(c_int), intent(out) :: ierr
    !! Error code

    call normalize_all_trajectories(trajectories, trajectories_norm, n_factors, n_samples, n_timepoints, ierr)

end subroutine normalize_all_trajectories_C