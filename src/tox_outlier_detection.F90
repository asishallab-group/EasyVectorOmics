!> Module for trajectory explanatory contribution analysis (TCA)
!! Provides functionality for detecting outliers in trajectory and spike contributions
module tox_trajectory_contribution_analysis
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_errors, only: ERR_INVALID_INPUT, ERR_EMPTY_INPUT, &
                         ERR_ALLOC_FAIL, set_ok, set_err, is_err
    use f42_utils, only: calc_percentile, calc_percentile_alloc, sort_real
    
    implicit none
    
    private
    public :: calc_spike_thresholds, calc_integrated_threshold, &
              detect_outliers_integrated, detect_outliers_spike, &
              calc_spike_thresholds_alloc, calc_integrated_threshold_alloc
    
contains

    !> Calculate empirical thresholds for spike contributions (using pre-sorted permutations)
    !!
    !! Performance version that uses pre-computed permutation arrays to avoid internal sorting.
    !! Each row of permutation contains sorted indices for that timepoint.
    pure subroutine calc_spike_thresholds(spike_contribs, percentile_val, thresholds, permutation, ierr)
        real(real64), intent(in) :: spike_contribs(:, :)
        !! 2D array of spike contributions [n_timepoints, n_samples]
        real(real64), intent(in) :: percentile_val
        !! Percentile value for threshold (0.0-100.0)
        real(real64), intent(out) :: thresholds(:)
        !! 1D array of thresholds for each timepoint [n_timepoints]
        integer(int32), intent(in) :: permutation(:, :)
        !! Pre-computed permutation indices for each timepoint [n_timepoints, n_samples]
        integer(int32), intent(out) :: ierr
        !! Error code
        
        integer(int32) :: n_timepoints, n_samples, i
        real(real64) :: percent
        
        ! Initialize error
        call set_ok(ierr)
        
        ! Input validation
        n_timepoints = size(spike_contribs, 1)
        n_samples = size(spike_contribs, 2)
        
        if (n_timepoints == 0 .or. n_samples == 0) then
            call set_err(ierr, ERR_EMPTY_INPUT)
            return
        end if
        
        if (size(thresholds) /= n_timepoints) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if
        
        if (size(permutation, 1) /= n_timepoints .or. size(permutation, 2) /= n_samples) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if
        
        if (percentile_val < 0.0_real64 .or. percentile_val > 100.0_real64) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if
        
        ! Convert percentile to [0,1] range for calc_percentile
        percent = percentile_val / 100.0_real64
        
        ! Calculate threshold for each timepoint using pre-sorted permutations
        do i = 1, n_timepoints
            ! Use pre-computed permutation for this timepoint
            call calc_percentile(spike_contribs(i, :), permutation(i, :), percent, thresholds(i), ierr)
            if (is_err(ierr)) then
                return
            end if
        end do
        
    end subroutine calc_spike_thresholds

    !> Calculate empirical thresholds for spike contributions (with internal allocations and sorting)
    !!
    !! Convenience version that handles workspace allocation and sorting internally.
    !! For performance-critical code, use calc_spike_thresholds with pre-computed permutations.
    subroutine calc_spike_thresholds_alloc(spike_contribs, percentile_val, thresholds, ierr)
        real(real64), intent(in) :: spike_contribs(:, :)
        !! 2D array of spike contributions [n_timepoints, n_samples]
        real(real64), intent(in) :: percentile_val
        !! Percentile value for threshold (0.0-100.0)
        real(real64), intent(out) :: thresholds(:)
        !! 1D array of thresholds for each timepoint [n_timepoints]
        integer(int32), intent(out) :: ierr
        !! Error code
        
        integer(int32) :: n_timepoints, n_samples, i
        integer(int32), allocatable :: permutation(:, :), stack_left(:), stack_right(:)
        real(real64), allocatable :: timepoint_data(:)
        
        ! Initialize error
        call set_ok(ierr)
        
        ! Input validation
        n_timepoints = size(spike_contribs, 1)
        n_samples = size(spike_contribs, 2)
        
        if (n_timepoints == 0 .or. n_samples == 0) then
            call set_err(ierr, ERR_EMPTY_INPUT)
            return
        end if
        
        if (size(thresholds) /= n_timepoints) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if
        
        ! Allocate 2D permutation array and workspace
        allocate(permutation(n_timepoints, n_samples), stat=ierr)
        allocate(stack_left(n_samples), stat=ierr)
        allocate(stack_right(n_samples), stat=ierr)
        if (is_err(ierr)) then
            call set_err(ierr, ERR_ALLOC_FAIL)
            return
        end if
        
        allocate(timepoint_data(n_samples), stat=ierr)
        if (is_err(ierr)) then
            call set_err(ierr, ERR_ALLOC_FAIL)
            return
        end if
        
        ! Pre-compute permutations by sorting each timepoint
        do i = 1, n_timepoints
            ! Extract and sort data for this timepoint
            timepoint_data = spike_contribs(i, :)
            call sort_real(timepoint_data, permutation(i, :), stack_left, stack_right)
        end do
        
        ! Now call the main subroutine with pre-computed permutations
        call calc_spike_thresholds(spike_contribs, percentile_val, thresholds, permutation, ierr)
        
    end subroutine calc_spike_thresholds_alloc

    !> Calculate empirical threshold for integrated contributions (using pre-sorted permutation)
    !!
    !! Performance version that uses pre-computed permutation array to avoid internal sorting.
    pure subroutine calc_integrated_threshold(contributions, percentile_val, threshold, permutation, ierr)
        real(real64), intent(in) :: contributions(:)
        !! 1D array of integrated contributions [n_samples]
        real(real64), intent(in) :: percentile_val
        !! Percentile value for threshold (0.0-100.0)
        real(real64), intent(out) :: threshold
        !! Scalar threshold value
        integer(int32), intent(in) :: permutation(:)
        !! Pre-computed permutation indices for sorted contributions [n_samples]
        integer(int32), intent(out) :: ierr
        !! Error code
        
        integer(int32) :: n_samples
        real(real64) :: percent
        
        ! Initialize error
        call set_ok(ierr)
        
        ! Input validation
        n_samples = size(contributions)
        
        if (n_samples == 0) then
            call set_err(ierr, ERR_EMPTY_INPUT)
            return
        end if
        
        if (size(permutation) /= n_samples) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if
        
        if (percentile_val < 0.0_real64 .or. percentile_val > 100.0_real64) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if
        
        ! Convert percentile to [0,1] range for calc_percentile
        percent = percentile_val / 100.0_real64
        
        ! Calculate percentile threshold using pre-computed permutation
        call calc_percentile(contributions, permutation, percent, threshold, ierr)
        
    end subroutine calc_integrated_threshold

    !> Calculate empirical threshold for integrated contributions (with internal allocation and sorting)
    !!
    !! Convenience version that handles workspace allocation and sorting internally.
    !! For performance-critical code, use calc_integrated_threshold with pre-computed permutation.
    subroutine calc_integrated_threshold_alloc(contributions, percentile_val, threshold, ierr)
        real(real64), intent(in) :: contributions(:)
        !! 1D array of integrated contributions [n_samples]
        real(real64), intent(in) :: percentile_val
        !! Percentile value for threshold (0.0-100.0)
        real(real64), intent(out) :: threshold
        !! Scalar threshold value
        integer(int32), intent(out) :: ierr
        !! Error code
        
        integer(int32) :: n_samples
        integer(int32), allocatable :: permutation(:), stack_left(:), stack_right(:)
        real(real64), allocatable :: sorted_contributions(:)
        
        ! Initialize error
        call set_ok(ierr)
        
        ! Input validation
        n_samples = size(contributions)
        
        if (n_samples == 0) then
            call set_err(ierr, ERR_EMPTY_INPUT)
            return
        end if
        
        ! Allocate workspace
        allocate(permutation(n_samples), stat=ierr)
        allocate(stack_left(n_samples), stat=ierr)
        allocate(stack_right(n_samples), stat=ierr)
        if (is_err(ierr)) then
            call set_err(ierr, ERR_ALLOC_FAIL)
            return
        end if
        
        allocate(sorted_contributions(n_samples), stat=ierr)
        if (is_err(ierr)) then
            call set_err(ierr, ERR_ALLOC_FAIL)
            return
        end if
        
        ! Pre-compute permutation by sorting
        sorted_contributions = contributions
        call sort_real(sorted_contributions, permutation, stack_left, stack_right)
        
        ! call the main subroutine with pre-computed permutation
        call calc_integrated_threshold(contributions, percentile_val, threshold, permutation, ierr)        
    end subroutine calc_integrated_threshold_alloc

    !> Detect outliers in integrated (trajectory-level) contributions
    !!
    !! Identifies samples where the integrated contribution exceeds the empirical threshold.
    !! Used to find trajectories with significant overall explanatory contributions.
    pure subroutine detect_outliers_integrated(contributions, threshold, outlier_mask, ierr)
        real(real64), intent(in) :: contributions(:)
        !! Array of integrated contributions [n_samples]
        real(real64), intent(in) :: threshold
        !! Scalar threshold value for outlier detection
        logical, intent(out) :: outlier_mask(:)
        !! 1D logical array indicating outliers [n_samples]
        integer(int32), intent(out) :: ierr
        !! Error code
        integer(int32) :: n_samples, i
        
        ! Initialize error
        call set_ok(ierr)
        
        ! Input validation
        n_samples = size(contributions)
        
        if (n_samples == 0) then
            call set_err(ierr, ERR_EMPTY_INPUT)
            return
        end if
        
        if (size(outlier_mask) /= n_samples) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if
        
        ! Detect outliers
        do i = 1, n_samples
            outlier_mask(i) = contributions(i) > threshold
        end do
        
    end subroutine detect_outliers_integrated

    !> Detect outliers in spike contributions
    !!
    !! Identifies timepoint-sample pairs where spike contributions exceed their
    !! timepoint-specific empirical thresholds. Used to find local significant events.
    pure subroutine detect_outliers_spike(spike_contribs, thresholds, outlier_mask, ierr)
        real(real64), intent(in) :: spike_contribs(:)
        !! Array of spike contributions [n_timepoints, n_samples]
        real(real64), intent(in) :: thresholds(:)
        !! Array of thresholds for each timepoint [n_timepoints]
        logical, intent(out) :: outlier_mask(:)
        !! Logical array indicating outliers [n_timepoints, n_samples]
        integer(int32), intent(out) :: ierr
        !! Error code
        
        integer(int32) :: n_timepoints, i
        
        ! Initialize error
        call set_ok(ierr)
        
        ! Input validation
        n_timepoints = size(spike_contribs)
        
        if (n_timepoints == 0) then
            call set_err(ierr, ERR_EMPTY_INPUT)
            return
        end if
        
        if (size(thresholds) /= n_timepoints) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if
        
        if (size(outlier_mask) /= n_timepoints) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if
        
        ! Detect outliers
        do i = 1, n_timepoints
            outlier_mask(i) = spike_contribs(i) > thresholds(i)
        end do
        
    end subroutine detect_outliers_spike
end module tox_trajectory_contribution_analysis

!> C wrapper for calc_spike_thresholds (with pre-computed 2D permutations)
subroutine calc_spike_thresholds_C(spike_contribs, n_timepoints, n_samples, &
                                    percentile_val, thresholds, permutation, ierr) &
                                    bind(C, name="calc_spike_thresholds_C")
    use iso_c_binding, only: c_double, c_int
    use tox_errors, only: set_ok, set_err, ERR_EMPTY_INPUT
    use tox_trajectory_contribution_analysis, only: calc_spike_thresholds
    real(c_double), intent(in) :: spike_contribs(n_timepoints, n_samples)
    !! 2D array of spike contributions [n_timepoints, n_samples]
    integer(c_int), intent(in), value :: n_timepoints
    !! Number of timepoints
    integer(c_int), intent(in), value :: n_samples
    !! Number of samples
    real(c_double), intent(in), value :: percentile_val
    !! Percentile value for threshold (0.0-100.0)
    real(c_double), intent(out) :: thresholds(n_timepoints)
    !! 1D array of thresholds for each timepoint [n_timepoints]
    integer(c_int), intent(in) :: permutation(n_timepoints, n_samples)
    !! Pre-computed permutation indices for each timepoint [n_timepoints, n_samples]
    integer(c_int), intent(out) :: ierr
    !! Error code

    
    ! Initialize error
    call set_ok(ierr)
    
    ! Check for valid dimensions
    if (n_timepoints <= 0 .or. n_samples <= 0) then
        call set_err(ierr, ERR_EMPTY_INPUT)  ! ERR_EMPTY_INPUT
        return
    end if
    
    ! Call Fortran subroutine
    call calc_spike_thresholds(spike_contribs, percentile_val, &
                                thresholds, permutation, ierr)
    
end subroutine calc_spike_thresholds_c

!> C wrapper for calc_spike_thresholds_alloc (with internal allocations)
subroutine calc_spike_thresholds_alloc_C(spike_contribs, n_timepoints, n_samples, &
                                        percentile_val, thresholds, ierr) &
                                        bind(C, name="calc_spike_thresholds_alloc_C")
    use iso_c_binding, only: c_double, c_int
    use tox_errors, only: set_ok, set_err, ERR_EMPTY_INPUT
    use tox_trajectory_contribution_analysis, only: calc_spike_thresholds_alloc
    real(c_double), intent(in) :: spike_contribs(n_timepoints, n_samples)
    !! Array of spike contributions [n_timepoints, n_samples]
    integer(c_int), intent(in), value :: n_timepoints, n_samples
    !! Number of timepoints and samples
    real(c_double), intent(in), value :: percentile_val
    !! Percentile value for threshold (0.0-100.0)
    real(c_double), intent(out) :: thresholds(n_timepoints)
    !! 1D array of thresholds for each timepoint [n_timepoints]
    integer(c_int), intent(out) :: ierr
    !! Error code
    
    ! Initialize error
    call set_ok(ierr)
    
    ! Check for valid dimensions
    if (n_timepoints <= 0 .or. n_samples <= 0) then
        call set_err(ierr, ERR_EMPTY_INPUT)
        return
    end if

    ! Call Fortran subroutine
    call calc_spike_thresholds_alloc(spike_contribs, percentile_val, &
                                    thresholds, ierr)

end subroutine calc_spike_thresholds_alloc_c

!> C wrapper for calc_integrated_threshold (with pre-computed permutation)
subroutine calc_integrated_threshold_C(contributions, n_samples, percentile_val, &
                                        threshold, permutation, ierr) &
                                        bind(C, name="calc_integrated_threshold_C")
    use iso_c_binding, only: c_double, c_int
    use tox_errors, only: set_ok, set_err, ERR_EMPTY_INPUT
    use tox_trajectory_contribution_analysis, only: calc_integrated_threshold
    real(c_double), intent(in) :: contributions(n_samples)
    !! 1D array of integrated contributions [n_samples]
    integer(c_int), intent(in), value :: n_samples
    !! Number of samples
    real(c_double), intent(in), value :: percentile_val
    !! Percentile value for threshold (0.0-100.0)
    real(c_double), intent(out) :: threshold
    !! Scalar threshold value
    integer(c_int), intent(in) :: permutation(n_samples)
    !! Pre-computed permutation indices for sorted contributions [n_samples]
    integer(c_int), intent(out) :: ierr
    !! Error code
    
    ! Initialize error
    call set_ok(ierr)
    
    ! Check for valid dimensions
    if (n_samples <= 0) then
        call set_err(ierr, ERR_EMPTY_INPUT)  ! ERR_EMPTY_INPUT
        return
    end if
    
    ! Call Fortran subroutine
    call calc_integrated_threshold(contributions, percentile_val, &
                                threshold, permutation, ierr)
    
end subroutine calc_integrated_threshold_c

!> C wrapper for calc_integrated_threshold_alloc (with internal allocations)
subroutine calc_integrated_threshold_alloc_c(contributions, n_samples, percentile_val, &
                                            threshold, ierr) &
                                            bind(c, name="calc_integrated_threshold_alloc")
    use iso_c_binding, only: c_double, c_int
    use tox_errors, only: set_ok, set_err, ERR_EMPTY_INPUT
    use tox_trajectory_contribution_analysis, only: calc_integrated_threshold_alloc
    real(c_double), intent(in) :: contributions(n_samples)
    integer(c_int), intent(in), value :: n_samples
    real(c_double), intent(in), value :: percentile_val
    real(c_double), intent(out) :: threshold
    integer(c_int), intent(out) :: ierr
    
    ! Initialize error
    call set_ok(ierr)
    
    ! Check for valid dimensions
    if (n_samples <= 0) then
        call set_err(ierr, ERR_EMPTY_INPUT)  ! ERR_EMPTY_INPUT
        return
    end if
    
    ! Call Fortran subroutine
    call calc_integrated_threshold_alloc(contributions, percentile_val, &
                                        threshold, ierr)
    
end subroutine calc_integrated_threshold_alloc_c

!> C wrapper for detect_outliers_integrated
subroutine detect_outliers_integrated_c(contributions, n_samples, threshold, &
                                        outlier_mask, ierr) &
                                        bind(c, name="detect_outliers_integrated")
    use tox_conversions, only: logical_as_c_int
    use iso_c_binding, only: c_double, c_int
    use tox_errors, only: set_ok, set_err, ERR_EMPTY_INPUT, ERR_ALLOC_FAIL, is_err
    use tox_trajectory_contribution_analysis, only: detect_outliers_integrated
    real(c_double), intent(in) :: contributions(n_samples)
    integer(c_int), intent(in), value :: n_samples
    real(c_double), intent(in), value :: threshold
    integer(c_int), intent(out) :: outlier_mask(n_samples)  ! c_int for logical
    integer(c_int), intent(out) :: ierr
    
    logical, allocatable :: mask_1d(:)
    integer :: i
    
    ! Initialize error
    call set_ok(ierr)
    
    ! Check for valid dimensions
    if (n_samples <= 0) then
        call set_err(ierr, ERR_EMPTY_INPUT)  ! ERR_EMPTY_INPUT
        return
    end if

    allocate(mask_1d(n_samples), stat=ierr)
    if(is_err(ierr)) then
        call set_err(ierr, ERR_ALLOC_FAIL)
        return
    end if
    
    ! Call Fortran subroutine
    call detect_outliers_integrated(contributions, threshold, &
                                    mask_1d, ierr)
    
    ! Convert Fortran logical to C int (0=false, 1=true)
    if (.not. is_err(ierr)) then
        call logical_as_c_int(mask_1d, outlier_mask)
    else
        ! Initialize output to false on error
        do i = 1, n_samples
            outlier_mask(i) = 0_c_int
        end do
    end if
    
end subroutine detect_outliers_integrated_c

!> C wrapper for detect_outliers_spike
subroutine detect_outliers_spike_c(spike_contribs, n_timepoints, thresholds, &
                                    outlier_mask, ierr) &
                                    bind(c, name="detect_outliers_spike")
    use tox_conversions, only: logical_as_c_int
    use iso_c_binding, only: c_double, c_int
    use tox_errors, only: set_ok, set_err, is_err, ERR_EMPTY_INPUT, ERR_ALLOC_FAIL
    use tox_trajectory_contribution_analysis, only: detect_outliers_spike
    real(c_double), intent(in) :: spike_contribs(n_timepoints)
    integer(c_int), intent(in), value :: n_timepoints
    real(c_double), intent(in) :: thresholds(n_timepoints)
    integer(c_int), intent(out) :: outlier_mask(n_timepoints)  ! c_int for logical
    integer(c_int), intent(out) :: ierr
    
    logical, allocatable :: mask(:)
    integer :: i
    
    ! Initialize error
    call set_ok(ierr)
    
    ! Check for valid dimensions
    if (n_timepoints <= 0) then
        call set_err(ierr, ERR_EMPTY_INPUT)  ! ERR_EMPTY_INPUT
        return
    end if

    allocate(mask(n_timepoints), stat=ierr)
    if(is_err(ierr)) then
        call set_err(ierr, ERR_ALLOC_FAIL)
        return
    end if
    
    ! Call Fortran subroutine
    call detect_outliers_spike(spike_contribs, thresholds, mask, ierr)
    
    if(.not. is_err(ierr)) then
        ! Convert Fortran logical to C int (0=false, 1=true)
        call logical_as_c_int(mask, outlier_mask)
    else
        ! Initialize output to false on error
        do i = 1, n_timepoints
            outlier_mask(i) = 0_c_int
        end do
    end if
    
end subroutine detect_outliers_spike_c