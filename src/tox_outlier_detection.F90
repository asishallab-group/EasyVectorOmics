!> Module for trajectory explanatory contribution analysis (TCA)
!! Provides functionality for detecting outliers in trajectory and spike contributions
module tox_trajectory_contribution_analysis
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_errors, only: ERROR_OK, ERROR_INVALID_INPUT, ERROR_EMPTY_INPUT, &
                         ERROR_ALLOC_FAIL, set_ok, set_err_once, is_err
    use f42_utils, only: calc_percentile
    
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
    !!
    !! @param[in] spike_contribs 2D array of spike contributions [n_timepoints, n_samples]
    !! @param[in] percentile_val Percentile value for threshold (0.0-100.0)
    !! @param[out] thresholds 1D array of thresholds for each timepoint [n_timepoints]
    !! @param[in] permutation 2D workspace array of pre-sorted permutation indices [n_timepoints, n_samples]
    !! @param[out] ierr Error code (0 = success, >0 = error)
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
            call set_err_once(ierr, ERR_EMPTY_INPUT)
            return
        end if
        
        if (size(thresholds) /= n_timepoints) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if
        
        if (size(permutation, 1) /= n_timepoints .or. size(permutation, 2) /= n_samples) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if
        
        if (percentile_val < 0.0_real64 .or. percentile_val > 100.0_real64) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
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
    pure subroutine calc_spike_thresholds_alloc(spike_contribs, percentile_val, thresholds, ierr)
        real(real64), intent(in) :: spike_contribs(:, :)
        !! 2D array of spike contributions [n_timepoints, n_samples]
        real(real64), intent(in) :: percentile_val
        !! Percentile value for threshold (0.0-100.0)
        real(real64), intent(out) :: thresholds(:)
        !! 1D array of thresholds for each timepoint [n_timepoints]
        integer(int32), intent(out) :: ierr
        !! Error code
        
        integer(int32) :: n_timepoints, n_samples, i
        integer(int32), allocatable :: permutation(:, :)
        real(real64), allocatable :: timepoint_data(:)
        
        ! Initialize error
        call set_ok(ierr)
        
        ! Input validation
        n_timepoints = size(spike_contribs, 1)
        n_samples = size(spike_contribs, 2)
        
        if (n_timepoints == 0 .or. n_samples == 0) then
            call set_err_once(ierr, ERR_EMPTY_INPUT)
            return
        end if
        
        if (size(thresholds) /= n_timepoints) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if
        
        ! Allocate 2D permutation array and workspace
        allocate(permutation(n_timepoints, n_samples), stat=ierr)
        if (is_err(ierr)) then
            call set_err(ierr, ERR_ALLOC_FAIL)
            return
        end if
        
        allocate(timepoint_data(n_samples), stat=ierr)
        if (is_err(ierr)) then
            call set_err(ierr, ERR_ALLOC_FAIL)
            deallocate(permutation)
            return
        end if
        
        ! Pre-compute permutations by sorting each timepoint
        do i = 1, n_timepoints
            ! Extract and sort data for this timepoint
            timepoint_data = spike_contribs(i, :)
            call sort_array(timepoint_data, permutation(i, :))
        end do
        
        ! Now call the main subroutine with pre-computed permutations
        call calc_spike_thresholds(spike_contribs, percentile_val, thresholds, permutation, ierr)
        
        ! Clean up
        deallocate(permutation, timepoint_data)
        
    end subroutine calc_spike_thresholds_alloc

    !> Calculate empirical threshold for integrated contributions (using pre-sorted permutation)
    !!
    !! Performance version that uses pre-computed permutation array to avoid internal sorting.
    !!
    !! @param[in] contributions 1D array of integrated contributions [n_samples]
    !! @param[in] percentile_val Percentile value for threshold (0.0-100.0)
    !! @param[out] threshold Scalar threshold value
    !! @param[in] permutation Pre-computed permutation indices for sorted contributions [n_samples]
    !! @param[out] ierr Error code (0 = success, >0 = error)
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
            call set_err_once(ierr, ERR_EMPTY_INPUT)
            return
        end if
        
        if (size(permutation) /= n_samples) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if
        
        if (percentile_val < 0.0_real64 .or. percentile_val > 100.0_real64) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
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
    pure subroutine calc_integrated_threshold_alloc(contributions, percentile_val, threshold, ierr)
        real(real64), intent(in) :: contributions(:)
        !! 1D array of integrated contributions [n_samples]
        real(real64), intent(in) :: percentile_val
        !! Percentile value for threshold (0.0-100.0)
        real(real64), intent(out) :: threshold
        !! Scalar threshold value
        integer(int32), intent(out) :: ierr
        !! Error code
        
        integer(int32) :: n_samples
        integer(int32), allocatable :: permutation(:)
        real(real64), allocatable :: sorted_contributions(:)
        
        ! Initialize error
        call set_ok(ierr)
        
        ! Input validation
        n_samples = size(contributions)
        
        if (n_samples == 0) then
            call set_err_once(ierr, ERR_EMPTY_INPUT)
            return
        end if
        
        ! Allocate workspace
        allocate(permutation(n_samples), stat=ierr)
        if (is_err(ierr)) then
            call set_err(ierr, ERR_ALLOC_FAIL)
            return
        end if
        
        allocate(sorted_contributions(n_samples), stat=ierr)
        if (is_err(ierr)) then
            call set_err(ierr, ERR_ALLOC_FAIL)
            deallocate(permutation)
            return
        end if
        
        ! Pre-compute permutation by sorting
        sorted_contributions = contributions
        call sort_array(sorted_contributions, permutation)
        
        ! Now call the main subroutine with pre-computed permutation
        call calc_integrated_threshold(contributions, percentile_val, threshold, permutation, ierr)
        
        ! Clean up
        deallocate(permutation, sorted_contributions)
        
    end subroutine calc_integrated_threshold_alloc

    !> Detect outliers in integrated (trajectory-level) contributions
    !!
    !! Identifies samples where the integrated contribution exceeds the empirical threshold.
    !! Used to find trajectories with significant overall explanatory contributions.
    !!
    !! @param[in] contributions 1D array of integrated contributions [n_samples]
    !! @param[in] threshold Scalar threshold value for outlier detection
    !! @param[out] outlier_mask 1D logical array indicating outliers [n_samples]
    !! @param[out] ierr Error code (0 = success, >0 = error)
    pure subroutine detect_outliers_integrated(contributions, threshold, outlier_mask, ierr)
        real(real64), intent(in) :: contributions(:)
        real(real64), intent(in) :: threshold
        logical, intent(out) :: outlier_mask(:)
        integer(int32), intent(out) :: ierr
        
        integer(int32) :: n_samples, i
        
        ! Initialize error
        call set_ok(ierr)
        
        ! Input validation
        n_samples = size(contributions)
        
        if (n_samples == 0) then
            call set_err_once(ierr, ERR_EMPTY_INPUT)
            return
        end if
        
        if (size(outlier_mask) /= n_samples) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
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
    !!
    !! @param[in] spike_contribs 2D array of spike contributions [n_timepoints, n_samples]
    !! @param[in] thresholds 1D array of thresholds for each timepoint [n_timepoints]
    !! @param[out] outlier_mask 2D logical array indicating outliers [n_timepoints, n_samples]
    !! @param[out] ierr Error code (0 = success, >0 = error)
    pure subroutine detect_outliers_spike(spike_contribs, thresholds, outlier_mask, ierr)
        real(real64), intent(in) :: spike_contribs(:, :)
        real(real64), intent(in) :: thresholds(:)
        logical, intent(out) :: outlier_mask(:, :)
        integer(int32), intent(out) :: ierr
        
        integer(int32) :: n_timepoints, n_samples, i, j
        
        ! Initialize error
        call set_ok(ierr)
        
        ! Input validation
        n_timepoints = size(spike_contribs, 1)
        n_samples = size(spike_contribs, 2)
        
        if (n_timepoints == 0 .or. n_samples == 0) then
            call set_err_once(ierr, ERR_EMPTY_INPUT)
            return
        end if
        
        if (size(thresholds) /= n_timepoints) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if
        
        if (size(outlier_mask, 1) /= n_timepoints .or. &
            size(outlier_mask, 2) /= n_samples) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if
        
        ! Detect outliers
        do j = 1, n_samples
            do i = 1, n_timepoints
                outlier_mask(i, j) = spike_contribs(i, j) > thresholds(i)
            end do
        end do
        
    end subroutine detect_outliers_spike

end module tox_trajectory_contribution_analysis