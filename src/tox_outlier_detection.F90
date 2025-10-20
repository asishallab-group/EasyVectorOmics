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
        real(real64), intent(in) :: spike_contribs(:, :)
        !! Array of spike contributions [n_timepoints, n_samples]
        real(real64), intent(in) :: thresholds(:)
        !! Array of thresholds for each timepoint [n_timepoints]
        logical, intent(out) :: outlier_mask(:, :)
        !! Logical array indicating outliers [n_timepoints, n_samples]
        integer(int32), intent(out) :: ierr
        !! Error code
        
        integer(int32) :: n_timepoints, n_samples, i, j
        
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
        
        if (size(outlier_mask, 1) /= n_timepoints .or. &
            size(outlier_mask, 2) /= n_samples) then
            call set_err(ierr, ERR_INVALID_INPUT)
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

!> C wrapper for calc_spike_thresholds (with pre-computed 2D permutations)
subroutine calc_spike_thresholds_C(spike_contribs, n_timepoints, n_samples, &
                                    percentile_val, thresholds, permutation, ierr) &
                                    bind(C, name="calc_spike_thresholds_C")
    use iso_c_binding, only: c_double, c_int
    use iso_fortran_env, only: real64, int32
    use tox_errors, only: set_ok, set_err, is_err, ERR_EMPTY_INPUT
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
    
    real(real64) :: contribs_2d(n_timepoints, n_samples)
    real(real64) :: thresh_1d(n_timepoints)
    integer(int32) :: perm_2d(n_timepoints, n_samples)
    integer :: i, j
    
    ! Initialize error
    call set_ok(ierr)
    
    ! Check for valid dimensions
    if (n_timepoints <= 0 .or. n_samples <= 0) then
        call set_err(ierr, ERR_EMPTY_INPUT)  ! ERR_EMPTY_INPUT
        return
    end if
    
    ! Convert input arrays (ORDER='F' in Python so same memory layout as Fortran)
    do j = 1, n_samples
        do i = 1, n_timepoints
            contribs_2d(i, j) = real(spike_contribs(i, j), real64)
            perm_2d(i, j) = int(permutation(i, j), int32)
        end do
    end do
    
    ! Call Fortran subroutine
    call calc_spike_thresholds(contribs_2d, real(percentile_val, real64), &
                                thresh_1d, perm_2d, ierr)
    
    ! Convert output back to C
    if (.not. is_err(ierr)) then
        do i = 1, n_timepoints
            thresholds(i) = real(thresh_1d(i), c_double)
        end do
    end if
    
end subroutine calc_spike_thresholds_c

!> C wrapper for calc_spike_thresholds_alloc (with internal allocations)
subroutine calc_spike_thresholds_alloc_C(spike_contribs, n_timepoints, n_samples, &
                                        percentile_val, thresholds, ierr) &
                                        bind(C, name="calc_spike_thresholds_alloc_C")
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
    
    real(real64) :: contribs_2d(n_timepoints, n_samples)
    !! 2D array for Fortran spike contributions
    real(real64) :: thresh_1d(n_timepoints)
    !! 1D array for Fortran thresholds

    integer :: i, j
    
    ! Initialize error
    call set_ok(ierr)
    
    ! Check for valid dimensions
    if (n_timepoints <= 0 .or. n_samples <= 0) then
        call set_err(fortran_ierr, ERR_EMPTY_INPUT)
        return
    end if
    
    ! Convert input array (ORDER='F' in Python so same memory layout as Fortran)
    do j = 1, n_samples
        do i = 1, n_timepoints
            contribs_2d(i, j) = real(spike_contribs(i, j), real64)
        end do
    end do
    
    ! Call Fortran subroutine
    call calc_spike_thresholds_alloc(contribs_2d, real(percentile_val, real64), &
                                    thresh_1d, fortran_ierr)
    
    ! Convert output back to C
    if (.not. is_err(fortran_ierr)) then
        do i = 1, n_timepoints
            thresholds(i) = real(thresh_1d(i), c_double)
        end do
    end if
    
    ! Set error code for C
    ierr = fortran_ierr
    
end subroutine calc_spike_thresholds_alloc_c

!> C wrapper for calc_integrated_threshold (with pre-computed permutation)
subroutine calc_integrated_threshold_C(contributions, n_samples, percentile_val, &
                                        threshold, permutation, ierr) &
                                        bind(C, name="calc_integrated_threshold_C")
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
    
    real(real64) :: contribs_1d(n_samples)
    integer(int32) :: perm_1d(n_samples)
    real(real64) :: threshold_f
    integer :: i
    
    ! Initialize error
    call set_ok(ierr)
    
    ! Check for valid dimensions
    if (n_samples <= 0) then
        call set_err(fortran_ierr, ERR_EMPTY_INPUT)  ! ERR_EMPTY_INPUT
        return
    end if
    
    ! Convert input arrays
    do i = 1, n_samples
        contribs_1d(i) = real(contributions(i), real64)
        perm_1d(i) = int(permutation(i), int32)
    end do
    
    ! Call Fortran subroutine
    call calc_integrated_threshold(contribs_1d, real(percentile_val, real64), &
                                threshold_f, perm_1d, ierr)
    
    ! Convert output back to C
    if (.not. is_err(ierr)) then
        threshold = real(threshold_f, c_double)
    else
        threshold = 0.0_c_double
    end if
    
end subroutine calc_integrated_threshold_c

!> C wrapper for calc_integrated_threshold_alloc (with internal allocations)
subroutine calc_integrated_threshold_alloc_c(contributions, n_samples, percentile_val, &
                                            threshold, ierr) &
    bind(c, name="calc_integrated_threshold_alloc")
    real(c_double), intent(in) :: contributions(n_samples)
    integer(c_int), intent(in), value :: n_samples
    real(c_double), intent(in), value :: percentile_val
    real(c_double), intent(out) :: threshold
    integer(c_int), intent(out) :: ierr
    
    real(real64) :: contribs_1d(n_samples)
    real(real64) :: threshold_f
    integer(int32) :: fortran_ierr
    integer :: i
    
    ! Initialize error
    call set_ok(fortran_ierr)
    ierr = ERR_OK
    
    ! Check for valid dimensions
    if (n_samples <= 0) then
        call set_err(fortran_ierr, 202)  ! ERR_EMPTY_INPUT
        ierr = fortran_ierr
        return
    end if
    
    ! Convert input array
    do i = 1, n_samples
        contribs_1d(i) = real(contributions(i), real64)
    end do
    
    ! Call Fortran subroutine
    call calc_integrated_threshold_alloc(contribs_1d, real(percentile_val, real64), &
                                        threshold_f, fortran_ierr)
    
    ! Convert output back to C
    if (.not. is_err(fortran_ierr)) then
        threshold = real(threshold_f, c_double)
    else
        threshold = 0.0_c_double
    end if
    
    ! Set error code for C
    ierr = fortran_ierr
    
end subroutine calc_integrated_threshold_alloc_c

!> C wrapper for detect_outliers_integrated
subroutine detect_outliers_integrated_c(contributions, n_samples, threshold, &
                                        outlier_mask, ierr) &
    bind(c, name="detect_outliers_integrated")
    real(c_double), intent(in) :: contributions(n_samples)
    integer(c_int), intent(in), value :: n_samples
    real(c_double), intent(in), value :: threshold
    integer(c_int), intent(out) :: outlier_mask(n_samples)  ! c_int for logical
    integer(c_int), intent(out) :: ierr
    
    real(real64) :: contribs_1d(n_samples)
    logical :: mask_1d(n_samples)
    integer(int32) :: fortran_ierr
    integer :: i
    
    ! Initialize error
    call set_ok(fortran_ierr)
    ierr = ERR_OK
    
    ! Check for valid dimensions
    if (n_samples <= 0) then
        call set_err(fortran_ierr, 202)  ! ERR_EMPTY_INPUT
        ierr = fortran_ierr
        return
    end if
    
    ! Convert input array
    do i = 1, n_samples
        contribs_1d(i) = real(contributions(i), real64)
    end do
    
    ! Call Fortran subroutine
    call detect_outliers_integrated(contribs_1d, real(threshold, real64), &
                                    mask_1d, fortran_ierr)
    
    ! Convert Fortran logical to C int (0=false, 1=true)
    if (.not. is_err(fortran_ierr)) then
        do i = 1, n_samples
            if (mask_1d(i)) then
                outlier_mask(i) = 1_c_int
            else
                outlier_mask(i) = 0_c_int
            end if
        end do
    else
        ! Initialize output to false on error
        do i = 1, n_samples
            outlier_mask(i) = 0_c_int
        end do
    end if
    
    ! Set error code for C
    ierr = fortran_ierr
    
end subroutine detect_outliers_integrated_c

!> C wrapper for detect_outliers_spike
subroutine detect_outliers_spike_c(spike_contribs, n_timepoints, n_samples, thresholds, &
                                    outlier_mask, ierr) &
    bind(c, name="detect_outliers_spike")
    real(c_double), intent(in) :: spike_contribs(n_timepoints, n_samples)
    integer(c_int), intent(in), value :: n_timepoints, n_samples
    real(c_double), intent(in) :: thresholds(n_timepoints)
    integer(c_int), intent(out) :: outlier_mask(n_timepoints, n_samples)  ! c_int for logical
    integer(c_int), intent(out) :: ierr
    
    real(real64) :: contribs_2d(n_timepoints, n_samples)
    real(real64) :: thresh_1d(n_timepoints)
    logical :: mask_2d(n_timepoints, n_samples)
    integer(int32) :: fortran_ierr
    integer :: i, j
    
    ! Initialize error
    call set_ok(fortran_ierr)
    ierr = ERR_OK
    
    ! Check for valid dimensions
    if (n_timepoints <= 0 .or. n_samples <= 0) then
        call set_err(fortran_ierr, 202)  ! ERR_EMPTY_INPUT
        ierr = fortran_ierr
        return
    end if
    
    ! Convert input arrays (ORDER='F' in Python so same memory layout as Fortran)
    do j = 1, n_samples
        do i = 1, n_timepoints
            contribs_2d(i, j) = real(spike_contribs(i, j), real64)
        end do
    end do
    
    ! Convert thresholds array
    do i = 1, n_timepoints
        thresh_1d(i) = real(thresholds(i), real64)
    end do
    
    ! Call Fortran subroutine
    call detect_outliers_spike(contribs_2d, thresh_1d, mask_2d, fortran_ierr)
    
    ! Convert Fortran logical to C int (0=false, 1=true)
    if (.not. is_err(fortran_ierr)) then
        do j = 1, n_samples
            do i = 1, n_timepoints
                if (mask_2d(i, j)) then
                    outlier_mask(i, j) = 1_c_int
                else
                    outlier_mask(i, j) = 0_c_int
                end if
            end do
        end do
    else
        ! Initialize output to false on error
        do j = 1, n_samples
            do i = 1, n_timepoints
                outlier_mask(i, j) = 0_c_int
            end do
        end do
    end if
    
    ! Set error code for C
    ierr = fortran_ierr
    
end subroutine detect_outliers_spike_c