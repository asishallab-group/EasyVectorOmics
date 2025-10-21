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
    pure subroutine calc_spike_thresholds(spike_contribs, percentile_val, thresholds, permutation, ierr)
        real(real64), intent(in) :: spike_contribs(:, :)
        !! 2D array of spike contributions [n_samples, n_genes] - rows=samples, columns=genes
        real(real64), intent(in) :: percentile_val
        !! Percentile value for threshold (0.0-100.0)
        real(real64), intent(out) :: thresholds(:)
        !! 1D array of thresholds for each GENE [n_genes]
        integer(int32), intent(in) :: permutation(:, :)
        !! Pre-computed permutation indices [n_samples, n_genes] - each COLUMN contains sorted indices for a gene
        integer(int32), intent(out) :: ierr
        !! Error code
        
        integer(int32) :: n_samples, n_genes, j
        
        ! Initialize error
        call set_ok(ierr)
        
        ! Input validation
        n_samples = size(spike_contribs, 1)
        n_genes = size(spike_contribs, 2)
        
        if (n_samples == 0 .or. n_genes == 0) then
            call set_err(ierr, ERR_EMPTY_INPUT)
            return
        end if
        
        if (size(thresholds) /= n_genes) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if
        
        ! CORRECTED: Check permutation dimensions match data
        if (size(permutation, 1) /= n_samples .or. size(permutation, 2) /= n_genes) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if
        
        if (percentile_val < 0.0_real64 .or. percentile_val > 100.0_real64) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if
        
        ! Calculate threshold for each GENE using pre-sorted permutations
        do j = 1, n_genes
            ! CORRECTED: Use the j-th COLUMN of permutation for the j-th gene
            call calc_percentile(spike_contribs(:, j), permutation(:, j), percentile_val, thresholds(j), ierr)
            if (is_err(ierr)) then
                return
            end if
        end do
        
    end subroutine calc_spike_thresholds

    !> Calculate empirical thresholds for spike contributions (with internal allocations and sorting)
    subroutine calc_spike_thresholds_alloc(spike_contribs, percentile_val, thresholds, ierr)
        real(real64), intent(in) :: spike_contribs(:, :)
        !! 2D array of spike contributions [n_samples, n_genes] - rows=samples, columns=genes
        real(real64), intent(in) :: percentile_val
        !! Percentile value for threshold (0.0-100.0)
        real(real64), intent(out) :: thresholds(:)
        !! 1D array of thresholds for each GENE [n_genes]
        integer(int32), intent(out) :: ierr
        !! Error code
        
        integer(int32) :: n_samples, n_genes, j, i, k
        integer(int32), allocatable :: permutation(:, :), stack_left(:), stack_right(:), temp_perm(:)
        real(real64), allocatable :: gene_data(:)
        
        ! Initialize error
        call set_ok(ierr)
        
        ! Input validation
        n_samples = size(spike_contribs, 1)
        n_genes = size(spike_contribs, 2)
        
        if (n_samples == 0 .or. n_genes == 0) then
            call set_err(ierr, ERR_EMPTY_INPUT)
            return
        end if
        
        if (size(thresholds) /= n_genes) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if
        
        ! CORRECTED: permutation(n_samples, n_genes) to match data layout
        allocate(permutation(n_samples, n_genes), stat=ierr)
        allocate(temp_perm(n_samples), stat=ierr)
        if (is_err(ierr)) then
            call set_err(ierr, ERR_ALLOC_FAIL)
            return
        end if

        do i = 1, n_genes
            permutation(:, i) = [(k, k = 1, n_samples)]
        end do
        
        allocate(stack_left(n_samples), stat=ierr)
        if (is_err(ierr)) then
            call set_err(ierr, ERR_ALLOC_FAIL)
            deallocate(permutation)
            return
        end if
        
        allocate(stack_right(n_samples), stat=ierr)
        if (is_err(ierr)) then
            call set_err(ierr, ERR_ALLOC_FAIL)
            deallocate(permutation, stack_left)
            return
        end if
        
        allocate(gene_data(n_samples), stat=ierr)
        if (is_err(ierr)) then
            call set_err(ierr, ERR_ALLOC_FAIL)
            deallocate(permutation, stack_left, stack_right)
            return
        end if
        
        ! Pre-compute permutations by sorting each GENE (column)
        do j = 1, n_genes
            ! Extract data for this gene (column) - this is efficient in column-major
            gene_data = spike_contribs(:, j)
            
            ! CORRECTED: Use the j-th COLUMN of permutation, not row
            call sort_real(gene_data, permutation(:, j), stack_left, stack_right)
        end do
        
        ! Now call the main subroutine with pre-computed permutations
        call calc_spike_thresholds(spike_contribs, percentile_val, thresholds, permutation, ierr)
        
        ! Clean up allocations
        deallocate(permutation, stack_left, stack_right, gene_data)
        
    end subroutine calc_spike_thresholds_alloc

    !> Calculate empirical threshold for integrated (trajectory-level) contributions
    !!
    !! Determines which SAMPLES are significant overall across all genes.
    !! Uses pre-computed permutation array to avoid internal sorting.
    pure subroutine calc_integrated_threshold(contributions, percentile_val, threshold, permutation, ierr)
        real(real64), intent(in) :: contributions(:)
        !! 1D array of integrated contributions [n_samples] - one value per SAMPLE
        real(real64), intent(in) :: percentile_val
        !! Percentile value for threshold (0.0-100.0)
        real(real64), intent(out) :: threshold
        !! Scalar threshold value
        integer(int32), intent(in) :: permutation(:)
        !! Pre-computed permutation indices for sorted contributions [n_samples]
        integer(int32), intent(out) :: ierr
        !! Error code
        
        ! Input validation is handled by calc_percentile
        call calc_percentile(contributions, permutation, percentile_val, threshold, ierr)
        
    end subroutine calc_integrated_threshold

    !> Calculate empirical threshold for integrated contributions (with internal allocation and sorting)
    !!
    !! Convenience version that handles workspace allocation and sorting internally.
    !! For performance-critical code, use calc_integrated_threshold with pre-computed permutation.
    subroutine calc_integrated_threshold_alloc(contributions, percentile_val, threshold, ierr)
        real(real64), intent(in) :: contributions(:)
        !! 1D array of integrated contributions [n_samples] - one value per SAMPLE
        real(real64), intent(in) :: percentile_val
        !! Percentile value for threshold (0.0-100.0)
        real(real64), intent(out) :: threshold
        !! Scalar threshold value
        integer(int32), intent(out) :: ierr
        !! Error code
        
        ! Use the existing calc_percentile_alloc from f42_utils
        call calc_percentile_alloc(contributions, percentile_val, threshold, ierr)
        
    end subroutine calc_integrated_threshold_alloc

    !> Detect outliers in integrated (trajectory-level) contributions
    !!
    !! Identifies SAMPLES where the integrated contribution exceeds the empirical threshold.
    !! Used to find trajectories with significant overall explanatory contributions.
    pure subroutine detect_outliers_integrated(contributions, threshold, outlier_mask, ierr)
        real(real64), intent(in) :: contributions(:)
        !! Array of integrated contributions [n_samples] - one value per SAMPLE
        real(real64), intent(in) :: threshold
        !! Scalar threshold value for outlier detection
        logical, intent(out) :: outlier_mask(:)
        !! 1D logical array indicating outlier SAMPLES [n_samples]
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
        
        ! Detect outlier SAMPLES
        do i = 1, n_samples
            outlier_mask(i) = contributions(i) > threshold
        end do
        
    end subroutine detect_outliers_integrated

    !> Detect outliers in spike contributions
    !!
    !! Identifies sample-gene pairs where spike contributions exceed their
    !! gene-specific empirical thresholds. Used to find local significant events.
    !! Data layout: spike_contribs(n_samples, n_genes) - rows=samples, columns=genes
    pure subroutine detect_outliers_spike(spike_contribs, thresholds, outlier_mask, ierr)
        real(real64), intent(in) :: spike_contribs(:, :)
        !! 2D array of spike contributions [n_samples, n_genes] - rows=samples, columns=genes
        real(real64), intent(in) :: thresholds(:)
        !! Array of thresholds for each GENE [n_genes]
        logical, intent(out) :: outlier_mask(:, :)
        !! 2D logical array indicating outliers [n_samples, n_genes]
        integer(int32), intent(out) :: ierr
        !! Error code
        
        integer(int32) :: n_samples, n_genes, i, j
        
        ! Initialize error
        call set_ok(ierr)
        
        ! Input validation
        n_samples = size(spike_contribs, 1)
        n_genes = size(spike_contribs, 2)
        
        if (n_samples == 0 .or. n_genes == 0) then
            call set_err(ierr, ERR_EMPTY_INPUT)
            return
        end if
        
        if (size(thresholds) /= n_genes) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if
        
        if (size(outlier_mask, 1) /= n_samples .or. size(outlier_mask, 2) /= n_genes) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if
        
        ! Outer loop over genes (columns), inner loop over samples (rows)
        do j = 1, n_genes
            do i = 1, n_samples
                outlier_mask(i, j) = spike_contribs(i, j) > thresholds(j)
            end do
        end do
        
    end subroutine detect_outliers_spike
end module tox_trajectory_contribution_analysis

!> C wrapper for calc_spike_thresholds (with pre-computed 2D permutations)
!! Data layout: spike_contribs(n_samples, n_genes) - rows=samples, columns=genes
subroutine calc_spike_thresholds_C(spike_contribs, n_samples, n_genes, &
                                    percentile_val, thresholds, permutation, ierr) &
                                    bind(C, name="calc_spike_thresholds_C")
    use iso_c_binding, only: c_double, c_int
    use tox_errors, only: set_ok, set_err, ERR_EMPTY_INPUT
    use tox_trajectory_contribution_analysis, only: calc_spike_thresholds
    real(c_double), intent(in) :: spike_contribs(n_samples, n_genes)
    !! 2D array of spike contributions [n_samples, n_genes] - rows=samples, columns=genes
    integer(c_int), intent(in), value :: n_samples
    !! Number of samples
    integer(c_int), intent(in), value :: n_genes
    !! Number of genes
    real(c_double), intent(in), value :: percentile_val
    !! Percentile value for threshold (0.0-100.0)
    real(c_double), intent(out) :: thresholds(n_genes)
    !! 1D array of thresholds for each gene [n_genes]
    integer(c_int), intent(in) :: permutation(n_samples, n_genes)  ! FIXED: [n_samples, n_genes]
    !! Pre-computed permutation indices [n_samples, n_genes] - each COLUMN for a gene
    integer(c_int), intent(out) :: ierr
    !! Error code

    
    ! Initialize error
    call set_ok(ierr)
    
    ! Check for valid dimensions
    if (n_samples <= 0 .or. n_genes <= 0) then
        call set_err(ierr, ERR_EMPTY_INPUT)
        return
    end if
    
    ! Call Fortran subroutine - FIXED: Now dimensions match
    call calc_spike_thresholds(spike_contribs, percentile_val, &
                                thresholds, permutation, ierr)
    
end subroutine calc_spike_thresholds_C

!> C wrapper for calc_spike_thresholds_alloc (with internal allocations)
!! Data layout: spike_contribs(n_samples, n_genes) - rows=samples, columns=genes
subroutine calc_spike_thresholds_alloc_C(spike_contribs, n_samples, n_genes, &
                                        percentile_val, thresholds, ierr) &
                                        bind(C, name="calc_spike_thresholds_alloc_C")
    use iso_c_binding, only: c_double, c_int
    use tox_errors, only: set_ok, set_err, ERR_EMPTY_INPUT
    use tox_trajectory_contribution_analysis, only: calc_spike_thresholds_alloc
    real(c_double), intent(in) :: spike_contribs(n_samples, n_genes)
    !! Array of spike contributions [n_samples, n_genes] - rows=samples, columns=genes
    integer(c_int), intent(in), value :: n_samples, n_genes
    !! Number of samples and genes
    real(c_double), intent(in), value :: percentile_val
    !! Percentile value for threshold (0.0-100.0)
    real(c_double), intent(out) :: thresholds(n_genes)
    !! 1D array of thresholds for each gene [n_genes]
    integer(c_int), intent(out) :: ierr
    !! Error code
    
    ! Initialize error
    call set_ok(ierr)
    
    ! Check for valid dimensions
    if (n_samples <= 0 .or. n_genes <= 0) then
        call set_err(ierr, ERR_EMPTY_INPUT)
        return
    end if

    ! Call Fortran subroutine
    call calc_spike_thresholds_alloc(spike_contribs, percentile_val, &
                                    thresholds, ierr)

end subroutine calc_spike_thresholds_alloc_C

!> C wrapper for calc_integrated_threshold (with pre-computed permutation)
!! Input: contributions(n_samples) - one integrated contribution per sample
subroutine calc_integrated_threshold_C(contributions, n_samples, percentile_val, &
                                        threshold, permutation, ierr) &
                                        bind(C, name="calc_integrated_threshold_C")
    use iso_c_binding, only: c_double, c_int
    use tox_errors, only: set_ok, set_err, ERR_EMPTY_INPUT
    use tox_trajectory_contribution_analysis, only: calc_integrated_threshold
    real(c_double), intent(in) :: contributions(n_samples)
    !! 1D array of integrated contributions [n_samples] - one value per sample
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
        call set_err(ierr, ERR_EMPTY_INPUT)
        return
    end if
    
    ! Call Fortran subroutine
    call calc_integrated_threshold(contributions, percentile_val, &
                                threshold, permutation, ierr)
    
end subroutine calc_integrated_threshold_C

!> C wrapper for calc_integrated_threshold_alloc (with internal allocations)
!! Input: contributions(n_samples) - one integrated contribution per sample
subroutine calc_integrated_threshold_alloc_C(contributions, n_samples, percentile_val, &
                                            threshold, ierr) &
                                            bind(C, name="calc_integrated_threshold_alloc_C")
    use iso_c_binding, only: c_double, c_int
    use tox_errors, only: set_ok, set_err, ERR_EMPTY_INPUT
    use tox_trajectory_contribution_analysis, only: calc_integrated_threshold_alloc
    real(c_double), intent(in) :: contributions(n_samples)
    !! 1D array of integrated contributions [n_samples] - one value per sample
    integer(c_int), intent(in), value :: n_samples
    !! Number of samples
    real(c_double), intent(in), value :: percentile_val
    !! Percentile value for threshold (0.0-100.0)
    real(c_double), intent(out) :: threshold
    !! Scalar threshold value
    integer(c_int), intent(out) :: ierr
    !! Error code
    
    ! Initialize error
    call set_ok(ierr)
    
    ! Check for valid dimensions
    if (n_samples <= 0) then
        call set_err(ierr, ERR_EMPTY_INPUT)
        return
    end if
    
    ! Call Fortran subroutine
    call calc_integrated_threshold_alloc(contributions, percentile_val, &
                                        threshold, ierr)
    
end subroutine calc_integrated_threshold_alloc_C

!> C wrapper for detect_outliers_integrated
!! Identifies outlier samples based on integrated contributions
subroutine detect_outliers_integrated_C(contributions, n_samples, threshold, &
                                        outlier_mask, ierr) &
                                        bind(C, name="detect_outliers_integrated_C")
    use tox_conversions, only: logical_as_c_int
    use iso_c_binding, only: c_double, c_int
    use tox_errors, only: set_ok, set_err, ERR_EMPTY_INPUT, ERR_ALLOC_FAIL, is_err
    use tox_trajectory_contribution_analysis, only: detect_outliers_integrated
    real(c_double), intent(in) :: contributions(n_samples)
    !! Array of integrated contributions [n_samples] - one value per sample
    integer(c_int), intent(in), value :: n_samples
    !! Number of samples
    real(c_double), intent(in), value :: threshold
    !! Scalar threshold value for outlier detection
    integer(c_int), intent(out) :: outlier_mask(n_samples)  ! c_int for logical
    !! 1D array indicating outlier samples [n_samples]
    integer(c_int), intent(out) :: ierr
    !! Error code
    
    logical, allocatable :: mask_1d(:)
    integer :: i
    
    ! Initialize error
    call set_ok(ierr)
    
    ! Check for valid dimensions
    if (n_samples <= 0) then
        call set_err(ierr, ERR_EMPTY_INPUT)
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
    
    deallocate(mask_1d)
    
end subroutine detect_outliers_integrated_C

!> C wrapper for detect_outliers_spike
!! Data layout: spike_contribs(n_samples, n_genes) - rows=samples, columns=genes
!! Identifies outlier gene-sample pairs using gene-specific thresholds
subroutine detect_outliers_spike_C(spike_contribs, n_samples, n_genes, thresholds, &
                                    outlier_mask, ierr) &
                                    bind(C, name="detect_outliers_spike_C")
    use tox_conversions, only: logical_as_c_int
    use iso_c_binding, only: c_double, c_int
    use tox_errors, only: set_ok, set_err, is_err, ERR_EMPTY_INPUT, ERR_ALLOC_FAIL
    use tox_trajectory_contribution_analysis, only: detect_outliers_spike
    real(c_double), intent(in) :: spike_contribs(n_samples, n_genes)
    !! 2D array of spike contributions [n_samples, n_genes] - rows=samples, columns=genes
    integer(c_int), intent(in), value :: n_samples
    !! Number of samples
    integer(c_int), intent(in), value :: n_genes
    !! Number of genes
    real(c_double), intent(in) :: thresholds(n_genes)
    !! Array of thresholds for each gene [n_genes]
    integer(c_int), intent(out) :: outlier_mask(n_samples, n_genes)  ! c_int for logical
    !! 2D array indicating outlier gene-sample pairs [n_samples, n_genes]
    integer(c_int), intent(out) :: ierr
    !! Error code
    
    logical, allocatable :: mask_2d(:,:)
    integer :: i, j
    
    ! Initialize error
    call set_ok(ierr)
    
    ! Check for valid dimensions
    if (n_samples <= 0 .or. n_genes <= 0) then
        call set_err(ierr, ERR_EMPTY_INPUT)
        return
    end if

    allocate(mask_2d(n_samples, n_genes), stat=ierr)
    if(is_err(ierr)) then
        call set_err(ierr, ERR_ALLOC_FAIL)
        return
    end if
    
    ! Call Fortran subroutine
    call detect_outliers_spike(spike_contribs, thresholds, mask_2d, ierr)
    
    if(.not. is_err(ierr)) then
        ! Convert Fortran logical to C int (0=false, 1=true)
        call logical_as_c_int(mask_2d, outlier_mask)
    else
        ! Initialize output to false on error
        do j = 1, n_genes
            do i = 1, n_samples
                outlier_mask(i, j) = 0_c_int
            end do
        end do
    end if
    
    deallocate(mask_2d)
    
end subroutine detect_outliers_spike_C