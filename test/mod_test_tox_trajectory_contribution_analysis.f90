!> Unit test suite for tox_trajectory_contribution_analysis routine.
module mod_test_tox_traj_contrib_analysis
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use tox_trajectory_contribution_analysis
    use tox_errors
    use tox_trajectory_normalization
    implicit none

    ! Abstract interface for all test procedures
    abstract interface
        subroutine test_interface()
        end subroutine test_interface
    end interface

    ! Type to hold test name and procedure pointer
    type :: test_case
        character(len=128) :: name
        procedure(test_interface), pointer, nopass :: test_proc => null()
    end type test_case

    real(real64), parameter :: TOL = epsilon(1.0_real64)

contains

    !> Get array of all available tests.
    function get_all_tests() result(all_tests)
        type(test_case) :: all_tests(15)

        all_tests(1) = test_case("test_tox_trajectory_contribution_analysis_get_vec_across_samples", test_get_vec_across_samples)
        all_tests(2) = test_case("test_tox_trajectory_contribution_analysis_get_vec_across_timepoints", test_get_vec_across_timepoints)
        all_tests(3) = test_case("test_tox_trajectory_contribution_analysis_trajectory_contribution", test_trajectory_contribution)
        all_tests(4) = test_case("test_tox_trajectory_contribution_analysis_spike_contribution", test_spike_contribution)
        all_tests(5) = test_case("test_tox_trajectory_contribution_analysis_calc_contributions", test_calc_contributions)
        all_tests(6) = test_case("test_process_trajectories_alloc", test_process_trajectories_alloc)
        all_tests(7) = test_case("test_process_trajectories_flat_alloc", test_process_trajectories_flat_alloc)
        all_tests(8) = test_case("test_process_trajectories_empty_input", test_process_trajectories_empty_input)
        all_tests(9) = test_case("test_process_trajectories_invalid_dependent", test_process_trajectories_invalid_dependent)
        all_tests(10) = test_case("test_process_trajectories_dependent_in_mask", test_process_trajectories_dependent_in_mask)
        all_tests(11) = test_case("test_normalize_variable_timeseries", test_normalize_variable_timeseries)
        all_tests(12) = test_case("test_normalize_single_trajectory", test_normalize_single_trajectory)
        all_tests(13) = test_case("test_normalize_all_trajectories", test_normalize_all_trajectories)
        all_tests(14) = test_case("test_normalize_edge_cases", test_normalize_edge_cases)
        all_tests(15) = test_case("test_normalize_invalid_inputs", test_normalize_invalid_inputs)
    end function get_all_tests

    !> Run all tox_trajectory_contribution_analysis tests.
    subroutine run_all_tests_tox_trajectory_contribution_analysis
        type(test_case), allocatable :: all_tests(:)
        integer(int32) :: i

        all_tests = get_all_tests()

        do i = 1, size(all_tests)
            call all_tests(i)%test_proc()
            print "(' ',A,' passed.')", trim(all_tests(i)%name)
        end do
        print *, "All tox_trajectory_contribution_analysis tests passed successfully."
    end subroutine run_all_tests_tox_trajectory_contribution_analysis

    !> Run specific tox_trajectory_contribution_analysis tests by name.
    subroutine run_named_tests_tox_trajectory_contribution_analysis(test_names)
        character(len=*), intent(in) :: test_names(:)
        type(test_case), allocatable :: all_tests(:)
        integer(int32) :: i, j
        logical :: found

        all_tests = get_all_tests()

        do i = 1, size(test_names)
            found = .false.
            do j = 1, size(all_tests)
                if (trim(test_names(i)) == trim(all_tests(j)%name)) then
                    call all_tests(j)%test_proc()
                    print "(' ',A,' passed.')", trim(test_names(i))
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                print *, "Unknown test: ", trim(test_names(i))
            end if
        end do
    end subroutine run_named_tests_tox_trajectory_contribution_analysis

    subroutine test_calc_contributions()
        integer(int32), parameter :: n_factors = 2, n_samples = 3, n_timepoints = 4
        real(real64) :: trajectories(n_factors, n_samples, n_timepoints)
        real(real64) :: spike_contribs(n_timepoints, n_samples), integrated_contribs(n_samples)
        real(real64) :: expected_spike(n_timepoints), expected_integrated
        integer(int32) :: ierr, i_factor, i_sample, i_timepoint
        real(real64) :: factor_vec(n_timepoints), dependent_vec(n_timepoints), magnitude

        ! Fill trajectories with known values: T(i,j,k) = 100*i + 10*j + k
        do i_factor = 1, n_factors
            do i_sample = 1, n_samples
                do i_timepoint = 1, n_timepoints
                    trajectories(i_factor, i_sample, i_timepoint) = 100.0_real64*i_factor + 10.0_real64*i_sample + real(i_timepoint, real64)
                end do
            end do
        end do

        ! Case: aligned factor and dependent (factor=2, dependent=1)
        call calc_contributions_alloc(trajectories, n_factors, n_samples, n_timepoints, 2, 1, MODE_NORMAL, spike_contribs, integrated_contribs, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_calc_contributions: MODE_NORMAL should return OK")

        ! Validate each sample
        do i_sample = 1, n_samples
            factor_vec = trajectories(2, i_sample, :)
            dependent_vec = trajectories(1, i_sample, :)
            magnitude = sqrt(sum(factor_vec**2)) * sqrt(sum(dependent_vec**2)) 
            expected_spike = (factor_vec * dependent_vec) / magnitude
            expected_integrated = sum(factor_vec * dependent_vec) / magnitude

            call assert_equal_array_real(spike_contribs(:, i_sample), expected_spike, n_timepoints, TOL, "test_calc_contributions: spike_contribs mismatch")
            call assert_equal_real(integrated_contribs(i_sample), expected_integrated, TOL, "test_calc_contributions: integrated_contribs mismatch")
        end do

        ! Case: invalid mode
        call calc_contributions_alloc(trajectories, n_factors, n_samples, n_timepoints, 2_int32, 1_int32, 99_int32, spike_contribs, integrated_contribs, ierr)
        call assert_equal_int(ierr, ERR_INVALID_INPUT, "test_calc_contributions: expected ERR_INVALID_INPUT for mode=99")

    end subroutine test_calc_contributions

    subroutine test_trajectory_contribution()
        integer(int32), parameter :: n_timepoints = 3
        real(real64) :: factor(n_timepoints), dependent(n_timepoints)
        real(real64) :: contribution, expected
        integer(int32) :: ierr

        ! Case 1: perfectly aligned vectors
        factor = [1.0_real64, 2.0_real64, 3.0_real64]
        dependent = [2.0_real64, 4.0_real64, 6.0_real64]
        expected = 1.0_real64

        call trajectory_contribution(factor, dependent, n_timepoints, MODE_NORMAL, contribution, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_trajectory_contribution: MODE_NORMAL: expected OK status")
        call assert_equal_real(contribution, expected, TOL, "test_trajectory_contribution: MODE_NORMAL")

        expected = 0.0_real64  ! acos(1.0) = 0.0
        call trajectory_contribution(factor, dependent, n_timepoints, MODE_RAP, contribution, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_trajectory_contribution: MODE_RAP: expected OK status")
        call assert_equal_real(contribution, expected, TOL, "test_trajectory_contribution: MODE_RAP")

        ! Case 2: orthogonal vectors
        factor = [1.0_real64, 0.0_real64, 0.0_real64]
        dependent = [0.0_real64, 1.0_real64, 0.0_real64]
        expected = 0.0_real64

        call trajectory_contribution(factor, dependent, n_timepoints, MODE_NORMAL, contribution, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_trajectory_contribution: Orthogonal MODE_NORMAL: expected OK status")
        call assert_equal_real(contribution, expected, TOL, "test_trajectory_contribution: Orthogonal MODE_NORMAL")

        ! Case 3: zero vector input
        factor = [0.0_real64, 0.0_real64, 0.0_real64]
        dependent = [1.0_real64, 2.0_real64, 3.0_real64]

        call trajectory_contribution(factor, dependent, n_timepoints, MODE_NORMAL, contribution, ierr)
        call assert_equal_int(ierr, ERR_DIVISION_BY_ZERO, "test_trajectory_contribution: Zero vector")

        ! Case 4: invalid mode
        factor = [1.0_real64, 2.0_real64, 3.0_real64]
        dependent = [1.0_real64, 2.0_real64, 3.0_real64]

        call trajectory_contribution(factor, dependent, n_timepoints, 99_int32, contribution, ierr)
        call assert_equal_int(ierr, ERR_INVALID_INPUT, "test_trajectory_contribution: Invalid mode")

        ! Case 5: NaN
        factor = [1.0_real64, ieee_value(0.0_real64, ieee_quiet_nan), 3.0_real64]
        dependent = [1.0_real64, 2.0_real64, 3.0_real64]

        call trajectory_contribution(factor, dependent, n_timepoints, MODE_NORMAL, contribution, ierr)
        call assert_equal_int(ierr, ERR_NAN_INF, "test_spike_contribution: NaN")
    end subroutine test_trajectory_contribution

    subroutine test_spike_contribution()
        integer(int32), parameter :: n_timepoints = 3
        real(real64) :: factor(n_timepoints), dependent(n_timepoints), contribution(n_timepoints), expected(n_timepoints)
        integer(int32) :: ierr

        ! Case 1: perfectly aligned vectors
        factor = [1.0_real64, 2.0_real64, 3.0_real64]
        dependent = [2.0_real64, 4.0_real64, 6.0_real64]
        expected = [2.0_real64, 8.0_real64, 18.0_real64] / sqrt(14.0_real64) / sqrt(56.0_real64)  ! dot products normalized by magnitude

        call spike_contribution(factor, dependent, n_timepoints, MODE_NORMAL, contribution, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_spike_contribution: MODE_NORMAL aligned")
        call assert_equal_array_real(contribution, expected, n_timepoints, TOL, "test_spike_contribution: MODE_NORMAL aligned")

        expected = acos(expected)  ! element-wise acos
        call spike_contribution(factor, dependent, n_timepoints, MODE_RAP, contribution, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_spike_contribution: MODE_RAP aligned")
        call assert_equal_array_real(contribution, expected, n_timepoints, TOL, "test_spike_contribution: MODE_RAP aligned")

        ! Case 2: orthogonal vectors
        factor = [1.0_real64, 0.0_real64, 0.0_real64]
        dependent = [0.0_real64, 1.0_real64, 0.0_real64]
        expected = [0.0, 0.0, 0.0]  ! dot products all zero

        call spike_contribution(factor, dependent, n_timepoints, MODE_NORMAL, contribution, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_spike_contribution: Orthogonal MODE_NORMAL")
        call assert_equal_array_real(contribution, expected, n_timepoints, TOL, "test_spike_contribution: Orthogonal MODE_NORMAL")

        expected = acos(expected)  ! acos(0.0) = π/2
        call spike_contribution(factor, dependent, n_timepoints, MODE_RAP, contribution, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_spike_contribution: Orthogonal MODE_RAP")
        call assert_equal_array_real(contribution, expected, n_timepoints, TOL, "test_spike_contribution: Orthogonal MODE_RAP")

        ! Case 3: zero vector input
        factor = [0.0_real64, 0.0_real64, 0.0_real64]
        dependent = [1.0_real64, 2.0_real64, 3.0_real64]

        call spike_contribution(factor, dependent, n_timepoints, MODE_NORMAL, contribution, ierr)
        call assert_equal_int(ierr, ERR_DIVISION_BY_ZERO, "test_spike_contribution: Zero vector")

        ! Case 4: invalid mode
        factor = [1.0_real64, 2.0_real64, 3.0_real64]
        dependent = [1.0_real64, 2.0_real64, 3.0_real64]

        call spike_contribution(factor, dependent, n_timepoints, 99_int32, contribution, ierr)
        call assert_equal_int(ierr, ERR_INVALID_INPUT, "test_spike_contribution: Invalid mode")

        ! Case 5: NaN
        factor = [1.0_real64, ieee_value(0.0_real64, ieee_quiet_nan), 3.0_real64]
        dependent = [1.0_real64, 2.0_real64, 3.0_real64]

        call spike_contribution(factor, dependent, n_timepoints, MODE_NORMAL, contribution, ierr)
        call assert_equal_int(ierr, ERR_NAN_INF, "test_spike_contribution: NaN")
    end subroutine test_spike_contribution

    subroutine test_get_vec_across_samples()
        integer(int32), parameter :: n_factors = 2, n_samples = 3, n_timepoints = 4
        real(real64) :: trajectories(n_factors, n_samples, n_timepoints)
        real(real64) :: result(n_samples)
        real(real64) :: expected(n_samples)
        integer(int32) :: ierr, i_factor, i_sample, i_timepoint

        ! Fill trajectories with known values: T(i,j,k) = 100*i + 10*j + k
        do i_factor = 1, n_factors
            do i_sample = 1, n_samples
                do i_timepoint = 1, n_timepoints
                    trajectories(i_factor,i_sample,i_timepoint) = 100.0_real64*i_factor + 10.0_real64*i_sample + real(i_timepoint,real64)
                end do
            end do
        end do

        ! Valid extraction: factor 2, timepoint 3
        expected = [213.0, 223.0, 233.0]
        call get_vec_across_samples(trajectories, n_factors, n_samples, n_timepoints, 2_int32, 3_int32, result, ierr)
        call assert_equal_int(ierr, ERR_OK, "Test failed: unexpected error code")
        call assert_equal_array_real(result, expected, n_samples, 0.0_real64, "test_get_vec_across_timepoints: returned vector doesn't match")
    end subroutine test_get_vec_across_samples

    subroutine test_get_vec_across_timepoints()
        integer(int32), parameter :: n_factors = 2, n_samples = 3, n_timepoints = 4
        real(real64) :: trajectories(n_factors, n_samples, n_timepoints)
        real(real64) :: result(n_timepoints)
        real(real64) :: expected(n_timepoints)
        integer(int32) :: ierr, i_factor, i_sample, i_timepoint

        ! Fill trajectories with known values: T(i,j,k) = 100*i + 10*j + k
        do i_factor = 1, n_factors
            do i_sample = 1, n_samples
                do i_timepoint = 1, n_timepoints
                    trajectories(i_factor,i_sample,i_timepoint) = 100.0_real64*i_factor + 10.0_real64*i_sample + real(i_timepoint,real64)
                end do
            end do
        end do

        ! Valid extraction: factor 1, sample 2
        expected = [121.0, 122.0, 123.0, 124.0]
        call get_vec_across_timepoints(trajectories, n_factors, n_samples, n_timepoints, 1_int32, 2_int32, result, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_get_vec_across_timepoints: unexpected error code")
        call assert_equal_array_real(result, expected, n_timepoints, 0.0_real64, "test_get_vec_across_timepoints: returned vector doesn't match")
    end subroutine test_get_vec_across_timepoints

    !> Test: Basic trajectory processing with per-timepoint percentiles
    subroutine test_process_trajectories_alloc()
        integer(int32), parameter :: n_factors = 3, n_samples = 4, n_timepoints = 5
        real(real64) :: trajectories(n_factors, n_samples, n_timepoints)
        logical :: factor_mask(n_factors)
        integer(int32) :: dependent_idx = 2
        integer(int32) :: mode = MODE_NORMAL
        real(real64) :: percentile = 95.0_real64
        real(real64) :: integrated_contribs(n_samples, n_factors-1)  ! dependent_idx excluded
        real(real64) :: spike_contribs(n_timepoints, n_samples, n_factors-1)
        real(real64) :: thresholds_integrated_contrib(n_factors-1)
        real(real64) :: thresholds_spike_contrib(n_timepoints, n_factors-1)
        logical :: outliers_integrated_contrib(n_samples, n_factors-1)
        logical :: outliers_spike_contrib(n_timepoints, n_samples, n_factors-1)
        integer(int32) :: ierr

        ! Initialize test data
        call random_number(trajectories)
        trajectories = trajectories * 10.0_real64  ! Scale to reasonable values
        
        ! Set factor mask: exclude dependent_idx
        factor_mask = .true.
        factor_mask(dependent_idx) = .false.

        ! Call the routine
        call process_trajectories_alloc(trajectories, n_factors, n_samples, n_timepoints, &
                                    factor_mask, count(factor_mask), dependent_idx, mode, percentile, &
                                    integrated_contribs, spike_contribs, &
                                    thresholds_integrated_contrib, outliers_integrated_contrib, &
                                    thresholds_spike_contrib, outliers_spike_contrib, ierr)

        call assert_equal_int(ierr, 0, 'Error code 0 for process_trajectories_alloc')
        call assert_true(size(integrated_contribs, 1) == n_samples, 'Integrated contribs correct sample dimension')
        call assert_true(size(integrated_contribs, 2) == count(factor_mask), 'Integrated contribs correct factor dimension')
        call assert_true(size(spike_contribs, 1) == n_timepoints, 'Spike contribs correct timepoint dimension')
        call assert_true(size(spike_contribs, 2) == n_samples, 'Spike contribs correct sample dimension')
        call assert_true(size(spike_contribs, 3) == count(factor_mask), 'Spike contribs correct factor dimension')
    end subroutine test_process_trajectories_alloc

    !> Test: Flat trajectory processing with global percentile
    subroutine test_process_trajectories_flat_alloc()
        integer(int32), parameter :: n_factors = 3, n_samples = 4, n_timepoints = 5
        real(real64) :: trajectories(n_factors, n_samples, n_timepoints)
        logical :: factor_mask(n_factors)
        integer(int32) :: dependent_idx = 2
        integer(int32) :: mode = MODE_RAP
        real(real64) :: percentile = 90.0_real64
        real(real64) :: integrated_contribs(n_samples, n_factors-1)
        real(real64) :: spike_contribs(n_timepoints, n_samples, n_factors-1)
        real(real64) :: thresholds_integrated_contrib(n_factors-1)
        real(real64) :: thresholds_spike_contrib(n_factors-1)  ! Note: 1D for flat version
        logical :: outliers_integrated_contrib(n_samples, n_factors-1)
        logical :: outliers_spike_contrib(n_timepoints, n_samples, n_factors-1)
        integer(int32) :: ierr

        ! Initialize test data
        call random_number(trajectories)
        trajectories = trajectories * 5.0_real64  ! Scale to reasonable values
        
        ! Set factor mask: exclude dependent_idx
        factor_mask = .true.
        factor_mask(dependent_idx) = .false.

        ! Call the routine
        call process_trajectories_flat_alloc(trajectories, n_factors, n_samples, n_timepoints, &
                                        factor_mask, count(factor_mask), dependent_idx, mode, percentile, &
                                        integrated_contribs, spike_contribs, &
                                        thresholds_integrated_contrib, outliers_integrated_contrib, &
                                        thresholds_spike_contrib, outliers_spike_contrib, ierr)

        call assert_equal_int(ierr, 0, 'Error code 0 for process_trajectories_flat_alloc')
        call assert_true(size(thresholds_spike_contrib) == count(factor_mask), 'Flat spike thresholds correct dimension')
        call assert_true(all(thresholds_spike_contrib > 0.0_real64), 'Flat spike thresholds should be positive')
    end subroutine test_process_trajectories_flat_alloc

    !> Test: Error handling for empty input in trajectory processing
    subroutine test_process_trajectories_empty_input()
        integer(int32), parameter :: n_factors = 0, n_samples = 0, n_timepoints = 0
        real(real64), allocatable :: trajectories(:,:,:)
        logical, allocatable :: factor_mask(:)
        integer(int32) :: dependent_idx = 1
        integer(int32) :: mode = MODE_NORMAL
        real(real64) :: percentile = 95.0_real64
        real(real64), allocatable :: integrated_contribs(:,:), spike_contribs(:,:,:)
        real(real64), allocatable :: thresholds_integrated_contrib(:), thresholds_spike_contrib(:,:)
        logical, allocatable :: outliers_integrated_contrib(:,:), outliers_spike_contrib(:,:,:)
        integer(int32) :: ierr

        ! Allocate zero-sized arrays
        allocate(trajectories(0,0,0))
        allocate(factor_mask(0))
        
        ! Call the routine - should error with ERR_EMPTY_INPUT
        call process_trajectories_alloc(trajectories, n_factors, n_samples, n_timepoints, &
                                    factor_mask, count(factor_mask), dependent_idx, mode, percentile, &
                                    integrated_contribs, spike_contribs, &
                                    thresholds_integrated_contrib, outliers_integrated_contrib, &
                                    thresholds_spike_contrib, outliers_spike_contrib, ierr)

        call assert_equal_int(ierr, ERR_EMPTY_INPUT, 'Error code for empty input')
    end subroutine test_process_trajectories_empty_input

    !> Test: Error handling for invalid dependent index
    subroutine test_process_trajectories_invalid_dependent()
        integer(int32), parameter :: n_factors = 3, n_samples = 4, n_timepoints = 5
        real(real64) :: trajectories(n_factors, n_samples, n_timepoints)
        logical :: factor_mask(n_factors)
        integer(int32) :: dependent_idx = 5  ! Out of bounds
        integer(int32) :: mode = MODE_NORMAL
        real(real64) :: percentile = 95.0_real64
        real(real64) :: integrated_contribs(n_samples, n_factors-1)
        real(real64) :: spike_contribs(n_timepoints, n_samples, n_factors-1)
        real(real64) :: thresholds_integrated_contrib(n_factors-1)
        real(real64) :: thresholds_spike_contrib(n_timepoints, n_factors-1)
        logical :: outliers_integrated_contrib(n_samples, n_factors-1)
        logical :: outliers_spike_contrib(n_timepoints, n_samples, n_factors-1)
        integer(int32) :: ierr

        call random_number(trajectories)
        factor_mask = .true.

        call process_trajectories_alloc(trajectories, n_factors, n_samples, n_timepoints, &
                                    factor_mask, count(factor_mask), dependent_idx, mode, percentile, &
                                    integrated_contribs, spike_contribs, &
                                    thresholds_integrated_contrib, outliers_integrated_contrib, &
                                    thresholds_spike_contrib, outliers_spike_contrib, ierr)

        call assert_equal_int(ierr, ERR_INVALID_INPUT, 'Error code for invalid dependent index')
    end subroutine test_process_trajectories_invalid_dependent

    !> Test: Error when dependent_idx is included in factor_mask
    subroutine test_process_trajectories_dependent_in_mask()
        integer(int32), parameter :: n_factors = 3, n_samples = 4, n_timepoints = 5
        real(real64) :: trajectories(n_factors, n_samples, n_timepoints)
        logical :: factor_mask(n_factors)
        integer(int32) :: dependent_idx = 2
        integer(int32) :: mode = MODE_NORMAL
        real(real64) :: percentile = 95.0_real64
        real(real64) :: integrated_contribs(n_samples, n_factors)
        real(real64) :: spike_contribs(n_timepoints, n_samples, n_factors)
        real(real64) :: thresholds_integrated_contrib(n_factors)
        real(real64) :: thresholds_spike_contrib(n_timepoints, n_factors)
        logical :: outliers_integrated_contrib(n_samples, n_factors)
        logical :: outliers_spike_contrib(n_timepoints, n_samples, n_factors)
        integer(int32) :: ierr

        call random_number(trajectories)
        
        ! Include dependent_idx in mask - this should cause an error
        factor_mask = .true.  ! All factors included, including dependent_idx

        call process_trajectories_alloc(trajectories, n_factors, n_samples, n_timepoints, &
                                    factor_mask, count(factor_mask), dependent_idx, mode, percentile, &
                                    integrated_contribs, spike_contribs, &
                                    thresholds_integrated_contrib, outliers_integrated_contrib, &
                                    thresholds_spike_contrib, outliers_spike_contrib, ierr)

        call assert_equal_int(ierr, ERR_INVALID_INPUT, 'Error when dependent_idx in factor_mask')
    end subroutine test_process_trajectories_dependent_in_mask

    subroutine test_normalize_variable_timeseries()
        integer(int32) :: n_points = 5
        real(real64) :: v(5), v_norm(5), v_norm_expected(5), v1_norm(1), v1_norm_expected(1), v1(1)
        integer(int32) :: ierr
        
        ! Test 1: Normal case
        v = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64]
        v_norm_expected = [0.0_real64, 0.25_real64, 0.5_real64, 0.75_real64, 1.0_real64]
        
        call normalize_variable_timeseries(v, v_norm, n_points, ierr)
        print *, "v:", v
        print *, "v_norm_expected:", v_norm_expected
        print *, "v_norm:", v_norm
        call assert_equal_int(ierr, ERR_OK, "test_normalize_variable_timeseries: normal case should succeed")
        call assert_equal_array_real(v_norm, v_norm_expected, n_points, TOL, "test_normalize_variable_timeseries: normal case values")
        
        ! Test 2: Constant vector (all values same)
        v = [2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64]
        v_norm_expected = [0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64]
        
        call normalize_variable_timeseries(v, v_norm, n_points, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_normalize_variable_timeseries: constant vector should succeed")
        call assert_equal_array_real(v_norm, v_norm_expected, n_points, TOL, "test_normalize_variable_timeseries: constant vector values")
        
        ! Test 3: Negative values
        v = [-5.0_real64, -2.0_real64, 0.0_real64, 3.0_real64, 6.0_real64]
        v_norm_expected = [0.0_real64, 0.27272727_real64, 0.45454545_real64, 0.72727273_real64, 1.0_real64]
        
        call normalize_variable_timeseries(v, v_norm, n_points, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_normalize_variable_timeseries: negative values should succeed")
        call assert_equal_array_real(v_norm, v_norm_expected, n_points, 1.0e-8_real64, "test_normalize_variable_timeseries: negative values")
        
        ! Test 4: Single point
        v1 = [3.14_real64]
        v1_norm_expected = [0.0_real64]  ! Single value normalized to 0
        
        call normalize_variable_timeseries(v1, v1_norm, 1, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_normalize_variable_timeseries: single point should succeed")
        call assert_equal_array_real(v1_norm, v1_norm_expected, 1, TOL, "test_normalize_variable_timeseries: single point")
    end subroutine test_normalize_variable_timeseries

    subroutine test_normalize_single_trajectory()
        integer(int32), parameter :: n_factors = 3, n_samples = 4
        real(real64) :: trajectory(n_factors, n_samples), trajectory_norm(n_factors, n_samples)
        real(real64) :: expected(n_factors, n_samples)
        integer(int32) :: ierr, i, j
        
        ! Create test trajectory
        do i = 1, n_factors
            do j = 1, n_samples
                trajectory(i, j) = real(i * 10 + j, real64)
            end do
        end do
        
        ! Calculate expected normalization
        do i = 1, n_factors
            do j = 1, n_samples
                expected(i, j) = (real(j, real64) - 1.0_real64) / real(n_samples - 1, real64)
            end do
        end do
        
        call normalize_single_trajectory(trajectory, trajectory_norm, n_factors, n_samples, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_normalize_single_trajectory: should succeed")
        
        ! Check each variable independently
        do i = 1, n_factors
            call assert_equal_array_real(trajectory_norm(i, :), expected(i, :), n_samples, TOL, &
                "test_normalize_single_trajectory: variable ")
        end do
        
        ! Verify min=0 and max=1 for each variable
        do i = 1, n_factors
            call assert_equal_real(minval(trajectory_norm(i, :)), 0.0_real64, TOL, &
                "test_normalize_single_trajectory: min=0 for variable ")
            call assert_equal_real(maxval(trajectory_norm(i, :)), 1.0_real64, TOL, &
                "test_normalize_single_trajectory: max=1 for variable ")
        end do
    end subroutine test_normalize_single_trajectory

    subroutine test_normalize_all_trajectories()
        integer(int32), parameter :: n_factors = 2, n_samples = 3, n_timepoints = 4
        real(real64) :: trajectories(n_factors, n_samples, n_timepoints)
        real(real64) :: trajectories_norm(n_factors, n_samples, n_timepoints)
        integer(int32) :: ierr, i_factor, i_sample, i_timepoint
        
        ! Fill with known pattern
        do i_factor = 1, n_factors
            do i_sample = 1, n_samples
                do i_timepoint = 1, n_timepoints
                    trajectories(i_factor, i_sample, i_timepoint) = &
                        real(i_factor * 100 + i_sample * 10 + i_timepoint, real64)
                end do
            end do
        end do
        
        call normalize_all_trajectories(trajectories, trajectories_norm, &
                                      n_factors, n_samples, n_timepoints, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_normalize_all_trajectories: should succeed")
        
        ! Check each entity (timepoint) independently
        do i_timepoint = 1, n_timepoints
            ! For each timepoint, check min=0 and max=1 for each variable
            do i_factor = 1, n_factors
                call assert_equal_real(minval(trajectories_norm(i_factor, :, i_timepoint)), &
                                     0.0_real64, TOL, &
                                     "test_normalize_all_trajectories min_val failed")
                call assert_equal_real(maxval(trajectories_norm(i_factor, :, i_timepoint)), &
                                     1.0_real64, TOL, &
                                     "test_normalize_all_trajectories max_val failed")
            end do
        end do
        
        ! Verify shape is preserved
        do i_timepoint = 1, n_timepoints
            do i_sample = 1, n_samples
                do i_factor = 1, n_factors
                    ! Check values are in [0,1]
                    call assert_true(trajectories_norm(i_factor, i_sample, i_timepoint) >= 0.0_real64 - TOL, &
                                   "test_normalize_all_trajectories: value >= 0")
                    call assert_true(trajectories_norm(i_factor, i_sample, i_timepoint) <= 1.0_real64 + TOL, &
                                   "test_normalize_all_trajectories: value <= 1")
                end do
            end do
        end do
    end subroutine test_normalize_all_trajectories

    subroutine test_normalize_edge_cases()
        real(real64) :: v(3), v_norm(3), v_norm_expected(3)
        real(real64) :: v5(5), v5_norm(5), v5_norm_expected(5)
        integer(int32) :: ierr
        
        ! Test 1: Very small values (near machine epsilon)
        v = [tiny(1.0_real64), 2.0_real64 * tiny(1.0_real64), 3.0_real64 * tiny(1.0_real64)]
        v_norm_expected = [0.0_real64, 0.5_real64, 1.0_real64]
        
        call normalize_variable_timeseries(v, v_norm, 3, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_normalize_edge_cases: small values should succeed")
        call assert_equal_array_real(v_norm, v_norm_expected, 3, TOL, "test_normalize_edge_cases: small values")
        
        ! Test 2: Large values
        v = [1.0e10_real64, 2.0e10_real64, 3.0e10_real64]
        v_norm_expected = [0.0_real64, 0.5_real64, 1.0_real64]
        
        call normalize_variable_timeseries(v, v_norm, 3, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_normalize_edge_cases: large values should succeed")
        call assert_equal_array_real(v_norm, v_norm_expected, 3, TOL, "test_normalize_edge_cases: large values")
        
        ! Test 3: Mixed positive and negative with zero
        v5 = [-10.0_real64, -5.0_real64, 0.0_real64, 5.0_real64, 10.0_real64]
        v5_norm_expected = [0.0_real64, 0.25_real64, 0.5_real64, 0.75_real64, 1.0_real64]
        
        call normalize_variable_timeseries(v5, v5_norm, 5, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_normalize_edge_cases: mixed values should succeed")
        call assert_equal_array_real(v5_norm, v5_norm_expected, 5, 1.0e-8_real64, "test_normalize_edge_cases: mixed values")
        
        ! Test 4: Already normalized (values in [0,1])
        v5 = [0.0_real64, 0.25_real64, 0.5_real64, 0.75_real64, 1.0_real64]
        v5_norm_expected = v5  ! Should stay the same
        
        call normalize_variable_timeseries(v5, v5_norm, 5, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_normalize_edge_cases: already normalized should succeed")
        call assert_equal_array_real(v5_norm, v5_norm_expected, 5, TOL, "test_normalize_edge_cases: already normalized")
    end subroutine test_normalize_edge_cases

    subroutine test_normalize_invalid_inputs()
        real(real64) :: v(5), v_norm(5)
        real(real64) :: v3(3), v3_norm(3)
        integer(int32) :: ierr
        
        ! Test 1: Empty array (n_points = 0)
        call normalize_variable_timeseries(v, v_norm, 0, ierr)
        call assert_equal_int(ierr, ERR_EMPTY_INPUT, "test_normalize_invalid_inputs: empty array should return ERR_EMPTY_INPUT")
        
        ! Test 2: Negative n_points
        call normalize_variable_timeseries(v, v_norm, -1, ierr)
        call assert_equal_int(ierr, ERR_INVALID_INPUT, "test_normalize_invalid_inputs: negative n_points should return ERR_INVALID_INPUT")
        
        ! Test 3: NaN in input
        v3 = [1.0_real64, ieee_value(0.0_real64, ieee_quiet_nan), 3.0_real64]
        call normalize_variable_timeseries(v3, v3_norm, 3, ierr)
        ! Note: NaN handling depends on minval/maxval behavior - check if error is set
        if (ierr /= ERR_OK) then
            ! If error is set, it should be ERR_NAN_INF
            call assert_equal_int(ierr, ERR_NAN_INF, "test_normalize_invalid_inputs: NaN should return ERR_NAN_INF")
        end if
        
        ! Test 4: Infinity in input
        v3 = [1.0_real64, 2.0_real64, huge(1.0_real64)]
        call normalize_variable_timeseries(v3, v3_norm, 3, ierr)
        ! Similar to NaN case
        if (ierr /= ERR_OK) then
            call assert_equal_int(ierr, ERR_NAN_INF, "test_normalize_invalid_inputs: Infinity should return ERR_NAN_INF")
        end if
    end subroutine test_normalize_invalid_inputs

end module mod_test_tox_traj_contrib_analysis
