!> Unit test suite for tox_trajectory_contribution_analysis routine.
module mod_test_tox_trajectory_contribution_analysis
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use tox_trajectory_contribution_analysis
    use tox_errors
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
        type(test_case) :: all_tests(5)

        all_tests(1) = test_case("test_tox_trajectory_contribution_analysis_get_vec_across_samples", test_get_vec_across_samples)
        all_tests(2) = test_case("test_tox_trajectory_contribution_analysis_get_vec_across_timepoints", test_get_vec_across_timepoints)
        all_tests(3) = test_case("test_tox_trajectory_contribution_analysis_trajectory_contribution", test_trajectory_contribution)
        all_tests(4) = test_case("test_tox_trajectory_contribution_analysis_spike_contribution", test_spike_contribution)
        all_tests(5) = test_case("test_tox_trajectory_contribution_analysis_calc_contributions", test_calc_contributions)
    end function get_all_tests

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
end module mod_test_tox_trajectory_contribution_analysis
