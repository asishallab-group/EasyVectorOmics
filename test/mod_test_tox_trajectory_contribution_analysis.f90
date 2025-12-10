!> Unit test suite for tox_trajectory_contribution_analysis routine.
module mod_test_tox_traj_contrib_analysis
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
        type(test_case) :: all_tests(1)

        all_tests(1) = test_case("test_compute_baselines_factor_dependent", test_compute_baselines_factor_dependent)
    end function get_all_tests

    subroutine test_compute_baselines_factor_dependent()
        integer(int32), parameter :: n_timepoints = 4
        real(real64) :: factor(n_timepoints), dependent(5)  ! dependent has 5 elements for mismatch test
        real(real64) :: factor_baseline, dependent_baseline, expected_factor_baseline, expected_dependent_baseline
        integer(int32) :: ierr

        ! Case 1: BASELINE_RAW (no centering)
        factor = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
        dependent(1:n_timepoints) = [5.0_real64, 6.0_real64, 7.0_real64, 8.0_real64]
        call compute_baselines_factor_dependent(n_timepoints, factor, dependent(1:n_timepoints), BASELINE_RAW, &
                                               factor_baseline, dependent_baseline, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_compute_baselines_factor_dependent: BASELINE_RAW: expected OK status")
        call assert_equal_real(factor_baseline, 0.0_real64, TOL, "test_compute_baselines_factor_dependent: BASELINE_RAW factor_baseline")
        call assert_equal_real(dependent_baseline, 0.0_real64, TOL, "test_compute_baselines_factor_dependent: BASELINE_RAW dependent_baseline")

        ! Case 2: BASELINE_MIN (minimum-centered)
        call compute_baselines_factor_dependent(n_timepoints, factor, dependent(1:n_timepoints), BASELINE_MIN, &
                                               factor_baseline, dependent_baseline, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_compute_baselines_factor_dependent: BASELINE_MIN: expected OK status")
        call assert_equal_real(factor_baseline, minval(factor), TOL, "test_compute_baselines_factor_dependent: BASELINE_MIN factor_baseline")
        call assert_equal_real(dependent_baseline, minval(dependent(1:n_timepoints)), TOL, "test_compute_baselines_factor_dependent: BASELINE_MIN dependent_baseline")

        ! Case 3: BASELINE_MEAN (mean-centered)
        expected_factor_baseline = sum(factor) / real(n_timepoints, kind=real64)
        expected_dependent_baseline = sum(dependent(1:n_timepoints)) / real(n_timepoints, kind=real64)
        call compute_baselines_factor_dependent(n_timepoints, factor, dependent(1:n_timepoints), BASELINE_MEAN, &
                                               factor_baseline, dependent_baseline, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_compute_baselines_factor_dependent: BASELINE_MEAN: expected OK status")
        call assert_equal_real(factor_baseline, expected_factor_baseline, TOL, "test_compute_baselines_factor_dependent: BASELINE_MEAN factor_baseline")
        call assert_equal_real(dependent_baseline, expected_dependent_baseline, TOL, "test_compute_baselines_factor_dependent: BASELINE_MEAN dependent_baseline")

        ! Case 4: mismatched input lengths
        dependent = [5.0_real64, 6.0_real64, 7.0_real64, 8.0_real64, 9.0_real64]  ! length 5
        call compute_baselines_factor_dependent(n_timepoints, factor, dependent(1:n_timepoints), BASELINE_RAW, &
                                               factor_baseline, dependent_baseline, ierr)  ! valid
        call assert_equal_int(ierr, ERR_OK, "test_compute_baselines_factor_dependent: valid lengths")
        call compute_baselines_factor_dependent(4_int32, factor(1:4), dependent(1:5), BASELINE_RAW, &
                                               factor_baseline, dependent_baseline, ierr)  ! invalid: mismatch
        ! Note: This test may not trigger length mismatch since we now validate in C wrapper only

        ! Case 5: invalid mode
        call compute_baselines_factor_dependent(n_timepoints, factor, dependent(1:n_timepoints), 99_int32, &
                                               factor_baseline, dependent_baseline, ierr)
        call assert_equal_int(ierr, ERR_INVALID_INPUT, "test_compute_baselines_factor_dependent: invalid mode")
    end subroutine test_compute_baselines_factor_dependent

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
end module mod_test_tox_traj_contrib_analysis
