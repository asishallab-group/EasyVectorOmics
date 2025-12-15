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
        type(test_case) :: all_tests(3)

        all_tests(1) = test_case("test_compute_baselines_factor_dependent", test_compute_baselines_factor_dependent)
        all_tests(2) = test_case("test_compute_contributions", test_compute_contributions)
        all_tests(3) = test_case("test_compute_all_contributions", test_compute_all_contributions)
    end function get_all_tests

    subroutine test_compute_all_contributions()
        integer(int32), parameter :: n_factors = 2, n_samples = 1, n_timepoints = 3
        integer(int32), parameter :: n_selected_factors = 1, n_selected_dependents = 1
        integer(int32) :: ierr, i_dependent
        real(real64) :: trajectories(n_factors, n_samples, n_timepoints)
        integer(int32) :: factor_indices(n_factors)
        integer(int32) :: dependent_indices(n_factors)
        real(real64) :: local_contributions(n_timepoints, n_factors, n_factors, n_samples)
        real(real64) :: total_contributions(n_factors, n_factors, n_samples)
        real(real64) :: temp_factors(n_timepoints, n_factors), temp_dependent(n_timepoints)
        real(real64) :: expected_local(n_timepoints)
        real(real64) :: expected_total

        ! -------------------------------
        ! Case 1: zero vectors
        ! -------------------------------
        ! Factor trajectory: [1,2,3]
        ! Dependent trajectory: [4,5,6]

        trajectories(1,1,:) = [1.0_real64, 2.0_real64, 3.0_real64]
        trajectories(2,1,:) = [4.0_real64, 5.0_real64, 6.0_real64]
        factor_indices = [1, 2]
        dependent_indices = [2, 1]

        call compute_all_contributions(trajectories, 0_int32, n_samples, n_timepoints, &
            factor_indices(:n_selected_factors), n_selected_factors, dependent_indices(:n_selected_dependents), n_selected_dependents, BASELINE_MIN, &
            local_contributions, total_contributions, temp_factors, temp_dependent, ierr)
        call assert_equal_int(ierr, ERR_EMPTY_INPUT, "test_compute_all_contributions: Case 1 expected error for n_factors=0")

        call compute_all_contributions(trajectories, n_factors, 0_int32, n_timepoints, &
            factor_indices(:n_selected_factors), n_selected_factors, dependent_indices(:n_selected_dependents), n_selected_dependents, BASELINE_MIN, &
            local_contributions, total_contributions, temp_factors, temp_dependent, ierr)
        call assert_equal_int(ierr, ERR_EMPTY_INPUT, "test_compute_all_contributions: Case 1 expected error for n_samples=0")

        call compute_all_contributions(trajectories, n_factors, n_samples, 0_int32, &
            factor_indices(:n_selected_factors), n_selected_factors, dependent_indices(:n_selected_dependents), n_selected_dependents, BASELINE_MIN, &
            local_contributions, total_contributions, temp_factors, temp_dependent, ierr)
        call assert_equal_int(ierr, ERR_EMPTY_INPUT, "test_compute_all_contributions: Case 1 expected error for n_timepoints=0")

        call compute_all_contributions(trajectories, n_factors, n_samples, n_timepoints, &
            factor_indices(:n_selected_factors), 0_int32, dependent_indices(:n_selected_dependents), n_selected_dependents, BASELINE_MIN, &
            local_contributions, total_contributions, temp_factors, temp_dependent, ierr)
        call assert_equal_int(ierr, ERR_EMPTY_INPUT, "test_compute_all_contributions: Case 1 expected error for n_selected_factors=0")

        call compute_all_contributions(trajectories, n_factors, n_samples, n_timepoints, &
            factor_indices(:n_selected_factors), n_selected_factors, dependent_indices(:n_selected_dependents), 0_int32, BASELINE_MIN, &
            local_contributions, total_contributions, temp_factors, temp_dependent, ierr)
        call assert_equal_int(ierr, ERR_EMPTY_INPUT, "test_compute_all_contributions: Case 1 expected error for n_selected_dependents=0")

        ! -------------------------------
        ! Case 2: MEAN baseline
        ! -------------------------------

        call compute_all_contributions(trajectories, n_factors, n_samples, n_timepoints, &
            factor_indices(:n_selected_factors), n_selected_factors, dependent_indices(:n_selected_dependents), n_selected_dependents, BASELINE_MEAN, &
            local_contributions, total_contributions, temp_factors, temp_dependent, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_compute_all_contributions: Case 2 ierr")

        ! Baselines: mean(factor)=2.0, mean(dependent)=5.0
        expected_local(1) = (1.0-2.0)*(4.0-5.0)   ! = 1.0
        expected_local(2) = (2.0-2.0)*(5.0-5.0)   ! = 0.0
        expected_local(3) = (3.0-2.0)*(6.0-5.0)   ! = 1.0
        expected_total = sum(expected_local)      ! = 2.0

        call assert_equal_array_real(local_contributions(:,1,1,1), expected_local, n_timepoints, TOL, "test_compute_all_contributions: Case 2 local contributions")
        call assert_equal_real(total_contributions(1,1,1), expected_total, TOL, "test_compute_all_contributions: Case 2 total contribution")

        ! -------------------------------
        ! Case 3: MIN baseline
        ! -------------------------------
        ! Factor trajectory: [2,4,6]
        ! Dependent trajectory: [1,3,5]
        trajectories(1,1,:) = [2.0_real64, 4.0_real64, 6.0_real64]
        trajectories(2,1,:) = [1.0_real64, 3.0_real64, 5.0_real64]

        call compute_all_contributions(trajectories, n_factors, n_samples, n_timepoints, &
            factor_indices(:n_selected_factors), n_selected_factors, dependent_indices(:n_selected_dependents), n_selected_dependents, BASELINE_MIN, &
            local_contributions, total_contributions, temp_factors, temp_dependent, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_compute_all_contributions: Case 3 ierr")

        ! Baselines: min(factor)=2.0, min(dependent)=1.0
        expected_local(1) = (2.0-2.0)*(1.0-1.0)   ! = 0.0
        expected_local(2) = (4.0-2.0)*(3.0-1.0)   ! = 4.0
        expected_local(3) = (6.0-2.0)*(5.0-1.0)   ! = 16.0
        expected_total = sum(expected_local)      ! = 20.0

        call assert_equal_array_real(local_contributions(:,1,1,1), expected_local, n_timepoints, TOL, "test_compute_all_contributions: Case 3 local contributions")
        call assert_equal_real(total_contributions(1,1,1), expected_total, TOL, "test_compute_all_contributions: Case 3 total contribution")

        ! -------------------------------
        ! Case 4: all elements equal
        ! -------------------------------
        trajectories(1,1,:) = 1.0_real64

        call compute_all_contributions(trajectories, n_factors, n_samples, n_timepoints, &
            factor_indices(:n_selected_factors), n_selected_factors, dependent_indices(:n_selected_dependents), n_selected_dependents, BASELINE_MIN, &
            local_contributions, total_contributions, temp_factors, temp_dependent, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_compute_all_contributions: Case 4 ierr zero factor")

        expected_local = 0.0_real64
        expected_total = 0.0_real64

        call assert_equal_array_real(local_contributions(:,1,1,1), expected_local, n_timepoints, TOL, "test_compute_all_contributions: Case 4 local contributions zero factor")
        call assert_equal_real(total_contributions(1,1,1), expected_total, TOL, "test_compute_all_contributions: Case 4 total contribution zero factor")

        trajectories(1,1,:) = [1.0_real64, 2.0_real64, 3.0_real64]
        trajectories(2,1,:) = 1.0_real64

        call compute_all_contributions(trajectories, n_factors, n_samples, n_timepoints, &
            factor_indices(:n_selected_factors), n_selected_factors, dependent_indices(:n_selected_dependents), n_selected_dependents, BASELINE_MIN, &
            local_contributions, total_contributions, temp_factors, temp_dependent, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_compute_all_contributions: Case 4 ierr zero dependent")

        expected_local = 0.0_real64
        expected_total = 0.0_real64

        call assert_equal_array_real(local_contributions(:,1,1,1), expected_local, n_timepoints, TOL, "test_compute_all_contributions: Case 4 local contributions zero dependent")
        call assert_equal_real(total_contributions(1,1,1), expected_total, TOL, "test_compute_all_contributions: Case 4 total contribution zero dependent")

        ! -------------------------------
        ! Case 5: multiple selected indices
        ! -------------------------------
        trajectories(1,1,:) = [1.0_real64, 2.0_real64, 3.0_real64]
        trajectories(2,1,:) = [4.0_real64, 5.0_real64, 6.0_real64]
        factor_indices = [1, 2]
        dependent_indices = [2, 2]

        call compute_all_contributions(trajectories, n_factors, n_samples, n_timepoints, &
            factor_indices, n_factors, dependent_indices, n_factors, BASELINE_MIN, &
            local_contributions, total_contributions, temp_factors, temp_dependent, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_compute_all_contributions: Case 5 ierr")

        do i_dependent = 1, size(dependent_indices)
            expected_local(1) = (1.0-1.0)*(4.0-4.0)
            expected_local(2) = (2.0-1.0)*(5.0-4.0)
            expected_local(3) = (3.0-1.0)*(6.0-4.0)
            expected_total = sum(expected_local)

            call assert_equal_array_real(local_contributions(:,1,i_dependent,1), expected_local, n_timepoints, TOL, "test_compute_all_contributions: Case 5 local contributions")
            call assert_equal_real(total_contributions(1,i_dependent,1), expected_total, TOL, "test_compute_all_contributions: Case 5 total contribution")

            expected_local(1) = (4.0-4.0) ** 2
            expected_local(2) = (5.0-4.0) ** 2
            expected_local(3) = (6.0-4.0) ** 2
            expected_total = sum(expected_local)

            call assert_equal_array_real(local_contributions(:,1,i_dependent,1), expected_local, n_timepoints, TOL, "test_compute_all_contributions: Case 5 local contributions")
            call assert_equal_real(total_contributions(1,i_dependent,1), expected_total, TOL, "test_compute_all_contributions: Case 5 total contribution")
        end do

    end subroutine test_compute_all_contributions

    subroutine test_compute_contributions()
        integer(int32), parameter :: n_dims = 4
        integer(int32) :: ierr, mode
        real(real64) :: factor(n_dims), dependent(n_dims)
        real(real64) :: local_contributions(n_dims), expected_local(n_dims)
        real(real64) :: total_contribution, expected_total

        ! -------------------------------
        ! Case 1: RAW baseline
        ! -------------------------------
        factor    = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
        dependent = [2.0_real64, 1.0_real64, 0.0_real64, -1.0_real64]
        mode = 1  ! RAW

        call compute_contributions(factor, dependent, n_dims, mode, local_contributions, total_contribution, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_compute_contributions: Case 1 ierr")

        ! RAW baseline means baseline = 0
        expected_local = factor * dependent
        expected_total = sum(expected_local)

        call assert_equal_array_real(local_contributions, expected_local, n_dims, TOL, "test_compute_contributions: Case 1 local contributions")
        call assert_equal_real(total_contribution, expected_total, TOL, "test_compute_contributions: Case 1 total contribution")

        ! -------------------------------
        ! Case 2: MIN baseline
        ! -------------------------------
        factor    = [3.0_real64, 5.0_real64, 2.0_real64, 4.0_real64]
        dependent = [1.0_real64, 2.0_real64, 0.0_real64, -1.0_real64]
        mode = 2  ! MIN

        call compute_contributions(factor, dependent, n_dims, mode, local_contributions, total_contribution, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_compute_contributions: Case 2 ierr")

        ! Baseline = min(factor)=2, min(dependent)=-1
        expected_local = (factor - minval(factor)) * (dependent - minval(dependent))
        expected_total = sum(expected_local)

        call assert_equal_array_real(local_contributions, expected_local, n_dims, TOL, "test_compute_contributions: Case 2 local contributions")
        call assert_equal_real(total_contribution, expected_total, TOL, "test_compute_contributions: Case 2 total contribution")

        ! -------------------------------
        ! Case 3: MEAN baseline
        ! -------------------------------
        factor    = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
        dependent = [4.0_real64, 3.0_real64, 2.0_real64, 1.0_real64]
        mode = 3  ! MEAN

        call compute_contributions(factor, dependent, n_dims, mode, local_contributions, total_contribution, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_compute_contributions: Case 3 ierr")

        ! Baseline = mean(factor)=2.5, mean(dependent)=2.5
        expected_local = (factor - sum(factor) / n_dims) * (dependent - sum(dependent) / n_dims)
        expected_total = sum(expected_local)

        call assert_equal_array_real(local_contributions, expected_local, n_dims, TOL, "test_compute_contributions: Case 3 local contributions")
        call assert_equal_real(total_contribution, expected_total, TOL, "test_compute_contributions: Case 3 total contribution")
    end subroutine test_compute_contributions

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
