! filepath: test/mod_test_trajectory_normalization.f90
!> Unit test suite for trajectory_normalization routines.
module mod_test_trajectory_normalization
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_trajectory_normalization
    use tox_errors
    implicit none

    abstract interface
        subroutine test_interface()
        end subroutine test_interface
    end interface

    type :: test_case
        character(len=128) :: name
        procedure(test_interface), pointer, nopass :: test_proc => null()
    end type test_case

    real(real64), parameter :: TOL = 1.0e-10_real64

contains

    function get_all_tests() result(all_tests)
        type(test_case) :: all_tests(3)

        all_tests(1) = test_case("test_normalize_variable_timeseries", test_normalize_variable_timeseries)
        all_tests(2) = test_case("test_normalize_single_trajectory", test_normalize_single_trajectory)
        all_tests(3) = test_case("test_normalize_all_trajectories", test_normalize_all_trajectories)
    end function get_all_tests

   !> Run all tox_trajectory_normalization tests.
    subroutine run_all_tests_trajectory_normalization
        type(test_case) :: all_tests(3)
        integer :: i


        all_tests = get_all_tests()
        do i = 1, size(all_tests)
            call all_tests(i)%test_proc()
            print *, trim(all_tests(i)%name), " passed."
        end do
        print *, "All trajectory_normalization tests passed successfully."
    end subroutine run_all_tests_trajectory_normalization

    !> Run specific tox_trajectory_normalization tests by name.
    subroutine run_named_tests_trajectory_normalization(test_names)
        character(len=*), intent(in) :: test_names(:)
        type(test_case) :: all_tests(3)
        integer :: i, j
        logical :: found

        all_tests = get_all_tests()

        do i = 1, size(test_names)
            found = .false.
            do j = 1, size(all_tests)
                if (trim(test_names(i)) == trim(all_tests(j)%name)) then
                    call all_tests(j)%test_proc()
                    print *, trim(all_tests(i)%name), " passed."
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                print *, "Unknown test: ", trim(test_names(i))
            end if
        end do
    end subroutine run_named_tests_trajectory_normalization

    !> Test the normalization of variable timeseries.
    subroutine test_normalize_variable_timeseries()
        integer(int32), parameter :: n_points = 5
        real(real64) :: v(n_points), v_norm(n_points), expected(n_points)
        integer(int32) :: ierr, status

        v = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64]
        expected = [0.0_real64, 0.25_real64, 0.5_real64, 0.75_real64, 1.0_real64]

        call normalize_variable_timeseries(v, v_norm, n_points, ierr, status)

        call assert_equal_int(ierr, ERR_OK, "test_normalize_variable_timeseries: ierr")
        call assert_equal_array_real(v_norm, expected, n_points, TOL, "test_normalize_variable_timeseries: values")
    end subroutine test_normalize_variable_timeseries

    !> Test the normalization of a single trajectory.
    subroutine test_normalize_single_trajectory()
        integer(int32), parameter :: n_factors = 2, n_timepoints = 4
        real(real64) :: trajectory(n_timepoints, n_factors), trajectory_norm(n_timepoints, n_factors)
        real(real64) :: expected(n_timepoints)
        integer(int32) :: ierr, status, i

        trajectory(:, 1) = [10.0_real64, 11.0_real64, 12.0_real64, 13.0_real64]
        trajectory(:, 2) = [20.0_real64, 21.0_real64, 22.0_real64, 23.0_real64]
        expected = [0.0_real64, 1.0_real64/3.0_real64, 2.0_real64/3.0_real64, 1.0_real64]

        call normalize_single_trajectory(trajectory, trajectory_norm, n_factors, n_timepoints, ierr, status)

        call assert_equal_int(ierr, ERR_OK, "test_normalize_single_trajectory: ierr")
        do i = 1, n_factors
            call assert_equal_array_real(trajectory_norm(:, i), expected, n_timepoints, TOL, "test_normalize_single_trajectory: factor")
        end do
    end subroutine test_normalize_single_trajectory

    !> Test the normalization of all trajectories in a 3D array.
    subroutine test_normalize_all_trajectories()
        integer(int32), parameter :: n_factors = 2, n_samples = 2, n_timepoints = 4
        real(real64) :: trajectories(n_factors, n_samples, n_timepoints)
        real(real64) :: trajectories_norm(n_factors, n_samples, n_timepoints)
        integer(int32) :: ierr, status, i_factor, i_sample

        trajectories(1, 1, :) = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
        trajectories(2, 1, :) = [10.0_real64, 20.0_real64, 30.0_real64, 40.0_real64]
        trajectories(1, 2, :) = [5.0_real64, 6.0_real64, 7.0_real64, 8.0_real64]
        trajectories(2, 2, :) = [2.0_real64, 4.0_real64, 6.0_real64, 8.0_real64]

        call normalize_all_trajectories(trajectories, trajectories_norm, n_factors, n_samples, n_timepoints, ierr, status)

        call assert_equal_int(ierr, ERR_OK, "test_normalize_all_trajectories: ierr")
        do i_sample = 1, n_samples
            do i_factor = 1, n_factors
                call assert_equal_real(minval(trajectories_norm(i_factor, i_sample, :)), 0.0_real64, TOL, "test_normalize_all_trajectories: min")
                call assert_equal_real(maxval(trajectories_norm(i_factor, i_sample, :)), 1.0_real64, TOL, "test_normalize_all_trajectories: max")
            end do
        end do
    end subroutine test_normalize_all_trajectories

end module mod_test_trajectory_normalization
