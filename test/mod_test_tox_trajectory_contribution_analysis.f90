!> Unit test suite for tox_trajectory_contribution_analysis routine.
module mod_test_tox_trajectory_contribution_analysis
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
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
        type(test_case) :: all_tests(2)

        all_tests(1) = test_case("test_tox_trajectory_contribution_analysis_get_vec_across_samples", test_get_vec_across_samples)
        all_tests(2) = test_case("test_tox_trajectory_contribution_analysis_get_vec_across_timepoints", test_get_vec_across_timepoints)
    end function get_all_tests

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
        call assert_true(is_ok(ierr), "Test failed: unexpected error code")
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
        call assert_true(is_ok(ierr), "test_get_vec_across_timepoints: unexpected error code")
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
