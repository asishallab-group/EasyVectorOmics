!> Unit test suite for tox_paralog_analysis routine.
module mod_test_tox_paralog_analysis
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding
    use tox_paralog_analysis
    use tox_errors
    implicit none

    ! Abstract interface for all test procedures
    abstract interface
        subroutine test_interface()
        end subroutine test_interface
    end interface

    ! Type to hold test name and procedure pointer
    type :: test_case
        character(len=64) :: name
        procedure(test_interface), pointer, nopass :: test_proc => null()
    end type test_case

    real(real64), parameter :: TOL = 1d-12

contains

    !> Get array of all available tests.
    function get_all_tests() result(all_tests)
        type(test_case) :: all_tests(4)

        all_tests(1) = test_case("test_tox_paralog_analysis_mask_set_state", test_mask_set_state)
        all_tests(2) = test_case("test_tox_paralog_analysis_mask_check_state", test_mask_check_state)
        all_tests(3) = test_case("test_tox_paralog_analysis_mask_get_first_successor_idx", test_mask_get_first_successor_idx)
        all_tests(4) = test_case("test_tox_paralog_analysis_calc_work_arr_paralog_subsets_size", test_calc_work_arr_paralog_subsets_size)
    end function get_all_tests

    subroutine test_calc_work_arr_paralog_subsets_size
        integer(int32), parameter :: n_paralogs = 100
    end subroutine test_calc_work_arr_paralog_subsets_size

    subroutine test_mask_set_state
        integer(int32), parameter :: n_paralogs = 32 + 27
        integer(int32), parameter :: mask_size = 2
        integer(int32), dimension(mask_size) :: expected_mask
        integer(int32), dimension(mask_size) :: actual_mask
        integer(int32) :: ierr, paralog

        expected_mask = 0
        actual_mask = 0

        ! set first paralog
        paralog = 1
        expected_mask(1) = ibset(expected_mask(1), paralog - 1)
        call mask_set_state(actual_mask, paralog, .true., ierr)
        call assert_true(is_ok(ierr), "test_tox_paralog_analysis_mask_set_state: could not set first paralog")
        call assert_equal_array_int(actual_mask, expected_mask, mask_size, "test_tox_paralog_analysis_mask_set_state: mismatched mask setting first paralog")

        ! set last paralog
        paralog = n_paralogs
        expected_mask(2) = ibset(expected_mask(2), paralog - 32 - 1)
        call mask_set_state(actual_mask, paralog, .true., ierr)
        call assert_true(is_ok(ierr), "test_tox_paralog_analysis_mask_set_state: could not set last paralog")
        call assert_equal_array_int(actual_mask, expected_mask, mask_size, "test_tox_paralog_analysis_mask_set_state: mismatched mask setting last paralog")

        ! set 32nd paralog
        paralog = 32
        expected_mask(1) = ibset(expected_mask(1), paralog - 1)
        call mask_set_state(actual_mask, paralog, .true., ierr)
        call assert_true(is_ok(ierr), "test_tox_paralog_analysis_mask_set_state: could not set 32nd paralog")
        call assert_equal_array_int(actual_mask, expected_mask, mask_size, "test_tox_paralog_analysis_mask_set_state: mismatched mask setting 32nd paralog")

        ! unset all
        call mask_set_state(actual_mask, 1, .false., ierr)
        call assert_true(is_ok(ierr), "test_tox_paralog_analysis_mask_set_state: could not unset first paralog")
        call mask_set_state(actual_mask, n_paralogs, .false., ierr)
        call assert_true(is_ok(ierr), "test_tox_paralog_analysis_mask_set_state: could not set last paralog")
        call mask_set_state(actual_mask, 32, .false., ierr)
        call assert_true(is_ok(ierr), "test_tox_paralog_analysis_mask_set_state: could not set 32nd paralog")

        call assert_true(all(actual_mask == 0), "test_tox_paralog_analysis_mask_set_state: not all unset")
    end subroutine test_mask_set_state

    subroutine test_mask_check_state
        integer(int32), parameter :: n_paralogs = 32 + 27
        integer(int32), parameter :: mask_size = 2
        integer(int32), dimension(mask_size) :: mask
        integer(int32) :: paralog, i

        mask = 0

        do i = 1, 32
            call assert_false(mask_check_state(mask, i), "test_tox_paralog_analysis_mask_check_state: state should be false")
        end do

        ! set first paralog
        paralog = 1
        mask(1) = ibset(mask(1), paralog - 1)
        call assert_true(mask_check_state(mask, paralog), "test_tox_paralog_analysis_mask_check_state: first paralog wrong state")

        ! set last paralog
        paralog = n_paralogs
        mask(2) = ibset(mask(2), paralog - 32 - 1)
        call assert_true(mask_check_state(mask, paralog), "test_tox_paralog_analysis_mask_check_state: last paralog wrong state")

        ! set 32nd paralog
        paralog = 32
        mask(1) = ibset(mask(1), paralog - 1)
        call assert_true(mask_check_state(mask, paralog), "test_tox_paralog_analysis_mask_check_state: 32nd paralog wrong state")
    end subroutine test_mask_check_state

    subroutine test_mask_get_first_successor_idx
        integer(int32), parameter :: n_paralogs = 32 + 27
        integer(int32), parameter :: mask_size = 2
        integer(int32), dimension(mask_size) :: mask
        integer(int32) :: paralog, ierr

        mask = 0

        call assert_equal_int(mask_get_first_successor_idx(mask, n_paralogs), 1, "test_tox_paralog_analysis_mask_get_first_successor_idx: wrong number of zeros")

        do paralog = 1, n_paralogs
            call mask_set_state(mask, paralog, .true., ierr)
            call assert_true(is_ok(ierr), "test_tox_paralog_analysis_mask_get_first_successor_idx: Unexpected error when setting paralog active")
            call assert_equal_int(mask_get_first_successor_idx(mask, n_paralogs), paralog + 1, "test_tox_paralog_analysis_mask_get_first_successor_idx: wrong number of zeros")
        end do
    end subroutine test_mask_get_first_successor_idx

    !> Run all tox_paralog_analysis tests.
    subroutine run_all_tests_tox_paralog_analysis
        type(test_case), allocatable :: all_tests(:)
        integer(int32) :: i

        all_tests = get_all_tests()

        do i = 1, size(all_tests)
            call all_tests(i)%test_proc()
            print *, trim(all_tests(i)%name), " passed."
        end do
        print *, "All tox_paralog_analysis tests passed successfully."
    end subroutine run_all_tests_tox_paralog_analysis

    !> Run specific tox_paralog_analysis tests by name.
    subroutine run_named_tests_tox_paralog_analysis(test_names)
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
                    print *, trim(test_names(i)), " passed."
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                print *, "Unknown test: ", trim(test_names(i))
            end if
        end do
    end subroutine run_named_tests_tox_paralog_analysis
end module mod_test_tox_paralog_analysis
