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
        type(test_case) :: all_tests(6)

        all_tests(1) = test_case("test_tox_paralog_analysis_mask_set_state", test_mask_set_state)
        all_tests(2) = test_case("test_tox_paralog_analysis_mask_check_state", test_mask_check_state)
        all_tests(3) = test_case("test_tox_paralog_analysis_mask_get_first_successor_idx", test_mask_get_first_successor_idx)
        all_tests(4) = test_case("test_tox_paralog_analysis_calc_work_arr_paralog_subsets_size", test_calc_work_arr_paralog_subsets_size)
        all_tests(5) = test_case("test_tox_paralog_analysis_test_filter_paralogs_by_pattern", test_filter_paralogs_by_pattern)
        all_tests(6) = test_case("test_tox_paralog_analysis_test_mask_chunk_count", test_mask_chunk_count)
    end function get_all_tests

    subroutine test_mask_chunk_count
        integer(int32) :: i, n_chunks, n_expected_chunks

        do n_expected_chunks = 1, 10
            do i = (n_expected_chunks - 1) * 32 + 1, n_expected_chunks * 32
                call mask_chunk_count(i, n_chunks)
                call assert_equal_int(n_chunks, n_expected_chunks, "mask_chunk_count: calculated chunk count differs from expected")
            end do
        end do
    end subroutine test_mask_chunk_count

    subroutine test_filter_paralogs_by_pattern
        integer(int32), parameter :: n_paralogs = 16
        integer(int32), dimension(1) :: mask
        real(real64), parameter :: threshold = 0.5
        real(real64), dimension(n_paralogs) :: paralog_angles
        integer(int32) :: ierr, i_paralog, n_in_filtered

        paralog_angles = 0.5 * threshold
        paralog_angles(1:n_paralogs:2) = 2 * threshold
        paralog_angles(1:4) = 2 * threshold

        call filter_paralogs_by_pattern(SUBFUNC_PATTERN, paralog_angles, threshold, n_paralogs, mask, size(mask), ierr)
        n_in_filtered = 0
        do i_paralog = 1, n_paralogs
            if (mask_check_state(mask, i_paralog)) then
                n_in_filtered = n_in_filtered + 1
            end if
        end do
        call assert_equal_int(n_in_filtered, count(paralog_angles >= threshold), "test_filter_paralogs_by_pattern: wrong filtering for subfunctionalization")
        call filter_paralogs_by_pattern(DOSAGE_PATTERN, paralog_angles, threshold, n_paralogs, mask, size(mask), ierr)
        n_in_filtered = 0
        do i_paralog = 1, n_paralogs
            if (mask_check_state(mask, i_paralog)) then
                n_in_filtered = n_in_filtered + 1
            end if
        end do
        call assert_equal_int(n_in_filtered, count(paralog_angles <= threshold), "test_filter_paralogs_by_pattern: wrong filtering for subfunctionalization")
    end subroutine test_filter_paralogs_by_pattern

    subroutine test_calc_work_arr_paralog_subsets_size
        integer(int32), parameter :: n_dims = 10
        integer(int32), parameter :: n_paralogs = 16, n_paralogs_overflow = 100
        integer(int32) :: i_paralog, max_subset_size_all_active, work_array_size, ierr, n_results, max_subset_size_overflown
        integer(int32) :: mask_all_active(1), active_mask(1), mask_all_active_overflow(4)
        integer(int32), allocatable :: work_arr_paralog_subsets(:, :)

        real(real64), dimension(n_dims) :: ancestor
        real(real64), dimension(n_dims, n_paralogs) :: paralogs
        real(real64), dimension(n_paralogs) :: temp_paralog_vector, subfunc_paralog_norms, subfunc_temp_work_array
        integer(int32), dimension(n_paralogs) :: subfunc_sorted_paralog_norms_perm
        real(real64), parameter :: rdi_threshold = 0.5


        call set_ok(ierr)

        ! stress the detect_patterns: Exploit an edge case where the whole working array is in use at some point to ensure correct size calculation

        do i_paralog = 1, n_paralogs
            call mask_set_state(mask_all_active, i_paralog, .true., ierr)
            call assert_true(is_ok(ierr), "test_calc_work_arr_paralog_subsets_size: unexpected error when enabling paralog in mask")

            subfunc_sorted_paralog_norms_perm(i_paralog) = n_paralogs - i_paralog + 1
        end do

        ancestor = 1.0
        paralogs = 0.0
        paralogs(:, n_paralogs) = 1.0 ! residual with last paralog active will produce a norm below rdi_threshold -> only subsets with last paralog included (cannot be extended) will be results
        subfunc_paralog_norms = 1.0
        subfunc_paralog_norms(n_paralogs) = 0.0 ! last paralog norm is lower residual -> no subset candidate will be pruned

        do i_paralog = 1, n_paralogs
            max_subset_size_all_active = i_paralog
            call calc_work_arr_paralog_subsets_size(max_subset_size_all_active, n_paralogs, work_array_size, mask_all_active, size(mask_all_active), ierr)

            allocate(work_arr_paralog_subsets(1, work_array_size + 1))
            work_arr_paralog_subsets = 0
            call detect_patterns(ancestor, paralogs, n_paralogs, n_dims, rdi_threshold, SUBFUNC_PATTERN, mask_all_active, size(mask_all_active), n_results, max_subset_size_all_active, work_arr_paralog_subsets, work_array_size + 1, active_mask, temp_paralog_vector, subfunc_paralog_norms=subfunc_paralog_norms, subfunc_sorted_paralog_norms_perm=subfunc_sorted_paralog_norms_perm, subfunc_temp_work_array=subfunc_temp_work_array, ierr=ierr)

            ! masks have at least one active bit -> non-zero
            ! masks also won't be reset to zero, as new added masks overwrite them anyway.
            ! Thus, all calculated needed space should be used during detection -> non-zero
            call assert_equal_int(count(work_arr_paralog_subsets /= 0), work_array_size, "test_calc_work_arr_paralog_subsets_size: less subsets used than expected")
            call assert_true(is_ok(ierr), "test_calc_work_arr_paralog_subsets_size: unexpected error when detecting patterns")
            deallocate(work_arr_paralog_subsets)
        end do

        do i_paralog = 1, n_paralogs_overflow
            call mask_set_state(mask_all_active_overflow, i_paralog, .true., ierr)
            call assert_true(is_ok(ierr), "test_calc_work_arr_paralog_subsets_size: unexpected error when enabling paralog in oerflow mask")
        end do
        max_subset_size_overflown = 16
        call calc_work_arr_paralog_subsets_size(max_subset_size_overflown, n_paralogs_overflow, work_array_size, mask_all_active_overflow, size(mask_all_active_overflow), ierr)
        call assert_not_equal_int(max_subset_size_overflown, 16_int32, "test_calc_work_arr_paralog_subsets_size: for overflow the max subset size should be different to input")
    end subroutine test_calc_work_arr_paralog_subsets_size

    subroutine test_mask_set_state
        integer(int32), parameter :: n_paralogs = 32 + 27
        integer(int32), parameter :: mask_size = 2
        integer(int32), dimension(mask_size) :: expected_mask
        integer(int32), dimension(mask_size) :: actual_mask
        integer(int32) :: ierr, paralog

        call set_ok(ierr)

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

        call set_ok(ierr)
        
        mask = 0

        call assert_equal_int(mask_get_first_successor_idx(mask), 1, "test_tox_paralog_analysis_mask_get_first_successor_idx: wrong number of zeros")

        do paralog = 1, n_paralogs
            call mask_set_state(mask, paralog, .true., ierr)
            call assert_true(is_ok(ierr), "test_tox_paralog_analysis_mask_get_first_successor_idx: Unexpected error when setting paralog active")
            call assert_equal_int(mask_get_first_successor_idx(mask), paralog + 1, "test_tox_paralog_analysis_mask_get_first_successor_idx: wrong number of zeros")
        end do

        do paralog = 1, n_paralogs - 1
            call mask_set_state(mask, paralog, .false., ierr)
            call assert_true(is_ok(ierr), "test_tox_paralog_analysis_mask_get_first_successor_idx: Unexpected error when setting paralog active")
            call assert_equal_int(mask_get_first_successor_idx(mask), n_paralogs + 1, "test_tox_paralog_analysis_mask_get_first_successor_idx: wrong number of zeros")
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
