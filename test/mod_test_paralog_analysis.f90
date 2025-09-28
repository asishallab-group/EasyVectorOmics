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
        character(len=128) :: name
        procedure(test_interface), pointer, nopass :: test_proc => null()
    end type test_case

    real(real64), parameter :: TOL = 1d-12

contains

    !> Get array of all available tests.
    function get_all_tests() result(all_tests)
        type(test_case) :: all_tests(14)

        all_tests(1) = test_case("test_tox_paralog_analysis_mask_set_state", test_mask_set_state)
        all_tests(2) = test_case("test_tox_paralog_analysis_mask_check_state", test_mask_check_state)
        all_tests(3) = test_case("test_tox_paralog_analysis_mask_get_first_successor_idx", test_mask_get_first_successor_idx)
        all_tests(4) = test_case("test_tox_paralog_analysis_calc_work_arr_paralog_subsets_size", test_calc_work_arr_paralog_subsets_size)
        all_tests(5) = test_case("test_tox_paralog_analysis_filter_paralogs_by_pattern", test_filter_paralogs_by_pattern)
        all_tests(6) = test_case("test_tox_paralog_analysis_mask_chunk_count", test_mask_chunk_count)
        all_tests(7) = test_case("test_tox_paralog_analysis_add_new_active_mask", test_add_new_active_mask)
        all_tests(8) = test_case("test_tox_paralog_analysis_add_to_results", test_add_to_results)
        all_tests(9) = test_case("test_tox_paralog_analysis_take_active_mask", test_take_active_mask)
        all_tests(10) = test_case("test_tox_paralog_analysis_fill_array_with_minvals_for_each_idx", test_fill_array_with_minvals_for_each_idx)
        all_tests(11) = test_case("test_tox_paralog_analysis_angle_between", test_angle_between)
        all_tests(12) = test_case("test_tox_paralog_analysis_generate_subsets", test_generate_subsets)
        all_tests(13) = test_case("test_tox_paralog_analysis_detect_patterns_perfect_subfunc_split", test_detect_patterns_perfect_subfunc_split)
        all_tests(14) = test_case("test_tox_paralog_analysis_detect_patterns_subfunc_at_angle_margin", test_detect_patterns_subfunc_at_angle_margin)
    end function get_all_tests

    subroutine test_detect_patterns_perfect_subfunc_split
        use f42_utils, only: radians

        integer(int32), parameter :: n_paralogs = 3, n_dims = 3, n_mask_chunks = 1
        real(real64), parameter ::  rdi_threshold = 1e-6

        integer(int32), dimension(n_mask_chunks) :: filtered_paralogs_mask, active_mask
        real(real64), dimension(n_dims) :: ancestor, temp_paralog_vector
        real(real64), dimension(n_paralogs) :: temp_work_array, paralog_norms, paralog_angles
        integer(int32), dimension(n_paralogs) :: sorted_paralog_norms_perm
        real(real64), dimension(n_dims, n_paralogs) :: paralogs
        integer(int32), dimension(:, :), allocatable :: work_arr_paralog_subsets

        integer(int32) :: n_results, max_subset_size, work_array_size, ierr, i_paralog
        real(real64) :: prefilter_threshold

        call set_ok(ierr)

        max_subset_size = n_paralogs
        ancestor = [ 1, 1, 0 ]
        paralogs(:, 1) = [ 1, 0, 0 ]
        paralogs(:, 2) = [ 0, 1, 0 ]
        paralogs(:, 3) = [ 0, 0, 1 ]
        paralog_norms = 1
        sorted_paralog_norms_perm = [ 1, 2, 3 ]
        prefilter_threshold = radians(30.0_real64)

        do i_paralog = 1, n_paralogs
            call angle_between(paralogs(:, i_paralog), ancestor, n_dims, paralog_angles(i_paralog))
        end do

        call filter_paralogs_by_pattern(SUBFUNC_PATTERN, paralog_angles, prefilter_threshold, n_paralogs, filtered_paralogs_mask, n_mask_chunks, ierr)
        call assert_true(is_ok(ierr), "test_detect_patterns_perfect_subfunc_split: unexpected error when filtering paralogs for subfunctionalization")
        call calc_work_arr_paralog_subsets_size(max_subset_size, n_paralogs, work_array_size, filtered_paralogs_mask, n_mask_chunks, ierr)
        call assert_true(is_ok(ierr), "test_detect_patterns_perfect_subfunc_split: unexpected error when calculating work array size")
        allocate(work_arr_paralog_subsets(n_mask_chunks, work_array_size))

        call detect_subfunctionalization(ancestor, paralogs, n_paralogs, n_dims, rdi_threshold, filtered_paralogs_mask, n_mask_chunks, n_results, max_subset_size, work_arr_paralog_subsets, work_array_size, active_mask, temp_paralog_vector, paralog_norms, sorted_paralog_norms_perm, temp_work_array, ierr)
        call assert_true(is_ok(ierr), "test_detect_patterns_perfect_subfunc_split: unexpected error when detecting subfunctionalization")
        call assert_equal_int(n_results, 1_int32, "test_detect_patterns_perfect_subfunc_split: expected only one result for subfunctionalization")
        call assert_equal_int(work_arr_paralog_subsets(1, 1), 3_int32, "test_detect_patterns_perfect_subfunc_split: expected result mask to be 3=0b011 for subfunctionalization")

        deallocate(work_arr_paralog_subsets)
        call filter_paralogs_by_pattern(DOSAGE_PATTERN, paralog_angles, prefilter_threshold, n_paralogs, filtered_paralogs_mask, n_mask_chunks, ierr)
        call assert_true(is_ok(ierr), "test_detect_patterns_perfect_subfunc_split: unexpected error when filtering paralogs for dosage effect")
        call calc_work_arr_paralog_subsets_size(max_subset_size, n_paralogs, work_array_size, filtered_paralogs_mask, n_mask_chunks, ierr)
        call assert_true(is_ok(ierr), "test_detect_patterns_perfect_subfunc_split: unexpected error when calculating work array size")
        allocate(work_arr_paralog_subsets(n_mask_chunks, work_array_size))

        call detect_dosage_effect(ancestor, paralogs, n_paralogs, n_dims, rdi_threshold, filtered_paralogs_mask, n_mask_chunks, n_results, max_subset_size, work_arr_paralog_subsets, work_array_size, active_mask, temp_paralog_vector, ierr)
        call assert_true(is_ok(ierr), "test_detect_patterns_perfect_subfunc_split: unexpected error when detecting dosage effect")
        call assert_equal_int(n_results, 0_int32, "test_detect_patterns_perfect_subfunc_split: expected no results for dosage effect")
    end subroutine test_detect_patterns_perfect_subfunc_split

    subroutine test_detect_patterns_subfunc_at_angle_margin
        use f42_utils, only: radians

        integer(int32), parameter :: n_paralogs = 2, n_dims = 3, n_mask_chunks = 1
        real(real64), parameter ::  rdi_threshold = 1e-3

        integer(int32), dimension(n_mask_chunks) :: filtered_paralogs_mask, active_mask
        real(real64), dimension(n_dims) :: ancestor, temp_paralog_vector
        real(real64), dimension(n_paralogs) :: temp_work_array, paralog_norms, paralog_angles
        integer(int32), dimension(n_paralogs) :: sorted_paralog_norms_perm
        real(real64), dimension(n_dims, n_paralogs) :: paralogs
        integer(int32), dimension(:, :), allocatable :: work_arr_paralog_subsets

        integer(int32) :: n_results, max_subset_size, work_array_size, ierr, i_paralog
        real(real64) :: prefilter_threshold

        call set_ok(ierr)

        ! should not pass
        max_subset_size = n_paralogs
        ancestor = [ 1, 1, 0 ]
        paralogs(:, 1) = [ 1.0_real64, 0.0_real64, 0.001_real64 ]
        paralogs(:, 2) = [ 0.0_real64, 1.0_real64, 0.001_real64 ]
        paralog_norms = sqrt(1.0_real64 + 0.001_real64 ** 2)
        sorted_paralog_norms_perm = [ 1, 2 ]
        prefilter_threshold = radians(44.9_real64)

        do i_paralog = 1, n_paralogs
            call angle_between(paralogs(:, i_paralog), ancestor, n_dims, paralog_angles(i_paralog))
        end do

        call filter_paralogs_by_pattern(SUBFUNC_PATTERN, paralog_angles, prefilter_threshold, n_paralogs, filtered_paralogs_mask, n_mask_chunks, ierr)
        call assert_true(is_ok(ierr), "test_detect_patterns_subfunc_at_angle_margin: unexpected error when filtering paralogs for subfunctionalization")
        call calc_work_arr_paralog_subsets_size(max_subset_size, n_paralogs, work_array_size, filtered_paralogs_mask, n_mask_chunks, ierr)
        call assert_true(is_ok(ierr), "test_detect_patterns_subfunc_at_angle_margin: unexpected error when calculating work array size")
        allocate(work_arr_paralog_subsets(n_mask_chunks, work_array_size))

        call detect_subfunctionalization(ancestor, paralogs, n_paralogs, n_dims, rdi_threshold, filtered_paralogs_mask, n_mask_chunks, n_results, max_subset_size, work_arr_paralog_subsets, work_array_size, active_mask, temp_paralog_vector, paralog_norms, sorted_paralog_norms_perm, temp_work_array, ierr)
        call assert_true(is_ok(ierr), "test_detect_patterns_subfunc_at_angle_margin: unexpected error when detecting subfunctionalization")
        call assert_equal_int(n_results, 0_int32, "test_detect_patterns_subfunc_at_angle_margin: expected only one result for subfunctionalization")

        deallocate(work_arr_paralog_subsets)
        call filter_paralogs_by_pattern(DOSAGE_PATTERN, paralog_angles, prefilter_threshold, n_paralogs, filtered_paralogs_mask, n_mask_chunks, ierr)
        call assert_true(is_ok(ierr), "test_detect_patterns_subfunc_at_angle_margin: unexpected error when filtering paralogs for dosage effect")
        call calc_work_arr_paralog_subsets_size(max_subset_size, n_paralogs, work_array_size, filtered_paralogs_mask, n_mask_chunks, ierr)
        call assert_true(is_ok(ierr), "test_detect_patterns_subfunc_at_angle_margin: unexpected error when calculating work array size")
        allocate(work_arr_paralog_subsets(n_mask_chunks, work_array_size))

        call detect_dosage_effect(ancestor, paralogs, n_paralogs, n_dims, rdi_threshold, filtered_paralogs_mask, n_mask_chunks, n_results, max_subset_size, work_arr_paralog_subsets, work_array_size, active_mask, temp_paralog_vector, ierr)
        call assert_true(is_ok(ierr), "test_detect_patterns_subfunc_at_angle_margin: unexpected error when detecting dosage effect")
        call assert_equal_int(n_results, 0_int32, "test_detect_patterns_subfunc_at_angle_margin: expected no results for dosage effect")
        deallocate(work_arr_paralog_subsets)

        ! should pass
        max_subset_size = n_paralogs
        ancestor = [ 1, 1, 0 ]
        paralogs(:, 1) = [ 1.0_real64, 0.0_real64, 0.0004_real64 ]
        paralogs(:, 2) = [ 0.0_real64, 1.0_real64, 0.0004_real64 ]
        paralog_norms = sqrt(1.0_real64 + 0.0004_real64 ** 2)
        sorted_paralog_norms_perm = [ 1, 2 ]
        prefilter_threshold = radians(44.9_real64)

        do i_paralog = 1, n_paralogs
            call angle_between(paralogs(:, i_paralog), ancestor, n_dims, paralog_angles(i_paralog))
        end do

        call filter_paralogs_by_pattern(SUBFUNC_PATTERN, paralog_angles, prefilter_threshold, n_paralogs, filtered_paralogs_mask, n_mask_chunks, ierr)
        call assert_true(is_ok(ierr), "test_detect_patterns_subfunc_at_angle_margin: unexpected error when filtering paralogs for subfunctionalization")
        call calc_work_arr_paralog_subsets_size(max_subset_size, n_paralogs, work_array_size, filtered_paralogs_mask, n_mask_chunks, ierr)
        call assert_true(is_ok(ierr), "test_detect_patterns_subfunc_at_angle_margin: unexpected error when calculating work array size")
        allocate(work_arr_paralog_subsets(n_mask_chunks, work_array_size))

        call detect_subfunctionalization(ancestor, paralogs, n_paralogs, n_dims, rdi_threshold, filtered_paralogs_mask, n_mask_chunks, n_results, max_subset_size, work_arr_paralog_subsets, work_array_size, active_mask, temp_paralog_vector, paralog_norms, sorted_paralog_norms_perm, temp_work_array, ierr)
        call assert_true(is_ok(ierr), "test_detect_patterns_subfunc_at_angle_margin: unexpected error when detecting subfunctionalization")
        call assert_equal_int(n_results, 1_int32, "test_detect_patterns_subfunc_at_angle_margin: expected only one result for subfunctionalization")
        call assert_equal_int(work_arr_paralog_subsets(1, 1), 3_int32, "test_detect_patterns_subfunc_at_angle_margin: expected result mask to be 3=0b011 for subfunctionalization")

        deallocate(work_arr_paralog_subsets)
        call filter_paralogs_by_pattern(DOSAGE_PATTERN, paralog_angles, prefilter_threshold, n_paralogs, filtered_paralogs_mask, n_mask_chunks, ierr)
        call assert_true(is_ok(ierr), "test_detect_patterns_subfunc_at_angle_margin: unexpected error when filtering paralogs for dosage effect")
        call calc_work_arr_paralog_subsets_size(max_subset_size, n_paralogs, work_array_size, filtered_paralogs_mask, n_mask_chunks, ierr)
        call assert_true(is_ok(ierr), "test_detect_patterns_subfunc_at_angle_margin: unexpected error when calculating work array size")
        allocate(work_arr_paralog_subsets(n_mask_chunks, work_array_size))

        call detect_dosage_effect(ancestor, paralogs, n_paralogs, n_dims, rdi_threshold, filtered_paralogs_mask, n_mask_chunks, n_results, max_subset_size, work_arr_paralog_subsets, work_array_size, active_mask, temp_paralog_vector, ierr)
        call assert_true(is_ok(ierr), "test_detect_patterns_subfunc_at_angle_margin: unexpected error when detecting dosage effect")
        call assert_equal_int(n_results, 0_int32, "test_detect_patterns_subfunc_at_angle_margin: expected no results for dosage effect")
    end subroutine test_detect_patterns_subfunc_at_angle_margin

    subroutine test_generate_subsets
    end subroutine test_generate_subsets

    subroutine test_angle_between
        use f42_utils, only: PI

        integer(int32), parameter :: n_dims = 5
        real(real64), dimension(n_dims) :: v1, v2
        real(real64) :: angle

        ! test v2 == v1
        v1 = [1, 2, 3, 4, 5]
        v2 = v1
        call angle_between(v1, v2, n_dims, angle)
        call assert_equal_real(angle, 0.0_real64, epsilon(angle), "test_angle_between: vector should have zero angle to itself")
        v2 = v1 / 2
        call angle_between(v1, v2, n_dims, angle)
        call assert_equal_real(angle, 0.0_real64, epsilon(angle), "test_angle_between: vector should have zero angle to itself")

        ! test v2 points in opposite direction
        v2 = -v1
        call angle_between(v1, v2, n_dims, angle)
        call assert_equal_real(angle, PI, epsilon(angle), "test_angle_between: opposite vector should have 180 deg angle")

        ! test v2 perpendicular to v1
        v1 = [1, 0, 0, 0, 0]
        v2 = [0, 1, 0, 0, 0]
        call angle_between(v1, v2, n_dims, angle)
        call assert_equal_real(angle, PI / 2, epsilon(angle), "test_angle_between: opposite vector should have 90 deg angle")

        ! test v2 45 deg to v1
        v1 = [1, 0, 0, 0, 0]
        v2 = [0.7071, 0.7071, 0.0, 0.0, 0.0]
        call angle_between(v1, v2, n_dims, angle)
        call assert_equal_real(angle, PI / 4, epsilon(angle), "test_angle_between: opposite vector should have 90 deg angle")
    end subroutine test_angle_between

    subroutine test_fill_array_with_minvals_for_each_idx
        integer(int32), parameter :: src_arr_len = 10
        real(real64), dimension(src_arr_len) :: out_arr, src_arr
        integer(int32), dimension(src_arr_len) :: sorted_src_arr_perm
        integer(int32) :: ierr, i_element

        call set_ok(ierr)

        src_arr = [ 0.85117158_real64, -0.35000591_real64,  0.0880246_real64 ,  1.2852401_real64 ,  1.47949047_real64, -0.28126471_real64,  0.82010833_real64,  0.03206986_real64,  0.24169488_real64, 0.39394542_real64 ]
        sorted_src_arr_perm = [ 2, 6, 8, 3, 9, 10, 7, 1, 4, 5 ]
        call fill_array_with_minvals_for_each_idx(out_arr, src_arr, sorted_src_arr_perm, src_arr_len, ierr)
        call assert_true(is_ok(ierr), "test_fill_array_with_minvals_for_each_idx: unexpected error when calling routine")

        do i_element = 1, src_arr_len
            call assert_equal_real(out_arr(i_element), minval(src_arr(i_element:src_arr_len)), 0.0_real64, "test_fill_array_with_minvals_for_each_idx: min value does not match")
        end do
    end subroutine test_fill_array_with_minvals_for_each_idx

    subroutine test_take_active_mask
        integer(int32), parameter :: n_mask_chunks = 1, n_subsets = 5
        integer(int32), dimension(n_mask_chunks, n_subsets) :: subsets
        integer(int32), dimension(n_mask_chunks, n_subsets) :: original_subsets
        integer(int32), dimension(n_mask_chunks) :: taken_active_mask
        integer(int32) :: n_results, n_active_masks, n_new_active_masks, n_actual_active_masks, ierr, i_subset

        call set_ok(ierr)

        original_subsets = reshape([1,2,3,4,5], [n_mask_chunks, n_subsets])
        
        do n_new_active_masks = 0, n_subsets
            do n_results = 0, n_subsets - n_new_active_masks
                subsets = original_subsets
                do n_active_masks = n_subsets - n_new_active_masks - n_results, 0, -1
                    n_actual_active_masks = n_active_masks
                    call take_active_mask(subsets, n_mask_chunks, n_subsets, n_results, n_actual_active_masks, n_new_active_masks, taken_active_mask, ierr)
                    if (n_active_masks == 0) then
                        call assert_false(is_ok(ierr), "test_take_active_mask: expected an error when taking an active mask from zero active masks")
                    else
                        call assert_true(is_ok(ierr), "test_take_active_mask: unexpected error when taking active mask")
                        call assert_equal_int(n_actual_active_masks, n_active_masks - 1, "test_take_active_mask: number of active masks did not change")
                        call assert_equal_int(taken_active_mask(1), original_subsets(1, n_results + n_active_masks), "test_take_active_mask: taken mask does not match")
                        call assert_equal_array_int(subsets(:, 1:n_results), original_subsets(:, 1:n_results), n_results, "test_take_active_mask: results have changed")
                        do i_subset = 1, n_new_active_masks
                            call assert_array_int_contains(subsets(:, n_results+n_actual_active_masks+1:n_results+n_actual_active_masks+n_new_active_masks), original_subsets(1, n_subsets-i_subset+1), n_new_active_masks, "test_take_active_mask: new active masks changed")
                        end do
                    end if
                end do
            end do
        end do
    end subroutine test_take_active_mask

    subroutine test_add_to_results
        integer(int32), parameter :: n_mask_chunks = 1, n_subsets = 5
        integer(int32), dimension(n_mask_chunks, n_subsets) :: subsets
        integer(int32), dimension(n_mask_chunks, n_subsets) :: original_subsets
        integer(int32), dimension(n_mask_chunks) :: new_active_mask
        integer(int32) :: n_results, n_active_masks, n_new_active_masks, n_actual_results, ierr, i_subset

        call set_ok(ierr)

        original_subsets = reshape([1,2,3,4,5], [n_mask_chunks, n_subsets])
        new_active_mask = 42
        
        do n_new_active_masks = 0, n_subsets
            do n_active_masks = 0, n_subsets - n_new_active_masks
                subsets = original_subsets
                do n_results = 0, n_subsets - n_new_active_masks - n_active_masks
                    n_actual_results = n_results
                    call add_to_results(subsets, n_mask_chunks, n_subsets, n_actual_results, n_active_masks, n_new_active_masks, new_active_mask, ierr)
                    if (n_results + n_active_masks + n_new_active_masks == n_subsets) then
                        call assert_false(is_ok(ierr), "test_add_to_results: expected an error when adding a result to full array")
                    else
                        call assert_true(is_ok(ierr), "test_add_to_results: unexpected error when adding result")
                        call assert_equal_int(n_actual_results, n_results + 1, "test_add_to_results: number of results did not change")
                        call assert_true(all(subsets(:, 1:n_actual_results) == new_active_mask(1)), "test_add_to_results: results don't match")
                        do i_subset = 1, n_active_masks
                            call assert_array_int_contains(subsets(:, n_actual_results+1:n_actual_results+n_active_masks), original_subsets(1, i_subset), n_active_masks, "test_add_to_results: active masks changed")
                        end do
                        do i_subset = 1, n_new_active_masks
                            call assert_array_int_contains(subsets(:, n_actual_results+n_active_masks+1:n_actual_results+n_active_masks+n_new_active_masks), original_subsets(1, n_active_masks+i_subset), n_new_active_masks, "test_add_to_results: new active masks changed")
                        end do
                    end if
                end do
            end do
        end do
    end subroutine test_add_to_results

    subroutine test_add_new_active_mask
        integer(int32), parameter :: n_mask_chunks = 1, n_subsets = 5
        integer(int32), dimension(n_mask_chunks, n_subsets) :: subsets
        integer(int32), dimension(n_mask_chunks, n_subsets) :: original_subsets
        integer(int32), dimension(n_mask_chunks) :: new_active_mask
        integer(int32) :: n_results, n_active_masks, n_new_active_masks, n_actual_new_active_masks, ierr

        call set_ok(ierr)

        original_subsets = reshape([1,2,3,4,5], [n_mask_chunks, n_subsets])
        new_active_mask = 42
        do n_results = 0, n_subsets
            do n_active_masks = 0, n_subsets - n_results
                subsets = original_subsets
                do n_new_active_masks = 0, n_subsets - n_results - n_active_masks
                    n_actual_new_active_masks = n_new_active_masks
                    call add_new_active_mask(subsets, n_mask_chunks, n_subsets, n_results, n_active_masks, n_actual_new_active_masks, new_active_mask, ierr)
                    if (n_results + n_active_masks + n_new_active_masks == n_subsets) then
                        call assert_false(is_ok(ierr), "test_add_new_active_mask: expected an error when adding a new active mask to full array")
                    else
                        call assert_true(is_ok(ierr), "test_add_new_active_mask: unexpected error when adding new active mask")
                        call assert_equal_int(n_actual_new_active_masks, n_new_active_masks + 1, "test_add_new_active_mask: number of new active masks did not change")
                        call assert_equal_array_int(subsets(:, 1:n_results), original_subsets(:, 1:n_results), n_results, "test_add_new_active_mask: results have changed")
                        call assert_equal_array_int(subsets(:, n_results+1:n_results+n_active_masks), original_subsets(:, n_results+1:n_results+n_active_masks), n_active_masks, "test_add_new_active_mask: active masks have changed")
                        call assert_true(all(subsets(:, n_results+n_active_masks+1:n_results+n_active_masks+n_actual_new_active_masks) == new_active_mask(1)), "test_add_new_active_mask: new active masks don't match")
                    end if
                end do
            end do
        end do
    end subroutine test_add_new_active_mask

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
        ! integer(int32), parameter :: n_dims = 10
        ! integer(int32), parameter :: n_paralogs = 16, n_paralogs_overflow = 100
        ! integer(int32) :: i_paralog, max_subset_size_all_active, work_array_size, ierr, n_results, max_subset_size_overflown
        ! integer(int32) :: mask_all_active(1), active_mask(1), mask_all_active_overflow(4)
        ! integer(int32), allocatable :: work_arr_paralog_subsets(:, :)

        ! real(real64), dimension(n_dims) :: ancestor
        ! real(real64), dimension(n_dims, n_paralogs) :: paralogs
        ! real(real64), dimension(n_paralogs) :: temp_paralog_vector, subfunc_paralog_norms, subfunc_temp_work_array
        ! integer(int32), dimension(n_paralogs) :: subfunc_sorted_paralog_norms_perm
        ! real(real64), parameter :: rdi_threshold = 0.5


        ! call set_ok(ierr)

        ! ! stress the detect_patterns: Exploit an edge case where the whole working array is in use at some point to ensure correct size calculation

        ! do i_paralog = 1, n_paralogs
        !     call mask_set_state(mask_all_active, i_paralog, .true., ierr)
        !     call assert_true(is_ok(ierr), "test_calc_work_arr_paralog_subsets_size: unexpected error when enabling paralog in mask")

        !     subfunc_sorted_paralog_norms_perm(i_paralog) = n_paralogs - i_paralog + 1
        ! end do

        ! ancestor = 1.0
        ! paralogs = 0.0
        ! paralogs(:, n_paralogs) = 1.0 ! residual with last paralog active will produce a norm below rdi_threshold -> only subsets with last paralog included (cannot be extended) will be results
        ! subfunc_paralog_norms = 1.0
        ! subfunc_paralog_norms(n_paralogs) = 0.0 ! last paralog norm is lower residual -> no subset candidate will be pruned

        ! do i_paralog = 1, n_paralogs
        !     max_subset_size_all_active = i_paralog
        !     call calc_work_arr_paralog_subsets_size(max_subset_size_all_active, n_paralogs, work_array_size, mask_all_active, size(mask_all_active), ierr)

        !     allocate(work_arr_paralog_subsets(1, work_array_size + 1))
        !     work_arr_paralog_subsets = 0
        !     call detect_patterns(ancestor, paralogs, n_paralogs, n_dims, rdi_threshold, SUBFUNC_PATTERN, mask_all_active, size(mask_all_active), n_results, max_subset_size_all_active, work_arr_paralog_subsets, work_array_size + 1, active_mask, temp_paralog_vector, subfunc_paralog_norms=subfunc_paralog_norms, subfunc_sorted_paralog_norms_perm=subfunc_sorted_paralog_norms_perm, subfunc_temp_work_array=subfunc_temp_work_array, ierr=ierr)

        !     ! masks have at least one active bit -> non-zero
        !     ! masks also won't be reset to zero, as new added masks overwrite them anyway.
        !     ! Thus, all calculated needed space should be used during detection -> non-zero
        !     call assert_equal_int(count(work_arr_paralog_subsets /= 0), work_array_size, "test_calc_work_arr_paralog_subsets_size: different count of subsets used than expected")
        !     call assert_true(is_ok(ierr), "test_calc_work_arr_paralog_subsets_size: unexpected error when detecting patterns")
        !     deallocate(work_arr_paralog_subsets)
        ! end do

        ! do i_paralog = 1, n_paralogs_overflow
        !     call mask_set_state(mask_all_active_overflow, i_paralog, .true., ierr)
        !     call assert_true(is_ok(ierr), "test_calc_work_arr_paralog_subsets_size: unexpected error when enabling paralog in oerflow mask")
        ! end do
        ! max_subset_size_overflown = 16
        ! call calc_work_arr_paralog_subsets_size(max_subset_size_overflown, n_paralogs_overflow, work_array_size, mask_all_active_overflow, size(mask_all_active_overflow), ierr)
        ! call assert_not_equal_int(max_subset_size_overflown, 16_int32, "test_calc_work_arr_paralog_subsets_size: for overflow the max subset size should be different to input")
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
