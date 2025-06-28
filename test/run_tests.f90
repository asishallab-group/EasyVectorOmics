program run_tests
  use mod_test_sorting
  use mod_test_normalize_by_std_dev
  use mod_test_quantile_normalization

  implicit none
  integer :: nargs, i
  character(len=32) :: arg

  nargs = command_argument_count()
  if (nargs == 0) then
    call run_all_tests_sorting()
    call run_all_tests_normalize_by_std_dev()
    call run_all_tests_quantile_normalization()

  else
    do i = 1, nargs
      call get_command_argument(i, arg)
      select case (trim(arg))
      ! --- Sorting tests ---
      case ("test_sort_real")
        call test_sort_real()
        print *, "test_sort_real passed."
      case ("test_sort_integer")
        call test_sort_integer()
        print *, "test_sort_integer passed."
      case ("test_sort_character")
        call test_sort_character()
        print *, "test_sort_character passed."
      case ("test_sort_integer_ascending")
        call test_sort_integer_ascending()
        print *, "test_sort_integer_ascending passed."
      case ("test_sort_real_descending")
        call test_sort_real_descending()
        print *, "test_sort_real_descending passed."
      case ("test_sort_char_random")
        call test_sort_char_random()
        print *, "test_sort_char_random passed."
      case ("test_sort_sorted_stability")
        call test_sort_sorted_stability()
        print *, "test_sort_sorted_stability passed."
      case ("test_sort_empty_array")
        call test_sort_empty_array()
        print *, "test_sort_empty_array passed."
      case ("test_sort_large_random")
        call test_sort_large_random()
        print *, "test_sort_large_random passed."

      ! --- Normalization tests ---
      case ("test_normalize_by_std_dev_basic")
        call test_normalize_by_std_dev_basic()
        print *, "test_normalize_by_std_dev_basic passed."
      case ("test_normalize_by_std_dev_constant_rows")
        call test_normalize_by_std_dev_constant_rows()
        print *, "test_normalize_by_std_dev_constant_rows passed."
      case ("test_normalize_by_std_dev_large_numbers")
        call test_normalize_by_std_dev_large_numbers()
        print *, "test_normalize_by_std_dev_large_numbers passed."
      case ("test_identity_matrix")
        call test_identity_matrix()
        print *, "test_identity_matrix passed."
      case ("test_zero_rows")
        call test_zero_rows()
        print *, "test_zero_rows passed."
      case ("test_negative_rows")
        call test_negative_rows()
        print *, "test_negative_rows passed."
      case ("test_large_random_matrix")
        call test_large_random_matrix()
        print *, "test_large_random_matrix passed."
      case ("test_single_nonzero")
        call test_single_nonzero()
        print *, "test_single_nonzero passed."
      case ("test_small_large_values")
        call test_small_large_values()
        print *, "test_small_large_values passed."
      case ("test_nan_inf_input")
        call test_nan_inf_input()
        print *, "test_nan_inf_input passed."
      case ("test_single_row_col")
        call test_single_row_col()
        print *, "test_single_row_col passed."
      case ("test_empty_matrix")
        call test_empty_matrix()
        print *, "test_empty_matrix passed."
      case ("test_symmetric_rows")
        call test_symmetric_rows()
        print *, "test_symmetric_rows passed."

      ! --- Quantile normalization tests ---
      case ("test_qn_preserves_dimensions")
        call test_qn_preserves_dimensions()
        print *, "test_qn_preserves_dimensions passed."
      case ("test_qn_identical_rows")
        call test_qn_identical_rows()
        print *, "test_qn_identical_rows passed."
      case ("test_qn_no_nans_and_standardizes")
        call test_qn_no_nans_and_standardizes()
        print *, "test_qn_no_nans_and_standardizes passed."
      case ("test_qn_single_row")
        call test_qn_single_row()
        print *, "test_qn_single_row passed."
      case ("test_qn_single_column")
        call test_qn_single_column()
        print *, "test_qn_single_column passed."
      case ("test_qn_all_equal")
        call test_qn_all_equal()
        print *, "test_qn_all_equal passed."
      case ("test_qn_large_random")
        call test_qn_large_random()
        print *, "test_qn_large_random passed."
      case ("test_qn_negative_values")
        call test_qn_negative_values()
        print *, "test_qn_negative_values passed."
      case ("test_qn_zero_matrix")
        call test_qn_zero_matrix()
        print *, "test_qn_zero_matrix passed."
      case ("test_qn_sorted_input")
        call test_qn_sorted_input()
        print *, "test_qn_sorted_input passed."
      case ("test_qn_reverse_sorted")
        call test_qn_reverse_sorted()
        print *, "test_qn_reverse_sorted passed."

      case default
        print *, "Unknown test: ", trim(arg)
      end select
    end do
  end if
end program run_tests