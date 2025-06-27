program run_tests
  use mod_test_sorting
  implicit none
  integer :: nargs, i
  character(len=32) :: arg

  nargs = command_argument_count()
  if (nargs == 0) then
    call run_all_tests_sorting()
  else
    do i = 1, nargs
      call get_command_argument(i, arg)
      select case (trim(arg))
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
      case default
        print *, "Unknown test: ", trim(arg)
      end select
    end do
  end if
end program run_tests