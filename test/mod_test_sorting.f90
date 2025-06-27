!> @brief Unit test suite for tox_sorting module using asserts.
!! @details
!! This module tests sorting routines (real, integer, character)
!! using fixed data and known expected outputs.
module mod_test_sorting
  use tox_sorting
  use asserts
  implicit none
  public 

contains

  !> @brief Run all sorting tests.
  !! Calls all test subroutines in this module.
  subroutine run_all_tests_sorting()
    call test_sort_real()
    call test_sort_integer()
    call test_sort_character()
    call test_sort_integer_ascending()
    call test_sort_real_descending()
    call test_sort_char_random()
    call test_sort_sorted_stability()
    call test_sort_empty_array()
    call test_sort_large_random()
    print *, "All sorting tests passed!"
  end subroutine run_all_tests_sorting

  !> @brief Test sorting of a real array.
  !! Checks that the permutation matches the expected result.
  subroutine test_sort_real()
    real(8), dimension(5) :: data = [3.0d0, 1.0d0, 5.0d0, 2.0d0, 4.0d0]
    integer, dimension(5) :: perm, expected = [2, 4, 1, 5, 3]
    integer :: stack_left(20), stack_right(20), i
    perm = [(i, i = 1, 5)]
    call sort_array(data, perm, stack_left, stack_right)
    call assert_equal_array_int(perm, expected, 5, "test_sort_real: permutation mismatch")
  end subroutine test_sort_real

  !> @brief Test sorting of an integer array.
  !! Checks that the permutation matches the expected result.
  subroutine test_sort_integer()
    integer, dimension(4) :: data = [10, 3, 7, 1]
    integer, dimension(4) :: perm, expected = [4, 2, 3, 1]
    integer :: stack_left(20), stack_right(20)
    integer :: i
    perm = [(i, i = 1, 4)]
    call sort_array(data, perm, stack_left, stack_right)
    call assert_equal_array_int(perm, expected, 4, "test_sort_integer: permutation mismatch")
  end subroutine test_sort_integer

  !> @brief Test sorting of a character array.
  !! Checks that the permutation matches the expected result.
  subroutine test_sort_character()
    character(len=6), dimension(3) :: data = ['delta ', 'alpha ', 'beta  ']
    integer, dimension(3) :: perm, expected = [2, 3, 1]
    integer :: stack_left(20), stack_right(20)
    integer :: i
    perm = [(i, i = 1, 3)]
    call sort_array(data, perm, stack_left, stack_right)
    call assert_equal_array_int(perm, expected, 3, "test_sort_character: permutation mismatch")
  end subroutine test_sort_character

  !> @brief Test sorting of an integer array in ascending order.
  !! Checks that the sorted values match the expected result.
  subroutine test_sort_integer_ascending()
    integer, dimension(5) :: data = [5, 2, 9, 1, 6]
    integer, dimension(5) :: perm
    integer :: stack_left(20), stack_right(20)
    integer :: i
    perm = [(i, i=1,5)]
    call sort_array(data, perm, stack_left, stack_right)
    call assert_equal_array_int(data(perm), [1,2,5,6,9], 5, "test_sort_integer_ascending: sorted values mismatch")
  end subroutine test_sort_integer_ascending

  !> @brief Test sorting of a real array in ascending order.
  !! Checks that the sorted values match the expected result.
  subroutine test_sort_real_descending()
    real(8), dimension(5) :: data = [3.5d0, 2.2d0, 8.8d0, 1.1d0, 7.7d0]
    integer, dimension(5) :: perm
    integer :: stack_left(20), stack_right(20)
    integer :: i
    perm = [(i, i=1,5)]
    call sort_array(data, perm, stack_left, stack_right)
    call assert_equal_array_real(data(perm), [1.1d0, 2.2d0, 3.5d0, 7.7d0, 8.8d0], 5, 1d-12, &
                            "test_sort_real_descending: sorted values mismatch")
  end subroutine test_sort_real_descending

  !> @brief Test sorting of a random character array.
  !! Checks that the sorted values match the expected result.
  subroutine test_sort_char_random()
    character(len=8), dimension(5) :: data = ['dog     ', 'apple   ', 'zebra   ', 'cat     ', 'bird    ']
    character(len=8), dimension(5) :: expected = ['apple   ', 'bird    ', 'cat     ', 'dog     ', 'zebra   ']
    integer, dimension(5) :: perm
    integer :: stack_left(20), stack_right(20)
    integer :: i
    perm = [(i, i=1,5)]
    call sort_array(data, perm, stack_left, stack_right)
    call assert_true(all(data(perm) == expected), "test_sort_char_random: sorted values mismatch")
  end subroutine test_sort_char_random

  !> @brief Test sorting of an already sorted integer array.
  !! Checks that the sorted values remain unchanged.
  subroutine test_sort_sorted_stability()
    integer, dimension(5) :: data = [1,2,3,4,5]
    integer, dimension(5) :: perm
    integer :: stack_left(20), stack_right(20)
    integer :: i
    perm = [(i, i=1,5)]
    call sort_array(data, perm, stack_left, stack_right)
    call assert_equal_array_int(data(perm), [1,2,3,4,5], 5, "test_sort_sorted_stability: sorted values mismatch")
  end subroutine test_sort_sorted_stability

  !> @brief Test sorting of an empty array.
  !! Checks that sorting does not crash on empty input.
  subroutine test_sort_empty_array()
    integer, dimension(0) :: data
    integer, dimension(0) :: perm
    integer :: stack_left(1), stack_right(1)
    call sort_array(data, perm, stack_left, stack_right)
    ! No assertion needed: just check no crash
  end subroutine test_sort_empty_array

  !> @brief Test sorting of a large random integer array.
  !! Checks that the sorted values match the expected result.
  subroutine test_sort_large_random()
    real(8), allocatable :: rdata(:)
    integer, allocatable :: data(:), perm(:), sorted(:)
    integer :: n
    integer, allocatable :: stack_left(:), stack_right(:), dummy_perm(:)
    integer :: i

    n = 1000
    allocate(rdata(n), data(n), perm(n), sorted(n))
    allocate(stack_left(64), stack_right(64))
    allocate(dummy_perm(n))

    call random_seed()
    call random_number(rdata)
    data = int(rdata * 10000)

    sorted = data
    dummy_perm = [(i, i=1,n)]
    call sort_array(sorted, dummy_perm, stack_left, stack_right)
    sorted = sorted(dummy_perm)

    perm = [(i, i=1,n)]
    call sort_array(data, perm, stack_left, stack_right)

    call assert_equal_array_int(data(perm), sorted, n, "test_sort_large_random: sorted values mismatch")

    deallocate(rdata, data, perm, sorted, stack_left, stack_right, dummy_perm)
  end subroutine test_sort_large_random

end module mod_test_sorting