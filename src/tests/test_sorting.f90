!> @brief Unit test suite for tox_sorting module.
!! @details
!! This module tests sorting routines (real, integer, character)
!! using fixed data and known expected outputs.
module test_sorting
  use tox_sorting
  use, intrinsic :: iso_fortran_env, only: error_unit
  implicit none
  private
  public :: run_all_tests

contains

  subroutine run_all_tests()
    call test_sort_real()
    call test_sort_integer()
    call test_sort_character()
    call test_sort_integer_ascending()
    call test_sort_real_descending()
    call test_sort_char_random()
    call test_sort_sorted_stability()
    call test_sort_empty_array()
    call test_sort_large_random()
  end subroutine run_all_tests

  subroutine test_sort_real()
    real(8), dimension(5) :: data = [3.0d0, 1.0d0, 5.0d0, 2.0d0, 4.0d0]
    integer, dimension(5) :: perm, expected = [2, 4, 1, 5, 3]
    integer :: stack_left(20), stack_right(20), i
    write(*,*) 'Test: sort_real'
    write(*,*) 'Input :', data
    perm = [(i, i = 1, 5)]
    call sort_array(data, perm, stack_left, stack_right)
    write(*,*) 'Perm  :', perm
    write(*,*) 'Sorted:', data(perm)
    do i = 1, 5
      if (perm(i) /= expected(i)) then
        write(error_unit, *) 'FAIL: test_sort_real at i=', i, ' got=', perm(i), ' expected=', expected(i)
        return
      end if
    end do
    write(*,*) 'PASS: test_sort_real'
  end subroutine test_sort_real

  subroutine test_sort_integer()
    integer, dimension(4) :: data = [10, 3, 7, 1]
    integer, dimension(4) :: perm, expected = [4, 2, 3, 1]
    integer :: stack_left(20), stack_right(20), i
    write(*,*) 'Test: sort_integer'
    write(*,*) 'Input :', data
    perm = [(i, i = 1, 4)]
    call sort_array(data, perm, stack_left, stack_right)
    write(*,*) 'Perm  :', perm
    write(*,*) 'Sorted:', data(perm)
    do i = 1, 4
      if (perm(i) /= expected(i)) then
        write(error_unit, *) 'FAIL: test_sort_integer at i=', i, ' got=', perm(i), ' expected=', expected(i)
        return
      end if
    end do
    write(*,*) 'PASS: test_sort_integer'
  end subroutine test_sort_integer

  subroutine test_sort_character()
    character(len=6), dimension(3) :: data = ['delta ', 'alpha ', 'beta  ']
    integer, dimension(3) :: perm, expected = [2, 3, 1]
    integer :: stack_left(20), stack_right(20), i
    write(*,*) 'Test: sort_character'
    do i = 1, 3
      write(*,'(A,I0,A,A)') 'Input [', i, ']: ', trim(data(i))
    end do
    perm = [(i, i = 1, 3)]
    call sort_array(data, perm, stack_left, stack_right)
    do i = 1, 3
      write(*,'(A,I0,A,A)') 'Sorted [', i, ']: ', trim(data(perm(i)))
    end do
    do i = 1, 3
      if (perm(i) /= expected(i)) then
        write(error_unit, *) 'FAIL: test_sort_character at i=', i, ' got=', perm(i), ' expected=', expected(i)
        return
      end if
    end do
    write(*,*) 'PASS: test_sort_character'
  end subroutine test_sort_character

  subroutine test_sort_integer_ascending()
    integer, dimension(5) :: data = [5, 2, 9, 1, 6]
    integer, dimension(5) :: perm
    integer :: i, stack_left(20), stack_right(20)
    write(*,*) 'Test: sort_integer_ascending'
    write(*,*) 'Input :', data
    perm = [(i, i=1,5)]
    call sort_array(data, perm, stack_left, stack_right)
    write(*,*) 'Perm  :', perm
    write(*,*) 'Sorted:', data(perm)
    if (any(data(perm) /= [1,2,5,6,9])) then
      write(error_unit,*) 'FAIL: test_sort_integer_ascending'
    else
      write(*,*) 'PASS: test_sort_integer_ascending'
    end if
  end subroutine test_sort_integer_ascending

  subroutine test_sort_real_descending()
    real(8), dimension(5) :: data = [3.5d0, 2.2d0, 8.8d0, 1.1d0, 7.7d0]
    integer, dimension(5) :: perm
    integer :: i, stack_left(20), stack_right(20)
    write(*,*) 'Test: sort_real_descending'
    write(*,*) 'Input :', data
    perm = [(i, i=1,5)]
    call sort_array(data, perm, stack_left, stack_right)
    write(*,*) 'Perm  :', perm
    write(*,*) 'Sorted:', data(perm)
    if (any(data(perm) /= [1.1d0, 2.2d0, 3.5d0, 7.7d0, 8.8d0])) then
      write(error_unit,*) 'FAIL: test_sort_real_descending'
    else
      write(*,*) 'PASS: test_sort_real_descending'
    end if
  end subroutine test_sort_real_descending

  subroutine test_sort_char_random()
    character(len=8), dimension(5) :: data = ['dog     ', 'apple   ', 'zebra   ', 'cat     ', 'bird    ']
    character(len=8), dimension(5) :: expected = ['apple   ', 'bird    ', 'cat     ', 'dog     ', 'zebra   ']
    integer, dimension(5) :: perm
    integer :: i, stack_left(20), stack_right(20)
    write(*,*) 'Test: sort_char_random'
    do i = 1, 5
      write(*,'(A,I0,A,A)') 'Input [', i, ']: ', trim(data(i))
    end do
    perm = [(i, i=1,5)]
    call sort_array(data, perm, stack_left, stack_right)
    do i = 1, 5
      write(*,'(A,I0,A,A)') 'Sorted [', i, ']: ', trim(data(perm(i)))
    end do
    if (any(data(perm) /= expected)) then
      write(error_unit,*) 'FAIL: test_sort_char_random'
    else
      write(*,*) 'PASS: test_sort_char_random'
    end if
  end subroutine test_sort_char_random

  subroutine test_sort_sorted_stability()
    integer, dimension(5) :: data = [1,2,3,4,5]
    integer, dimension(5) :: perm
    integer :: i, stack_left(20), stack_right(20)
    write(*,*) 'Test: sort_sorted_stability'
    write(*,*) 'Input :', data
    perm = [(i, i=1,5)]
    call sort_array(data, perm, stack_left, stack_right)
    write(*,*) 'Perm  :', perm
    write(*,*) 'Sorted:', data(perm)
    if (any(data(perm) /= [1,2,3,4,5])) then
      write(error_unit,*) 'FAIL: test_sort_sorted_stability'
    else
      write(*,*) 'PASS: test_sort_sorted_stability'
    end if
  end subroutine test_sort_sorted_stability

  subroutine test_sort_empty_array()
    integer, dimension(0) :: data
    integer, dimension(0) :: perm
    integer :: stack_left(1), stack_right(1)
    write(*,*) 'Test: sort_empty_array'
    call sort_array(data, perm, stack_left, stack_right)
    write(*,*) 'PASS: test_sort_empty_array (no crash)'
  end subroutine test_sort_empty_array

  subroutine test_sort_large_random()
    real(8), allocatable :: rdata(:)
    integer, allocatable :: data(:), perm(:), sorted(:)
    integer :: i, n
    integer, allocatable :: stack_left(:), stack_right(:)
    integer, allocatable :: dummy_perm(:)

    n = 1000
    allocate(rdata(n), data(n), perm(n), sorted(n))
    allocate(stack_left(64), stack_right(64))
    allocate(dummy_perm(n))

    call random_seed()
    call random_number(rdata)
    data = int(rdata * 10000)

    write(*,*) 'Test: sort_large_random'
    write(*,*) 'Input (first 10):', data(1:10)

    sorted = data
    dummy_perm = [(i, i=1,n)]
    call sort_array(sorted, dummy_perm, stack_left, stack_right)
    sorted = sorted(dummy_perm)

    perm = [(i, i=1,n)]
    call sort_array(data, perm, stack_left, stack_right)

    write(*,*) 'Sorted (first 10):', data(perm(1:10))

    if (any(data(perm) /= sorted)) then
      write(error_unit,*) 'FAIL: test_sort_large_random'
    else
      write(*,*) 'PASS: test_sort_large_random'
    end if

    deallocate(rdata, data, perm, sorted, stack_left, stack_right, dummy_perm)
  end subroutine test_sort_large_random

end module test_sorting