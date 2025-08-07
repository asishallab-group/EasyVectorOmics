!> @brief Utility module for data analysis.
!> @details This module provides general-purpose utility functions for data analysis, to be used as needed.

module f42_utils
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none

  ! Expose specific functions C/R
  public :: sort_real, sort_integer, sort_character
  ! Only internal use since R and Python can't resolve interfaces (can't use sort_array)
  public :: sort_array

  !> Generic interface for internal Fortran use only
  interface sort_array
    module procedure sort_real, sort_integer, sort_character
  end interface sort_array
contains

    !> Sort a real(real64) array indirectly using quicksort.
  !>  Creates a sorted version of the array by reordering the `perm` vector.
  !> The original data in `array` remains unchanged.
  !> 
  !> @param array        Real input array to sort<br>
  !> @param perm         Permutation vector that will be sorted<br>
  !> @param stack_left   Manual stack of left indices for quicksort recursion<br>
  !> @param stack_right  Manual stack of right indices for quicksort recursion<br>
  pure subroutine sort_real(array, perm, stack_left, stack_right)
    real(real64), intent(in) :: array(:)
    integer, intent(inout) :: perm(:)
    integer, intent(inout) :: stack_left(:), stack_right(:)
    call quicksort_real(array, perm, size(array), stack_left, stack_right)
  end subroutine sort_real

  !> Sort an integer array indirectly using quicksort.
  !>  Similar to `sort_real`, but for integer input.
  pure subroutine sort_integer(array, perm, stack_left, stack_right)
    integer, intent(in) :: array(:)
    integer, intent(inout) :: perm(:)
    integer, intent(inout) :: stack_left(:), stack_right(:)
    call quicksort_int(array, perm, size(array), stack_left, stack_right)
  end subroutine sort_integer

  !> Sort a character array indirectly using quicksort.
  !>  Uses lexicographic ordering and permutation vector sorting.
  pure subroutine sort_character(array, perm, stack_left, stack_right)
    character(len=*), intent(in) :: array(:)
    integer, intent(inout) :: perm(:)
    integer, intent(inout) :: stack_left(:), stack_right(:)
    call quicksort_char(array, perm, size(array), stack_left, stack_right)
  end subroutine sort_character

  !> Internal quicksort implementation for real arrays.
  !>  Sorts indirectly using the permutation vector `perm`. Manual stack replaces recursion.
  pure subroutine quicksort_real(array, perm, n, stack_left, stack_right)
    real(real64), intent(in) :: array(:)
    integer, intent(inout) :: perm(:)
    integer, intent(in) :: n
    integer, intent(inout) :: stack_left(:), stack_right(:)
    integer :: left, right, i, j, top, pivot_idx
    real(real64) :: pivot_val

    top = 1
    stack_left(top) = 1
    stack_right(top) = n

    ! Iterative quicksort using explicit stack
    do while (top > 0)
      left = stack_left(top)
      right = stack_right(top)
      top = top - 1

      if (left >= right) cycle

      ! Select pivot and initialize pointers
      pivot_idx = (left + right) / 2
      pivot_val = array(perm(pivot_idx))
      i = left
      j = right

      ! Partitioning loop
      do
        do while (array(perm(i)) < pivot_val)
          i = i + 1
        end do
        do while (array(perm(j)) > pivot_val)
          j = j - 1
        end do
        if (i <= j) then
          call swap_int(perm(i), perm(j))
          i = i + 1
          j = j - 1
        end if
        if (i > j) exit
      end do

      ! Push new ranges onto stack
      if (left < j) then
        top = top + 1
        stack_left(top) = left
        stack_right(top) = j
      end if
      if (i < right) then
        top = top + 1
        stack_left(top) = i
        stack_right(top) = right
      end if
    end do
  end subroutine quicksort_real

  !> Internal quicksort implementation for integer arrays.
  !>  Indirectly sorts `array` using `perm`, same algorithm as `quicksort_real`.
  pure subroutine quicksort_int(array, perm, n, stack_left, stack_right)
    integer, intent(in) :: array(:)
    integer, intent(inout) :: perm(:)
    integer, intent(in) :: n
    integer, intent(inout) :: stack_left(:), stack_right(:)
    integer :: left, right, i, j, top, pivot_idx
    integer :: pivot_val

    top = 1
    stack_left(top) = 1
    stack_right(top) = n

    do while (top > 0)
      left = stack_left(top)
      right = stack_right(top)
      top = top - 1

      if (left >= right) cycle

      pivot_idx = (left + right) / 2
      pivot_val = array(perm(pivot_idx))
      i = left
      j = right

      do
        do while (array(perm(i)) < pivot_val)
          i = i + 1
        end do
        do while (array(perm(j)) > pivot_val)
          j = j - 1
        end do
        if (i <= j) then
          call swap_int(perm(i), perm(j))
          i = i + 1
          j = j - 1
        end if
        if (i > j) exit
      end do

      if (left < j) then
        top = top + 1
        stack_left(top) = left
        stack_right(top) = j
      end if
      if (i < right) then
        top = top + 1
        stack_left(top) = i
        stack_right(top) = right
      end if
    end do
  end subroutine quicksort_int

  !> Internal quicksort implementation for character arrays.
  !>  Lexicographic quicksort using string comparison, indirect via `perm`.
  pure subroutine quicksort_char(array, perm, n, stack_left, stack_right)
    character(len=*), intent(in) :: array(:)
    integer, intent(inout) :: perm(:)
    integer, intent(in) :: n
    integer, intent(inout) :: stack_left(:), stack_right(:)
    integer :: left, right, i, j, top, pivot_idx
    character(len=len(array)) :: pivot_val

    top = 1
    stack_left(top) = 1
    stack_right(top) = n

    do while (top > 0)
      left = stack_left(top)
      right = stack_right(top)
      top = top - 1

      if (left >= right) cycle

      pivot_idx = (left + right) / 2
      pivot_val = array(perm(pivot_idx))
      i = left
      j = right

      do
        do while (array(perm(i)) < pivot_val)
          i = i + 1
        end do
        do while (array(perm(j)) > pivot_val)
          j = j - 1
        end do
        if (i <= j) then
          call swap_int(perm(i), perm(j))
          i = i + 1
          j = j - 1
        end if
        if (i > j) exit
      end do

      if (left < j) then
        top = top + 1
        stack_left(top) = left
        stack_right(top) = j
      end if
      if (i < right) then
        top = top + 1
        stack_left(top) = i
        stack_right(top) = right
      end if
    end do
  end subroutine quicksort_char

  !> Swap two integer values in-place.
  !> @param a First integer to swap<br>
  !> @param b Second integer to swap<br>
  pure subroutine swap_int(a, b)
    integer, intent(inout) :: a, b
    integer :: temp
    temp = a; a = b; b = temp
  end subroutine swap_int


end module f42_utils



! === R WRAPPERS ===
subroutine sort_real_r(array, perm, stack_left, stack_right, n)
  use, intrinsic :: iso_fortran_env, only: real64
  use f42_utils, only: sort_real
  implicit none
  integer, intent(in) :: n
  real(real64), intent(in) :: array(n)
  integer, intent(inout) :: perm(n), stack_left(n), stack_right(n)

  call sort_real(array, perm, stack_left, stack_right)
end subroutine sort_real_r

subroutine sort_integer_r(array, perm, stack_left, stack_right, n)
  use f42_utils, only: sort_integer
  implicit none
  integer, intent(in) :: n
  integer, intent(in) :: array(n)
  integer, intent(inout) :: perm(n), stack_left(n), stack_right(n)

  call sort_integer(array, perm, stack_left, stack_right)
end subroutine sort_integer_r

subroutine sort_character_r(char_matrix, perm, stack_left, stack_right, n, strlen)
  use f42_utils, only: sort_character
  implicit none
  integer, intent(in) :: n, strlen
  character(len=1), intent(in) :: char_matrix(strlen, n)
  integer, intent(inout) :: perm(n), stack_left(n), stack_right(n)

  character(len=strlen) :: array(n)
  integer :: i, j

  ! Reconstruct each full string from the character matrix
  do i = 1, n
    do j = 1, strlen
      array(i)(j:j) = char_matrix(j, i)
    end do
  end do

  ! Call the original sorting routine
  call sort_character(array, perm, stack_left, stack_right)
end subroutine


! === C WRAPPERS ===
subroutine sort_real_c(array, perm, stack_left, stack_right, n) bind(C, name="sort_real_c")
  use iso_c_binding
  use f42_utils, only: sort_real
  implicit none
  integer(c_int), intent(in), value :: n
  real(c_double), intent(in) :: array(n)
  integer(c_int), intent(inout) :: perm(n)
  integer(c_int), intent(inout) :: stack_left(n)
  integer(c_int), intent(inout) :: stack_right(n)

  call sort_real(array, perm, stack_left, stack_right)
end subroutine sort_real_c

subroutine sort_integer_c(array, perm, stack_left, stack_right, n) bind(C, name="sort_integer_c")
  use iso_c_binding
  use f42_utils, only: sort_integer
  implicit none
  integer(c_int), intent(in), value :: n
  integer(c_int), intent(in) :: array(n)
  integer(c_int), intent(inout) :: perm(n)
  integer(c_int), intent(inout) :: stack_left(n)
  integer(c_int), intent(inout) :: stack_right(n)

  call sort_integer(array, perm, stack_left, stack_right)
end subroutine sort_integer_c

subroutine sort_character_c(char_matrix, perm, stack_left, stack_right, n, strlen) bind(C, name="sort_character_c")
  use iso_c_binding
  use f42_utils, only: sort_character
  implicit none
  integer(c_int), intent(in), value :: n, strlen
  character(kind=c_char,len=1), intent(in) :: char_matrix(strlen, n)
  integer(c_int), intent(inout) :: perm(n)
  integer(c_int), intent(inout) :: stack_left(n)
  integer(c_int), intent(inout) :: stack_right(n)

  character(len=strlen) :: array(n)
  integer :: i, j

  ! Reconstruct Fortran strings from the C-style char matrix
  do i = 1, n
    do j = 1, strlen
      array(i)(j:j) = char_matrix(j, i)
    end do
  end do

  call sort_character(array, perm, stack_left, stack_right)
end subroutine sort_character_c