!> TOX sorting module for TensorOmics.
!> 
!> Provides non-recursive quicksort implementations for real, integer, and character arrays.
!> Sorting is performed indirectly via a permutation vector, preserving the original arrays.
!> Suitable for expression data analysis workflows requiring stable and indirect sorting.
module tox_sorting
  implicit none
  private

  public :: sort_array

  !> Generic interface to sort real, integer, or character arrays.
  !>  Calls one of the module procedures `sort_real`, `sort_integer`, or `sort_character`
  !> depending on the input array type. Sorting is done indirectly through a permutation vector.
  interface sort_array
    module procedure sort_real, sort_integer, sort_character
  end interface sort_array

contains

  !> Sort a real(8) array indirectly using quicksort.
  !>  Creates a sorted version of the array by reordering the `perm` vector.
  !> The original data in `array` remains unchanged.
  !> 
  !> @param array        Real input array to sort<br>
  !> @param perm         Permutation vector that will be sorted<br>
  !> @param stack_left   Manual stack of left indices for quicksort recursion<br>
  !> @param stack_right  Manual stack of right indices for quicksort recursion<br>
  subroutine sort_real(array, perm, stack_left, stack_right)
    real(8), intent(in) :: array(:)
    integer, intent(inout) :: perm(:)
    integer, intent(inout) :: stack_left(:), stack_right(:)
    call quicksort_real(array, perm, size(array), stack_left, stack_right)
  end subroutine sort_real

  !> Sort an integer array indirectly using quicksort.
  !>  Similar to `sort_real`, but for integer input.
  subroutine sort_integer(array, perm, stack_left, stack_right)
    integer, intent(in) :: array(:)
    integer, intent(inout) :: perm(:)
    integer, intent(inout) :: stack_left(:), stack_right(:)
    call quicksort_int(array, perm, size(array), stack_left, stack_right)
  end subroutine sort_integer

  !> Sort a character array indirectly using quicksort.
  !>  Uses lexicographic ordering and permutation vector sorting.
  subroutine sort_character(array, perm, stack_left, stack_right)
    character(len=*), intent(in) :: array(:)
    integer, intent(inout) :: perm(:)
    integer, intent(inout) :: stack_left(:), stack_right(:)
    call quicksort_char(array, perm, size(array), stack_left, stack_right)
  end subroutine sort_character

  !> Internal quicksort implementation for real arrays.
  !>  Sorts indirectly using the permutation vector `perm`. Manual stack replaces recursion.
  subroutine quicksort_real(array, perm, n, stack_left, stack_right)
    real(8), intent(in) :: array(:)
    integer, intent(inout) :: perm(:)
    integer, intent(in) :: n
    integer, intent(inout) :: stack_left(:), stack_right(:)
    integer :: left, right, i, j, top, pivot_idx
    real(8) :: pivot_val

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
  subroutine quicksort_int(array, perm, n, stack_left, stack_right)
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
  subroutine quicksort_char(array, perm, n, stack_left, stack_right)
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
  subroutine swap_int(a, b)
    integer, intent(inout) :: a, b
    integer :: temp
    temp = a; a = b; b = temp
  end subroutine swap_int

end module tox_sorting
