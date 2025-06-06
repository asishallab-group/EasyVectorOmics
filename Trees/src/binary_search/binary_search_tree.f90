! filepath: /home/aaron/EasyVectorOmics/tensor-omics/BST/src/binary_search_tree.f90
!> \file
!! \brief Flat-index-based Binary Search Tree (BST) utilities for 1D range queries.
!!
!! This module provides routines to build a BST index (via sorting), access sorted values,
!! and perform efficient range queries over a real-valued array.
module binary_search_tree
  use tox_sorting
  implicit none
  public :: build_bst_index, get_sorted_value, bst_range_query
contains
  
  !> \brief Build the BST index by sorting indices using values in x.
  !! \param x Input real array to be indexed.
  !! \param n Number of elements in x.
  !! \param ix Output integer array holding the permutation index such that x(ix) is sorted.
  !! \param stack_left Manual stack of left indices for quicksort recursion (preallocated)
  !! \param stack_right Manual stack of right indices for quicksort recursion (preallocated)
  subroutine build_bst_index(x, n, ix, stack_left, stack_right)
    real(8), intent(in) :: x(:)
    integer, intent(in) :: n
    integer, intent(out) :: ix(:)
    integer, intent(inout) :: stack_left(:), stack_right(:)
    integer :: i

    print *, 'build_bst_index: size(x) =', size(x)
    print *, 'build_bst_index: size(ix) =', size(ix)
    print *, 'build_bst_index: size(stack_left) =', size(stack_left)
    print *, 'build_bst_index: size(stack_right) =', size(stack_right)

    do i = 1, n
      ix(i) = i
    end do
    call sort_array(x, ix, stack_left, stack_right)
  end subroutine build_bst_index

  !> \brief Get the value at the sorted position.
  !! \param x Input real array.
  !! \param ix Permutation index array.
  !! \param i Sorted position (1-based).
  !! \return Value x(ix(i)), i.e., the i-th smallest value in x.
  function get_sorted_value(x, ix, i) result(val)
    real(8), intent(in) :: x(:)
    integer, intent(in) :: ix(:)
    integer, intent(in) :: i
    real(8) :: val
    val = x(ix(i))
  end function get_sorted_value

  !> \brief Perform a 1D range query over the sorted index.
  !! \param x Input real array.
  !! \param ix Permutation index array (sorted).
  !! \param n Number of elements.
  !! \param lo Lower bound of range (inclusive).
  !! \param hi Upper bound of range (inclusive).
  !! \param out_ix Output array of indices in x that match the range.
  !! \param out_n Number of matches found.
  subroutine bst_range_query(x, ix, n, lo, hi, out_ix, out_n)
    real(8), intent(in) :: x(:)
    integer, intent(in) :: ix(:)
    integer, intent(in) :: n
    real(8), intent(in) :: lo, hi
    integer, intent(out) :: out_ix(:)
    integer, intent(out) :: out_n
    integer :: i
    out_n = 0
    do i = 1, n
      if (x(ix(i)) >= lo .and. x(ix(i)) <= hi) then
        out_n = out_n + 1
        out_ix(out_n) = ix(i)
      else if (x(ix(i)) > hi) then
        exit
      end if
    end do
  end subroutine bst_range_query

end module binary_search_tree

!> \brief Wrapper for getting range query usable by R
subroutine bst_range_query_r(x, ix, n, lo, hi, out_ix, out_n)
  use binary_search_tree
  implicit none
  real(8), intent(in) :: x(n)
  integer, intent(in) :: ix(n)
  integer, intent(in) :: n
  real(8), intent(in) :: lo, hi
  integer, intent(out) :: out_ix(n)
  integer, intent(out) :: out_n

  call bst_range_query(x, ix, n, lo, hi, out_ix, out_n)
end subroutine bst_range_query_r

!> \brief Wrapper for building BST index usable by R
subroutine build_bst_index_r(x, n, ix, stack_left, stack_right)
  use binary_search_tree
  implicit none
  integer, intent(in) :: n
  real(8), intent(in) :: x(n)

  integer, intent(out) :: ix(n)
  integer, intent(inout) :: stack_left(n), stack_right(n)

  print *, 'Size n:', n
  print *, 'Size x:', size(x)
  print *, 'Size ix:', size(ix)
  print *, 'Size stack_left:', size(stack_left)
  print *, 'Size stack_right:', size(stack_right)
  call build_bst_index(x, n, ix, stack_left, stack_right)
end subroutine build_bst_index_r