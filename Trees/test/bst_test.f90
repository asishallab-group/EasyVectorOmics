! filepath: ./tensor-omics/BST/src/bst_test.f90
!> \file
!! \brief Unit test for the binary_search_tree module.
!!
!! Tests BST index construction and verifies monotonicity of sorted access.
program test_bst
  use binary_search_tree
  implicit none
  integer, parameter :: n = 1000
  real(8) :: x(n)
  integer :: ix(n)
  integer :: stack_left(n), stack_right(n)
  integer :: i, res_ix(n), res_n
  logical :: is_sorted
  real(8) :: val

  call random_array(x, n)
  call build_bst_index(x, n, ix, stack_left, stack_right)

  ! Test 1: Check monotonicity of x(ix)
  is_sorted = .true.
  do i = 2, n
    if (x(ix(i)) < x(ix(i-1))) then
      is_sorted = .false.
      exit
    end if
  end do
  if (is_sorted) then
    print *, 'BST index test PASSED: x(ix) is monotonic non-decreasing.'
  else
    print *, 'BST index test FAILED: x(ix) is not monotonic.'
  end if

  ! Test 2: get_sorted_value
  val = get_sorted_value(x, ix, n/2)
  print *, 'get_sorted_value(x, ix, n/2) = ', val

  ! Test 3: bst_range_query
  call bst_range_query(x, ix, n, 0.2d0, 0.8d0, res_ix, res_n)
  print *, 'bst_range_query: Number of values in [0.2, 0.8] =', res_n
  if (res_n > 0) then
    print *, 'First value in range: ', x(res_ix(1))
    print *, 'Last value in range: ', x(res_ix(res_n))
  end if

contains

  !> \brief Fill an array with random real values in [0,1).
  !! \param arr Output real array.
  !! \param n Number of elements.
  subroutine random_array(arr, nval)
    real(8), intent(out) :: arr(:)
    integer, intent(in) :: nval
    integer :: j
    call random_seed()
    do j = 1, nval
      call random_number(arr(j))
    end do
  end subroutine random_array

end program test_bst