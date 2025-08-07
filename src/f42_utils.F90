!> Utility module for data analysis.
!| This module provides general-purpose utility functions for data analysis, to be used as needed.

module f42_utils
  use, intrinsic :: iso_fortran_env, only: real64, int32
  implicit none

  public :: sort_real, sort_integer, sort_character
  public :: sort_array

  interface sort_array
    module procedure sort_real, sort_integer, sort_character
  end interface sort_array
contains

  !> Sort a real array indirectly using quicksort.
  !| Creates a sorted version of the array by reordering the `perm` vector. The original data in `array` remains unchanged.
  pure subroutine sort_real(array, perm, stack_left, stack_right)
    !| Real input array to sort
    real(real64), intent(in) :: array(:)
    !| Permutation vector that will be sorted
    integer(int32), intent(inout) :: perm(:)
    !| Manual stack of left indices for quicksort recursion
    integer(int32), intent(inout) :: stack_left(:)
    !| Manual stack of right indices for quicksort recursion
    integer(int32), intent(inout) :: stack_right(:)
    call quicksort_real(array, perm, size(array), stack_left, stack_right)
  end subroutine sort_real

  !> Sort an integer array indirectly using quicksort.
  !| Similar to `sort_real`, but for integer input.
  pure subroutine sort_integer(array, perm, stack_left, stack_right)
    !| Integer input array to sort
    integer(int32), intent(in) :: array(:)
    !| Permutation vector that will be sorted
    integer(int32), intent(inout) :: perm(:)
    !| Manual stack of left indices for quicksort recursion
    integer(int32), intent(inout) :: stack_left(:)
    !| Manual stack of right indices for quicksort recursion
    integer(int32), intent(inout) :: stack_right(:)
    call quicksort_int(array, perm, size(array), stack_left, stack_right)
  end subroutine sort_integer

  !> Sort a character array indirectly using quicksort.
  !| Uses lexicographic ordering and permutation vector sorting.
  pure subroutine sort_character(array, perm, stack_left, stack_right)
    !| Character input array to sort
    character(len=*), intent(in) :: array(:)
    !| Permutation vector that will be sorted
    integer(int32), intent(inout) :: perm(:)
    !| Manual stack of left indices for quicksort recursion
    integer(int32), intent(inout) :: stack_left(:)
    !| Manual stack of right indices for quicksort recursion
    integer(int32), intent(inout) :: stack_right(:)
    call quicksort_char(array, perm, size(array), stack_left, stack_right)
  end subroutine sort_character

  !> Internal quicksort implementation for real arrays.
  !| Sorts indirectly using the permutation vector `perm`. Manual stack replaces recursion.
  pure subroutine quicksort_real(array, perm, n, stack_left, stack_right)
    !| Real input array to sort
    real(real64), intent(in) :: array(:)
    !| Permutation vector that will be sorted
    integer(int32), intent(inout) :: perm(:)
    !| Size of the array
    integer(int32), intent(in) :: n
    !| Manual stack of left indices for quicksort recursion
    integer(int32), intent(inout) :: stack_left(:)
    !| Manual stack of right indices for quicksort recursion
    integer(int32), intent(inout) :: stack_right(:)

    integer(int32) :: left, right, i, j, top, pivot_idx
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
  !| Indirectly sorts `array` using `perm`, same algorithm as `quicksort_real`.
  pure subroutine quicksort_int(array, perm, n, stack_left, stack_right)
    !| Integer input array to sort
    integer(int32), intent(in) :: array(:)
    !| Permutation vector that will be sorted
    integer(int32), intent(inout) :: perm(:)
    !| Size of the array
    integer(int32), intent(in) :: n
    !| Manual stack of left indices for quicksort recursion
    integer(int32), intent(inout) :: stack_left(:)
    !| Manual stack of right indices for quicksort recursion
    integer(int32), intent(inout) :: stack_right(:)

    integer(int32) :: left, right, i, j, top, pivot_idx
    integer(int32) :: pivot_val

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
  !| Lexicographic quicksort using string comparison, indirect via `perm`.
  pure subroutine quicksort_char(array, perm, n, stack_left, stack_right)
    !| Character input array to sort
    character(len=*), intent(in) :: array(:)
    !| Permutation vector that will be sorted
    integer(int32), intent(inout) :: perm(:)
    !| Size of the array
    integer(int32), intent(in) :: n
    !| Manual stack of left indices for quicksort recursion
    integer(int32), intent(inout) :: stack_left(:)
    !| Manual stack of right indices for quicksort recursion
    integer(int32), intent(inout) :: stack_right(:)

    integer(int32) :: left, right, i, j, top, pivot_idx
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
  pure subroutine swap_int(a, b)
    !| First integer to swap
    integer(int32), intent(inout) :: a
    !| Second integer to swap
    integer(int32), intent(inout) :: b
    integer(int32) :: temp
    temp = a; a = b; b = temp
  end subroutine swap_int

  !> Finds the indices of the true values in a logical mask.
  pure subroutine which(mask, n, idx_out, m_max, m_out)
    !| Logical array of size n.
    logical, intent(in) :: mask(:)
    !| Size of the mask.
    integer(int32), intent(in) :: n
    !| Integer array to store the indices of true values.
    integer(int32), intent(out) :: idx_out(:)
    !| Maximum size of idx_out.
    integer(int32), intent(in) :: m_max
    !| Actual size of idx_out (number of true values found).
    integer(int32), intent(out) :: m_out
    integer(int32) :: i, count
    count = 0
    idx_out = 0  ! Initialize to avoid garbage values
    do i = 1, n
      if (mask(i)) then
        count = count + 1
        if (count <= m_max) then
          idx_out(count) = i
        end if
      end if
    end do
    m_out = count
  end subroutine which

  !> Performs LOESS smoothing on a set of data points.
  !| Smooths y_ref at x_query using reference points x_ref, y_ref, and kernel parameters. Optionally excludes points using mask_in.
  pure subroutine loess_smooth_2d(n_total, n_target, x_ref, y_ref, indices_used, x_query, &
                          kernel_sigma, kernel_cutoff, y_out, workspace_weights, &
                          workspace_values, mask_in)
    !| Total number of reference points.
    integer(int32), intent(in) :: n_total
    !| Number of target points to smooth.
    integer(int32), intent(in) :: n_target
    !| Reference x-coordinates.
    real(real64), intent(in) :: x_ref(n_total)
    !| Reference y-coordinates (length n_total).
    real(real64), intent(in) :: y_ref(n_total)
    !| Target x-coordinates to smooth.
    real(real64), intent(in) :: x_query(n_target)
    !| Indices of reference points used for smoothing.
    integer(int32), intent(in) :: indices_used(n_total)
    !| Bandwidth parameter for the kernel.
    real(real64), intent(in) :: kernel_sigma
    !| Cutoff for the kernel.
    real(real64), intent(in) :: kernel_cutoff
    !| Output smoothed values (length n_target).
    real(real64), intent(out) :: y_out(n_target)
    !| Temporary array for weights.
    real(real64), intent(inout) :: workspace_weights(n_total)
    !| Temporary array for values.
    real(real64), intent(inout) :: workspace_values(n_total)
    !| Optional logical mask for explicit exclusion.
    logical, intent(in), optional :: mask_in(n_total)

    integer(int32) :: q, i, idx, m_out, valid_indices(n_total)
    real(real64) :: query_x, ref_x, delta, sum_weights, weight
    logical :: mask(n_total)
    logical :: found_exact
    integer(int32) :: exact_idx
    real(real64) :: min_dist
    integer(int32) :: min_idx, j
    integer(int32) :: n_mask_true, last_true_idx
    logical :: in_range(n_total)

    do q = 1, n_target
      query_x = x_query(q)
      sum_weights = 0.0
      y_out(q) = 0.0
      if (present(mask_in)) then
        do i = 1, n_total
          idx = indices_used(i)
          mask(i) = mask_in(idx)
        end do
        ! If the mask leaves only one allowed point, return it directly and cycle
        n_mask_true = 0
        last_true_idx = -1
        do i = 1, n_total
          if (mask(i)) then
            n_mask_true = n_mask_true + 1
            last_true_idx = indices_used(i)
          end if
        end do
        if (n_mask_true == 1) then
          y_out(q) = y_ref(last_true_idx)
          cycle
        end if
      else
        mask = .true.
      end if
      ! --- Check for exact match ---
      found_exact = .false.
      exact_idx = -1
      do i = 1, n_total
        idx = indices_used(i)
        if (query_x == x_ref(idx)) then
          found_exact = .true.
          exact_idx = idx
          exit
        end if
      end do
      if (found_exact) then
        y_out(q) = y_ref(exact_idx)
        cycle
      end if
      ! --- Usual LOESS smoothing ---
      
      do i = 1, n_total
        idx = indices_used(i)
        ref_x = x_ref(idx)
        delta = abs(query_x - ref_x)
        in_range(i) = (delta <= kernel_cutoff * kernel_sigma)
        mask(i) = mask(i) .and. in_range(i)
      end do
      call which(mask, n_total, valid_indices, n_total, m_out)
      if (m_out == 1) then
        idx = indices_used(valid_indices(1))
        y_out(q) = y_ref(idx)
        cycle
      end if
      do i = 1, m_out
        idx = indices_used(valid_indices(i))
        ref_x = x_ref(idx)
        delta = abs(query_x - ref_x)
        weight = exp(-(delta / kernel_sigma)**2)
        sum_weights = sum_weights + weight
        y_out(q) = y_out(q) + weight * y_ref(idx)
      end do
      if (sum_weights > 0.0) then
        y_out(q) = y_out(q) / sum_weights
      else
        ! Extrapolation: return y_ref at the closest x_ref (among indices_used)
        min_dist = abs(query_x - x_ref(indices_used(1)))
        min_idx = indices_used(1)
        do j = 2, n_total
          if (abs(query_x - x_ref(indices_used(j))) < min_dist) then
            min_dist = abs(query_x - x_ref(indices_used(j)))
            min_idx = indices_used(j)
          end if
        end do
        y_out(q) = y_ref(min_idx)
      end if
    end do
  end subroutine loess_smooth_2d

end module f42_utils



! === R WRAPPERS ===

!> R wrapper for loess_smooth_2d.
!| Calls loess_smooth_2d with standard Fortran types for R interface.
subroutine loess_smooth_2d_r(n_total, n_target, x_ref, y_ref, indices_used, x_query, &
    kernel_sigma, kernel_cutoff, y_out, workspace_weights, workspace_values, mask_in)
  use f42_utils, only: loess_smooth_2d
  use, intrinsic :: iso_fortran_env, only: real64, int32
  implicit none
  !| Total number of reference points.
  integer(int32), intent(in) :: n_total
  !| Number of target points to smooth.
  integer(int32), intent(in) :: n_target
  !| Reference x-coordinates.
  real(real64), intent(in) :: x_ref(n_total)
  !| Reference y-coordinates (length n_total).
  real(real64), intent(in) :: y_ref(n_total)
  !| Target x-coordinates to smooth.
  real(real64), intent(in) :: x_query(n_target)
  !| Indices of reference points used for smoothing.
  integer(int32), intent(in) :: indices_used(n_total)
  !| Bandwidth parameter for the kernel.
  real(real64), intent(in) :: kernel_sigma
  !| Cutoff for the kernel.
  real(real64), intent(in) :: kernel_cutoff
  !| Output smoothed values (length n_target).
  real(real64), intent(out) :: y_out(n_target)
  !| Temporary array for weights.
  real(real64), intent(inout) :: workspace_weights(n_total)
  !| Temporary array for values.
  real(real64), intent(inout) :: workspace_values(n_total)
  !| Logical mask for explicit exclusion.
  logical, intent(in) :: mask_in(n_total)
  integer(int32) :: i
  logical :: any_true
  any_true = .false.
  do i = 1, n_total
    if (mask_in(i)) then
      any_true = .true.
      exit
    end if
  end do
  if (any_true) then
    call loess_smooth_2d(n_total, n_target, x_ref, y_ref, indices_used, x_query, &
      kernel_sigma, kernel_cutoff, y_out, workspace_weights, workspace_values, mask_in)
  else
    call loess_smooth_2d(n_total, n_target, x_ref, y_ref, indices_used, x_query, &
      kernel_sigma, kernel_cutoff, y_out, workspace_weights, workspace_values)
  end if
end subroutine loess_smooth_2d_r

! === C WRAPPERS ===

!> C wrapper for which.
!| Converts integer mask to logical and calls which.
subroutine which_c(mask, n, idx_out, m_max, m_out) bind(C, name="which_c")
  use iso_c_binding
  use, intrinsic :: iso_fortran_env, only: int32
  use f42_utils, only: which
  implicit none
  !| Size of the mask.
  integer(c_int), intent(in), value :: n
  !| Maximum size of idx_out.
  integer(c_int), intent(in), value :: m_max
  !| Integer mask array (0/1 values).
  integer(c_int), intent(in) :: mask(n)
  !| Output array for indices of true values.
  integer(c_int), intent(out) :: idx_out(m_max)
  !| Actual size of idx_out (number of true values found).
  integer(c_int), intent(out) :: m_out
  logical :: mask_f(n)
  integer(int32) :: i
  do i = 1, n
    mask_f(i) = (mask(i) /= 0)
  end do
  call which(mask_f, n, idx_out, m_max, m_out)
end subroutine which_c

!> C wrapper for loess_smooth_2d.
!| Converts integer mask to logical and calls loess_smooth_2d for C interface.
subroutine loess_smooth_2d_c(n_total, n_target, x_ref, y_ref, indices_used, x_query, &
    kernel_sigma, kernel_cutoff, y_out, workspace_weights, workspace_values, mask_in) bind(C, name="loess_smooth_2d_c")
  use iso_c_binding
  use, intrinsic :: iso_fortran_env, only: int32
  use f42_utils, only: loess_smooth_2d
  implicit none
  !| Total number of reference points.
  integer(c_int), intent(in), value :: n_total
  !| Number of target points to smooth.
  integer(c_int), intent(in), value :: n_target
  !| Reference x-coordinates.
  real(c_double), intent(in) :: x_ref(n_total)
  !| Reference y-coordinates (length n_total).
  real(c_double), intent(in) :: y_ref(n_total)
  !| Target x-coordinates to smooth.
  real(c_double), intent(in) :: x_query(n_target)
  !| Indices of reference points used for smoothing.
  integer(c_int), intent(in) :: indices_used(n_total)
  !| Bandwidth parameter for the kernel.
  real(c_double), intent(in), value :: kernel_sigma
  !| Cutoff for the kernel.
  real(c_double), intent(in), value :: kernel_cutoff
  !| Output smoothed values (length n_target).
  real(c_double), intent(out) :: y_out(n_target)
  !| Temporary array for weights.
  real(c_double), intent(inout) :: workspace_weights(n_total)
  !| Temporary array for values.
  real(c_double), intent(inout) :: workspace_values(n_total)
  !| Integer mask array (0/1 values).
  integer(c_int), intent(in) :: mask_in(n_total)
  logical :: mask_f(n_total)
  integer(int32) :: i
  logical :: all_true

  do i = 1, size(mask_in)
    mask_f(i) = (mask_in(i) /= 0)
  end do

  call loess_smooth_2d(n_total, n_target, x_ref, y_ref, indices_used, x_query, &
    kernel_sigma, kernel_cutoff, y_out, workspace_weights, workspace_values, mask_f)

end subroutine loess_smooth_2d_c
