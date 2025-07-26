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

  !> Finds the indices of the true values in a logical mask.
  !! @param mask Logical array of size n.
  !! @param n Size of the mask.
  !! @param idx_out Integer array to store the indices of true values.
  !! @param m_max Maximum size of idx_out.
  !! @param m_out Actual size of idx_out (number of true values found).
  pure subroutine which(mask, n, idx_out, m_max, m_out)
    logical, intent(in) :: mask(:)
    integer, intent(in) :: n
    integer, intent(out) :: idx_out(:)
    integer, intent(in) :: m_max
    integer, intent(out) :: m_out
    integer :: i, count
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
  !! @param n_total Total number of reference points.
  !! @param n_target Number of target points to smooth.
  !! @param x_ref Reference x-coordinates.
  !! @param y_ref Reference y-coordinates (length n_total).
  !! @param indices_used Indices of reference points used for smoothing.
  !! @param x_query Target x-coordinates to smooth.
  !! @param kernel_sigma Bandwidth parameter for the kernel.
  !! @param kernel_cutoff Cutoff for the kernel.
  !! @param y_out Output smoothed values (length n_target).
  !! @param workspace_weights Temporary array for weights.
  !! @param workspace_values Temporary array for values.
  !! @param mask_in Optional logical mask for explicit exclusion.
  pure subroutine loess_smooth_2d(n_total, n_target, x_ref, y_ref, indices_used, x_query, &
                          kernel_sigma, kernel_cutoff, y_out, workspace_weights, &
                          workspace_values, mask_in)

    integer, intent(in) :: n_total, n_target
    real(real64), intent(in) :: x_ref(n_total), y_ref(n_total), x_query(n_target)
    integer, intent(in) :: indices_used(n_total)
    real(real64), intent(in) :: kernel_sigma, kernel_cutoff
    real(real64), intent(out) :: y_out(n_target)
    real(real64), intent(inout) :: workspace_weights(n_total), workspace_values(n_total)
    logical, intent(in), optional :: mask_in(n_total)

    integer :: q, i, idx, m_out, valid_indices(n_total)
    real(real64) :: query_x, ref_x, delta, sum_weights, weight
    logical :: mask(n_total)
    logical :: found_exact
    integer :: exact_idx
    real(real64) :: min_dist
    integer :: min_idx, j
    integer :: n_mask_true, last_true_idx
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
        ! Si la máscara deja solo un punto permitido, devolverlo directamente y cycle
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
      end do
      do i = 1, n_total
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

subroutine which_r(mask, n, idx_out, m_max, m_out)
  use f42_utils, only: which
  implicit none
  integer, intent(in) :: n, m_max
  logical, intent(in) :: mask(n)
  integer, intent(out) :: idx_out(m_max), m_out
  call which(mask, n, idx_out, m_max, m_out)
end subroutine which_r

subroutine loess_smooth_2d_r(n_total, n_target, x_ref, y_ref, indices_used, x_query, &
    kernel_sigma, kernel_cutoff, y_out, workspace_weights, workspace_values, mask_in)
  use f42_utils, only: loess_smooth_2d
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  integer, intent(in) :: n_total, n_target
  real(real64), intent(in) :: x_ref(n_total), y_ref(n_total), x_query(n_target)
  integer, intent(in) :: indices_used(n_total)
  real(real64), intent(in) :: kernel_sigma, kernel_cutoff
  real(real64), intent(out) :: y_out(n_target)
  real(real64), intent(inout) :: workspace_weights(n_total), workspace_values(n_total)
  logical, intent(in) :: mask_in(n_total)
  integer :: i
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

subroutine sort_character_c(array, perm, stack_left, stack_right, n, strlen) bind(C, name="sort_character_c")
  use iso_c_binding
  use f42_utils, only: sort_character
  implicit none
  integer(c_int), intent(in), value :: n, strlen
  character(len=strlen), intent(in) :: array(n)
  integer(c_int), intent(inout) :: perm(n)
  integer(c_int), intent(inout) :: stack_left(n)
  integer(c_int), intent(inout) :: stack_right(n)

  call sort_character(array, perm, stack_left, stack_right)
end subroutine sort_character_c


subroutine which_c(mask, n, idx_out, m_max, m_out) bind(C, name="which_c")
  use iso_c_binding
  use f42_utils, only: which
  implicit none
  integer(c_int), intent(in), value :: n, m_max
  integer(c_int), intent(in) :: mask(n)
  integer(c_int), intent(out) :: idx_out(m_max), m_out
  logical :: mask_f(n)
  integer :: i
  do i = 1, n
    mask_f(i) = (mask(i) /= 0)
  end do
  call which(mask_f, n, idx_out, m_max, m_out)
end subroutine which_c


subroutine loess_smooth_2d_c(n_total, n_target, x_ref, y_ref, indices_used, x_query, &
    kernel_sigma, kernel_cutoff, y_out, workspace_weights, workspace_values, mask_in) bind(C, name="loess_smooth_2d_c")
  use iso_c_binding
  use f42_utils, only: loess_smooth_2d
  implicit none
  integer(c_int), intent(in), value :: n_total, n_target
  real(c_double), intent(in) :: x_ref(n_total), y_ref(n_total), x_query(n_target)
  integer(c_int), intent(in) :: indices_used(n_total)
  real(c_double), intent(in), value :: kernel_sigma, kernel_cutoff
  real(c_double), intent(out) :: y_out(n_target)
  real(c_double), intent(inout) :: workspace_weights(n_total), workspace_values(n_total)
  integer(c_int), intent(in) :: mask_in(n_total)
  logical :: mask_f(n_total)
  integer :: i
  logical :: all_true


  do i = 1, size(mask_in)
    mask_f(i) = (mask_in(i) /= 0)
  end do

  call loess_smooth_2d(n_total, n_target, x_ref, y_ref, indices_used, x_query, &
    kernel_sigma, kernel_cutoff, y_out, workspace_weights, workspace_values, mask_f)

end subroutine loess_smooth_2d_c
