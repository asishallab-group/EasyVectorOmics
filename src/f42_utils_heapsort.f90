module f42_utils_heapsort
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use tox_errors, only: ERR_OK, ERR_INVALID_INPUT, ERR_EMPTY_INPUT, set_ok, set_err_once
  implicit none

  public :: sort_real, sort_integer, sort_character
  public :: sort_array

  interface sort_array
    module procedure sort_real, sort_integer, sort_character
  end interface sort_array

contains
  !> Sort a real array indirectly using heapsort.  
  !> The sorted array is returned via the perm array, which contains the indices of the original array in sorted order.
  pure subroutine sort_real(array, perm)
    real(real64), intent(in) ::array(:)
    integer(int32), intent(inout) ::perm(size(array))
    call heapsort_real(array, perm)
  end subroutine sort_real  
  !> Sort an integer array indirectly using heapsort.  
   !| Similar to `sort_real`, but for integer input.
   pure subroutine sort_integer(array, perm)
   integer(int32), intent(in) ::array(:)
   integer(int32), intent(inout) ::perm(size(array))
   call heapsort_integer(array, perm)
  end subroutine sort_integer  
  !> Sort a character array indirectly using heapsort.  
  !| Similar to `sort_real`, but for character input.     
  pure subroutine sort_character(array, perm)
    character(len=*), intent(in) ::array(:)
    integer(int32), intent(inout) ::perm(size(array))
    call heapsort_character(array, perm)
  end subroutine sort_character  

  !> Heapsort implementation for real arrays.
  pure subroutine heapsort_real(array, perm)
  implicit none
    real(real64), intent(in) ::array(:)
    integer(int32), intent(inout) ::perm(:)
    integer :: n, i, temp
      n = size(array)

      ! Initialize permutation vector
      do i = 1, n
        perm(i) = i
      end do

    ! Build max-heap
    do i = n / 2, 1, -1
      call heapify_real(array, perm, n, i)
    end do
    ! Heap sort
    do i = n, 2, -1
      temp = perm(1)
      perm(1) = perm(i)
      perm(i) = temp
      call heapify_real(array, perm, i - 1, 1)
    end do  


    contains
    !> Helper function to maintain the heap property.
     pure subroutine heapify(array, perm, heap_size, root)
      real(real64), intent(in) ::array(:)
      integer(int32), intent(inout) ::perm(:)
      integer, intent(in) ::heap_size, root
      integer :: largest, left, right, temp

      largest = i
      left = 2 * i
      right = 2 * i + 1
      if (left <= heap_size .and. array(perm(left)) > array(perm(largest))) then
        largest = left
      end if    
      if (right <= heap_size .and. array(perm(right)) > array(perm(largest))) then
        largest = right
      end if
`      if (largest /= root) then
        temp = perm(root)
        perm(root) = perm(largest)
        perm(largest) = temp
        call heapify(array, perm, heap_size, largest)
      end if
    end subroutine heapify
  end subroutine heapsort_real

  !> Heapsort implementation for integer arrays. 
  pure subroutine heapsort_integer(array, perm)
  implicit none
    integer(int32), intent(in) ::array(:)
    integer(int32), intent(inout) ::perm(:)
    integer :: n, i, temp
      n = size(array)   
      ! Initialize permutation vector
      do i = 1, n
        perm(i) = i
      end do  
    ! Build max-heap
    do i = n / 2, 1, -1
      call heapify_integer(array, perm, n, i)
    end do
    ! Heap sort
    do i = n, 2, -1
      temp = perm(1)
      perm(1) = perm(i)
      perm(i) = temp
      call heapify_integer(array, perm, i - 1, 1)
    end do
    contains
    !> Helper function to maintain the heap property.
     pure subroutine heapify_integer(array, perm, heap_size, root)
      integer(int32), intent(in) ::array(:)
      integer(int32), intent(inout) ::perm(:)     
      integer, intent(in) ::heap_size, root
      integer :: largest, left, right, temp
      largest = i
      left = 2 * i
      right = 2 * i + 1
      if (left <= heap_size .and. array(perm(left)) > array(perm(largest))) then
        largest = left
      end if    
      if (right <= heap_size .and. array(perm(right)) > array(perm(largest))) then
        largest = right   
      end if
      if (largest /= root) then
        temp = perm(root)
        perm(root) = perm(largest)
        perm(largest) = temp
        call heapify_integer(array, perm, heap_size, largest)
      end if
    end subroutine heapify_integer
  end subroutine heapsort_integer 

!> Heapsort implementation for character arrays. 
  pure subroutine heapsort_character(array, perm)
  implicit none
    character(len=*), intent(in) ::array(:)
    integer(int32), intent(inout) ::perm(:)
    integer :: n, i, temp
      n = size(array)   
      ! Initialize permutation vector
      do i = 1, n
        perm(i) = i
      end do  
    ! Build max-heap
    do i = n / 2, 1, -1
      call heapify_character(array, perm, n, i)
    end do
    ! Heap sort
    do i = n, 2, -1
      temp = perm(1)
      perm(1) = perm(i)
      perm(i) = temp
      call heapify_character(array, perm, i - 1, 1)
    end do
    contains
    !> Helper function to maintain the heap property.
     pure subroutine heapify_character(array, perm, heap_size, root)
      character(len=*), intent(in) ::array(:)
      integer(int32), intent(inout) ::perm(:)     
      integer, intent(in) ::heap_size, root
      integer :: largest, left, right, temp
      largest = i
      left = 2 * i
      right = 2 * i + 1
      if (left <= heap_size .and. array(perm(left)) > array(perm(largest))) then
        largest = left
      end if    
      if (right <= heap_size .and. array(perm(right)) > array(perm(largest))) then
        largest = right   
      end if
      if (largest /= root) then
        temp = perm(root)
        perm(root) = perm(largest)
        perm(largest) = temp
        call heapify_character(array, perm, heap_size, largest)
      end if
    end subroutine heapify_character
  end subroutine heapsort_character 

  !> Swap two integer values in-place.
    pure subroutine swap_int(a, b)
        integer(int32), intent(inout) ::a, b
        integer(int32) ::temp
        temp = a
        a = b
        b = temp
    end subroutine swap_int

      !> Finds the indices of the true values in a logical mask.
    pure subroutine which(mask, n, idx_out, m_max, m_out, ierr)
        logical, intent(in) ::mask(:)
        integer(int32), intent(in) ::n
        integer(int32), intent(out) ::idx_out(:)
        integer(int32), intent(in) ::m_max
        integer(int32), intent(out) ::m_out
        integer(int32), intent(out) ::ierr
        integer(int32) :: i, count
    
        ! Initialize output
        m_out = 0
        ierr = ERR_OK
        if (n <= 0 .or. size(mask) < n .or. m_max <= 0 .or. size(idx_out) < m_max) then
            ierr = ERR_INVALID_INPUT
            return
        end if
    
        count = 0
        do i = 1, n
            if (mask(i)) then
            count = count + 1
            if (count <= m_max) then
                idx_out(count) = i
            end if
            end if
        end do
    
        m_out = count
        if (m_out > m_max) then
            ierr = ERR_EMPTY_INPUT
        end if
    
        end subroutine which

        !> Performs LOESS smoothing on a set of data points.
  !| Smooths y_ref at x_query using reference points x_ref, y_ref, and kernel parameters.
  !| The user must pre-filter data and provide only valid indices in indices_used.
  pure subroutine loess_smooth_2d(n_total, n_target, x_ref, y_ref, indices_used, n_used, x_query, &
                          kernel_sigma, kernel_cutoff, y_out, ierr)
    !| Total number of reference points.
    integer(int32), intent(in) :: n_total
    !| Number of target points to smooth.
    integer(int32), intent(in) :: n_target
    !| Reference x-coordinates.
    real(real64), intent(in) :: x_ref(n_total)
    !| Reference y-coordinates (length n_total).
    real(real64), intent(in) :: y_ref(n_total)
    !| Indices of reference points used for smoothing (only valid indices).
    integer(int32), intent(in) :: indices_used(n_used)
    !| Number of indices actually used for smoothing.
    integer(int32), intent(in) :: n_used
    !| Target x-coordinates to smooth.
    real(real64), intent(in) :: x_query(n_target)
    !| Bandwidth parameter for the kernel.
    real(real64), intent(in) :: kernel_sigma
    !| Cutoff for the kernel.
    real(real64), intent(in) :: kernel_cutoff
    !| Output smoothed values (length n_target).
    real(real64), intent(out) :: y_out(n_target)
    !| Error code: 0=ok, 201=invalid input, 202=empty input
    integer(int32), intent(out) :: ierr

    integer(int32) :: q, i, idx
    real(real64) :: query_x, ref_x, delta, sum_weights, weight
    real(real64) :: min_dist
    integer(int32) :: min_idx
    logical :: exact_match_found, use_kernel

    ! Initialize error code
    call set_ok(ierr)

    ! Input validation
    if (n_total <= 0 .or. n_target <= 0) then
      call set_err_once(ierr, ERR_EMPTY_INPUT)
      y_out = 0.0_real64
      return
    end if
    
    if (n_used <= 0) then
      call set_err_once(ierr, ERR_EMPTY_INPUT)
      y_out = 0.0_real64
      return
    end if
    
    if (kernel_sigma < 0.0_real64 .or. kernel_cutoff < 0.0_real64) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      y_out = 0.0_real64
      return
    end if
    
    ! Validate array sizes
    if (size(x_ref) < n_total .or. size(y_ref) < n_total .or. &
        size(indices_used) < n_used .or. size(x_query) < n_target .or. &
        size(y_out) < n_target) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      y_out = 0.0_real64
      return
    end if
    
    ! Validate indices are within bounds
    do i = 1, n_used
      if (indices_used(i) < 1 .or. indices_used(i) > n_total) then
        call set_err_once(ierr, ERR_INVALID_INPUT)
        y_out = 0.0_real64
        return
      end if
    end do

    ! Check if we should use kernel smoothing
    use_kernel = (kernel_sigma > 0.0_real64)

    do q = 1, n_target
      query_x = x_query(q)
      sum_weights = 0.0_real64
      y_out(q) = 0.0_real64
      min_dist = huge(1.0_real64)
      min_idx = indices_used(1)
      exact_match_found = .false.

      ! Process all reference points
      do i = 1, n_used
        idx = indices_used(i)
        ref_x = x_ref(idx)
        delta = abs(query_x - ref_x)
        
        ! Check for exact match
        if (delta == 0.0_real64) then
          y_out(q) = y_ref(idx)
          exact_match_found = .true.
          exit
        end if
        
        ! Track closest point for potential fallback
        if (delta < min_dist) then
          min_dist = delta
          min_idx = idx
        end if
        
        ! Apply kernel smoothing if enabled and within cutoff
        if (use_kernel .and. delta <= kernel_cutoff * kernel_sigma) then
          weight = exp(-(delta / kernel_sigma)**2)
          sum_weights = sum_weights + weight
          y_out(q) = y_out(q) + weight * y_ref(idx)
        end if
      end do

      ! Finalize result if no exact match was found
      if (.not. exact_match_found) then
        if (sum_weights > 0.0_real64) then
          ! We have weighted average from kernel smoothing
          y_out(q) = y_out(q) / sum_weights
        else
          ! Fallback: use nearest neighbor
          y_out(q) = y_ref(min_idx)
        end if
      end if
    end do
  end subroutine loess_smooth_2d

end module f42_utils



! === R WRAPPERS ===

!> R wrapper for loess_smooth_2d.
!| Direct wrapper - user must pre-filter indices in R before calling.
subroutine loess_smooth_2d_r(n_total, n_target, x_ref, y_ref, indices_used, n_used, x_query, &
    kernel_sigma, kernel_cutoff, y_out, ierr)
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
  !| Indices of reference points used for smoothing (pre-filtered).
  integer(int32), intent(in) :: indices_used(n_used)
  !| Number of indices actually used for smoothing.
  integer(int32), intent(in) :: n_used
  !| Target x-coordinates to smooth.
  real(real64), intent(in) :: x_query(n_target)
  !| Bandwidth parameter for the kernel.
  real(real64), intent(in) :: kernel_sigma
  !| Cutoff for the kernel.
  real(real64), intent(in) :: kernel_cutoff
  !| Output smoothed values (length n_target).
  real(real64), intent(out) :: y_out(n_target)
  !| Error code: 0=ok, 201=invalid input, 202=empty input
  integer(int32), intent(out) :: ierr
  
  call loess_smooth_2d(n_total, n_target, x_ref, y_ref, indices_used, n_used, x_query, &
    kernel_sigma, kernel_cutoff, y_out, ierr)
end subroutine loess_smooth_2d_r

! === C WRAPPERS ===

!> C wrapper for which.
!| Converts integer mask to logical and calls which.
subroutine which_c(mask, n, idx_out, m_max, m_out, ierr) bind(C, name="which_c")
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
  !| Error code: 0=ok, 201=invalid input, 202=empty input
  integer(c_int), intent(out) :: ierr
  logical :: mask_f(n)
  integer(int32) :: i, ierr_f
  do i = 1, n
    mask_f(i) = (mask(i) /= 0)
  end do
  call which(mask_f, n, idx_out, m_max, m_out, ierr_f)
  ierr = ierr_f
end subroutine which_c

!> C wrapper for loess_smooth_2d.
!| Direct wrapper - user must pre-filter indices in C before calling.
subroutine loess_smooth_2d_c(n_total, n_target, x_ref, y_ref, indices_used, n_used, x_query, &
    kernel_sigma, kernel_cutoff, y_out, ierr) bind(C, name="loess_smooth_2d_c")
  use iso_c_binding, only : c_int, c_double
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
  !| Indices of reference points used for smoothing (pre-filtered).
  integer(c_int), intent(in) :: indices_used(n_used)
  !| Number of indices actually used for smoothing.
  integer(c_int), intent(in), value :: n_used
  !| Target x-coordinates to smooth.
  real(c_double), intent(in) :: x_query(n_target)
  !| Bandwidth parameter for the kernel.
  real(c_double), intent(in), value :: kernel_sigma
  !| Cutoff for the kernel.
  real(c_double), intent(in), value :: kernel_cutoff
  !| Output smoothed values (length n_target).
  real(c_double), intent(out) :: y_out(n_target)
  !| Error code: 0=ok, 201=invalid input, 202=empty input
  integer(c_int), intent(out) :: ierr

  integer(int32) :: ierr_f
  call loess_smooth_2d(n_total, n_target, x_ref, y_ref, indices_used, n_used, x_query, &
    kernel_sigma, kernel_cutoff, y_out, ierr_f)
  ierr = ierr_f

end subroutine loess_smooth_2d_c


end module f42_utils_heapsort



  
       



