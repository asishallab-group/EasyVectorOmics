!> \file
!! \brief Flat-index-based Binary Search Tree (BST) utilities for 1D range queries.
!!
!! This module provides routines to build a BST index (via sorting), access sorted values,
!! and perform efficient range queries over a real-valued array.
module binary_search_tree
  use f42_utils, only: sort_array
  use iso_fortran_env, only: int32, real64
  use tox_errors, only: ERR_OK, ERR_INVALID_INPUT, ERR_EMPTY_INPUT, ERR_DIM_MISMATCH, ERR_SIZE_MISMATCH, is_ok, set_err_once, set_ok
  implicit none
  public :: build_bst_index, get_sorted_value, bst_range_query
contains
  
  !> \brief Build the BST index by sorting indices using values in x.
  subroutine build_bst_index(values, num_values, sorted_indices, left_stack, right_stack, ierr)
    real(real64), intent(in) :: values(num_values)      
    !! Input real array to be indexed
    integer(int32), intent(in) :: num_values           
    !! Number of elements in values array
    integer(int32), intent(out) :: sorted_indices(num_values)  
    !! Output permutation index
    integer(int32), intent(inout) :: left_stack(num_values)    
    !! Manual stack for left indices
    integer(int32), intent(inout) :: right_stack(num_values)   
    !! Manual stack for right indices
    integer(int32), intent(out) :: ierr                   
    !! Error code
    integer(int32) :: idx

    call set_ok(ierr)
    
    ! Input validation
    if (num_values < 0) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      return
    end if
    
    if (num_values == 0) then
      call set_err_once(ierr, ERR_EMPTY_INPUT)
      return
    end if
    
    if (size(values) < num_values .or. size(sorted_indices) < num_values .or. &
        size(left_stack) < num_values .or. size(right_stack) < num_values) then
      call set_err_once(ierr, ERR_DIM_MISMATCH)
      return
    end if

    do idx = 1, num_values
      sorted_indices(idx) = idx
    end do
    call sort_array(values, sorted_indices, left_stack, right_stack)
  end subroutine build_bst_index

  !> \brief Get the value at the sorted position.
  function get_sorted_value(values, sorted_indices, position, ierr) result(sorted_value)
    real(real64), intent(in) :: values(:)              
    !! Input real array
    integer(int32), intent(in) :: sorted_indices(:)    
    !! Permutation index array
    integer(int32), intent(in) :: position             
    !! Sorted position (1-based)
    integer(int32), intent(out) :: ierr                
    !! Error code
    real(real64) :: sorted_value

    call set_ok(ierr)
    
    ! Input validation
    if (position < 1 .or. position > size(sorted_indices)) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      return
    end if
    
    if (sorted_indices(position) < 1 .or. sorted_indices(position) > size(values)) then
      call set_err_once(ierr, ERR_DIM_MISMATCH)
      return
    end if
    
    sorted_value = values(sorted_indices(position))
  end function get_sorted_value

  !> \brief Perform a 1D range query over the sorted index.
  subroutine bst_range_query(values, sorted_indices, num_values, lower_bound, upper_bound, &
                            output_indices, num_matches, ierr)
    real(real64), intent(in) :: values(num_values)      
    !! Input real array
    integer(int32), intent(in) :: sorted_indices(num_values)  
    !! Permutation index array (sorted)
    integer(int32), intent(in) :: num_values           
    !! Number of elements
    real(real64), intent(in) :: lower_bound            
    !! Lower bound of range (inclusive)
    real(real64), intent(in) :: upper_bound            
    !! Upper bound of range (inclusive)
    integer(int32), intent(out) :: output_indices(num_values)  
    !! Output array of matching indices
    integer(int32), intent(out) :: num_matches         
    !! Number of matches found
    integer(int32), intent(out) :: ierr                
    !! Error code
    integer(int32) :: idx

    call set_ok(ierr)
    
    ! Input validation
    if (num_values < 0) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      return
    end if
    
    if (num_values == 0) then
      call set_err_once(ierr, ERR_EMPTY_INPUT)
      return
    end if
    
    if (lower_bound > upper_bound) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      return
    end if
    
    if (size(values) < num_values .or. size(sorted_indices) < num_values .or. &
        size(output_indices) < num_values) then
      call set_err_once(ierr, ERR_DIM_MISMATCH)
      return
    end if

    num_matches = 0
    do idx = 1, num_values
      if (values(sorted_indices(idx)) >= lower_bound .and. &
          values(sorted_indices(idx)) <= upper_bound) then
        num_matches = num_matches + 1
        output_indices(num_matches) = sorted_indices(idx)
      else if (values(sorted_indices(idx)) > upper_bound) then
        exit
      end if
    end do
  end subroutine bst_range_query

end module binary_search_tree

!> \brief Wrapper for getting range query usable by R
subroutine bst_range_query_r(values, sorted_indices, num_values, lower_bound, upper_bound, &
                            output_indices, num_matches, ierr)

  use binary_search_tree, only: bst_range_query
  use iso_fortran_env, only: int32, real64
  use tox_errors, only: ERR_OK

  implicit none
  integer(int32), intent(in) :: num_values           
  !! Number of elements
  real(real64), intent(in) :: values(num_values)      
  !! Input real array
  integer(int32), intent(in) :: sorted_indices(num_values)  
  !! Permutation index array
  real(real64), intent(in) :: lower_bound            
  !! Lower bound of range
  real(real64), intent(in) :: upper_bound            
  !! Upper bound of range
  integer(int32), intent(out) :: output_indices(num_values)  
  !! Output array of matching indices
  integer(int32), intent(out) :: num_matches         
  !! Number of matches found
  integer(int32), intent(out) :: ierr                
  !! Error code

  call bst_range_query(values, sorted_indices, num_values, lower_bound, upper_bound, &
                      output_indices, num_matches, ierr)
end subroutine bst_range_query_r

!> \brief Wrapper for building BST index usable by R
subroutine build_bst_index_r(values, num_values, sorted_indices, left_stack, right_stack, ierr)
  use binary_search_tree, only: build_bst_index
  use iso_fortran_env, only: int32, real64
  use tox_errors, only: ERR_OK
  implicit none
  integer(int32), intent(in) :: num_values           
  !! Number of elements
  real(real64), intent(in) :: values(num_values)     
  !! Input real array
  integer(int32), intent(out) :: sorted_indices(num_values)  
  !! Output permutation index
  integer(int32), intent(inout) :: left_stack(num_values)    
  !! Manual stack for left indices
  integer(int32), intent(inout) :: right_stack(num_values)   
  !! Manual stack for right indices
  integer(int32), intent(out) :: ierr                
  !! Error code

  call build_bst_index(values, num_values, sorted_indices, left_stack, right_stack, ierr)
end subroutine build_bst_index_r

!> \brief Wrapper using C for getting range query usable by python
subroutine bst_range_query_C(values, sorted_indices, num_values, lower_bound, upper_bound, &
                            output_indices, num_matches, ierr) bind(C, name='bst_range_query_C')
  use iso_c_binding, only: c_int, c_double, c_f_pointer, c_loc
  use binary_search_tree, only: bst_range_query
  use tox_errors, only: ERR_OK
  implicit none
  real(c_double), intent(in), target :: values(*)    
  !! Input real array (C-style)
  integer(c_int), intent(in), target :: sorted_indices(*)  
  !! Permutation index array (C-style)
  integer(c_int), value :: num_values               
  !! Number of elements
  real(c_double), value :: lower_bound              
  !! Lower bound of range
  real(c_double), value :: upper_bound              
  !! Upper bound of range
  integer(c_int), intent(out), target :: output_indices(*)  
  !! Output array (C-style)
  integer(c_int), intent(out) :: num_matches        
  !! Number of matches found
  integer(c_int), intent(out) :: ierr                
  !! Error code

  ! Zero-copy pointer association
  real(c_double), pointer :: values_shape(:)
  integer(c_int), pointer :: sorted_indices_shape(:)
  integer(c_int), pointer :: output_indices_shape(:)

  call c_f_pointer(c_loc(values(1)), values_shape, [num_values])
  call c_f_pointer(c_loc(sorted_indices(1)), sorted_indices_shape, [num_values])
  call c_f_pointer(c_loc(output_indices(1)), output_indices_shape, [num_values])

  call bst_range_query(values_shape, sorted_indices_shape, num_values, lower_bound, upper_bound, &
                      output_indices_shape, num_matches, ierr)
end subroutine bst_range_query_C

!> \brief Wrapper using C for building BST index usable by python
subroutine build_bst_index_C(values, num_values, sorted_indices, left_stack, right_stack, ierr) &
                            bind(C, name='build_bst_index_C')
  use iso_c_binding, only: c_int, c_double, c_f_pointer, c_loc
  use binary_search_tree
  use tox_errors, only: ERR_OK
  implicit none
  integer(c_int), value :: num_values               
  !! Number of elements
  real(c_double), intent(in), target :: values(num_values)    
  !! Input real array (C-style)
  integer(c_int), intent(out), target :: sorted_indices(num_values)  
  !! Output permutation index (C-style)
  integer(c_int), intent(inout), target :: left_stack(*)    
  !! Manual stack for left indices
  integer(c_int), intent(inout), target :: right_stack(*)   
  !! Manual stack for right indices
  integer(c_int), intent(out) :: ierr                
  !! Error code

  ! Zero-copy pointer association
  real(c_double), pointer :: values_shape(:)
  integer(c_int), pointer :: sorted_indices_shape(:)
  integer(c_int), pointer :: left_stack_shape(:)
  integer(c_int), pointer :: right_stack_shape(:)

  call c_f_pointer(c_loc(values(1)), values_shape, [num_values])
  call c_f_pointer(c_loc(sorted_indices(1)), sorted_indices_shape, [num_values])
  call c_f_pointer(c_loc(left_stack(1)), left_stack_shape, [num_values])
  call c_f_pointer(c_loc(right_stack(1)), right_stack_shape, [num_values])

  call build_bst_index(values_shape, num_values, sorted_indices_shape, left_stack_shape, right_stack_shape, ierr)
end subroutine build_bst_index_C