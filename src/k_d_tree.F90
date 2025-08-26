module kd_tree
    use f42_utils, only: sort_array
    use iso_fortran_env, only: int32, real64
    use tox_errors, only: ERR_OK, ERR_INVALID_INPUT, ERR_EMPTY_INPUT, ERR_DIM_MISMATCH, ERR_SIZE_MISMATCH, set_ok, set_err_once, is_ok
    implicit none
    private
    public :: build_kd_index, build_spherical_kd, get_kd_point

contains

    !> \brief Build a k-d tree index using a stack-based, non-recursive approach.
    subroutine build_kd_index(points, num_dimensions, num_points, kd_indices, dimension_order, &
                            workspace, value_buffer, permutation, left_stack, right_stack, ierr)
        integer(int32), intent(in) :: num_dimensions      
        !! Number of dimensions
        integer(int32), intent(in) :: num_points          
        !! Number of points
        real(real64), intent(in) :: points(num_dimensions, num_points)  
        !! Data points
        integer(int32), intent(in) :: dimension_order(num_dimensions)   
        !! Dimension order (by variance)
        integer(int32), intent(out) :: kd_indices(num_points)           
        !! Output index array (k-d tree order)
        integer(int32), intent(inout) :: workspace(num_points)          
        !! Workspace array
        real(real64), intent(inout) :: value_buffer(num_points)         
        !! Workspace for sorting
        integer(int32), intent(inout) :: permutation(num_points)        
        !! Workspace for sorting
        integer(int32), intent(inout) :: left_stack(:)                  
        !! Workspace for sorting
        integer(int32), intent(inout) :: right_stack(:)                 
        !! Workspace for sorting
        integer(int32), intent(out) :: ierr                             
        !! Error code
        
        integer(int32), parameter :: max_depth = 64
        integer(int32) :: recursion_stack(3, max_depth)  !! Stack for l, r, depth
        integer(int32) :: stack_top
        integer(int32) :: left_idx, right_idx, mid_idx, current_dim, current_depth
        integer(int32) :: i

        call set_ok(ierr)
        
        ! Input validation
        if (num_dimensions <= 0 .or. num_points < 0) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if
        
        if (num_points == 0) then
            ! Empty array is a valid case, just return success
            call set_ok(ierr)
            return
        end if
        
        if (size(points, 1) < num_dimensions .or. size(points, 2) < num_points .or. &
            size(kd_indices) < num_points .or. size(workspace) < num_points .or. &
            size(value_buffer) < num_points .or. size(permutation) < num_points .or. &
            size(left_stack) < num_points .or. size(right_stack) < num_points) then
            call set_err_once(ierr, ERR_DIM_MISMATCH)
            return
        end if
        
        if (any(dimension_order < 1) .or. any(dimension_order > num_dimensions)) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if

        !! Initialize kd_indices to 1:num_points (original indices)
        do i = 1, num_points
            kd_indices(i) = i
        end do

        stack_top = 1
        recursion_stack(1, 1) = 1
        recursion_stack(2, 1) = num_points
        recursion_stack(3, 1) = 0

        do while (stack_top > 0)
            left_idx = recursion_stack(1, stack_top)
            right_idx = recursion_stack(2, stack_top)
            current_depth = recursion_stack(3, stack_top)
            stack_top = stack_top - 1

            if (right_idx <= left_idx) cycle

            !! Choose split dimension by cycling through dimension_order
            current_dim = dimension_order(mod(current_depth, num_dimensions) + 1)

            !! Find median index
            mid_idx = left_idx + (right_idx - left_idx) / 2

            !! Partition kd_indices(left_idx:right_idx) by points(current_dim, kd_indices(:))
            call partial_sort_by_dimension(points, num_dimensions, kd_indices, left_idx, right_idx, &
                                        current_dim, mid_idx, workspace, value_buffer, permutation, &
                                        left_stack, right_stack, ierr)
            if (.not. is_ok(ierr)) return

            !! Push right and left intervals onto stack
            if (mid_idx < right_idx) then
                stack_top = stack_top + 1
                recursion_stack(1, stack_top) = mid_idx + 1
                recursion_stack(2, stack_top) = right_idx
                recursion_stack(3, stack_top) = current_depth + 1
            end if
            if (left_idx < mid_idx) then
                stack_top = stack_top + 1
                recursion_stack(1, stack_top) = left_idx
                recursion_stack(2, stack_top) = mid_idx - 1
                recursion_stack(3, stack_top) = current_depth + 1
            end if
        end do
    end subroutine build_kd_index

    !> \brief Helper: sorts kd_indices(left_idx:right_idx) by points(dimension, kd_indices(:))
    subroutine partial_sort_by_dimension(points, num_dimensions, kd_indices, left_idx, right_idx, &
                                        dim, mid_idx, workspace, value_buffer, permutation, &
                                        left_stack, right_stack, ierr)
        use f42_utils, only: sort_array
        integer(int32), intent(in) :: num_dimensions      
        !! Number of dimensions
        integer(int32), intent(in) :: left_idx            
        !! Left index of subarray
        integer(int32), intent(in) :: right_idx           
        !! Right index of subarray
        integer(int32), intent(in) :: dim         
        !! Dimension to sort by
        integer(int32), intent(in) :: mid_idx             
        !! Target median index
        real(real64), intent(in) :: points(num_dimensions, *)  
        !! Input points array
        integer(int32), intent(inout) :: kd_indices(:)         
        !! Index array to modify
        integer(int32), intent(inout) :: workspace(:)          
        !! Workspace array
        real(real64), intent(inout) :: value_buffer(:)         
        !! Buffer for dimension values
        integer(int32), intent(inout) :: permutation(:)        
        !! Permutation array
        integer(int32), intent(inout) :: left_stack(:)         
        !! Stack for left indices
        integer(int32), intent(inout) :: right_stack(:)        
        !! Stack for right indices
        integer(int32), intent(out) :: ierr                    
        !! Error code
        
        integer(int32) :: subarray_size, i

        call set_ok(ierr)
        
        ! Input validation
        if (left_idx < 1 .or. right_idx > size(kd_indices) .or. left_idx > right_idx) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if
        
        if (dim < 1 .or. dim > num_dimensions) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if
        
        if (mid_idx < left_idx .or. mid_idx > right_idx) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if

        subarray_size = right_idx - left_idx + 1
        if (subarray_size <= 1) return
        
        !! Fill value_buffer with the values of points(dimension, kd_indices(left_idx:right_idx))
        do i = 1, subarray_size
            value_buffer(i) = points(dim, kd_indices(left_idx + i - 1))
            permutation(i) = i
        end do

        call sort_array(value_buffer(1:subarray_size), permutation(1:subarray_size), left_stack, right_stack)

        !! Reorder kd_indices(left_idx:right_idx) according to permutation
        do i = 1, subarray_size
            if (permutation(i) < 1 .or. permutation(i) > subarray_size) then
                ierr = ERR_INVALID_INPUT
                return
            end if
            workspace(i) = kd_indices(left_idx + permutation(i) - 1)
        end do
        do i = 1, subarray_size
            kd_indices(left_idx + i - 1) = workspace(i)
        end do
    end subroutine partial_sort_by_dimension

    !> \brief Build spherical k-d tree index
    subroutine build_spherical_kd(vectors, num_dimensions, num_vectors, sphere_indices, &
                                dimension_order, workspace, value_buffer, permutation, &
                                left_stack, right_stack, ierr)

        integer(int32), intent(in) :: num_dimensions      
        !! Number of dimensions
        integer(int32), intent(in) :: num_vectors         
        !! Number of vectors
        real(real64), intent(in) :: vectors(num_dimensions, num_vectors)  
        !! Input unit vectors
        integer(int32), intent(out) :: sphere_indices(num_vectors)  
        !! Output index array
        integer(int32), intent(inout) :: dimension_order(num_dimensions)  
        !! Dimension order
        integer(int32), intent(inout) :: workspace(num_vectors)     
        !! Workspace array
        real(real64), intent(inout) :: value_buffer(num_vectors)    
        !! Value buffer
        integer(int32), intent(inout) :: permutation(num_vectors)   
        !! Permutation array
        integer(int32), intent(inout) :: left_stack(:)              
        !! Left stack
        integer(int32), intent(inout) :: right_stack(:)             
        !! Right stack
        integer(int32), intent(out) :: ierr                         
        !! Error code

        call build_kd_index(vectors, num_dimensions, num_vectors, sphere_indices, dimension_order, &
                          workspace, value_buffer, permutation, left_stack, right_stack, ierr)
    end subroutine build_spherical_kd

    !> \brief Get value from sorted array
    function get_value_sorted(values, indices, position, ierr) result(sorted_value)
        real(real64), intent(in) :: values(:)            
        !! Input array
        integer(int32), intent(in) :: indices(:)         
        !! Index array
        integer(int32), intent(in) :: position           
        !! Position in sorted array
        integer(int32), intent(out) :: ierr              
        !! Error code
        real(real64) :: sorted_value

        call set_ok(ierr)
        
        ! Input validation
        if (position < 1 .or. position > size(indices)) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if
        
        if (size(indices) == 0) then
            call set_err_once(ierr, ERR_EMPTY_INPUT)
            return
        end if
        
        if (indices(position) < 1 .or. indices(position) > size(values)) then
            call set_err_once(ierr, ERR_DIM_MISMATCH)
            return
        end if
        
        sorted_value = values(indices(position))
    end function get_value_sorted

    !> \brief Get point from KD index
    subroutine get_kd_point(points, kd_indices, position, point_values, ierr)
        real(real64), intent(in) :: points(:, :)         
        !! Input points
        integer(int32), intent(in) :: kd_indices(:)      
        !! KD index array
        integer(int32), intent(in) :: position           
        !! Position in index
        real(real64), intent(out) :: point_values(:)     
        !! Output point values
        integer(int32), intent(out) :: ierr              
        !! Error code

        call set_ok(ierr)
        
        ! Input validation
        if (position < 1 .or. position > size(kd_indices)) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if
        
        if (kd_indices(position) < 1 .or. kd_indices(position) > size(points, 2)) then
            call set_err_once(ierr, ERR_DIM_MISMATCH)
            return
        end if
        
        if (size(point_values) < size(points, 1)) then
            call set_err_once(ierr, ERR_DIM_MISMATCH)
            return
        end if
        
        point_values = points(:, kd_indices(position))
    end subroutine get_kd_point

end module kd_tree

!> \brief R interface for building KD index
subroutine build_kd_index_r(points, num_dimensions, num_points, kd_indices, dimension_order, &
                          workspace, value_buffer, permutation, left_stack, right_stack, ierr)
    use kd_tree, only: build_kd_index
    use iso_fortran_env, only: int32, real64
    implicit none
    integer(int32), intent(in) :: num_dimensions      
    !! Number of dimensions
    integer(int32), intent(in) :: num_points          
    !! Number of points
    real(real64), intent(in) :: points(num_dimensions, num_points)  
    !! Input points
    integer(int32), intent(in) :: dimension_order(num_dimensions)   
    !! Dimension order
    integer(int32), intent(out) :: kd_indices(num_points)           
    !! Output indices
    integer(int32), intent(inout) :: workspace(num_points)          
    !! Workspace
    real(real64), intent(inout) :: value_buffer(num_points)         
    !! Value buffer
    integer(int32), intent(inout) :: permutation(num_points)        
    !! Permutation array
    integer(int32), intent(inout) :: left_stack(num_points)         
    !! Left stack
    integer(int32), intent(inout) :: right_stack(num_points)        
    !! Right stack
    integer(int32), intent(out) :: ierr                
    !! Error code

    call build_kd_index(points, num_dimensions, num_points, kd_indices, dimension_order, &
                      workspace, value_buffer, permutation, left_stack, right_stack, ierr)
end subroutine build_kd_index_r

!> \brief R interface for building spherical KD index
subroutine build_spherical_kd_r(vectors, num_dimensions, num_vectors, sphere_indices, &
                              dimension_order, workspace, value_buffer, permutation, &
                              left_stack, right_stack, ierr)
    use kd_tree, only: build_spherical_kd
    use iso_fortran_env, only: int32, real64
    implicit none
    integer(int32), intent(in) :: num_dimensions      
    !! Number of dimensions
    integer(int32), intent(in) :: num_vectors         
    !! Number of vectors
    real(real64), intent(in) :: vectors(num_dimensions, num_vectors)  
    !! Input vectors
    integer(int32), intent(out) :: sphere_indices(num_vectors)  
    !! Output indices
    integer(int32), intent(inout) :: dimension_order(num_dimensions)  
    !! Dimension order
    integer(int32), intent(inout) :: workspace(num_vectors)     
    !! Workspace
    real(real64), intent(inout) :: value_buffer(num_vectors)    
    !! Value buffer
    integer(int32), intent(inout) :: permutation(num_vectors)   
    !! Permutation array
    integer(int32), intent(inout) :: left_stack(num_vectors)    
    !! Left stack
    integer(int32), intent(inout) :: right_stack(num_vectors)   
    !! Right stack
    integer(int32), intent(out) :: ierr                
    !! Error code

    call build_spherical_kd(vectors, num_dimensions, num_vectors, sphere_indices, dimension_order, &
                          workspace, value_buffer, permutation, left_stack, right_stack, ierr)
end subroutine build_spherical_kd_r

!> \brief C interface for building KD index
subroutine build_kd_index_C(points_flat, num_dimensions, num_points, kd_indices, dimension_order, &
                          workspace, value_buffer, permutation, left_stack, right_stack, ierr) &
                          bind(C, name="build_kd_index_C")
    use iso_c_binding, only: c_int, c_double, c_f_pointer, c_loc
    use kd_tree, only: build_kd_index
    implicit none
    integer(c_int), value :: num_dimensions          
    !! Number of dimensions
    integer(c_int), value :: num_points              
    !! Number of points
    real(c_double), intent(in), target :: points_flat(num_points)  
    !! Flattened points array
    integer(c_int), intent(in), target :: dimension_order(num_dimensions)  
    !! Dimension order
    integer(c_int), intent(out), target :: kd_indices(num_points)      
    !! Output indices
    integer(c_int), intent(inout), target :: workspace(*)     
    !! Workspace
    real(c_double), intent(inout), target :: value_buffer(*)  
    !! Value buffer
    integer(c_int), intent(inout), target :: permutation(*)   
    !! Permutation array
    integer(c_int), intent(inout), target :: left_stack(*)    
    !! Left stack
    integer(c_int), intent(inout), target :: right_stack(*)   
    !! Right stack
    integer(c_int), intent(out) :: ierr                
    !! Error code

    real(c_double), pointer :: points(:, :)
    integer(c_int), pointer :: dimension_order_p(:)
    integer(c_int), pointer :: kd_indices_p(:)
    integer(c_int), pointer :: workspace_p(:)
    real(c_double), pointer :: value_buffer_p(:)
    integer(c_int), pointer :: permutation_p(:)
    integer(c_int), pointer :: left_stack_p(:), right_stack_p(:)

    call c_f_pointer(c_loc(points_flat(1)), points, [num_dimensions, num_points])
    call c_f_pointer(c_loc(dimension_order(1)), dimension_order_p, [num_dimensions])
    call c_f_pointer(c_loc(kd_indices(1)), kd_indices_p, [num_points])
    call c_f_pointer(c_loc(workspace(1)), workspace_p, [num_points])
    call c_f_pointer(c_loc(value_buffer(1)), value_buffer_p, [num_points])
    call c_f_pointer(c_loc(permutation(1)), permutation_p, [num_points])
    call c_f_pointer(c_loc(left_stack(1)), left_stack_p, [num_points])
    call c_f_pointer(c_loc(right_stack(1)), right_stack_p, [num_points])

    call build_kd_index(points, num_dimensions, num_points, kd_indices_p, dimension_order_p, &
                      workspace_p, value_buffer_p, permutation_p, left_stack_p, right_stack_p, ierr)
end subroutine build_kd_index_C