module kd_tree
    use tox_sorting
    implicit none
    private
    public :: build_kd_index, build_spherical_kd 

contains

    !> Build a k-d tree index using a stack-based, non-recursive approach.
    !! Arguments:
    !!   X         - real(kind=8), dimension(d, n): data points
    !!   d         - integer: number of dimensions
    !!   n         - integer: number of points
    !!   kd_ix     - integer, dimension(n): output index array (k-d tree order)
    !!   dim_order - integer, dimension(d): dimension order (by variance)
    !!   work      - integer, dimension(n): workspace array
    !!   subarray  - real(8), dimension(n): workspace for sorting
    !!   perm      - integer, dimension(n): workspace for sorting
    !!   stack_left, stack_right - integer, dimension(max_depth): workspace for sorting
    subroutine build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, stack_left, stack_right)
        implicit none
        integer, intent(in) :: d, n
        real(kind=8), intent(in) :: X(d, n)
        integer, intent(in) :: dim_order(d)
        integer, intent(out) :: kd_ix(n)
        integer, intent(inout) :: work(n)
        real(8), intent(inout) :: subarray(n)
        integer, intent(inout) :: perm(n)
        integer, intent(inout) :: stack_left(:), stack_right(:)
        integer, parameter :: max_depth = 64
        integer :: stack(3, max_depth) !! l, r, depth
        integer :: stack_top
        integer :: l, r, mid, dim, depth
        integer :: i

        !> Initialize kd_ix to 1:n (original indices)
        do i = 1, n
            kd_ix(i) = i
        end do

        stack_top = 1
        stack(1, 1) = 1
        stack(2, 1) = n
        stack(3, 1) = 0

        do while (stack_top > 0)
            l = stack(1, stack_top)
            r = stack(2, stack_top)
            depth = stack(3, stack_top)
            stack_top = stack_top - 1

            if (r <= l) cycle

            !> Choose split dimension by cycling through dim_order
            dim = dim_order(mod(depth, d) + 1)

            !> Find median index
            mid = l + (r - l) / 2

            !> Partition kd_ix(l:r) by X(dim, kd_ix(:)), so that kd_ix(mid) is median in dim
            call partial_sort_by_dimension(X, d, kd_ix, l, r, dim, mid, work, subarray, perm, stack_left, stack_right)

            !print *, "kd_ix after partition:", kd_ix(l:r)

            !> Push right and left intervals onto stack
            if (mid < r) then
                stack_top = stack_top + 1
                stack(1, stack_top) = mid + 1
                stack(2, stack_top) = r
                stack(3, stack_top) = depth + 1
            end if
            if (l < mid) then
                stack_top = stack_top + 1
                stack(1, stack_top) = l
                stack(2, stack_top) = mid - 1
                stack(3, stack_top) = depth + 1
            end if
        end do
    end subroutine build_kd_index

    !> Helper: sorts kd_ix(l:r) by X(dim, kd_ix(:)) using tox_sorting sort_array
    subroutine partial_sort_by_dimension(X, d, kd_ix, l, r, dim, mid, work, subarray, perm, stack_left, stack_right)
        implicit none
        integer, intent(in) :: d, l, r, dim, mid
        real(kind=8), intent(in) :: X(d, *)
        integer, intent(inout) :: kd_ix(:)
        integer, intent(inout) :: work(:)
        real(8), intent(inout) :: subarray(:)
        integer, intent(inout) :: perm(:)
        integer, intent(inout) :: stack_left(:), stack_right(:)
        integer :: n_sub, i

        n_sub = r - l + 1
        if (n_sub <= 1) return
        !> Fill subarray with the values of X(dim, kd_ix(l:r))
        do i = 1, n_sub
            subarray(i) = X(dim, kd_ix(l + i - 1))
            perm(i) = i
        end do

        call sort_array(subarray(1:n_sub), perm(1:n_sub), stack_left, stack_right)

        !> Reorder kd_ix(l:r) according to perm
        do i = 1, n_sub
            if (perm(i) < 1 .or. perm(i) > n_sub) then
                print *, "ERROR: perm(", i, ") out of bounds: ", perm(i)
                stop 1
            end if
            work(i) = kd_ix(l + perm(i) - 1)
        end do
        do i = 1, n_sub
            kd_ix(l + i - 1) = work(i)
        end do

        !print *, "kd_ix after reorder:", kd_ix(l:r)

        ! No allocation or deallocation here; handled by parent
    end subroutine partial_sort_by_dimension

    subroutine build_spherical_kd(V, d, n, sphere_ix, dim_order, work, subarray, perm, stack_left, stack_right)
        ! V: real(kind=8), dimension(d, n), input, each column is unit length
        ! d: integer, input, number of dimensions
        ! n: integer, input, number of vectors
        ! sphere_ix: integer, dimension(n), output, index array for spherical k-d tree
        ! dim_order: integer, dimension(d), input/output, order of dimensions for splitting

        real(kind=8), intent(in) :: V(:,:)
        integer, intent(in) :: d, n
        integer, intent(out) :: sphere_ix(:)
        integer, intent(inout) :: dim_order(:)
        integer, intent(inout) :: work(n)
        real(8), intent(inout) :: subarray(n)
        integer, intent(inout) :: perm(n)
        integer, intent(inout) :: stack_left(:), stack_right(:)

        call build_kd_index(V, d, n, sphere_ix, dim_order, work, subarray , perm, stack_left, stack_right)
    end subroutine build_spherical_kd

    pure function get_value_sorted(x, ix, i) result(val)
        real(8), intent(in) :: x(:)
        integer, intent(in) :: ix(:)
        integer, intent(in) :: i
        real(8) :: val
        val = x(ix(i))
    end function get_value_sorted

    !> Get point from KD index
    subroutine get_kd_point(X, kd_ix, i, point)
        real(8), intent(in) :: X(:, :)
        integer, intent(in) :: kd_ix(:)
        integer, intent(in) :: i
        real(8), intent(out) :: point(:)
        point = X(:, kd_ix(i))
    end subroutine get_kd_point

end module kd_tree

subroutine build_kd_index_r(X, d, n, kd_ix, dim_order, work, subarray, perm, stack_left, stack_right)
    use kd_tree
    implicit none
    integer, intent(in) :: d, n
    real(kind=8), intent(in) :: X(d, n)
    integer, intent(in) :: dim_order(d)
    integer, intent(out) :: kd_ix(n)
    integer, intent(inout) :: work(n)
    real(8), intent(inout) :: subarray(n)
    integer, intent(inout) :: perm(n)
    integer, intent(inout) :: stack_left(n), stack_right(n)

    ! Call the build_kd_index from the kd_tree module
    call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, stack_left, stack_right)
end subroutine build_kd_index_r

subroutine build_spherical_kd_r(V, d, n, sphere_ix, dim_order, work, subarray, perm, stack_left, stack_right)
    use kd_tree
    implicit none
    real(kind=8), intent(in) :: V(d,n)
    integer, intent(in) :: d, n
    integer, intent(out) :: sphere_ix(n)
    integer, intent(inout) :: dim_order(n)
    integer, intent(inout) :: work(n)
    real(8), intent(inout) :: subarray(n)
    integer, intent(inout) :: perm(n)
    integer, intent(inout) :: stack_left(n), stack_right(n)

    ! Call the build_spherical_kd from the kd_tree module
    call build_spherical_kd(V, d, n, sphere_ix, dim_order, work, subarray, perm, stack_left, stack_right)
end subroutine build_spherical_kd_r

!> C interface: Build k-d tree index (zero-copy, short lines, no duplicate declarations)
subroutine build_kd_index_C(X_flat, d, n, kd_ix, dim_order, work, subarray, perm, stack_left, stack_right) &
    bind(C, name="build_kd_index_C")
    use iso_c_binding
    use kd_tree
    implicit none
    integer(c_int), value :: d, n
    real(c_double), intent(in), target :: X_flat(*)
    integer(c_int), intent(in), target :: dim_order(*)
    integer(c_int), intent(out), target :: kd_ix(*)
    integer(c_int), intent(inout), target :: work(*)
    real(c_double), intent(inout), target :: subarray(*)
    integer(c_int), intent(inout), target :: perm(*)
    integer(c_int), intent(inout), target :: stack_left(*), stack_right(*)

    real(c_double), pointer :: X(:,:)
    integer(c_int), pointer :: dim_order_p(:)
    integer(c_int), pointer :: kd_ix_p(:)
    integer(c_int), pointer :: work_p(:)
    real(c_double), pointer :: subarray_p(:)
    integer(c_int), pointer :: perm_p(:)
    integer(c_int), pointer :: stack_left_p(:), stack_right_p(:)

    call c_f_pointer(c_loc(X_flat(1)), X, [d, n])
    call c_f_pointer(c_loc(dim_order(1)), dim_order_p, [d])
    call c_f_pointer(c_loc(kd_ix(1)), kd_ix_p, [n])
    call c_f_pointer(c_loc(work(1)), work_p, [n])
    call c_f_pointer(c_loc(subarray(1)), subarray_p, [n])
    call c_f_pointer(c_loc(perm(1)), perm_p, [n])
    call c_f_pointer(c_loc(stack_left(1)), stack_left_p, [n])
    call c_f_pointer(c_loc(stack_right(1)), stack_right_p, [n])

    call build_kd_index(X, d, n, kd_ix_p, dim_order_p, work_p, subarray_p, perm_p, stack_left_p, stack_right_p)
end subroutine build_kd_index_C

!> C interface: Build spherical k-d tree (zero-copy, short lines, no duplicate declarations)
subroutine build_spherical_kd_C(V_flat, d, n, sphere_ix, dim_order, work, subarray, perm, stack_left, stack_right) &
    bind(C, name="build_spherical_kd_C")
    use iso_c_binding
    use kd_tree
    implicit none
    integer(c_int), value :: d, n
    real(c_double), intent(in), target :: V_flat(*)
    integer(c_int), intent(out), target :: sphere_ix(*)
    integer(c_int), intent(inout), target :: dim_order(*), work(*), perm(*)
    integer(c_int), intent(inout), target :: stack_left(*), stack_right(*)
    real(c_double), intent(inout), target :: subarray(*)

    real(c_double), pointer :: V(:,:)
    integer(c_int), pointer :: sphere_ix_p(:)
    integer(c_int), pointer :: dim_order_p(:)
    integer(c_int), pointer :: work_p(:)
    real(c_double), pointer :: subarray_p(:)
    integer(c_int), pointer :: perm_p(:)
    integer(c_int), pointer :: stack_left_p(:), stack_right_p(:)

    call c_f_pointer(c_loc(V_flat(1)), V, [d, n])
    call c_f_pointer(c_loc(sphere_ix(1)), sphere_ix_p, [n])
    call c_f_pointer(c_loc(dim_order(1)), dim_order_p, [n])
    call c_f_pointer(c_loc(work(1)), work_p, [n])
    call c_f_pointer(c_loc(subarray(1)), subarray_p, [n])
    call c_f_pointer(c_loc(perm(1)), perm_p, [n])
    call c_f_pointer(c_loc(stack_left(1)), stack_left_p, [n])
    call c_f_pointer(c_loc(stack_right(1)), stack_right_p, [n])

    call build_spherical_kd(V, d, n, sphere_ix_p, dim_order_p, work_p, subarray_p, perm_p, stack_left_p, stack_right_p)
end subroutine build_spherical_kd_C
