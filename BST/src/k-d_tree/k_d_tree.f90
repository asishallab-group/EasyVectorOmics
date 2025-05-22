module cartesian_kd_tree
    use tox_sorting
    implicit none
    private
    public :: build_kd_index

contains

    !> Build a k-d tree index using a stack-based, non-recursive approach.
    !! Arguments:
    !!   X         - real(kind=8), dimension(d, n): data points
    !!   d         - integer: number of dimensions
    !!   n         - integer: number of points
    !!   kd_ix     - integer, dimension(n): output index array (k-d tree order)
    !!   dim_order - integer, dimension(d): dimension order (by variance)
    !!   work      - integer, dimension(n): workspace array
    subroutine build_kd_index(X, d, n, kd_ix, dim_order, work)
        implicit none
        integer, intent(in) :: d, n
        real(kind=8), intent(in) :: X(d, n)
        integer, intent(in) :: dim_order(d)
        integer, intent(out) :: kd_ix(n)
        integer, intent(inout) :: work(n)
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
            mid = (l + r) / 2

            !> Partition kd_ix(l:r) by X(dim, kd_ix(:)), so that kd_ix(mid) is median in dim
            call partial_sort_by_dimension(X, d, kd_ix, l, r, dim, mid, work)

            !> Push right and left intervals onto stack
            if (mid + 1 < r) then
                stack_top = stack_top + 1
                stack(1, stack_top) = mid + 1
                stack(2, stack_top) = r
                stack(3, stack_top) = depth + 1
            end if
            if (l < mid - 1) then
                stack_top = stack_top + 1
                stack(1, stack_top) = l
                stack(2, stack_top) = mid - 1
                stack(3, stack_top) = depth + 1
            end if
        end do
    end subroutine build_kd_index

    !> Helper: sorts kd_ix(l:r) by X(dim, kd_ix(:)) using tox_sorting sort_array
    subroutine partial_sort_by_dimension(X, d, kd_ix, l, r, dim, mid, work)
        implicit none
        integer, intent(in) :: d, l, r, dim, mid
        real(kind=8), intent(in) :: X(d, *)
        integer, intent(inout) :: kd_ix(*)
        integer, intent(inout) :: work(*)
        integer :: n_sub, i
        ! The following arrays must be allocated by the caller and passed via 'work'
        real(8), pointer :: subarray(:)
        integer, pointer :: perm(:)
        integer, pointer :: stack_left(:), stack_right(:)

        n_sub = r - l + 1
        if (n_sub <= 1) return

        !> Assume subarray, perm, stack_left, stack_right are already associated with sufficient size
        !> Fill subarray with the values of X(dim, kd_ix(l:r))
        do i = 1, n_sub
            subarray(i) = X(dim, kd_ix(l + i - 1))
            perm(i) = i
        end do

        !> Sort subarray indirectly, perm will hold the sorted order
        call sort_array(subarray, perm, stack_left, stack_right)

        !> Reorder kd_ix(l:r) according to perm
        do i = 1, n_sub
            work(i) = kd_ix(l + perm(i) - 1)
        end do
        do i = 1, n_sub
            kd_ix(l + i - 1) = work(i)
        end do

        ! No allocation or deallocation here; handled by parent
    end subroutine partial_sort_by_dimension

end module cartesian_kd_tree