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

            !print *, "STACK POP: l=", l, " r=", r, " depth=", depth, " dim=", dim, " mid=", mid
            !print *, "kd_ix before partition:", kd_ix(l:r)

            !> Partition kd_ix(l:r) by X(dim, kd_ix(:)), so that kd_ix(mid) is median in dim
            call partial_sort_by_dimension(X, d, kd_ix, l, r, dim, mid, work, subarray, perm, stack_left, stack_right)

            !print *, "kd_ix after partition:", kd_ix(l:r)

            !> Push right and left intervals onto stack
            if (mid < r) then
                stack_top = stack_top + 1
                stack(1, stack_top) = mid + 1
                stack(2, stack_top) = r
                stack(3, stack_top) = depth + 1
                !print *, "STACK PUSH RIGHT: l=", mid+1, " r=", r, " depth=", depth+1
            end if
            if (l < mid) then
                stack_top = stack_top + 1
                stack(1, stack_top) = l
                stack(2, stack_top) = mid - 1
                stack(3, stack_top) = depth + 1
                !print *, "STACK PUSH LEFT: l=", l, " r=", mid-1, " depth=", depth+1
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

        !print *, "partial_sort_by_dimension: l=", l, " r=", r, " dim=", dim, " n_sub=", n_sub

        !> Fill subarray with the values of X(dim, kd_ix(l:r))
        do i = 1, n_sub
            subarray(i) = X(dim, kd_ix(l + i - 1))
            perm(i) = i
        end do

        !print *, "subarray before sort:", subarray(1:n_sub)
        !print *, "perm before sort:", perm(1:n_sub)
        call sort_array(subarray(1:n_sub), perm(1:n_sub), stack_left, stack_right)
        !print *, "subarray after sort:", subarray(perm(1:n_sub))
        !print *, "perm after sort:", perm(1:n_sub)

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