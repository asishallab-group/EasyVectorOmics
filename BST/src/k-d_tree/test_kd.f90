program test_kd_tree
    use cartesian_kd_tree
    use tox_sorting
    implicit none

    integer, parameter :: d = 2, n = 6
    real(8) :: X(d, n)
    integer :: kd_ix(n)
    integer :: dim_order(d)
    integer :: work(n)
    real(8) :: var(d)
    ! Workspace for partial_sort_by_dimension
    real(8), target :: subarray(n)
    integer, target :: perm(n), stack_left(64), stack_right(64)
    integer :: i

    ! Example 2D points
    X(:,1) = [1.0d0, 2.0d0]
    X(:,2) = [2.0d0, 3.0d0]
    X(:,3) = [3.0d0, 1.0d0]
    X(:,4) = [4.0d0, 0.0d0]
    X(:,5) = [0.0d0, 4.0d0]
    X(:,6) = [5.0d0, 2.0d0]

    ! Compute dimension order by variance (descending)
    do i = 1, d
        var(i) = sum((X(i,:) - sum(X(i,:))/n)**2) / n
    end do
    if (var(1) >= var(2)) then
        dim_order = [1,2]
    else
        dim_order = [2,1]
    end if

    ! Associate workspace pointers for partial_sort_by_dimension
    call associate_workspace(subarray, perm, stack_left, stack_right)

    ! Build the k-d tree index
    call build_kd_index(X, d, n, kd_ix, dim_order, work)

    print *, 'KD Tree index order:'
    print *, kd_ix

    call test_kd_tree_3d()
    call test_kd_tree_4d()

contains

    !> Dummy workspace association for test (Fortran pointers not used in this test)
    subroutine associate_workspace(sa, p, sl, sr)
        real(8), target :: sa(:)
        integer, target :: p(:), sl(:), sr(:)
        ! No-op for this test, as pointers are not used in the current implementation
    end subroutine associate_workspace

    subroutine test_kd_tree_3d()
        integer, parameter :: d = 3, n = 8
        real(8) :: X(d, n)
        integer :: kd_ix(n)
        integer :: dim_order(d)
        integer :: work(n)
        real(8) :: var(d)
        integer :: i
        ! Workspace for partial_sort_by_dimension
        real(8), target :: subarray(n)
        integer, target :: perm(n), stack_left(64), stack_right(64)

        ! Example 3D points
        X(:,1) = [1.0d0, 2.0d0, 3.0d0]
        X(:,2) = [2.0d0, 3.0d0, 1.0d0]
        X(:,3) = [3.0d0, 1.0d0, 2.0d0]
        X(:,4) = [4.0d0, 0.0d0, 4.0d0]
        X(:,5) = [0.0d0, 4.0d0, 0.0d0]
        X(:,6) = [5.0d0, 2.0d0, 5.0d0]
        X(:,7) = [1.5d0, 2.5d0, 3.5d0]
        X(:,8) = [2.5d0, 3.5d0, 1.5d0]

        ! Compute dimension order by variance (descending)
        do i = 1, d
            var(i) = sum((X(i,:) - sum(X(i,:))/n)**2) / n
        end do
        call sort_dim_order(var, dim_order, d)

        call build_kd_index(X, d, n, kd_ix, dim_order, work)
        print *, '3D KD Tree index order:'
        print *, kd_ix
    end subroutine test_kd_tree_3d

    subroutine test_kd_tree_4d()
        integer, parameter :: d = 4, n = 10
        real(8) :: X(d, n)
        integer :: kd_ix(n)
        integer :: dim_order(d)
        integer :: work(n)
        real(8) :: var(d)
        integer :: i
        ! Workspace for partial_sort_by_dimension
        real(8), target :: subarray(n)
        integer, target :: perm(n), stack_left(64), stack_right(64)

        ! Example 4D points
        X(:,1) = [1.0d0, 2.0d0, 3.0d0, 4.0d0]
        X(:,2) = [2.0d0, 3.0d0, 1.0d0, 0.0d0]
        X(:,3) = [3.0d0, 1.0d0, 2.0d0, 1.0d0]
        X(:,4) = [4.0d0, 0.0d0, 4.0d0, 2.0d0]
        X(:,5) = [0.0d0, 4.0d0, 0.0d0, 3.0d0]
        X(:,6) = [5.0d0, 2.0d0, 5.0d0, 4.0d0]
        X(:,7) = [1.5d0, 2.5d0, 3.5d0, 0.5d0]
        X(:,8) = [2.5d0, 3.5d0, 1.5d0, 1.5d0]
        X(:,9) = [3.5d0, 0.5d0, 2.5d0, 2.5d0]
        X(:,10)= [4.5d0, 1.5d0, 4.5d0, 3.5d0]

        ! Compute dimension order by variance (descending)
        do i = 1, d
            var(i) = sum((X(i,:) - sum(X(i,:))/n)**2) / n
        end do
        call sort_dim_order(var, dim_order, d)

        call build_kd_index(X, d, n, kd_ix, dim_order, work)
        print *, '4D KD Tree index order:'
        print *, kd_ix
    end subroutine test_kd_tree_4d

    !> Sorts dimension indices by descending variance using tox_sorting
    subroutine sort_dim_order(var, dim_order, d)
        real(8), intent(in) :: var(d)
        integer, intent(out) :: dim_order(d)
        integer, intent(in) :: d
        integer :: perm(d)
        integer :: stack_left(64), stack_right(64)
        integer :: i

        do i = 1, d
            perm(i) = i
        end do
        call sort_array(-var, perm, stack_left, stack_right) ! Sort descending
        do i = 1, d
            dim_order(i) = perm(i)
        end do
    end subroutine sort_dim_order

end program test_kd_tree
