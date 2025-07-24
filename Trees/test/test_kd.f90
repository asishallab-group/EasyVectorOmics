program test_kd_tree
    use kd_tree
    use tox_sorting
    implicit none

    integer, parameter :: max_d = 10, max_n = 100
    integer :: d, n
    real(8), allocatable :: X(:,:), V(:,:)
    integer, allocatable :: kd_ix(:), sphere_ix(:), dim_order(:), work(:)
    real(8), allocatable :: var(:)
    real(8), allocatable :: subarray(:)
    integer, allocatable :: perm(:), stack_left(:), stack_right(:)
    integer :: i, j, fail_count
    real(8), allocatable :: val(:)

    ! --- Test 1: 2D, n=6, cartesian ---
    d = 2; n = 6
    if (allocated(X)) deallocate(X)
    if (allocated(kd_ix)) deallocate(kd_ix)
    if (allocated(dim_order)) deallocate(dim_order)
    if (allocated(work)) deallocate(work)
    if (allocated(var)) deallocate(var)
    if (allocated(subarray)) deallocate(subarray)
    if (allocated(perm)) deallocate(perm)
    if (allocated(stack_left)) deallocate(stack_left)
    if (allocated(stack_right)) deallocate(stack_right)
    allocate(X(d,n), kd_ix(n), dim_order(d), work(n), var(d), subarray(n), perm(n), stack_left(64), stack_right(64))
    X(:,1) = [1.0d0, 2.0d0]
    X(:,2) = [2.0d0, 3.0d0]
    X(:,3) = [3.0d0, 1.0d0]
    X(:,4) = [4.0d0, 0.0d0]
    X(:,5) = [0.0d0, 4.0d0]
    X(:,6) = [5.0d0, 2.0d0]
    do i = 1, d
        var(i) = sum((X(i,:) - sum(X(i,:))/n)**2) / n
    end do
    if (var(1) >= var(2)) then
        dim_order = [1,2]
    else
        dim_order = [2,1]
    end if
    call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, stack_left, stack_right)
    print *, '2D Cartesian KD Tree index:', kd_ix
    call assert_permutation(kd_ix, n, '2D cartesian')

    ! --- Test 2: 3D, n=8, spherical ---
    d = 3; n = 8
    if (allocated(V)) deallocate(V)
    if (allocated(sphere_ix)) deallocate(sphere_ix)
    if (allocated(dim_order)) deallocate(dim_order)
    if (allocated(work)) deallocate(work)
    if (allocated(var)) deallocate(var)
    if (allocated(subarray)) deallocate(subarray)
    if (allocated(perm)) deallocate(perm)
    if (allocated(stack_left)) deallocate(stack_left)
    if (allocated(stack_right)) deallocate(stack_right)
    allocate(V(d,n), sphere_ix(n), dim_order(d), work(n), var(d), subarray(n), perm(n), stack_left(64), stack_right(64))
    call random_unit_vectors(V, d, n)
    do i = 1, d
        var(i) = sum((V(i,:) - sum(V(i,:))/n)**2) / n
    end do
    call sort_dim_order(var, dim_order, d)
    call build_spherical_kd(V, d, n, sphere_ix, dim_order, work, subarray, perm, stack_left, stack_right)
    print *, '3D Spherical KD Tree index:', sphere_ix
    call assert_permutation(sphere_ix, n, '3D spherical')

    ! --- Edge Case: n=0 ---
    d = 2; n = 0
    if (allocated(X)) deallocate(X)
    if (allocated(kd_ix)) deallocate(kd_ix)
    if (allocated(dim_order)) deallocate(dim_order)
    if (allocated(work)) deallocate(work)
    if (allocated(var)) deallocate(var)
    if (allocated(subarray)) deallocate(subarray)
    if (allocated(perm)) deallocate(perm)
    if (allocated(stack_left)) deallocate(stack_left)
    if (allocated(stack_right)) deallocate(stack_right)
    allocate(X(d,0), kd_ix(0), dim_order(d), work(0), var(d), subarray(0), perm(0), stack_left(64), stack_right(64))
    call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, stack_left, stack_right)
    print *, 'Edge n=0 passed.'

    ! --- Edge Case: n=1 ---
    d = 2; n = 1
    if (allocated(X)) deallocate(X)
    if (allocated(kd_ix)) deallocate(kd_ix)
    if (allocated(dim_order)) deallocate(dim_order)
    if (allocated(work)) deallocate(work)
    if (allocated(var)) deallocate(var)
    if (allocated(subarray)) deallocate(subarray)
    if (allocated(perm)) deallocate(perm)
    if (allocated(stack_left)) deallocate(stack_left)
    if (allocated(stack_right)) deallocate(stack_right)
    allocate(X(d,1), kd_ix(1), dim_order(d), work(1), var(d), subarray(1), perm(1), stack_left(64), stack_right(64))
    X(:,1) = [1.0d0, 2.0d0]
    dim_order = [1,2]
    call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, stack_left, stack_right)
    if (kd_ix(1) == 1) then
        print *, 'Edge n=1 passed.'
    else
        print *, 'Edge n=1 failed.'
    end if

    ! --- Edge Case: identische Punkte ---
    d = 3; n = 5
    if (allocated(X)) deallocate(X)
    if (allocated(kd_ix)) deallocate(kd_ix)
    if (allocated(dim_order)) deallocate(dim_order)
    if (allocated(work)) deallocate(work)
    if (allocated(var)) deallocate(var)
    if (allocated(subarray)) deallocate(subarray)
    if (allocated(perm)) deallocate(perm)
    if (allocated(stack_left)) deallocate(stack_left)
    if (allocated(stack_right)) deallocate(stack_right)
    allocate(X(d,n), kd_ix(n), dim_order(d), work(n), var(d), subarray(n), perm(n), stack_left(64), stack_right(64))
    X = 1.0d0
    dim_order = [1,2,3]
    call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, stack_left, stack_right)
    call assert_permutation(kd_ix, n, 'identical points')

    ! --- Edge Case: Einheitsvektoren (spherical) ---
    d = 4; n = 4
    if (allocated(V)) deallocate(V)
    if (allocated(sphere_ix)) deallocate(sphere_ix)
    if (allocated(dim_order)) deallocate(dim_order)
    if (allocated(work)) deallocate(work)
    if (allocated(var)) deallocate(var)
    if (allocated(subarray)) deallocate(subarray)
    if (allocated(perm)) deallocate(perm)
    if (allocated(stack_left)) deallocate(stack_left)
    if (allocated(stack_right)) deallocate(stack_right)
    allocate(V(d,n), sphere_ix(n), dim_order(d), work(n), var(d), subarray(n), perm(n), stack_left(64), stack_right(64))
    V = 0.0d0
    do i = 1, d
        V(i,i) = 1.0d0
    end do
    dim_order = [1,2,3,4]
    call build_spherical_kd(V, d, n, sphere_ix, dim_order, work, subarray, perm, stack_left, stack_right)
    call assert_permutation(sphere_ix, n, 'unit vectors')

    ! --- Edge Case: hohe Dimension, wenig Punkte ---
    d = 10; n = 3
    if (allocated(X)) deallocate(X)
    if (allocated(kd_ix)) deallocate(kd_ix)
    if (allocated(dim_order)) deallocate(dim_order)
    if (allocated(work)) deallocate(work)
    if (allocated(var)) deallocate(var)
    if (allocated(subarray)) deallocate(subarray)
    if (allocated(perm)) deallocate(perm)
    if (allocated(stack_left)) deallocate(stack_left)
    if (allocated(stack_right)) deallocate(stack_right)
    allocate(X(d,n), kd_ix(n), dim_order(d), work(n), var(d), subarray(n), perm(n), stack_left(64), stack_right(64))
    call random_matrix(X, d, n)
    dim_order = [(i, i=1,d)]
    call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, stack_left, stack_right)
    call assert_permutation(kd_ix, n, 'high-d low-n')

    ! --- Edge Case: 1D, n=10, sortiert ---
    d = 1; n = 10
    if (allocated(X)) deallocate(X)
    if (allocated(kd_ix)) deallocate(kd_ix)
    if (allocated(dim_order)) deallocate(dim_order)
    if (allocated(work)) deallocate(work)
    if (allocated(var)) deallocate(var)
    if (allocated(subarray)) deallocate(subarray)
    if (allocated(perm)) deallocate(perm)
    if (allocated(stack_left)) deallocate(stack_left)
    if (allocated(stack_right)) deallocate(stack_right)
    allocate(X(d,n), kd_ix(n), dim_order(d), work(n), var(d), subarray(n), perm(n), stack_left(64), stack_right(64))
    do i = 1, n
        X(1,i) = i
    end do
    dim_order = [1]
    call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, stack_left, stack_right)
    call assert_permutation(kd_ix, n, '1D sorted')

    ! --- Additional Test: 2D, n=2 (minimal nontrivial) ---
    d = 2; n = 2
    if (allocated(X)) deallocate(X)
    if (allocated(kd_ix)) deallocate(kd_ix)
    if (allocated(dim_order)) deallocate(dim_order)
    if (allocated(work)) deallocate(work)
    if (allocated(var)) deallocate(var)
    if (allocated(subarray)) deallocate(subarray)
    if (allocated(perm)) deallocate(perm)
    if (allocated(stack_left)) deallocate(stack_left)
    if (allocated(stack_right)) deallocate(stack_right)
    allocate(X(d,n), kd_ix(n), dim_order(d), work(n), var(d), subarray(n), perm(n), stack_left(64), stack_right(64))
    X(:,1) = [1.0d0, 2.0d0]
    X(:,2) = [2.0d0, 1.0d0]
    dim_order = [1,2]
    call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, stack_left, stack_right)
    call assert_permutation(kd_ix, n, '2D n=2')

    ! --- Additional Test: 1D, n=1 (trivial) ---
    d = 1; n = 1
    if (allocated(X)) deallocate(X)
    if (allocated(kd_ix)) deallocate(kd_ix)
    if (allocated(dim_order)) deallocate(dim_order)
    if (allocated(work)) deallocate(work)
    if (allocated(var)) deallocate(var)
    if (allocated(subarray)) deallocate(subarray)
    if (allocated(perm)) deallocate(perm)
    if (allocated(stack_left)) deallocate(stack_left)
    if (allocated(stack_right)) deallocate(stack_right)
    allocate(X(d,n), kd_ix(n), dim_order(d), work(n), var(d), subarray(n), perm(n), stack_left(64), stack_right(64))
    X(1,1) = 42.0d0
    dim_order = [1]
    call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, stack_left, stack_right)
    call assert_permutation(kd_ix, n, '1D n=1')

    ! --- Additional Test: 1D, n=2 (minimal nontrivial) ---
    d = 1; n = 2
    if (allocated(X)) deallocate(X)
    if (allocated(kd_ix)) deallocate(kd_ix)
    if (allocated(dim_order)) deallocate(dim_order)
    if (allocated(work)) deallocate(work)
    if (allocated(var)) deallocate(var)
    if (allocated(subarray)) deallocate(subarray)
    if (allocated(perm)) deallocate(perm)
    if (allocated(stack_left)) deallocate(stack_left)
    if (allocated(stack_right)) deallocate(stack_right)
    allocate(X(d,n), kd_ix(n), dim_order(d), work(n), var(d), subarray(n), perm(n), stack_left(64), stack_right(64))
    X(1,1) = 1.0d0
    X(1,2) = 2.0d0
    dim_order = [1]
    call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, stack_left, stack_right)
    call assert_permutation(kd_ix, n, '1D n=2')

    ! --- Additional Test: 3D, n=100 (larger random) ---
    d = 3; n = 100
    if (allocated(X)) deallocate(X)
    if (allocated(kd_ix)) deallocate(kd_ix)
    if (allocated(dim_order)) deallocate(dim_order)
    if (allocated(work)) deallocate(work)
    if (allocated(var)) deallocate(var)
    if (allocated(subarray)) deallocate(subarray)
    if (allocated(perm)) deallocate(perm)
    if (allocated(stack_left)) deallocate(stack_left)
    if (allocated(stack_right)) deallocate(stack_right)
    allocate(X(d,n), kd_ix(n), dim_order(d), work(n), var(d), subarray(n), perm(n), stack_left(128), stack_right(128))
    call random_matrix(X, d, n)
    dim_order = [1,2,3]
    call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, stack_left, stack_right)
    allocate(val(d))
    call get_kd_point(X, kd_ix, 4, val)
    print *, "Point: ", val
    call assert_permutation(kd_ix, n, '3D n=100')

    ! --- Additional Test: 5D, n=10 (medium random) ---
    d = 5; n = 10
    if (allocated(X)) deallocate(X)
    if (allocated(kd_ix)) deallocate(kd_ix)
    if (allocated(dim_order)) deallocate(dim_order)
    if (allocated(work)) deallocate(work)
    if (allocated(var)) deallocate(var)
    if (allocated(subarray)) deallocate(subarray)
    if (allocated(perm)) deallocate(perm)
    if (allocated(stack_left)) deallocate(stack_left)
    if (allocated(stack_right)) deallocate(stack_right)
    allocate(X(d,n), kd_ix(n), dim_order(d), work(n), var(d), subarray(n), perm(n), stack_left(64), stack_right(64))
    call random_matrix(X, d, n)
    dim_order = [1,2,3,4,5]
    call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, stack_left, stack_right)
    call assert_permutation(kd_ix, n, '5D n=10')

    ! --- Additional Test: 2D, n=6, reversed order ---
    d = 2; n = 6
    if (allocated(X)) deallocate(X)
    if (allocated(kd_ix)) deallocate(kd_ix)
    if (allocated(dim_order)) deallocate(dim_order)
    if (allocated(work)) deallocate(work)
    if (allocated(var)) deallocate(var)
    if (allocated(subarray)) deallocate(subarray)
    if (allocated(perm)) deallocate(perm)
    if (allocated(stack_left)) deallocate(stack_left)
    if (allocated(stack_right)) deallocate(stack_right)
    allocate(X(d,n), kd_ix(n), dim_order(d), work(n), var(d), subarray(n), perm(n), stack_left(64), stack_right(64))
    do i = 1, n
        X(:,i) = [real(n-i+1,8), real(i,8)]
    end do
    dim_order = [1,2]
    call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, stack_left, stack_right)
    call assert_permutation(kd_ix, n, '2D reversed')

    ! --- Additional Test: 2D, n=6, all zeros ---
    d = 2; n = 6
    if (allocated(X)) deallocate(X)
    if (allocated(kd_ix)) deallocate(kd_ix)
    if (allocated(dim_order)) deallocate(dim_order)
    if (allocated(work)) deallocate(work)
    if (allocated(var)) deallocate(var)
    if (allocated(subarray)) deallocate(subarray)
    if (allocated(perm)) deallocate(perm)
    if (allocated(stack_left)) deallocate(stack_left)
    if (allocated(stack_right)) deallocate(stack_right)
    allocate(X(d,n), kd_ix(n), dim_order(d), work(n), var(d), subarray(n), perm(n), stack_left(64), stack_right(64))
    X = 0.0d0
    dim_order = [1,2]
    call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, stack_left, stack_right)
    call assert_permutation(kd_ix, n, '2D all zeros')

    print *, 'All edge-case tests completed.'

contains

    subroutine random_unit_vectors(V, d, n)
        real(8), intent(out) :: V(d,n)
        integer, intent(in) :: d, n
        integer :: i
        real(8) :: norm
        call random_seed()
        do i = 1, n
            call random_number(V(:,i))
            V(:,i) = V(:,i) - 0.5d0
            norm = sqrt(sum(V(:,i)**2))
            if (norm > 0) V(:,i) = V(:,i) / norm
        end do
    end subroutine random_unit_vectors

    subroutine random_matrix(X, d, n)
        real(8), intent(out) :: X(d,n)
        integer, intent(in) :: d, n
        integer :: i, j
        call random_seed()
        do j = 1, n
            call random_number(X(:,j))
        end do
    end subroutine random_matrix

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
        call sort_array(-var, perm, stack_left, stack_right)
        do i = 1, d
            dim_order(i) = perm(i)
        end do
    end subroutine sort_dim_order

    subroutine assert_permutation(ix, n, label)
        integer, intent(in) :: ix(:), n
        character(*), intent(in) :: label
        integer :: seen(max_n), i
        seen = 0
        do i = 1, n
            if (ix(i) < 1 .or. ix(i) > n) then
                print *, 'FAIL:', label, 'index out of bounds:', ix(i)
                return
            end if
            seen(ix(i)) = seen(ix(i)) + 1
        end do
        if (any(seen(1:n) /= 1)) then
            print *, 'FAIL:', label, 'not a permutation:', ix(1:min(n,10))
        else
            print *, 'PASS:', label
        end if
    end subroutine assert_permutation

end program test_kd_tree
