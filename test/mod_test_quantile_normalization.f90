! filepath: test/mod_test_quantile_normalization.f90
module mod_test_quantile_normalization
  use asserts
  implicit none
  public
contains

  subroutine run_all_tests_quantile_normalization()
    call test_qn_preserves_dimensions()
    call test_qn_identical_rows()
    call test_qn_no_nans_and_standardizes()
    call test_qn_single_row()
    call test_qn_single_column()
    call test_qn_all_equal()
    call test_qn_large_random()
    call test_qn_negative_values()
    call test_qn_zero_matrix()
    call test_qn_sorted_input()
    call test_qn_reverse_sorted()
    print *, "All quantile_normalization tests passed successfully."
  end subroutine

  subroutine test_qn_preserves_dimensions()
    integer :: nrow, ncol, max_stack
    real(8), dimension(2,3) :: mat, result
    real(8) :: temp_col(2), rank_means(2)
    integer :: perm(2), stack_left(10), stack_right(10)
    nrow = 2; ncol = 3; max_stack = 10
    mat(:,1) = [0.1d0, 0.2d0]
    mat(:,2) = [0.3d0, 0.4d0]
    mat(:,3) = [0.5d0, 0.6d0]
    call quantile_normalization_r(nrow, ncol, mat, result, temp_col, rank_means, perm, stack_left, stack_right, max_stack)
    call assert_equal_int(size(result,1), nrow, "qn: row count not preserved")
    call assert_equal_int(size(result,2), ncol, "qn: col count not preserved")
    end subroutine

    subroutine test_qn_identical_rows()
    integer :: nrow, ncol, max_stack
    real(8), dimension(2,3) :: mat, result, expected
    real(8) :: temp_col(2), rank_means(2)
    integer :: perm(2), stack_left(10), stack_right(10)
    nrow = 2; ncol = 3; max_stack = 10
    mat(1,:) = [5.0d0, 5.0d0, 5.0d0]
    mat(2,:) = [5.0d0, 5.0d0, 5.0d0]
    expected = mat
    call quantile_normalization_r(nrow, ncol, mat, result, temp_col, rank_means, perm, stack_left, stack_right, max_stack)
    call assert_equal_array_real(result, expected, 6, 1d-12, "qn: identical rows not preserved")
    call assert_true(all(reshape(isfinite_mat(result), [size(result)])), "qn: result has non-finite values")
    end subroutine

    subroutine test_qn_no_nans_and_standardizes()
    integer :: nrow, ncol, i, max_stack
    real(8), dimension(2,3) :: mat, result
    real(8), allocatable :: sorted_cols(:,:)
    real(8) :: temp_col(2), rank_means(2)
    integer :: perm(2), stack_left(10), stack_right(10)
    nrow = 2; ncol = 3; max_stack = 10
    mat(:,1) = [2.0d0, 0.0d0]
    mat(:,2) = [5.0d0, 3.0d0]
    mat(:,3) = [7.0d0, 1.0d0]
    call quantile_normalization_r(nrow, ncol, mat, result, temp_col, rank_means, perm, stack_left, stack_right, max_stack)
    call assert_true(.not. any(reshape(isnan_mat(result), [size(result)])), "qn: result has NaN")
    call assert_true(all(reshape(isfinite_mat(result), [size(result)])), "qn: result has non-finite values")
    call assert_equal_int(size(result,1), nrow, "qn: row count not preserved")
    call assert_equal_int(size(result,2), ncol, "qn: col count not preserved")
    allocate(sorted_cols(nrow, ncol))
    do i = 1, ncol
        sorted_cols(:,i) = sort_vec(result(:,i))
    end do
    do i = 2, ncol
        call assert_equal_array_real(sorted_cols(:,i), sorted_cols(:,1), nrow, 1d-6, "qn: column distributions differ")
    end do
    deallocate(sorted_cols)
    end subroutine

    subroutine test_qn_single_row()
    integer :: nrow, ncol, max_stack
    real(8), dimension(1,4) :: mat, result, expected
    real(8) :: temp_col(1), rank_means(1)
    integer :: perm(1), stack_left(10), stack_right(10)
    nrow = 1; ncol = 4; max_stack = 10
    mat(1,:) = [1.0d0, 2.0d0, 3.0d0, 4.0d0]
    expected = mat
    call quantile_normalization_r(nrow, ncol, mat, result, temp_col, rank_means, perm, stack_left, stack_right, max_stack)
    call assert_equal_array_real(result, expected, 4, 1d-12, "qn: single row not preserved")
    end subroutine

    subroutine test_qn_single_column()
    integer :: nrow, ncol, max_stack
    real(8), dimension(3,1) :: mat, result, expected
    real(8) :: temp_col(3), rank_means(3)
    integer :: perm(3), stack_left(10), stack_right(10)
    nrow = 3; ncol = 1; max_stack = 10
    mat(:,1) = [7.0d0, 2.0d0, 5.0d0]
    expected = mat
    call quantile_normalization_r(nrow, ncol, mat, result, temp_col, rank_means, perm, stack_left, stack_right, max_stack)
    call assert_equal_array_real(result, expected, 3, 1d-12, "qn: single column not preserved")
    end subroutine

    subroutine test_qn_all_equal()
    integer :: nrow, ncol, max_stack
    real(8), dimension(4,4) :: mat, result, expected
    real(8) :: temp_col(4), rank_means(4)
    integer :: perm(4), stack_left(10), stack_right(10)
    nrow = 4; ncol = 4; max_stack = 10
    mat = 3.14d0
    expected = mat
    call quantile_normalization_r(nrow, ncol, mat, result, temp_col, rank_means, perm, stack_left, stack_right, max_stack)
    call assert_equal_array_real(result, expected, 16, 1d-12, "qn: all equal values not preserved")
    end subroutine

    subroutine test_qn_large_random()
    integer, parameter :: nrow=10, ncol=10
    integer :: i, max_stack
    real(8), dimension(nrow,ncol) :: mat, result
    real(8), allocatable :: sorted_cols(:,:)
    real(8) :: temp_col(nrow), rank_means(nrow)
    integer :: perm(nrow), stack_left(100), stack_right(100)
    max_stack = 100
    call random_seed()
    call random_number(mat)
    call quantile_normalization_r(nrow, ncol, mat, result, temp_col, rank_means, perm, stack_left, stack_right, max_stack)
    allocate(sorted_cols(nrow, ncol))
    do i = 1, ncol
        sorted_cols(:,i) = sort_vec(result(:,i))
    end do
    do i = 2, ncol
        call assert_equal_array_real(sorted_cols(:,i), sorted_cols(:,1), nrow, 1d-6, &
                                "qn: random matrix column distributions differ")
    end do
    deallocate(sorted_cols)
    end subroutine

    subroutine test_qn_negative_values()
    integer :: nrow, ncol, max_stack
    real(8), dimension(2,3) :: mat, result
    real(8) :: temp_col(2), rank_means(2)
    integer :: perm(2), stack_left(10), stack_right(10)
    nrow = 2; ncol = 3; max_stack = 10
    mat(:,1) = [-1.0d0, -2.0d0]
    mat(:,2) = [-3.0d0, -4.0d0]
    mat(:,3) = [-5.0d0, -6.0d0]
    call quantile_normalization_r(nrow, ncol, mat, result, temp_col, rank_means, perm, stack_left, stack_right, max_stack)
    ! You can add further asserts if you want to check the expected result
    end subroutine

    subroutine test_qn_zero_matrix()
    integer :: nrow, ncol, max_stack
    real(8), dimension(3,3) :: mat, result, expected
    real(8) :: temp_col(3), rank_means(3)
    integer :: perm(3), stack_left(10), stack_right(10)
    nrow = 3; ncol = 3; max_stack = 10
    mat = 0.0d0
    expected = 0.0d0
    call quantile_normalization_r(nrow, ncol, mat, result, temp_col, rank_means, perm, stack_left, stack_right, max_stack)
    call assert_equal_array_real(result, expected, 9, 1d-12, "qn: zero matrix not preserved")
    end subroutine

    subroutine test_qn_sorted_input()
    integer :: nrow, ncol, max_stack
    real(8), dimension(3,3) :: mat, result
    real(8) :: temp_col(3), rank_means(3)
    integer :: perm(3), stack_left(10), stack_right(10)
    nrow = 3; ncol = 3; max_stack = 10
    mat(:,1) = [1.0d0, 2.0d0, 3.0d0]
    mat(:,2) = [1.0d0, 2.0d0, 3.0d0]
    mat(:,3) = [1.0d0, 2.0d0, 3.0d0]
    call quantile_normalization_r(nrow, ncol, mat, result, temp_col, rank_means, perm, stack_left, stack_right, max_stack)
    call assert_equal_array_real(result, mat, 9, 1d-12, "qn: sorted input not preserved")
    end subroutine

    subroutine test_qn_reverse_sorted()
    integer :: nrow, ncol, max_stack
    real(8), dimension(3,3) :: mat, result, expected
    real(8) :: temp_col(3), rank_means(3)
    integer :: perm(3), stack_left(10), stack_right(10)
    nrow = 3; ncol = 3; max_stack = 10
    mat(:,1) = [3.0d0, 3.0d0, 3.0d0]
    mat(:,2) = [2.0d0, 2.0d0, 2.0d0]
    mat(:,3) = [1.0d0, 1.0d0, 1.0d0]
    expected = 2.0d0
    call quantile_normalization_r(nrow, ncol, mat, result, temp_col, rank_means, perm, stack_left, stack_right, max_stack)
    call assert_equal_array_real(result, expected, 9, 1d-12, "qn: reverse sorted input not normalized to mean")
    end subroutine

  function isfinite_mat(arr) result(mask)
    real(8), intent(in) :: arr(:,:)
    logical :: mask(size(arr,1), size(arr,2))
    integer :: i, j
    do i = 1, size(arr,1)
      do j = 1, size(arr,2)
        mask(i,j) = abs(arr(i,j)) < huge(arr(i,j)) .and. arr(i,j) == arr(i,j)
      end do
    end do
  end function isfinite_mat

  function isnan_mat(arr) result(mask)
    real(8), intent(in) :: arr(:,:)
    logical :: mask(size(arr,1), size(arr,2))
    integer :: i, j
    do i = 1, size(arr,1)
      do j = 1, size(arr,2)
        mask(i,j) = arr(i,j) /= arr(i,j)
      end do
    end do
  end function isnan_mat

  function sort_vec(vec) result(sorted)
    real(8), intent(in) :: vec(:)
    real(8) :: sorted(size(vec))
    integer :: i, j
    sorted = vec
    do i = 1, size(vec)-1
      do j = i+1, size(vec)
        if (sorted(j) < sorted(i)) then
          call swap(sorted(i), sorted(j))
        end if
      end do
    end do
  end function sort_vec

  subroutine swap(a, b)
    real(8), intent(inout) :: a, b
    real(8) :: tmp
    tmp = a
    a = b
    b = tmp
  end subroutine swap

end module mod_test_quantile_normalization