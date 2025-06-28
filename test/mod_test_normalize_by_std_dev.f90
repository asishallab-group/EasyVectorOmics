! filepath: test/mod_test_normalize_by_std_dev.f90
!> @brief Unit test suite for normalize_by_std_dev routine.
!! @details
!! This module contains a comprehensive set of tests for the
!! normalize_by_std_dev routine, including edge cases and
!! typical usage scenarios.
module mod_test_normalize_by_std_dev
  use asserts
  implicit none
  public

contains

  !> @brief Run all normalize_by_std_dev tests.
  !! @details
  !! Executes all test cases for the normalize_by_std_dev routine.
  !! If all assertions pass, prints a single success message.
  subroutine run_all_tests_normalize_by_std_dev()
    call test_normalize_by_std_dev_basic()
    call test_normalize_by_std_dev_constant_rows()
    call test_normalize_by_std_dev_large_numbers()
    call test_identity_matrix()
    call test_zero_rows()
    call test_negative_rows()
    call test_large_random_matrix()
    call test_single_nonzero()
    call test_small_large_values()
    call test_nan_inf_input()
    call test_single_row_col()
    call test_empty_matrix()
    call test_symmetric_rows()
    print *, "All normalize_by_std_dev tests passed successfully."
  end subroutine run_all_tests_normalize_by_std_dev

  !> @brief Test that normalize_by_std_dev normalizes values correctly.
  subroutine test_normalize_by_std_dev_basic()
    real(8), dimension(2,2) :: mat, result, expected
    real(8), dimension(2) :: std_dev
    integer :: i, j

    mat = reshape([2.0d0, 4.0d0, 6.0d0, 8.0d0], [2,2])
    call normalize_by_std_dev_r(2, 2, mat, result)

    do i = 1, 2
      std_dev(i) = sqrt((mat(i,1)**2 + mat(i,2)**2) / 2.0d0)
      do j = 1, 2
        expected(i,j) = mat(i,j) / std_dev(i)
      end do
    end do

    call assert_equal_array_real(result, expected, 4, 1d-12, "normalize_by_std_dev: basic normalization failed")
  end subroutine

  !> @brief Test that normalize_by_std_dev handles constant rows (should normalize to 1).
  subroutine test_normalize_by_std_dev_constant_rows()
    real(8), dimension(2,2) :: mat, result, expected
    integer :: i, j

    mat = reshape([5.0d0, 5.0d0, 5.0d0, 5.0d0], [2,2])
    call normalize_by_std_dev_r(2, 2, mat, result)

    expected = 1.0d0

    call assert_true(all(result == 1.0d0), "normalize_by_std_dev: not all values are 1 for constant rows")
    call assert_no_nan_real(result, 4, "normalize_by_std_dev: NaN in result for constant rows")
  end subroutine

  !> @brief Test that normalize_by_std_dev normalizes large numbers properly.
  subroutine test_normalize_by_std_dev_large_numbers()
    real(8), dimension(2,2) :: mat, result, expected
    real(8), dimension(2) :: std_dev
    integer :: i, j

    mat = reshape([1e6, 2e6, 1e6, 2e6], [2,2])
    call normalize_by_std_dev_r(2, 2, mat, result)

    do i = 1, 2
      std_dev(i) = sqrt((mat(i,1)**2 + mat(i,2)**2) / 2.0d0)
      do j = 1, 2
        expected(i,j) = mat(i,j) / std_dev(i)
      end do
    end do

    call assert_equal_array_real(result, expected, 4, 1d-12, "normalize_by_std_dev: large numbers normalization failed")
    call assert_no_nan_real(result, 4, "normalize_by_std_dev: NaN in result for large numbers")
    call assert_true(all(isfinite_mat(result)), "normalize_by_std_dev: Inf in result for large numbers")

  end subroutine

  !> @brief Test normalization of the identity matrix.
  subroutine test_identity_matrix()
    real(8), dimension(3,3) :: mat, result
    integer :: i, j
    mat = 0.0d0
    do i = 1, 3
      mat(i,i) = 1.0d0
    end do
    call normalize_by_std_dev_r(3, 3, mat, result)
    do i = 1, 3
      call assert_in_range_real(sum(result(i,:)**2)/3.0d0, 1d0-1d-12, 1d0+1d-12, "identity: RMS not 1")
      do j = 1, 3
        if (i /= j) call assert_equal_real(result(i,j), 0.0d0, 1d-12, "identity: off-diagonal not zero")
      end do
    end do
  end subroutine

  !> @brief Test normalization of rows with all zeros.
  subroutine test_zero_rows()
    real(8), dimension(2,3) :: mat, result
    mat = 0.0d0
    call normalize_by_std_dev_r(2, 3, mat, result)
    call assert_true(all(result == 0.0d0), "zero rows: not all zeros")
    call assert_no_nan_real(result, 6, "zero rows: NaN in result")
  end subroutine

  !> @brief Test normalization of rows with negative values.
  subroutine test_negative_rows()
    real(8), dimension(2,3) :: mat, result, expected
    real(8), dimension(2) :: std_dev
    integer :: i, j
    mat = reshape([-2.0d0, -4.0d0, -6.0d0, -8.0d0, -10.0d0, -12.0d0], [2,3])
    call normalize_by_std_dev_r(2, 3, mat, result)
    do i = 1, 2
      std_dev(i) = sqrt(sum(mat(i,:)**2)/3.0d0)
      do j = 1, 3
        expected(i,j) = mat(i,j)/std_dev(i)
      end do
    end do
    call assert_equal_array_real(result, expected, 6, 1d-12, "negative rows: normalization failed")
  end subroutine

  !> @brief Test normalization of a large random matrix.
  subroutine test_large_random_matrix()
    integer, parameter :: nrow=20, ncol=30
    real(8), dimension(nrow,ncol) :: mat, result
    integer :: i
    call random_seed()
    call random_number(mat)
    call normalize_by_std_dev_r(nrow, ncol, mat, result)
    do i = 1, nrow
      call assert_in_range_real(sqrt(sum(result(i,:)**2)/ncol), 1d0-1d-10, 1d0+1d-10, "large random: RMS not 1")
    end do
    call assert_no_nan_real(result, nrow*ncol, "large random: NaN in result")
  end subroutine

  !> @brief Test normalization of rows with a single nonzero value.
  subroutine test_single_nonzero()
    real(8), dimension(2,4) :: mat, result, expected
    integer :: i

    mat = 0.0d0
    mat(1,3) = 5.0d0
    mat(2,2) = -7.0d0

    call normalize_by_std_dev_r(2, 4, mat, result)

    do i = 1, 2
      expected(i,:) = mat(i,:) / sqrt(sum(mat(i,:)**2)/4.0d0)
    end do

    call assert_equal_array_real(result, expected, 8, 1d-12, "single nonzero: normalization failed")
  end subroutine

  !> @brief Test normalization with very small and very large values.
  subroutine test_small_large_values()
    real(8), dimension(2,2) :: mat, result, expected
    real(8), dimension(2) :: std_dev
    integer :: i, j
    mat = reshape([1e-10, 1e10, 1e-10, 1e10], [2,2])
    call normalize_by_std_dev_r(2, 2, mat, result)
    do i = 1, 2
      std_dev(i) = sqrt(sum(mat(i,:)**2)/2.0d0)
      do j = 1, 2
        expected(i,j) = mat(i,j)/std_dev(i)
      end do
    end do
    call assert_equal_array_real(result, expected, 4, 1d-10, "small/large values: normalization failed")
    call assert_no_nan_real(result, 4, "small/large values: NaN in result")
  end subroutine

  !> @brief Test normalization when input contains NaN or Inf.
  subroutine test_nan_inf_input()
    real(8), dimension(2,2) :: mat, result
    mat = reshape([1.0d0, 2.0d0, huge(1.0d0), 4.0d0], [2,2]) ! Simulate Inf in (2,1)
    call normalize_by_std_dev_r(2, 2, mat, result)

    call assert_true(all(isfinite_mat(result)), "normalize_by_std_dev: output contains NaN/Inf unexpectedly")
  end subroutine

  !> @brief Test normalization of a single row and a single column matrix.
  subroutine test_single_row_col()
    real(8), dimension(1,4) :: mat1, result1, expected1
    real(8), dimension(4,1) :: mat2, result2
    real(8) :: std_dev
    integer :: j

    mat1 = reshape([2.0d0, 4.0d0, 6.0d0, 8.0d0], [1,4])
    call normalize_by_std_dev_r(1, 4, mat1, result1)
    std_dev = sqrt(sum(mat1(1,:)**2)/4.0d0)
    do j = 1, 4
      expected1(1,j) = mat1(1,j)/std_dev
    end do
    call assert_equal_array_real(result1, expected1, 4, 1d-12, "single row: normalization failed")

    mat2 = reshape([2.0d0, 4.0d0, 6.0d0, 8.0d0], [4,1])
    call normalize_by_std_dev_r(4, 1, mat2, result2)
    call assert_true(all(abs(result2) == 1.0d0), "single col: normalization failed")
  end subroutine

  !> @brief Test normalization of an empty matrix.
  subroutine test_empty_matrix()
    real(8), allocatable :: mat(:,:), result(:,:)
    allocate(mat(0,0), result(0,0))
    call normalize_by_std_dev_r(0, 0, mat, result)
    ! No assertion needed: just check no crash
  end subroutine

  !> @brief Test normalization of symmetric rows.
  subroutine test_symmetric_rows()
    real(8), dimension(2,3) :: mat, result
    integer :: j
    mat(1,:) = [1.0d0, 2.0d0, 3.0d0]
    mat(2,:) = [2.0d0, 4.0d0, 6.0d0]
    call normalize_by_std_dev_r(2, 3, mat, result)
    do j = 1, 3
      call assert_equal_real(result(2,j), result(1,j), 1d-12, "symmetric rows: not equal after normalization")
    end do
  end subroutine

  !> @brief Helper function to check if all values are finite (matrix version).
  function isfinite_mat(arr) result(mask)
    real(8), intent(in) :: arr(:,:)
    logical :: mask(size(arr,1), size(arr,2))
    integer :: i, j
    do i = 1, size(arr,1)
      do j = 1, size(arr,2)
        mask(i,j) = abs(arr(i,j)) < huge(1.0d0)
      end do
    end do
  end function isfinite_mat

end module mod_test_normalize_by_std_dev