!> @file mod_test_loess_smoothing.f90
!> Unit test suite for LOESS smoothing and compute_edf from f42_utils
!> @details Unit tests for LOESS smoothing, compute_edf, including masking and edge cases.

module mod_test_loess_smoothing
  use asserts
  use f42_utils
  use tox_errors
  use, intrinsic :: iso_fortran_env, only: real64, int32
  implicit none
  public

  ! Abstract interface for all test procedures
  abstract interface
    subroutine test_interface()
    end subroutine test_interface
  end interface

  ! Type to hold test name and procedure pointer
  type :: test_case
    character(len=64) :: name
    procedure(test_interface), pointer, nopass :: test_proc => null()
  end type test_case

contains

  !> Get array of all available LOESS tests and EDF tests.
  subroutine get_all_tests(all_tests)
    type(test_case), intent(out) :: all_tests(22)
    all_tests(1) = test_case("test_loess_constant_input", test_loess_constant_input)
    all_tests(2) = test_case("test_loess_linear_trend", test_loess_linear_trend)
    all_tests(3) = test_case("test_loess_outlier_suppression", test_loess_outlier_suppression)
    all_tests(4) = test_case("test_loess_sparse_fallback", test_loess_sparse_fallback)
    all_tests(5) = test_case("test_loess_single_point", test_loess_single_point)
    all_tests(6) = test_case("test_loess_identical_points", test_loess_identical_points)
    all_tests(7) = test_case("test_loess_linear_interp", test_loess_linear_interp)
    all_tests(8) = test_case("test_loess_weight_decay", test_loess_weight_decay)
    all_tests(9) = test_case("test_loess_mask_exclusion", test_loess_mask_exclusion)
    all_tests(10) = test_case("test_loess_fallback", test_loess_fallback)
    all_tests(11) = test_case("test_loess_edge_query", test_loess_edge_query)
    all_tests(12) = test_case("test_loess_invalid_dimensions", test_loess_invalid_dimensions)
    all_tests(13) = test_case("test_loess_invalid_parameters", test_loess_invalid_parameters)
    all_tests(14) = test_case("test_loess_invalid_indices", test_loess_invalid_indices)
    all_tests(15) = test_case("test_edf_simple", test_edf_simple)
    all_tests(16) = test_case("test_edf_all_unique", test_edf_all_unique)
    all_tests(17) = test_case("test_edf_all_same", test_edf_all_same)
    all_tests(18) = test_case("test_edf_duplicates", test_edf_duplicates)
    all_tests(19) = test_case("test_edf_single_value", test_edf_single_value)
    all_tests(20) = test_case("test_edf_empty_input", test_edf_empty_input)
    all_tests(21) = test_case("test_edf_large_dataset", test_edf_large_dataset)
    all_tests(22) = test_case("test_edf_negative_values", test_edf_negative_values)
  end subroutine get_all_tests

  !> Run all LOESS smoothing tests.
  subroutine run_all_tests_loess_smoothing()
    type(test_case) :: all_tests(22)
    integer(int32) :: i
    call get_all_tests(all_tests)
    do i = 1, size(all_tests)
      call all_tests(i)%test_proc()
      print *, trim(all_tests(i)%name), " passed."
    end do
    print *, "All LOESS smoothing tests passed successfully."
  end subroutine run_all_tests_loess_smoothing

  !> Run specific LOESS smoothing tests by name.
  subroutine run_named_tests_loess_smoothing(test_names)
    character(len=*), intent(in) :: test_names(:)
    type(test_case) :: all_tests(22)
    integer(int32) :: i, j
    logical :: found
    call get_all_tests(all_tests)
    do i = 1, size(test_names)
      found = .false.
      do j = 1, size(all_tests)
        if (trim(test_names(i)) == trim(all_tests(j)%name)) then
          call all_tests(j)%test_proc()
          print *, trim(test_names(i)), " passed."
          found = .true.
          exit
        end if
      end do
      if (.not. found) then
        print *, "Unknown test: ", trim(test_names(i))
      end if
    end do
  end subroutine run_named_tests_loess_smoothing

  subroutine test_loess_constant_input()
    integer(int32), parameter :: n = 100
    real(real64) :: x_ref(n), y_ref(n), x_query(n/2), y_out(n/2)
    integer(int32) :: indices_used(n), i, ierr
    x_ref = 5.0_real64; y_ref = 10.0_real64; x_query = 5.0_real64; indices_used = [(i, i = 1, n)]
    call loess_smooth_2d(n, n/2, x_ref, y_ref, indices_used, n, x_query, 1.0_real64, 3.0_real64, y_out, ierr)
    call assert_equal_int(ierr, 0, 'Constant input error check')
    call assert_true(all(abs(y_out - 10.0_real64) < 1.0e-6_real64), 'Constant input test')
  end subroutine

  subroutine test_loess_linear_trend()
    integer(int32), parameter :: n = 100
    real(real64) :: x_ref(n), y_ref(n), x_query(n/2), y_out(n/2)
    integer(int32) :: indices_used(n), i, ierr
    x_ref = [(real(i, kind=real64), i = 1, n)]
    y_ref = 0.5_real64 * x_ref
    x_query = [(real(i, kind=real64) + 0.5_real64, i = 1, n/2)]
    indices_used = [(i, i = 1, n)]
    call loess_smooth_2d(n, n/2, x_ref, y_ref, indices_used, n, x_query, 1.0_real64, 3.0_real64, y_out, ierr)
    call assert_equal_int(ierr, 0, 'Linear trend error check')
    call assert_true(all(abs(y_out - 0.5_real64 * x_query) < 0.05_real64), 'Linear trend test')
  end subroutine

  subroutine test_loess_outlier_suppression()
    integer(int32), parameter :: n = 100
    real(real64) :: x_ref(n), y_ref(n), x_query(50), y_out(50)
    integer(int32) :: indices_used(n), i, ierr
    x_ref(1:n-1) = 10.0_real64
    x_ref(n) = 100.0_real64
    y_ref(1:n-1) = 5.0_real64
    y_ref(n) = 99.0_real64
    x_query = 10.0_real64
    indices_used = [(i, i = 1, n)]
    call loess_smooth_2d(n, 50, x_ref, y_ref, indices_used, n, x_query, 1.0_real64, 3.0_real64, y_out, ierr)
    call assert_equal_int(ierr, 0, 'Outlier suppression error check')
    call assert_true(all(abs(y_out - 5.0_real64) < 0.01_real64), 'Outlier suppression test')
  end subroutine

  subroutine test_loess_sparse_fallback()
    integer(int32), parameter :: n = 100
    real(real64) :: x_ref(n), y_ref(n), x_query(50), y_out(50)
    integer(int32) :: indices_used(n), i, ierr
    x_ref = [(real(i, kind=real64) * 100.0_real64, i = 1, n)]
    y_ref = [(real(i, kind=real64), i = 1, n)]
    x_query = [(real(i, kind=real64) * 100.0_real64 + 50.0_real64, i = 1, 50)]
    indices_used = [(i, i = 1, n)]
    call loess_smooth_2d(n, 50, x_ref, y_ref, indices_used, n, x_query, 1.0_real64, 3.0_real64, y_out, ierr)
    call assert_equal_int(ierr, 0, 'Sparse fallback error check')
    call assert_true(all(abs(y_out - y_ref(indices_used(1:50))) < 1.0e-6_real64), 'Sparse fallback test')
  end subroutine

  subroutine test_loess_single_point()
    real(real64) :: x1(1), y1(1), xq1(1), yout1(1)
    integer(int32) :: ierr
    x1 = 0.0_real64; y1(1) = 42.0_real64; xq1 = 0.0_real64
    call loess_smooth_2d(1, 1, x1, y1, (/1/), 1, xq1, 0.1_real64, 1.0_real64, yout1, ierr)
    call assert_equal_int(ierr, 0, 'Single point error check')
    call assert_equal_real(yout1(1), y1(1), 1e-6_real64, 'Single point test')
  end subroutine

  subroutine test_loess_identical_points()
    real(real64) :: x2(2), y2(2), xq2(1), yout2(1)
    integer(int32) :: idx2(2), ierr
    x2 = (/0.0_real64, 1.0_real64/); y2 = (/1.0_real64, 1.0_real64/); xq2 = 0.0_real64; idx2 = (/1,2/)
    call loess_smooth_2d(2, 1, x2, y2, idx2, 2, xq2, 0.1_real64, 1.0_real64, yout2, ierr)
    call assert_equal_int(ierr, 0, 'Identical points error check')
    call assert_equal_real(yout2(1), y2(1), 1e-6_real64, 'Identical points test')
  end subroutine

  subroutine test_loess_linear_interp()
    real(real64) :: x3(2), y3(2), xq3(1), yout3(1)
    integer(int32) :: idx3(2), ierr
    x3 = (/0.0_real64, 2.0_real64/); y3 = (/0.0_real64, 2.0_real64/); xq3 = 1.0_real64; idx3 = (/1,2/)
    call loess_smooth_2d(2, 1, x3, y3, idx3, 2, xq3, 0.5_real64, 3.0_real64, yout3, ierr)
    call assert_equal_int(ierr, 0, 'Linear interpolation error check')
    call assert_true(abs(yout3(1) - 1.0_real64) < 0.1_real64, 'Linear interpolation test')
  end subroutine

  subroutine test_loess_weight_decay()
    real(real64) :: x2(2), y2(2), xq2(1), yout2(1)
    integer(int32) :: idx2(2), ierr
    x2 = (/0.0_real64, 10.0_real64/); y2 = (/0.0_real64, 10.0_real64/); xq2 = 0.0_real64; idx2 = (/1,2/)
    call loess_smooth_2d(2, 1, x2, y2, idx2, 2, xq2, 1.0_real64, 3.0_real64, yout2, ierr)
    call assert_equal_int(ierr, 0, 'Weight decay error check')
    call assert_true(yout2(1) < 5.0_real64, 'Weight decay test')
  end subroutine

  subroutine test_loess_mask_exclusion()
    real(real64) :: x3(3), y3(3), xq3(1), yout3(1)
    integer(int32) :: indices_used(1), ierr
    x3 = (/0.0_real64, 10.0_real64, 20.0_real64/)
    y3 = (/0.0_real64, 10.0_real64, 20.0_real64/)
    xq3 = 0.0_real64
    ! Pre-filter: only use index 2 (the middle point)
    indices_used = (/2/)  ! Use only the middle point
    call loess_smooth_2d(3, 1, x3, y3, indices_used, 1, xq3, 10.0_real64, 3.0_real64, yout3, ierr)
    call assert_equal_int(ierr, 0, 'Mask exclusion error check')
    call assert_true(abs(yout3(1) - 10.0_real64) < 1.0_real64, 'Mask exclusion test')
  end subroutine

  subroutine test_loess_fallback()
    real(real64) :: x1(1), y1(1), xq1(1), yout1(1)
    integer(int32) :: ierr
    x1 = 0.0_real64; y1(1) = 123.0_real64; xq1 = 100.0_real64
    call loess_smooth_2d(1, 1, x1, y1, (/1/), 1, xq1, 0.1_real64, 1.0_real64, yout1, ierr)
    call assert_equal_int(ierr, 0, 'Fallback error check')
    call assert_equal_real(yout1(1), y1(1), 1e-6_real64, 'Fallback test')
  end subroutine

  subroutine test_loess_edge_query()
    ! Test edge query points (extrapolation)
    real(real64) :: x(5), y(5), xq(2), yout(2)
    integer(int32) :: idx(5), ierr
    x = (/1.0, 2.0, 3.0, 4.0, 5.0/)
    y = (/10.0, 20.0, 30.0, 40.0, 50.0/)
    idx = (/1,2,3,4,5/)
    xq = (/0.0, 6.0/)
    call loess_smooth_2d(5, 2, x, y, idx, 5, xq, 1.0_real64, 3.0_real64, yout, ierr)
    call assert_equal_int(ierr, 0, 'Edge query error check')
    call assert_true(abs(yout(1) - 10.0_real64) < 1.0 .and. abs(yout(2) - 50.0_real64) < 1.0, 'Edge query test')
  end subroutine

  ! Error handling tests
  subroutine test_loess_invalid_dimensions()
    real(real64) :: x(5), y(5), xq(2), yout(2)
    integer(int32) :: idx(5), ierr
    x = (/1.0, 2.0, 3.0, 4.0, 5.0/)
    y = (/10.0, 20.0, 30.0, 40.0, 50.0/)
    idx = (/1,2,3,4,5/)
    xq = (/0.0, 6.0/)
    
    ! Test with zero dimensions
    call loess_smooth_2d(0, 2, x, y, idx, 5, xq, 1.0_real64, 3.0_real64, yout, ierr)
    call assert_equal_int(ierr, 202, 'Zero n_total error check')
    
    call loess_smooth_2d(5, 0, x, y, idx, 5, xq, 1.0_real64, 3.0_real64, yout, ierr)
    call assert_equal_int(ierr, 202, 'Zero n_target error check')
    
    call loess_smooth_2d(5, 2, x, y, idx, 0, xq, 1.0_real64, 3.0_real64, yout, ierr)
    call assert_equal_int(ierr, 202, 'Zero n_used error check')
  end subroutine

  subroutine test_loess_invalid_parameters()
    real(real64) :: x(5), y(5), xq(2), yout(2)
    integer(int32) :: idx(5), ierr
    x = (/1.0, 2.0, 3.0, 4.0, 5.0/)
    y = (/10.0, 20.0, 30.0, 40.0, 50.0/)
    idx = (/1,2,3,4,5/)
    xq = (/0.0, 6.0/)
    
    ! Test with negative kernel_sigma
    call loess_smooth_2d(5, 2, x, y, idx, 5, xq, -1.0_real64, 3.0_real64, yout, ierr)
    call assert_equal_int(ierr, 201, 'Negative sigma error check')
    
    ! Test with negative kernel_cutoff
    call loess_smooth_2d(5, 2, x, y, idx, 5, xq, 1.0_real64, -3.0_real64, yout, ierr)
    call assert_equal_int(ierr, 201, 'Negative cutoff error check')
  end subroutine

  subroutine test_loess_invalid_indices()
    real(real64) :: x(5), y(5), xq(2), yout(2)
    integer(int32) :: idx(5), ierr
    x = (/1.0, 2.0, 3.0, 4.0, 5.0/)
    y = (/10.0, 20.0, 30.0, 40.0, 50.0/)
    xq = (/0.0, 6.0/)
    
    ! Test with index out of bounds (too high)
    idx = (/1, 2, 3, 4, 6/)  ! index 6 is out of bounds for n_total=5
    call loess_smooth_2d(5, 2, x, y, idx, 5, xq, 1.0_real64, 3.0_real64, yout, ierr)
    call assert_equal_int(ierr, 201, 'Index too high error check')
    
    ! Test with index out of bounds (too low)
    idx = (/0, 2, 3, 4, 5/)  ! index 0 is out of bounds
    call loess_smooth_2d(5, 2, x, y, idx, 5, xq, 1.0_real64, 3.0_real64, yout, ierr)
    call assert_equal_int(ierr, 201, 'Index too low error check')
  end subroutine

  ! ============================================================================
  ! EDF Tests
  ! ============================================================================

  !> @brief Test EDF computation with a simple dataset.
  !> Input: [1.0, 3.0, 1.0, 2.0, 3.0, 3.0]
  !> Expected unique values: [1.0, 2.0, 3.0]
  !> Expected CDF: [0.333..., 0.5, 1.0]
  subroutine test_edf_simple()
    real(real64), dimension(6) :: values = [1.0_real64, 3.0_real64, 1.0_real64, &
                                             2.0_real64, 3.0_real64, 3.0_real64]
    integer(int32) :: n_values = 6
    integer(int32) :: perm(6), stack_left(6), stack_right(6)
    real(real64) :: unique_values(6), cdf_values(6)
    real(real64) :: expected_unique(3), expected_cdf(3)
    integer(int32) :: n_unique, ierr
    
    expected_unique = [1.0_real64, 2.0_real64, 3.0_real64]
    expected_cdf = [2.0_real64/6.0_real64, 3.0_real64/6.0_real64, 6.0_real64/6.0_real64]
    
    call compute_edf(values, n_values, perm, stack_left, stack_right, &
                     unique_values, cdf_values, n_unique, ierr)
    
    call assert_equal_int(ierr, ERR_OK, "test_edf_simple: error code should be ERR_OK")
    call assert_equal_int(n_unique, 3, "test_edf_simple: should find 3 unique values")
    call assert_equal_array_real(unique_values(1:n_unique), expected_unique, n_unique, &
                                  1d-12, "test_edf_simple: unique values mismatch")
    call assert_equal_array_real(cdf_values(1:n_unique), expected_cdf, n_unique, &
                                  1d-12, "test_edf_simple: CDF values mismatch")
  end subroutine

  !> @brief Test EDF with all unique values.
  !> Each value appears exactly once, so CDF should be [1/n, 2/n, ..., n/n]
  subroutine test_edf_all_unique()
    real(real64), dimension(5) :: values = [5.0_real64, 1.0_real64, 3.0_real64, &
                                             2.0_real64, 4.0_real64]
    integer(int32) :: n_values = 5
    integer(int32) :: perm(5), stack_left(5), stack_right(5)
    real(real64) :: unique_values(5), cdf_values(5)
    real(real64) :: expected_unique(5), expected_cdf(5)
    integer(int32) :: n_unique, ierr, i
    
    expected_unique = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64]
    do i = 1, 5
      expected_cdf(i) = real(i, real64) / 5.0_real64
    end do
    
    call compute_edf(values, n_values, perm, stack_left, stack_right, &
                     unique_values, cdf_values, n_unique, ierr)
    
    call assert_equal_int(ierr, ERR_OK, "test_edf_all_unique: error code should be ERR_OK")
    call assert_equal_int(n_unique, 5, "test_edf_all_unique: should find 5 unique values")
    call assert_equal_array_real(unique_values(1:n_unique), expected_unique, n_unique, &
                                  1d-12, "test_edf_all_unique: unique values mismatch")
    call assert_equal_array_real(cdf_values(1:n_unique), expected_cdf, n_unique, &
                                  1d-12, "test_edf_all_unique: CDF values mismatch")
  end subroutine

  !> @brief Test EDF with all identical values.
  !> Should return 1 unique value with CDF = 1.0
  subroutine test_edf_all_same()
    real(real64), dimension(7) :: values = 2.5_real64
    integer(int32) :: n_values = 7
    integer(int32) :: perm(7), stack_left(7), stack_right(7)
    real(real64) :: unique_values(7), cdf_values(7)
    integer(int32) :: n_unique, ierr
    
    call compute_edf(values, n_values, perm, stack_left, stack_right, &
                     unique_values, cdf_values, n_unique, ierr)
    
    call assert_equal_int(ierr, ERR_OK, "test_edf_all_same: error code should be ERR_OK")
    call assert_equal_int(n_unique, 1, "test_edf_all_same: should find 1 unique value")
    call assert_equal_real(unique_values(1), 2.5_real64, 1d-12, "test_edf_all_same: unique value should be 2.5")
    call assert_equal_real(cdf_values(1), 1.0_real64, 1d-12, "test_edf_all_same: CDF should be 1.0")
  end subroutine

  !> @brief Test EDF with specific duplicate pattern.
  !> Input: [1.0, 1.0, 2.0, 2.0, 2.0, 3.0]
  !> Expected CDF: [2/6, 5/6, 6/6]
  subroutine test_edf_duplicates()
    real(real64), dimension(6) :: values = [1.0_real64, 1.0_real64, 2.0_real64, &
                                             2.0_real64, 2.0_real64, 3.0_real64]
    integer(int32) :: n_values = 6
    integer(int32) :: perm(6), stack_left(6), stack_right(6)
    real(real64) :: unique_values(6), cdf_values(6)
    real(real64) :: expected_cdf(3)
    integer(int32) :: n_unique, ierr
    
    expected_cdf = [2.0_real64/6.0_real64, 5.0_real64/6.0_real64, 1.0_real64]
    
    call compute_edf(values, n_values, perm, stack_left, stack_right, &
                     unique_values, cdf_values, n_unique, ierr)
    
    call assert_equal_int(ierr, ERR_OK, "test_edf_duplicates: error code should be ERR_OK")
    call assert_equal_int(n_unique, 3, "test_edf_duplicates: should find 3 unique values")
    call assert_equal_array_real(cdf_values(1:n_unique), expected_cdf, n_unique, &
                                  1d-12, "test_edf_duplicates: CDF values mismatch")
  end subroutine

  !> @brief Test EDF with a single value.
  subroutine test_edf_single_value()
    real(real64), dimension(1) :: values = [42.0_real64]
    integer(int32) :: n_values = 1
    integer(int32) :: perm(1), stack_left(1), stack_right(1)
    real(real64) :: unique_values(1), cdf_values(1)
    integer(int32) :: n_unique, ierr
    
    call compute_edf(values, n_values, perm, stack_left, stack_right, &
                     unique_values, cdf_values, n_unique, ierr)
    
    call assert_equal_int(ierr, ERR_OK, "test_edf_single_value: error code should be ERR_OK")
    call assert_equal_int(n_unique, 1, "test_edf_single_value: should find 1 unique value")
    call assert_equal_real(unique_values(1), 42.0_real64, 1d-12, "test_edf_single_value: unique value should be 42.0")
    call assert_equal_real(cdf_values(1), 1.0_real64, 1d-12, "test_edf_single_value: CDF should be 1.0")
  end subroutine

  !> @brief Test EDF with empty input (should return error).
  subroutine test_edf_empty_input()
    real(real64), dimension(5) :: values = 0.0_real64
    integer(int32) :: n_values = 0  ! Empty input
    integer(int32) :: perm(5), stack_left(5), stack_right(5)
    real(real64) :: unique_values(5), cdf_values(5)
    integer(int32) :: n_unique, ierr
    
    call compute_edf(values, n_values, perm, stack_left, stack_right, &
                     unique_values, cdf_values, n_unique, ierr)
    
    call assert_equal_int(ierr, ERR_EMPTY_INPUT, "test_edf_empty_input: should return ERR_EMPTY_INPUT")
    call assert_equal_int(n_unique, 0, "test_edf_empty_input: n_unique should be 0")
  end subroutine

  !> @brief Test EDF with a larger dataset (100 values).
  subroutine test_edf_large_dataset()
    integer(int32), parameter :: n = 100
    real(real64), dimension(n) :: values
    integer(int32) :: perm(n), stack_left(n), stack_right(n)
    real(real64) :: unique_values(n), cdf_values(n)
    integer(int32) :: n_unique, ierr, i
    
    ! Create values: [1, 1, 2, 2, 3, 3, ..., 50, 50] (each appears twice)
    do i = 1, n
      values(i) = real((i + 1) / 2, real64)
    end do
    
    call compute_edf(values, n, perm, stack_left, stack_right, &
                     unique_values, cdf_values, n_unique, ierr)
    
    call assert_equal_int(ierr, ERR_OK, "test_edf_large_dataset: error code should be ERR_OK")
    call assert_equal_int(n_unique, 50, "test_edf_large_dataset: should find 50 unique values")
    call assert_equal_real(cdf_values(n_unique), 1.0_real64, 1d-12, &
                           "test_edf_large_dataset: final CDF should be 1.0")
    call assert_equal_real(unique_values(1), 1.0_real64, 1d-12, &
                           "test_edf_large_dataset: first unique value should be 1.0")
    call assert_equal_real(unique_values(n_unique), 50.0_real64, 1d-12, &
                           "test_edf_large_dataset: last unique value should be 50.0")
  end subroutine

  !> @brief Test EDF with negative values and zero.
  subroutine test_edf_negative_values()
    real(real64), dimension(7) :: values = [-2.0_real64, 0.0_real64, -1.0_real64, &
                                             1.0_real64, -2.0_real64, 0.0_real64, 2.0_real64]
    integer(int32) :: n_values = 7
    integer(int32) :: perm(7), stack_left(7), stack_right(7)
    real(real64) :: unique_values(7), cdf_values(7)
    real(real64) :: expected_unique(5)
    integer(int32) :: n_unique, ierr
    
    expected_unique = [-2.0_real64, -1.0_real64, 0.0_real64, 1.0_real64, 2.0_real64]
    
    call compute_edf(values, n_values, perm, stack_left, stack_right, &
                     unique_values, cdf_values, n_unique, ierr)
    
    call assert_equal_int(ierr, ERR_OK, "test_edf_negative_values: error code should be ERR_OK")
    call assert_equal_int(n_unique, 5, "test_edf_negative_values: should find 5 unique values")
    call assert_equal_array_real(unique_values(1:n_unique), expected_unique, n_unique, &
                                  1d-12, "test_edf_negative_values: unique values mismatch")
    call assert_equal_real(cdf_values(n_unique), 1.0_real64, 1d-12, &
                           "test_edf_negative_values: final CDF should be 1.0")
  end subroutine

end module mod_test_loess_smoothing
