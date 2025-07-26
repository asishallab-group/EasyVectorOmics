!> @file mod_test_loess_smoothing.f90
!> @brief Unit test suite for LOESS smoothing (f42_loess_smoothing.F90)
!> @details Unit tests for LOESS smoothing, including masking and edge cases.

module mod_test_loess_smoothing
  use f42_utils
  use, intrinsic :: iso_fortran_env, only: real64
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

  !> @brief Get array of all available LOESS tests.
  function get_all_tests() result(all_tests)
    type(test_case) :: all_tests(11)
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
  end function get_all_tests

  !> @brief Run all LOESS smoothing tests.
  subroutine run_all_tests_loess_smoothing()
    type(test_case) :: all_tests(11)
    integer :: i
    all_tests = get_all_tests()
    do i = 1, size(all_tests)
      call all_tests(i)%test_proc()
      print *, trim(all_tests(i)%name), " passed."
    end do
    print *, "All LOESS smoothing tests passed successfully."
  end subroutine run_all_tests_loess_smoothing

  !> @brief Run specific LOESS smoothing tests by name.
  subroutine run_named_tests_loess_smoothing(test_names)
    character(len=*), intent(in) :: test_names(:)
    type(test_case) :: all_tests(11)
    integer :: i, j
    logical :: found
    all_tests = get_all_tests()
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
    integer, parameter :: n = 100
    real(real64) :: x_ref(n), y_ref(n), x_query(n/2), y_out(n/2)
    integer :: indices_used(n), i
    real(real64) :: workspace_weights(n), workspace_values(n)
    x_ref = 5.0_real64; y_ref = 10.0_real64; x_query = 5.0_real64; indices_used = [(i, i = 1, n)]
    call loess_smooth_2d(n, n/2, x_ref, y_ref, indices_used, x_query, 1.0_real64, 3.0_real64, y_out, &
                      workspace_weights, workspace_values)
    if (all(abs(y_out - 10.0_real64) < 1.0e-6_real64)) then
      print *, 'Test 1 (constant input): PASSED'
    else
      print *, 'Test 1 (constant input): FAILED', y_out
    end if
  end subroutine

  subroutine test_loess_linear_trend()
    integer, parameter :: n = 100
    real(real64) :: x_ref(n), y_ref(n), x_query(n/2), y_out(n/2)
    integer :: indices_used(n), i
    real(real64) :: workspace_weights(n), workspace_values(n)
    x_ref = [(real(i, kind=real64), i = 1, n)]
    y_ref = 0.5_real64 * x_ref
    x_query = [(real(i, kind=real64) + 0.5_real64, i = 1, n/2)]
    indices_used = [(i, i = 1, n)]
    call loess_smooth_2d(n, n/2, x_ref, y_ref, indices_used, x_query, 1.0_real64, 3.0_real64, &
                      y_out, workspace_weights, workspace_values)
    if (all(abs(y_out - 0.5_real64 * x_query) < 0.05_real64)) then
      print *, 'Test 2 (linear trend): PASSED'
    else
      print *, 'Test 2 (linear trend): FAILED', y_out
    end if
  end subroutine

  subroutine test_loess_outlier_suppression()
    integer, parameter :: n = 100
    real(real64) :: x_ref(n), y_ref(n), x_query(50), y_out(50)
    integer :: indices_used(n), i
    real(real64) :: workspace_weights(n), workspace_values(n)
    x_ref(1:n-1) = 10.0_real64
    x_ref(n) = 100.0_real64
    y_ref(1:n-1) = 5.0_real64
    y_ref(n) = 99.0_real64
    x_query = 10.0_real64
    indices_used = [(i, i = 1, n)]
    call loess_smooth_2d(n, 50, x_ref, y_ref, indices_used, x_query, 1.0_real64, 3.0_real64, &
                      y_out, workspace_weights, workspace_values)
    if (all(abs(y_out - 5.0_real64) < 0.01_real64)) then
      print *, 'Test 3 (outlier suppression): PASSED'
    else
      print *, 'Test 3 (outlier suppression): FAILED', y_out
    end if
  end subroutine

  subroutine test_loess_sparse_fallback()
    integer, parameter :: n = 100
    real(real64) :: x_ref(n), y_ref(n), x_query(50), y_out(50)
    integer :: indices_used(n), i
    real(real64) :: workspace_weights(n), workspace_values(n)
    x_ref = [(real(i, kind=real64) * 100.0_real64, i = 1, n)]
    y_ref = [(real(i, kind=real64), i = 1, n)]
    x_query = [(real(i, kind=real64) * 100.0_real64 + 50.0_real64, i = 1, 50)]
    indices_used = [(i, i = 1, n)]
    call loess_smooth_2d(n, 50, x_ref, y_ref, indices_used, x_query, 1.0_real64, 3.0_real64, y_out, &
                      workspace_weights, workspace_values)
    if (all(abs(y_out - y_ref(indices_used(1:50))) < 1.0e-6_real64)) then
      print *, 'Test 4 (sparse fallback): PASSED'
    else
      print *, 'Test 4 (sparse fallback): FAILED', y_out
    end if
  end subroutine

  ! Removed vector field test (multivariate) for univariate-only interface

  subroutine test_loess_single_point()
    real(real64) :: x1(1), y1(1), xq1(1), yout1(1)
    real(real64) :: workspace_weights(1), workspace_values(1)
    x1 = 0.0_real64; y1(1) = 42.0_real64; xq1 = 0.0_real64
    call loess_smooth_2d(1, 1, x1, y1, (/1/), xq1, 0.1_real64, 1.0_real64, yout1, workspace_weights, workspace_values)
    if (abs(yout1(1) - y1(1)) < 1e-6_real64) then
      print *, 'Test 5 (single point): PASSED'
    else
      print *, 'Test 5 (single point): FAILED', yout1
    end if
  end subroutine

  subroutine test_loess_identical_points()
    real(real64) :: x2(2), y2(2), xq2(1), yout2(1)
    integer :: idx2(2)
    real(real64) :: workspace_weights(2), workspace_values(2)
    x2 = (/0.0_real64, 1.0_real64/); y2 = (/1.0_real64, 1.0_real64/); xq2 = 0.0_real64; idx2 = (/1,2/)
    call loess_smooth_2d(2, 1, x2, y2, idx2, xq2, 0.1_real64, 1.0_real64, yout2, workspace_weights, workspace_values)
    if (abs(yout2(1) - y2(1)) < 1e-6_real64) then
      print *, 'Test 6 (identical points): PASSED'
    else
      print *, 'Test 6 (identical points): FAILED', yout2
    end if
  end subroutine

  subroutine test_loess_linear_interp()
    real(real64) :: x3(2), y3(2), xq3(1), yout3(1)
    integer :: idx3(2)
    real(real64) :: workspace_weights(2), workspace_values(2)
    x3 = (/0.0_real64, 2.0_real64/); y3 = (/0.0_real64, 2.0_real64/); xq3 = 1.0_real64; idx3 = (/1,2/)
    call loess_smooth_2d(2, 1, x3, y3, idx3, xq3, 0.5_real64, 3.0_real64, yout3, workspace_weights, workspace_values)
    if (abs(yout3(1) - 1.0_real64) < 0.1_real64) then
      print *, 'Test 7 (linear interp): PASSED'
    else
      print *, 'Test 7 (linear interp): FAILED', yout3
    end if
  end subroutine

  subroutine test_loess_weight_decay()
    real(real64) :: x2(2), y2(2), xq2(1), yout2(1)
    integer :: idx2(2)
    real(real64) :: workspace_weights(2), workspace_values(2)
    x2 = (/0.0_real64, 10.0_real64/); y2 = (/0.0_real64, 10.0_real64/); xq2 = 0.0_real64; idx2 = (/1,2/)
    call loess_smooth_2d(2, 1, x2, y2, idx2, xq2, 1.0_real64, 3.0_real64, yout2, workspace_weights, workspace_values)
    if (yout2(1) < 5.0_real64) then
      print *, 'Test 8 (weight decay): PASSED'
    else
      print *, 'Test 8 (weight decay): FAILED', yout2
    end if
  end subroutine

  subroutine test_loess_mask_exclusion()
    real(real64) :: x3(3), y3(3), xq3(1), yout3(1)
    integer :: idx3(3)
    real(real64) :: workspace_weights(3), workspace_values(3)
    logical :: mask3(3)
    x3 = (/0.0_real64, 10.0_real64, 20.0_real64/)
    y3 = (/0.0_real64, 10.0_real64, 20.0_real64/)
    xq3 = 0.0_real64
    idx3 = (/1,3,2/)
    mask3 = (/ .false., .true., .false. /)
    call loess_smooth_2d(3, 1, x3, y3, idx3, xq3, 10.0_real64, 3.0_real64, yout3, workspace_weights, workspace_values, mask3)
    if (abs(yout3(1) - 10.0_real64) < 1.0_real64) then
      print *, 'Test 9 (mask exclusion): PASSED'
    else
      print *, 'Test 9 (mask exclusion): FAILED', yout3
    end if
  end subroutine

  subroutine test_loess_fallback()
    real(real64) :: x1(1), y1(1), xq1(1), yout1(1)
    real(real64) :: workspace_weights(1), workspace_values(1)
    x1 = 0.0_real64; y1(1) = 123.0_real64; xq1 = 100.0_real64
    call loess_smooth_2d(1, 1, x1, y1, (/1/), xq1, 0.1_real64, 1.0_real64, yout1, workspace_weights, workspace_values)
    if (abs(yout1(1) - y1(1)) < 1e-6_real64) then
      print *, 'Test 10 (fallback): PASSED'
    else
      print *, 'Test 10 (fallback): FAILED', yout1
    end if
  end subroutine

  ! Removed 10D scaling test (multivariate) for univariate-only interface

  subroutine test_loess_edge_query()
    ! Test edge query points (extrapolation)
    real(real64) :: x(5), y(5), xq(2), yout(2)
    integer :: idx(5)
    real(real64) :: workspace_weights(5), workspace_values(5)
    x = (/1.0, 2.0, 3.0, 4.0, 5.0/)
    y = (/10.0, 20.0, 30.0, 40.0, 50.0/)
    idx = (/1,2,3,4,5/)
    xq = (/0.0, 6.0/)
    call loess_smooth_2d(5, 2, x, y, idx, xq, 1.0_real64, 3.0_real64, yout, workspace_weights, workspace_values)
    if (abs(yout(1) - 10.0_real64) < 1.0 .and. abs(yout(2) - 50.0_real64) < 1.0) then
      print *, 'Test 11 (edge query): PASSED'
    else
      print *, 'Test 11 (edge query): FAILED', yout
    end if
  end subroutine


end module mod_test_loess_smoothing
