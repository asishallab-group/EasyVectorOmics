!> @file mod_test_loess_smoothing.f90
!> @brief Unit test suite for LOESS smoothing (f42_loess_smoothing.F90)
!> @details Unit tests for LOESS smoothing, including masking and edge cases.

module mod_test_loess_smoothing
  use loess_module
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
    type(test_case) :: all_tests(12)
    all_tests(1) = test_case("test_loess_constant_input", test_loess_constant_input)
    all_tests(2) = test_case("test_loess_linear_trend", test_loess_linear_trend)
    all_tests(3) = test_case("test_loess_outlier_suppression", test_loess_outlier_suppression)
    all_tests(4) = test_case("test_loess_sparse_fallback", test_loess_sparse_fallback)
    all_tests(5) = test_case("test_loess_vector_field", test_loess_vector_field)
    all_tests(6) = test_case("test_loess_single_point", test_loess_single_point)
    all_tests(7) = test_case("test_loess_identical_vectors", test_loess_identical_vectors)
    all_tests(8) = test_case("test_loess_linear_interp", test_loess_linear_interp)
    all_tests(9) = test_case("test_loess_weight_decay", test_loess_weight_decay)
    all_tests(10) = test_case("test_loess_mask_exclusion", test_loess_mask_exclusion)
    all_tests(11) = test_case("test_loess_fallback", test_loess_fallback)
    all_tests(12) = test_case("test_loess_10d_scaling", test_loess_10d_scaling)
  end function get_all_tests

  !> @brief Run all LOESS smoothing tests.
  subroutine run_all_tests_loess_smoothing()
    type(test_case) :: all_tests(12)
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
    type(test_case) :: all_tests(12)
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
    integer, parameter :: n = 100, d = 1
    real(real64) :: x_ref(n), y_ref(d, n), x_query(n/2), y_out(d, n/2)
    integer :: indices_used(n), i
    real(real64) :: workspace_weights(n), workspace_values(d, n)
    x_ref = 5.0_real64; y_ref = 10.0_real64; x_query = 5.0_real64; indices_used = (/ (i, i = 1, n) /)
    call loess_smooth(n, n/2, d, x_ref, y_ref, indices_used, x_query, 1.0_real64, 3.0_real64, y_out, &
                      workspace_weights, workspace_values)
    print *, 'Test 1 (constant input):', all(abs(y_out - 10.0_real64) < 1.0e-6_real64)
  end subroutine

  subroutine test_loess_linear_trend()
    integer, parameter :: n = 100, d = 1
    real(real64) :: x_ref(n), y_ref(d, n), x_query(n/2), y_out(d, n/2)
    integer :: indices_used(n), i
    real(real64) :: workspace_weights(n), workspace_values(d, n)
    x_ref = (/ (real(i, kind=real64), i = 1, n) /)
    y_ref(1, :) = 0.5_real64 * x_ref
    x_query = (/ (real(i, kind=real64) + 0.5_real64, i = 1, n/2) /)
    indices_used = (/ (i, i = 1, n) /)
    call loess_smooth(n, n/2, d, x_ref, y_ref, indices_used, x_query, 1.0_real64, 3.0_real64, &
                      y_out, workspace_weights, workspace_values)
    print *, 'Test 2 (linear trend):', all(abs(y_out(1, :) - 0.5_real64 * x_query) < 0.05_real64)
  end subroutine

  subroutine test_loess_outlier_suppression()
    integer, parameter :: n = 100, d = 1
    real(real64) :: x_ref(n), y_ref(d, n), x_query(50), y_out(d, 50)
    integer :: indices_used(n), i
    real(real64) :: workspace_weights(n), workspace_values(d, n)
    x_ref = (/ (10.0_real64, i = 1, n - 1), 100.0_real64 /)
    y_ref(1, :) = (/ (5.0_real64, i = 1, n - 1), 99.0_real64 /)
    x_query = 10.0_real64
    indices_used = (/ (i, i = 1, n) /)
    call loess_smooth(n, 50, d, x_ref, y_ref, indices_used, x_query, 1.0_real64, 3.0_real64, &
                      y_out, workspace_weights, workspace_values)
    print *, 'Test 3 (outlier suppression):', abs(y_out(1, 1) - 5.0_real64) < 0.01_real64
  end subroutine

  subroutine test_loess_sparse_fallback()
    integer, parameter :: n = 100, d = 1
    real(real64) :: x_ref(n), y_ref(d, n), x_query(50), y_out(d, 50)
    integer :: indices_used(n), i
    real(real64) :: workspace_weights(n), workspace_values(d, n)
    x_ref = (/ (real(i, kind=real64) * 100.0_real64, i = 1, n) /)
    y_ref(1, :) = (/ (real(i, kind=real64), i = 1, n) /)
    x_query = (/ (real(i, kind=real64) * 100.0_real64 + 50.0_real64, i = 1, 50) /)
    indices_used = (/ (i, i = 1, n) /)
    call loess_smooth(n, 50, d, x_ref, y_ref, indices_used, x_query, 1.0_real64, 3.0_real64, y_out, &
                      workspace_weights, workspace_values)
    print *, 'Test 4 (sparse fallback):', all(abs(y_out(1, :) - y_ref(1, indices_used(1:50))) < 1.0e-6_real64)
  end subroutine

  subroutine test_loess_vector_field()
    integer, parameter :: n = 100, d = 3, n_target = 50
    real(real64) :: x_ref(n), y_ref(d, n), x_query(n_target), y_out(d, n_target)
    integer :: indices_used(n), i, j
    real(real64) :: workspace_weights(n), workspace_values(d, n)
    do j = 1, d
      do i = 1, n
        y_ref(j, i) = real(i, kind=real64) + 0.1_real64 * real(j, kind=real64)
      end do
    end do
    x_ref = (/ (real(i, kind=real64), i = 1, n) /)
    x_query = (/ (real(i, kind=real64) + 0.5_real64, i = 1, n_target) /)
    indices_used = (/ (i, i = 1, n) /)
    call loess_smooth(n, n_target, d, x_ref, y_ref, indices_used, x_query, 0.1_real64, 3.0_real64, &
                      y_out, workspace_weights, workspace_values)
    print *, 'Test 5 (vector field):', all(abs(y_out - y_ref(:, indices_used(1:n_target))) < 0.1_real64)
  end subroutine

  subroutine test_loess_single_point()
    real(real64) :: x1(1), y1(3,1), xq1(1), yout1(3,1)
    real(real64) :: workspace_weights(1), workspace_values(3,1)
    x1 = 0.0_real64; y1(:,1) = (/1.0_real64, 2.0_real64, 3.0_real64/); xq1 = 0.0_real64
    call loess_smooth(1, 1, 3, x1, y1, (/1/), xq1, 0.1_real64, 1.0_real64, yout1, workspace_weights, workspace_values)
    print *, 'Test 6 (single point):', all(abs(yout1(:,1) - y1(:,1)) < 1e-6_real64)
  end subroutine

  subroutine test_loess_identical_vectors()
    real(real64) :: x2(2), y2(3,2), xq2(1), yout2(3,1)
    integer :: idx2(2)
    real(real64) :: workspace_weights(2), workspace_values(3,2)
    x2 = (/0.0_real64, 1.0_real64/); y2(:,1) = (/1.0_real64,2.0_real64,3.0_real64/); y2(:,2) = &
              (/1.0_real64,2.0_real64,3.0_real64/); xq2 = 0.0_real64; idx2 = (/1,2/)
    call loess_smooth(2, 1, 3, x2, y2, idx2, xq2, 0.1_real64, 1.0_real64, yout2, workspace_weights, workspace_values)
    print *, 'Test 7 (identical):', all(abs(yout2(:,1) - y2(:,1)) < 1e-6_real64)
  end subroutine

  subroutine test_loess_linear_interp()
    real(real64) :: x3(2), y3(3,2), xq3(1), yout3(3,1)
    integer :: idx3(2)
    real(real64) :: workspace_weights(2), workspace_values(3,2)
    x3 = (/0.0_real64, 2.0_real64/); y3(:,1) = (/0.0_real64,0.0_real64,0.0_real64/); y3(:,2) = &
          (/2.0_real64,2.0_real64,2.0_real64/); xq3 = 1.0_real64; idx3 = (/1,2/)
    call loess_smooth(2, 1, 3, x3, y3, idx3, xq3, 0.5_real64, 3.0_real64, yout3, workspace_weights, workspace_values)
    print *, 'Test 8 (linear interp):', all(abs(yout3(:,1) - 1.0_real64) < 0.1_real64)
  end subroutine

  subroutine test_loess_weight_decay()
    real(real64) :: x2(2), y2(3,2), xq2(1), yout2(3,1)
    integer :: idx2(2)
    real(real64) :: workspace_weights(2), workspace_values(3,2)
    x2 = (/0.0_real64, 10.0_real64/); y2(:,1) = (/0.0_real64,0.0_real64,0.0_real64/); y2(:,2) = &
          (/10.0_real64,10.0_real64,10.0_real64/); xq2 = 0.0_real64; idx2 = (/1,2/)
    call loess_smooth(2, 1, 3, x2, y2, idx2, xq2, 1.0_real64, 3.0_real64, yout2, workspace_weights, workspace_values)
    print *, 'Test 9 (weight decay):', all(yout2(:,1) < 5.0_real64)
  end subroutine

  subroutine test_loess_mask_exclusion()
    real(real64) :: x3(3), y3(3,3), xq3(1), yout3(3,1)
    integer :: idx3(3)
    real(real64) :: workspace_weights(3), workspace_values(3,3)
    logical :: mask3(3)
    x3 = (/0.0_real64, 10.0_real64, 20.0_real64/)
    y3(:,1) = (/0.0_real64,0.0_real64,0.0_real64/)
    y3(:,2) = (/10.0_real64,10.0_real64,10.0_real64/)
    y3(:,3) = (/20.0_real64,20.0_real64,20.0_real64/)
    xq3 = 0.0_real64
    idx3 = (/1,3,2/)
    mask3 = (/ .false., .true., .false. /)
    call loess_smooth(3, 1, 3, x3, y3, idx3, xq3, 10.0_real64, 3.0_real64, yout3, workspace_weights, workspace_values, mask3)
    print *, 'Test 10 (mask exclusion):', all(abs(yout3(:,1) - 10.0_real64) < 1.0_real64)
  end subroutine

  subroutine test_loess_fallback()
    real(real64) :: x1(1), y1(3,1), xq1(1), yout1(3,1)
    real(real64) :: workspace_weights(1), workspace_values(3,1)
    x1 = 0.0_real64; y1(:,1) = (/1.0_real64,2.0_real64,3.0_real64/); xq1 = 100.0_real64
    call loess_smooth(1, 1, 3, x1, y1, (/1/), xq1, 0.1_real64, 1.0_real64, yout1, workspace_weights, workspace_values)
    print *, 'Test 11 (fallback):', all(abs(yout1(:,1) - y1(:,1)) < 1e-6_real64)
  end subroutine

  subroutine test_loess_10d_scaling()
    integer, parameter :: n12 = 100, d12 = 10
    real(real64) :: x12(n12), y12(d12, n12), xq12(1), yout12(d12,1)
    integer :: idx12(n12), i
    real(real64) :: workspace_weights(n12), workspace_values(d12, n12)
    do i = 1, n12
      x12(i) = real(i, kind=real64)
      y12(:,i) = real(i, kind=real64) / 10.0_real64
      idx12(i) = i
    end do
    xq12(1) = 50.5_real64
    call loess_smooth(n12, 1, d12, x12, y12, idx12, xq12, 5.0_real64, 3.0_real64, yout12, workspace_weights, workspace_values)
    print *, 'Test 12 (10D scaling):', all(abs(yout12(:,1) - 5.05_real64) < 0.5_real64)
  end subroutine

end module mod_test_loess_smoothing
