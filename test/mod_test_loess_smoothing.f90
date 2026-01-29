!> @file mod_test_loess_smoothing.f90
!> Unit test suite for LOESS smoothing from tox_loess
!> @details Unit tests for LOESS smoothing, including edge cases and input validation.
!
! IMPORTANT NOTE (API change):
! - loess_fit_plain / loess_fit_robust expect:
!     w(n)        (weights array, rank-1)
!     z(n,1)      (evaluation points, rank-2)
!     yhat(n)     (output, rank-1)
!     ierr        (output)
! - These routines do NOT support "m query points" directly. We therefore run queries by
!   placing desired query x values into the first entries of z(:,1), and filling the rest
!   with something valid (e.g., x_ref). Then we assert on yhat(1:m).

module mod_test_loess_smoothing
  use asserts
  use tox_errors, only: ERR_INVALID_INPUT
  use tox_loess,  only: loess_fit_plain, loess_fit_robust
  use, intrinsic :: iso_fortran_env, only: real64, int32
  implicit none
  public

  abstract interface
    subroutine test_interface()
    end subroutine test_interface
  end interface

  type :: test_case
    character(len=64) :: name
    procedure(test_interface), pointer, nopass :: test_proc => null()
  end type test_case

contains

  !> Get array of all available LOESS tests.
  subroutine get_all_tests(all_tests)
    type(test_case), intent(out) :: all_tests(14)

    all_tests(1)  = test_case("test_loess_constant_input",       test_loess_constant_input)
    all_tests(2)  = test_case("test_loess_linear_trend",         test_loess_linear_trend)
    all_tests(3)  = test_case("test_loess_outlier_suppression",  test_loess_outlier_suppression)
    all_tests(4)  = test_case("test_loess_sparse_like_spacing",  test_loess_sparse_like_spacing)
    all_tests(5)  = test_case("test_loess_single_point",         test_loess_single_point)
    all_tests(6)  = test_case("test_loess_identical_points",     test_loess_identical_points)
    all_tests(7)  = test_case("test_loess_linear_interp",        test_loess_linear_interp)
    all_tests(8)  = test_case("test_loess_weight_effect",        test_loess_weight_effect)
    all_tests(9)  = test_case("test_loess_edge_query",           test_loess_edge_query)
    all_tests(10) = test_case("test_loess_invalid_span",         test_loess_invalid_span)
    all_tests(11) = test_case("test_loess_invalid_indices",      test_loess_invalid_indices)
    all_tests(12) = test_case("test_loess_workspace_too_small",  test_loess_workspace_too_small)
    all_tests(13) = test_case("test_loess_degree_0_constant",    test_loess_degree_0_constant)
    all_tests(14) = test_case("test_loess_robust_same_as_plain", test_loess_robust_same_as_plain)
  end subroutine get_all_tests

  !> Run all LOESS smoothing tests.
  subroutine run_all_tests_loess_smoothing()
    type(test_case) :: all_tests(14)
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
    type(test_case) :: all_tests(14)
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

  ! ---------------------------------------------------------------------------
  ! TESTS
  ! ---------------------------------------------------------------------------

  !> Constant input: y should remain constant everywhere.
  subroutine test_loess_constant_input()
    integer(int32), parameter :: n = 100
    real(real64) :: x_ref(n), y_ref(n), w(n), z(n,1), yhat(n)
    integer(int32) :: iv(50000), ierr
    real(real64) :: wv(500000), diagl(n)

    x_ref = 5.0_real64
    y_ref = 10.0_real64
    w     = 1.0_real64
    z(:,1)= 5.0_real64

    iv = 1_int32

    call loess_fit_plain(n, x_ref, y_ref, w, z, 0.5_real64, 1_int32, n, .false., .false., &
                         iv, size(iv), wv, size(wv), diagl, yhat, ierr)

    call assert_equal_int(ierr, 0, "Constant input ierr==0")
    call assert_true(all(abs(yhat - 10.0_real64) < 1.0e-6_real64), "Constant input yhat constant")
  end subroutine test_loess_constant_input

  !> Linear trend: y = 0.5 x, should be approximated well.
  subroutine test_loess_linear_trend()
    integer(int32), parameter :: n = 100
    real(real64) :: x_ref(n), y_ref(n), w(n), z(n,1), yhat(n)
    integer(int32) :: iv(50000), ierr, i
    real(real64) :: wv(500000), diagl(n)

    do i = 1, n
      x_ref(i) = real(i, real64)
    end do
    y_ref = 0.5_real64 * x_ref
    w     = 1.0_real64
    z(:,1)= x_ref

    iv = 1_int32

    call loess_fit_plain(n, x_ref, y_ref, w, z, 0.7_real64, 1_int32, n, .false., .false., &
                         iv, size(iv), wv, size(wv), diagl, yhat, ierr)

    call assert_equal_int(ierr, 0, "Linear trend ierr==0")
    call assert_true(all(abs(yhat - y_ref) < 0.05_real64), "Linear trend approx ok")
  end subroutine test_loess_linear_trend

  !> Robust outlier suppression: one strong outlier should be down-weighted.
  subroutine test_loess_outlier_suppression()
    integer(int32), parameter :: n = 100
    real(real64) :: x_ref(n), y_ref(n), w(n), z(n,1), yhat(n)
    real(real64) :: diagl(n), rw(n), ww(n), res(n)
    integer(int32) :: iv(50000), pi(n), ierr, i
    real(real64) :: wv(500000)

    x_ref = 10.0_real64
    y_ref = 5.0_real64

    ! One outlier
    x_ref(n) = 100.0_real64
    y_ref(n) = 99.0_real64

    w     = 1.0_real64
    z(:,1)= 10.0_real64

    iv = 1_int32
    pi = [(i, i=1,n)]

    call loess_fit_robust(n, x_ref, y_ref, w, z, 0.7_real64, 1_int32, n, .false., .false., 3_int32, &
                          iv, size(iv), wv, size(wv), diagl, rw, ww, res, pi, yhat, ierr)

    call assert_equal_int(ierr, 0, "Robust outlier suppression ierr==0")
    call assert_true(all(abs(yhat - 5.0_real64) < 0.05_real64), "Outlier suppressed (near 5)")
  end subroutine test_loess_outlier_suppression

  !> Sparse-like spacing: large gaps in x should still produce finite outputs.
  subroutine test_loess_sparse_like_spacing()
    integer(int32), parameter :: n = 80
    real(real64) :: x_ref(n), y_ref(n), w(n), z(n,1), yhat(n)
    integer(int32) :: iv(50000), ierr, i
    real(real64) :: wv(500000), diagl(n)

    do i = 1, n
      x_ref(i) = 100.0_real64 * real(i, real64)
      y_ref(i) = real(i, real64)
    end do
    w     = 1.0_real64
    z(:,1)= x_ref

    iv = 1_int32

    call loess_fit_plain(n, x_ref, y_ref, w, z, 0.7_real64, 1_int32, n, .false., .false., &
                         iv, size(iv), wv, size(wv), diagl, yhat, ierr)

    call assert_equal_int(ierr, 0, "Sparse-like spacing ierr==0")
    call assert_true(all(.not. (yhat /= yhat)), "Sparse-like spacing: no NaNs")  ! NaN check: yhat==yhat
  end subroutine test_loess_sparse_like_spacing

  !> Single point: yhat should equal that y for any z (LOESS degenerates).
  subroutine test_loess_single_point()
    integer(int32), parameter :: n = 1
    real(real64) :: x_ref(n), y_ref(n), w(n), z(n,1), yhat(n)
    integer(int32) :: iv(50000), ierr
    real(real64) :: wv(500000), diagl(n)

    x_ref(1) = 0.0_real64
    y_ref(1) = 42.0_real64
    w(1)     = 1.0_real64
    z(1,1)   = 100.0_real64  ! query outside, still should return constant

    iv = 1_int32

    call loess_fit_plain(n, x_ref, y_ref, w, z, 1.0_real64, 1_int32, 1_int32, .false., .false., &
                         iv, size(iv), wv, size(wv), diagl, yhat, ierr)

    call assert_equal_int(ierr, 0, "Single point ierr==0")
    call assert_equal_real(yhat(1), 42.0_real64, 1.0e-6_real64, "Single point yhat==y")
  end subroutine test_loess_single_point

  !> Identical points: y constant, output should be constant.
  subroutine test_loess_identical_points()
    integer(int32), parameter :: n = 2
    real(real64) :: x_ref(n), y_ref(n), w(n), z(n,1), yhat(n)
    integer(int32) :: iv(50000), ierr
    real(real64) :: wv(500000), diagl(n)

    x_ref = [0.0_real64, 0.0_real64]
    y_ref = [1.0_real64, 1.0_real64]
    w     = 1.0_real64
    z(:,1)= [0.0_real64, 0.0_real64]

    iv = 1_int32

    call loess_fit_plain(n, x_ref, y_ref, w, z, 1.0_real64, 1_int32, 2_int32, .false., .false., &
                         iv, size(iv), wv, size(wv), diagl, yhat, ierr)

    call assert_equal_int(ierr, 0, "Identical points ierr==0")
    call assert_true(all(abs(yhat - 1.0_real64) < 1.0e-6_real64), "Identical points yhat constant")
  end subroutine test_loess_identical_points

  !> Linear interpolation-ish sanity: with two points (0,0) and (2,2), query at 1 should be near 1.
  subroutine test_loess_linear_interp()
    integer(int32), parameter :: n = 2
    real(real64) :: x_ref(n), y_ref(n), w(n), z(n,1), yhat(n)
    integer(int32) :: iv(50000), ierr
    real(real64) :: wv(500000), diagl(n)

    x_ref = [0.0_real64, 2.0_real64]
    y_ref = [0.0_real64, 2.0_real64]
    w     = 1.0_real64

    ! Put query at z(1)=1 and fill z(2) with something valid
    z(1,1) = 1.0_real64
    z(2,1) = 1.0_real64

    iv = 1_int32

    call loess_fit_plain(n, x_ref, y_ref, w, z, 1.0_real64, 1_int32, 2_int32, .false., .false., &
                         iv, size(iv), wv, size(wv), diagl, yhat, ierr)

    call assert_equal_int(ierr, 0, "Linear interp ierr==0")
    call assert_true(abs(yhat(1) - 1.0_real64) < 0.25_real64, "Linear interp yhat(1) near 1")
  end subroutine test_loess_linear_interp

  !> Weight effect: down-weight the second point strongly; near x=0 we should be closer to y=0.
  subroutine test_loess_weight_effect()
    integer(int32), parameter :: n = 2
    real(real64) :: x_ref(n), y_ref(n), w(n), z(n,1), yhat(n)
    integer(int32) :: iv(50000), ierr
    real(real64) :: wv(500000), diagl(n)

    x_ref = [0.0_real64, 10.0_real64]
    y_ref = [0.0_real64, 10.0_real64]
    w     = [1.0_real64, 1.0e-6_real64]   ! almost ignore second point

    z(:,1) = [0.0_real64, 0.0_real64]

    iv = 1_int32

    call loess_fit_plain(n, x_ref, y_ref, w, z, 1.0_real64, 1_int32, 2_int32, .false., .false., &
                         iv, size(iv), wv, size(wv), diagl, yhat, ierr)

    call assert_equal_int(ierr, 0, "Weight effect ierr==0")
    call assert_true(yhat(1) < 1.0_real64, "Weight effect: yhat near 0 when querying at 0")
  end subroutine test_loess_weight_effect

  !> Edge queries: test z outside min/max of x. Expect near boundary values (loess extrapolation is limited).
  subroutine test_loess_edge_query()
    integer(int32), parameter :: n = 5
    real(real64) :: x_ref(n), y_ref(n), w(n), z(n,1), yhat(n)
    integer(int32) :: iv(50000), ierr
    real(real64) :: wv(500000), diagl(n)

    x_ref = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64]
    y_ref = [10.0_real64, 20.0_real64, 30.0_real64, 40.0_real64, 50.0_real64]
    w     = 1.0_real64

    ! Queries in first two rows, fill rest with something valid
    z(1,1) = 0.0_real64
    z(2,1) = 6.0_real64
    z(3,1) = 3.0_real64
    z(4,1) = 3.0_real64
    z(5,1) = 3.0_real64

    iv = 1_int32

    call loess_fit_plain(n, x_ref, y_ref, w, z, 0.7_real64, 1_int32, 5_int32, .false., .false., &
                         iv, size(iv), wv, size(wv), diagl, yhat, ierr)

    call assert_equal_int(ierr, 0, "Edge query ierr==0")
    call assert_true(abs(yhat(1) - 10.0_real64) < 5.0_real64, "Edge low query near y(min)")
    call assert_true(abs(yhat(2) - 50.0_real64) < 5.0_real64, "Edge high query near y(max)")
  end subroutine test_loess_edge_query

  !> Invalid span: negative span should trigger ERR_INVALID_INPUT.
  subroutine test_loess_invalid_span()
    integer(int32), parameter :: n = 5
    real(real64) :: x_ref(n), y_ref(n), w(n), z(n,1), yhat(n)
    integer(int32) :: iv(50000), ierr
    real(real64) :: wv(500000), diagl(n)

    x_ref = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64]
    y_ref = [10.0_real64, 20.0_real64, 30.0_real64, 40.0_real64, 50.0_real64]
    w     = 1.0_real64
    z(:,1)= x_ref

    iv = 1_int32

    call loess_fit_plain(n, x_ref, y_ref, w, z, -1.0_real64, 1_int32, 5_int32, .false., .false., &
                         iv, size(iv), wv, size(wv), diagl, yhat, ierr)

    call assert_equal_int(ierr, ERR_INVALID_INPUT, "Invalid span -> ERR_INVALID_INPUT")
  end subroutine test_loess_invalid_span

  !> Invalid iv indices: iv contains values outside [1..n] should trigger ERR_INVALID_INPUT.
  subroutine test_loess_invalid_indices()
    integer(int32), parameter :: n = 5
    real(real64) :: x_ref(n), y_ref(n), w(n), z(n,1), yhat(n)
    integer(int32) :: iv(50000), ierr
    real(real64) :: wv(500000), diagl(n)

    x_ref = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64]
    y_ref = [10.0_real64, 20.0_real64, 30.0_real64, 40.0_real64, 50.0_real64]
    w     = 1.0_real64
    z(:,1)= x_ref

    iv = 1_int32
    iv(1) = 0_int32   ! invalid (too low)

    call loess_fit_plain(n, x_ref, y_ref, w, z, 0.7_real64, 1_int32, 5_int32, .false., .false., &
                         iv, size(iv), wv, size(wv), diagl, yhat, ierr)

    call assert_equal_int(ierr, ERR_INVALID_INPUT, "Invalid iv indices -> ERR_INVALID_INPUT")
  end subroutine test_loess_invalid_indices

  !> Workspace too small: should trigger ERR_INVALID_INPUT.
  subroutine test_loess_workspace_too_small()
    integer(int32), parameter :: n = 5
    real(real64) :: x_ref(n), y_ref(n), w(n), z(n,1), yhat(n)
    integer(int32) :: iv_small(10), ierr
    real(real64) :: wv_small(10), diagl(n)

    x_ref = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64]
    y_ref = [10.0_real64, 20.0_real64, 30.0_real64, 40.0_real64, 50.0_real64]
    w     = 1.0_real64
    z(:,1)= x_ref

    iv_small = 1_int32

    call loess_fit_plain(n, x_ref, y_ref, w, z, 0.7_real64, 1_int32, 5_int32, .false., .false., &
                         iv_small, size(iv_small), wv_small, size(wv_small), diagl, yhat, ierr)

    call assert_equal_int(ierr, ERR_INVALID_INPUT, "Workspace too small -> ERR_INVALID_INPUT")
  end subroutine test_loess_workspace_too_small

  !> Degree 0 with constant y: should still return constant.
  subroutine test_loess_degree_0_constant()
    integer(int32), parameter :: n = 50
    real(real64) :: x_ref(n), y_ref(n), w(n), z(n,1), yhat(n)
    integer(int32) :: iv(50000), ierr, i
    real(real64) :: wv(500000), diagl(n)

    do i = 1, n
      x_ref(i) = real(i, real64)
    end do
    y_ref = 7.0_real64
    w     = 1.0_real64
    z(:,1)= x_ref

    iv = 1_int32

    call loess_fit_plain(n, x_ref, y_ref, w, z, 0.7_real64, 0_int32, n, .false., .false., &
                         iv, size(iv), wv, size(wv), diagl, yhat, ierr)

    call assert_equal_int(ierr, 0, "Degree 0 constant ierr==0")
    call assert_true(all(abs(yhat - 7.0_real64) < 1.0e-6_real64), "Degree 0 constant preserved")
  end subroutine test_loess_degree_0_constant

  !> Robust equals plain in a clean dataset (no outliers): results should be close.
  subroutine test_loess_robust_same_as_plain()
    integer(int32), parameter :: n = 60
    real(real64) :: x_ref(n), y_ref(n), w(n), z(n,1), yhat_plain(n), yhat_robust(n)
    real(real64) :: diagl(n), rw(n), ww(n), res(n)
    integer(int32) :: iv(50000), pi(n), ierr1, ierr2, i
    real(real64) :: wv(500000)

    do i = 1, n
      x_ref(i) = real(i, real64)
    end do
    y_ref = 2.0_real64 * x_ref + 1.0_real64
    w     = 1.0_real64
    z(:,1)= x_ref

    iv = 1_int32
    pi = [(i, i=1,n)]

    call loess_fit_plain(n, x_ref, y_ref, w, z, 0.7_real64, 1_int32, n, .false., .false., &
                         iv, size(iv), wv, size(wv), diagl, yhat_plain, ierr1)

    call loess_fit_robust(n, x_ref, y_ref, w, z, 0.7_real64, 1_int32, n, .false., .false., 3_int32, &
                          iv, size(iv), wv, size(wv), diagl, rw, ww, res, pi, yhat_robust, ierr2)

    call assert_equal_int(ierr1, 0, "Plain ierr==0")
    call assert_equal_int(ierr2, 0, "Robust ierr==0")
    call assert_true(all(abs(yhat_plain - yhat_robust) < 0.1_real64), "Robust ~ Plain (clean data)")
  end subroutine test_loess_robust_same_as_plain

end module mod_test_loess_smoothing
