!> @brief Unit test suite for tissue versatility routines.
module mod_test_tissue_versatility
  use asserts
  use avmod
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

  !> @brief Get array of all available tests.
  function get_all_tests() result(all_tests)
    type(test_case) :: all_tests(11)
    all_tests(1) = test_case("test_uniform_expression", test_uniform_expression)
    all_tests(2) = test_case("test_single_axis_expression", test_single_axis_expression)
    all_tests(3) = test_case("test_null_vector", test_null_vector)
    all_tests(4) = test_case("test_partial_axis_selection", test_partial_axis_selection)
    all_tests(5) = test_case("test_mixed_vectors", test_mixed_vectors)
    all_tests(6) = test_case("test_angle_degrees", test_angle_degrees)
    all_tests(7) = test_case("test_multiple_vectors_selection", test_multiple_vectors_selection)
    all_tests(8) = test_case("test_high_dimensional_vectors", test_high_dimensional_vectors)
    all_tests(9) = test_case("test_randomized_vectors_axes", test_randomized_vectors_axes)
    all_tests(10) = test_case("test_numerical_stability", test_numerical_stability)
    all_tests(11) = test_case("test_invalid_input_no_axes", test_invalid_input_no_axes)
  end function get_all_tests

  !> @brief Run all tissue versatility tests.
  subroutine run_all_tests_tissue_versatility()
    type(test_case) :: all_tests(11)
    integer(int32) :: i
    all_tests = get_all_tests()
    do i = 1, size(all_tests)
      call all_tests(i)%test_proc()
      print *, trim(all_tests(i)%name), " passed."
    end do
    print *, "All tissue versatility tests passed successfully."
  end subroutine run_all_tests_tissue_versatility

  !> @brief Run specific tissue versatility tests by name.
  subroutine run_named_tests_tissue_versatility(test_names)
    character(len=*), intent(in) :: test_names(:)
    type(test_case) :: all_tests(11)
    integer(int32) :: i, j
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
  end subroutine run_named_tests_tissue_versatility


  !> @brief Test tissue versatility in higher dimensions (4D, 5D).
  subroutine test_high_dimensional_vectors()
    real(real64) :: expr4(4,1), tv4(1), angle4(1)
    real(real64) :: expr5(5,1), tv5(1), angle5(1)
    logical :: select_vec(1), select_axes4(4), select_axes5(5)
    expr4(:,1) = [1.0_real64, 1.0_real64, 1.0_real64, 1.0_real64]
    expr5(:,1) = [2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64]
    select_vec = [.true.]
    select_axes4 = [.true., .true., .true., .true.]
    select_axes5 = [.true., .true., .true., .true., .true.]
    call compute_tissue_versatility(4, 1, expr4, select_vec, 1, select_axes4, 4, tv4, angle4)
    call compute_tissue_versatility(5, 1, expr5, select_vec, 1, select_axes5, 5, tv5, angle5)
    call assert_equal_real(tv4(1), 0.0_real64, 1e-12_real64, "4D uniform TV")
    call assert_equal_real(angle4(1), 0.0_real64, 1e-12_real64, "4D uniform angle")
    call assert_equal_real(tv5(1), 0.0_real64, 1e-12_real64, "5D uniform TV")
    call assert_equal_real(angle5(1), 0.0_real64, 1e-12_real64, "5D uniform angle")
  end subroutine test_high_dimensional_vectors

  !> @brief Test tissue versatility with randomized vectors and axis selections.
  subroutine test_randomized_vectors_axes()
    integer(int32), parameter :: n_axes = 5, n_vecs = 4
    real(real64) :: expr(n_axes, n_vecs), tv(n_vecs), angle(n_vecs)
    logical :: select_vec(n_vecs), select_axes(n_axes)
    integer(int32) :: i
    call random_seed()
    call random_number(expr)
    select_vec = [.true., .true., .true., .true.]
    select_axes = [.true., .false., .true., .false., .true.]
    call compute_tissue_versatility(n_axes, n_vecs, expr, select_vec, 4, select_axes, 3, tv, angle)
    do i = 1, n_vecs
      call assert_true(tv(i) >= 0.0_real64 .and. tv(i) <= 1.0_real64, "Randomized TV in [0,1]")
      call assert_true(angle(i) >= 0.0_real64 .and. angle(i) <= 90.0_real64, "Randomized angle in [0,90]")
    end do
  end subroutine test_randomized_vectors_axes

  !> @brief Test numerical stability with very large and very small values.
  subroutine test_numerical_stability()
    real(real64) :: expr(3,2), tv(2), angle(2)
    logical :: select_vec(2), select_axes(3)
    expr(:,1) = [1e-15_real64, 1e-15_real64, 1e-15_real64]
    expr(:,2) = [1e15_real64, 1e15_real64, 1e15_real64]
    select_vec = [.true., .true.]
    select_axes = [.true., .true., .true.]
    call compute_tissue_versatility(3, 2, expr, select_vec, 2, select_axes, 3, tv, angle)
    call assert_equal_real(tv(1), 0.0_real64, 1e-12_real64, "Small numbers TV")
    call assert_equal_real(angle(1), 0.0_real64, 1e-12_real64, "Small numbers angle")
    call assert_equal_real(tv(2), 0.0_real64, 1e-12_real64, "Large numbers TV")
    call assert_equal_real(angle(2), 0.0_real64, 1e-12_real64, "Large numbers angle")
  end subroutine test_numerical_stability

  !> @brief Test invalid input: no axes selected (should return error indicator).
  subroutine test_invalid_input_no_axes()
    real(real64) :: expr(3,1), tv(1), angle(1)
    logical :: select_vec(1), select_axes(3)
    expr(:,1) = [1.0_real64, 2.0_real64, 3.0_real64]
    select_vec = [.true.]
    select_axes = [.false., .false., .false.]
    call compute_tissue_versatility(3, 1, expr, select_vec, 1, select_axes, 0, tv, angle)
    call assert_equal_real(tv(1), -1.0_real64, 0.0_real64, "No axes selected TV error")
  end subroutine test_invalid_input_no_axes

  !> @brief Test perfectly uniform expression (should yield TV=0).
  subroutine test_uniform_expression()
    real(real64) :: expr(3,1), tv(1), angle(1)
    logical :: select_vec(1), select_axes(3)
    expr(:,1) = [2.0_real64, 2.0_real64, 2.0_real64]
    select_vec = [.true.]
    select_axes = [.true., .true., .true.]
    call compute_tissue_versatility(3, 1, expr, select_vec, 1, select_axes, 3, tv, angle)
    call assert_equal_real(tv(1), 0.0_real64, 1e-12_real64, "Uniform expression TV")
    call assert_equal_real(angle(1), 0.0_real64, 1e-12_real64, "Uniform expression angle")
  end subroutine test_uniform_expression

  !> @brief Test expression in only one axis (should yield TV=1).
  subroutine test_single_axis_expression()
    real(real64) :: expr(3,1), tv(1), angle(1)
    logical :: select_vec(1), select_axes(3)
    expr(:,1) = [0.0_real64, 0.0_real64, 5.0_real64]
    select_vec = [.true.]
    select_axes = [.true., .true., .true.]
    call compute_tissue_versatility(3, 1, expr, select_vec, 1, select_axes, 3, tv, angle)
    call assert_equal_real(tv(1), 1.0_real64, 1e-12_real64, "Single axis TV")
    call assert_true(angle(1) > 0.0_real64, "Single axis angle > 0")
  end subroutine test_single_axis_expression

  !> @brief Test null vector (should yield TV=1, angle=90).
  subroutine test_null_vector()
    real(real64) :: expr(3,1), tv(1), angle(1)
    logical :: select_vec(1), select_axes(3)
    expr(:,1) = [0.0_real64, 0.0_real64, 0.0_real64]
    select_vec = [.true.]
    select_axes = [.true., .true., .true.]
    call compute_tissue_versatility(3, 1, expr, select_vec, 1, select_axes, 3, tv, angle)
    call assert_equal_real(tv(1), 1.0_real64, 1e-12_real64, "Null vector TV")
    call assert_equal_real(angle(1), 90.0_real64, 1e-12_real64, "Null vector angle")
  end subroutine test_null_vector

  !> @brief Test axis selection (subspace).
  subroutine test_partial_axis_selection()
    real(real64) :: expr(3,1), tv(1), angle(1)
    logical :: select_vec(1), select_axes(3)
    expr(:,1) = [1.0_real64, 2.0_real64, 3.0_real64]
    select_vec = [.true.]
    select_axes = [.true., .false., .true.]
    call compute_tissue_versatility(3, 1, expr, select_vec, 1, select_axes, 2, tv, angle)
    ! Should behave as a 2D vector [1,3]
    call assert_true(tv(1) >= 0.0_real64 .and. tv(1) <= 1.0_real64, "Partial axis TV in [0,1]")
    call assert_true(angle(1) >= 0.0_real64 .and. angle(1) <= 90.0_real64, "Partial axis angle in [0,90]")
  end subroutine test_partial_axis_selection

  !> @brief Test mixed vectors (one uniform, one single axis, one null).
  subroutine test_mixed_vectors()
    real(real64) :: expr(3,3), tv(3), angle(3)
    logical :: select_vec(3), select_axes(3)
    expr(:,1) = [1.0_real64, 1.0_real64, 1.0_real64] ! uniform
    expr(:,2) = [0.0_real64, 0.0_real64, 2.0_real64] ! single axis
    expr(:,3) = [0.0_real64, 0.0_real64, 0.0_real64] ! null
    select_vec = [.true., .true., .true.]
    select_axes = [.true., .true., .true.]
    call compute_tissue_versatility(3, 3, expr, select_vec, 3, select_axes, 3, tv, angle)
    call assert_equal_real(tv(1), 0.0_real64, 1e-12_real64, "Mixed: uniform TV")
    call assert_equal_real(tv(2), 1.0_real64, 1e-12_real64, "Mixed: single axis TV")
    call assert_equal_real(tv(3), 1.0_real64, 1e-12_real64, "Mixed: null TV")
    call assert_equal_real(angle(1), 0.0_real64, 1e-12_real64, "Mixed: uniform angle")
    call assert_true(angle(2) > 0.0_real64, "Mixed: single axis angle > 0")
    call assert_equal_real(angle(3), 90.0_real64, 1e-12_real64, "Mixed: null angle")
  end subroutine test_mixed_vectors

  !> @brief Test angle output in degrees for a known case.
  subroutine test_angle_degrees()
    real(real64) :: expr(2,1), tv(1), angle(1)
    logical :: select_vec(1), select_axes(2)
    expr(:,1) = [1.0_real64, 0.0_real64]
    select_vec = [.true.]
    select_axes = [.true., .true.]
    call compute_tissue_versatility(2, 1, expr, select_vec, 1, select_axes, 2, tv, angle)
    ! Angle should be 45 degrees (between [1,0] and [1,1])
    call assert_true(abs(angle(1) - 45.0_real64) < 1e-12_real64, "Angle output is 45 degrees")
  end subroutine test_angle_degrees

  !> @brief Test selection of multiple vectors.
  subroutine test_multiple_vectors_selection()
    real(real64) :: expr(2,3), tv(2), angle(2)
    logical :: select_vec(3), select_axes(2)
    expr(:,1) = [1.0_real64, 1.0_real64] ! uniform
    expr(:,2) = [0.0_real64, 2.0_real64] ! single axis
    expr(:,3) = [0.0_real64, 0.0_real64] ! null
    select_vec = [.true., .false., .true.]
    select_axes = [.true., .true.]
    call compute_tissue_versatility(2, 3, expr, select_vec, 2, select_axes, 2, tv, angle)
    call assert_equal_real(tv(1), 0.0_real64, 1e-12_real64, "Multiple: uniform TV")
    call assert_equal_real(tv(2), 1.0_real64, 1e-12_real64, "Multiple: null TV")
    call assert_equal_real(angle(1), 0.0_real64, 1e-12_real64, "Multiple: uniform angle")
    call assert_equal_real(angle(2), 90.0_real64, 1e-12_real64, "Multiple: null angle")
  end subroutine test_multiple_vectors_selection

end module mod_test_tissue_versatility