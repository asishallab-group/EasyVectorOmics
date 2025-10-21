!> @file mod_test_edf.f90
!> Unit test suite for EDF (Empirical Distribution Function)
!> Contains dedicated tests for compute_edf from f42_utils.

module mod_test_edf
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

  !> Get array of all available EDF tests.
  function get_all_tests_edf() result(all_tests)
    type(test_case) :: all_tests(8)
    all_tests(1) = test_case("test_edf_simple", test_edf_simple)
    all_tests(2) = test_case("test_edf_all_unique", test_edf_all_unique)
    all_tests(3) = test_case("test_edf_all_same", test_edf_all_same)
    all_tests(4) = test_case("test_edf_duplicates", test_edf_duplicates)
    all_tests(5) = test_case("test_edf_single_value", test_edf_single_value)
    all_tests(6) = test_case("test_edf_empty_input", test_edf_empty_input)
    all_tests(7) = test_case("test_edf_large_dataset", test_edf_large_dataset)
    all_tests(8) = test_case("test_edf_negative_values", test_edf_negative_values)
  end function get_all_tests_edf

  !> Run all EDF tests.
  subroutine run_all_tests_edf()
  type(test_case) :: all_tests(8)
  integer(int32) :: i
  all_tests = get_all_tests_edf()
    do i = 1, size(all_tests)
      call all_tests(i)%test_proc()
      print *, trim(all_tests(i)%name), " passed."
    end do
    print *, "All EDF tests passed successfully."
  end subroutine run_all_tests_edf

  !> Run specific EDF tests by name.
  subroutine run_named_tests_edf(test_names)
    character(len=*), intent(in) :: test_names(:)
  type(test_case) :: all_tests(8)
  integer(int32) :: i, j
  logical :: found
  all_tests = get_all_tests_edf()
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
  end subroutine run_named_tests_edf

  ! ============================================================================
  ! EDF Tests
  ! ============================================================================

  !> Test EDF computation with a simple dataset.
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

  !> Test EDF with all unique values.
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

  !> Test EDF with all identical values.
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

  !> Test EDF with specific duplicate pattern.
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

  !> Test EDF with a single value.
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

  !> Test EDF with empty input (should return error).
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

  !> Test EDF with a larger dataset (100 values).
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

  !> Test EDF with negative values and zero.
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

end module mod_test_edf
