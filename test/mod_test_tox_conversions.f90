!> Unit test suite for tox_conversions routine.
module mod_test_tox_conversions
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use iso_c_binding
    use tox_conversions
    implicit none

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

    integer(int32), parameter :: TEST_COUNT = 5
    real(real64), parameter :: TOL = 1d-12

contains

    !> Get array of all available tests.
    function get_all_tests() result(all_tests)
        type(test_case) :: all_tests(TEST_COUNT)

        all_tests(1) = test_case("test_tox_conversions_c_double_as_real64", test_tox_conversions_c_double_as_real64)
        all_tests(2) = test_case("test_tox_conversions_c_int_as_int32", test_tox_conversions_c_int_as_int32)
        all_tests(3) = test_case("test_tox_conversions_c_int_as_char", test_tox_conversions_c_int_as_char)
        all_tests(4) = test_case("test_tox_conversions_c_int_1d_as_string", test_tox_conversions_c_int_1d_as_string)
        all_tests(5) = test_case("test_tox_conversions_c_int_2d_as_string", test_tox_conversions_c_int_2d_as_string)
    end function get_all_tests

    subroutine test_tox_conversions_c_double_as_real64
        real(c_double) :: c_val
        real(real64) :: casted_scalar
        real(real64) :: expected_scalar
        real(real64) :: casted_array(1)

        c_val = 3.141592653589793_c_double
        expected_scalar = 3.141592653589793_real64
        call c_double_as_real64(c_val, casted_scalar)
        call assert_equal_real(casted_scalar, expected_scalar, TOL, "test_tox_conversions_c_double_as_real64: value mismatch")

        call c_double_as_real64([c_val], casted_array)
        call assert_equal_array_real(casted_array, [expected_scalar], 1, TOL, "test_tox_conversions_c_double_as_real64: value mismatch")
    end subroutine test_tox_conversions_c_double_as_real64

    subroutine test_tox_conversions_c_int_as_int32
        integer(c_int) :: c_val
        integer(int32) :: casted_scalar
        integer(int32) :: casted_array(1)
        integer(int32) :: expected_scalar

        c_val = 42_c_int
        expected_scalar = 42_int32
        call c_int_as_int32(c_val, casted_scalar)
        call assert_equal_int(casted_scalar, expected_scalar, "test_tox_conversions_c_int_as_int32: value mismatch")

        call c_int_as_int32([c_val], casted_array)
        call assert_equal_array_int(casted_array, [expected_scalar], 1, "test_tox_conversions_c_int_as_int32: value mismatch")
    end subroutine test_tox_conversions_c_int_as_int32

    subroutine test_tox_conversions_c_int_as_char
        integer(c_int) :: c_val
        character(len=1) :: casted_scalar
        character(len=1) :: casted_array(1)
        character(len=1) :: expected_scalar

        expected_scalar = "H"
        c_val = iachar(expected_scalar)
        call c_int_as_char(c_val, casted_scalar)
        call assert_true(casted_scalar == expected_scalar, "test_tox_conversions_c_int_as_char: value mismatch for scalar")

        call c_int_as_char([c_val], casted_array)
        call assert_true(casted_array(1) == expected_scalar, "test_tox_conversions_c_int_as_char: value mismatch for char array")

        ! non ascii
        expected_scalar = "?"
        call c_int_as_char(128_c_int, casted_scalar)
        call assert_true(casted_scalar == expected_scalar, "test_tox_conversions_c_int_as_char: value mismatch for non ascii (too large)")

        call c_int_as_char(-1_c_int, casted_scalar)
        call assert_true(casted_scalar == expected_scalar, "test_tox_conversions_c_int_as_char: value mismatch for non ascii (negative)")
    end subroutine test_tox_conversions_c_int_as_char

    subroutine test_tox_conversions_c_int_1d_as_string
        integer(c_int) :: c_int_array(5)
        character(len=:), allocatable :: f_char

        c_int_array = [72_c_int, 101_c_int, -11_c_int, 108_c_int, 111_c_int] ! 'H', 'e', non-ASCII, 'l', 'o'
        call c_int_1d_as_string(c_int_array, f_char)
        call assert_true(f_char == "He?lo", "test_tox_conversions_c_int_1d_as_string: value mismatch")

        ! test empty string
        c_int_array = [0_c_int, 101_c_int, -11_c_int, 108_c_int, 111_c_int]
        call c_int_1d_as_string(c_int_array, f_char)
        call assert_true(f_char == "", "test_tox_conversions_c_int_1d_as_string: value mismatch")
    end subroutine test_tox_conversions_c_int_1d_as_string

    subroutine test_tox_conversions_c_int_2d_as_string
        integer(c_int) :: c_int_array(5, 2)
        character(len=:), allocatable :: f_char(:)

        c_int_array(:, 1) = [72_c_int, 101_c_int, -11_c_int, 108_c_int, 111_c_int] ! 'H', 'e', non-ASCII, 'l', 'o'
        c_int_array(:, 2) = [0_c_int, 101_c_int, -11_c_int, 108_c_int, 111_c_int] ! string ends at first element

        call c_int_2d_as_string(c_int_array, f_char)
        call assert_equal_int(size(f_char, 1), 2, "test_tox_conversions_c_int_2d_as_string: Did not get two strings")
        call assert_true(f_char(1) == "He?lo" .and. f_char(2) == "", "test_tox_conversions_c_int_2d_as_string: value mismatch")
    end subroutine test_tox_conversions_c_int_2d_as_string

    !> Run all tox_conversions tests.
    subroutine run_all_tests_tox_conversions
        type(test_case) :: all_tests(TEST_COUNT)
        integer(int32) :: i

        all_tests = get_all_tests()

        do i = 1, size(all_tests)
            call all_tests(i)%test_proc()
            print *, trim(all_tests(i)%name), " passed."
        end do
        print *, "All tox_conversions tests passed successfully."
    end subroutine run_all_tests_tox_conversions

    !> Run specific tox_conversions tests by name.
    subroutine run_named_tests_tox_conversions(test_names)
        character(len=*), intent(in) :: test_names(:)
        type(test_case) :: all_tests(TEST_COUNT)
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
    end subroutine run_named_tests_tox_conversions
end module mod_test_tox_conversions
