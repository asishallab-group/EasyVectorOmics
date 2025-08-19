module tox_conversions
    use iso_fortran_env, only: int32, real64
    use iso_c_binding, only: c_int, c_double, c_null_char, c_double_complex
    implicit none

contains

    !> Casts a c_double value into real64, elemental -> any shape
    elemental subroutine c_double_as_real64(c_val, f_val)
        real(c_double), intent(in) :: c_val
         !! the element containing the c variant of the number
        real(real64), intent(out) :: f_val
         !! the element that will hold the real64 representation of the c_double

        f_val = real(c_val, kind=real64)
    end subroutine c_double_as_real64

    !> Casts a real64 value into c_double, elemental -> any shape
    elemental subroutine real64_as_c_double(f_val, c_val)
        real(real64), intent(in) :: f_val
         !! the element containing the fortran variant of the number
        real(c_double), intent(out) :: c_val
         !! the element that will hold the c_double representation of the real64

        c_val = real(f_val, kind=c_double)
    end subroutine real64_as_c_double

    !> Casts a c_double_complex value into fortran complex(real64), elemental -> any shape
    elemental subroutine c_complex_as_complex(c_val, f_val)
        complex(c_double_complex), intent(in) :: c_val
         !! the element containing the c variant of the number
        complex(real64), intent(out) :: f_val
         !! the element that will hold the complex(real64) representation of the c_double_complex

        f_val = cmplx(real(c_val), aimag(c_val), kind=real64)
    end subroutine c_complex_as_complex

    !> Casts a fortran complex(real64) value into c_double_complex, elemental -> any shape
    elemental subroutine complex_as_c_complex(f_val, c_val)
        complex(real64), intent(in) :: f_val
         !! the element containing the fortran variant of the number
        complex(c_double_complex), intent(out) :: c_val
         !! the element that will hold the c_double_complex representation of the complex(real64)

        c_val = cmplx(real(f_val), aimag(f_val), kind=c_double_complex)
    end subroutine complex_as_c_complex

    !> Casts a c_int value into int32, elemental -> any shape
    elemental subroutine c_int_as_int32(c_val, f_val)
        integer(c_int), intent(in) :: c_val
         !! the element containing the c variant of the number
        integer(int32), intent(out) :: f_val
         !! the element that will hold the int32 representation of the c_int

        f_val = int(c_val, kind=int32)
    end subroutine c_int_as_int32

    !> Casts a int32 value into c_int, elemental -> any shape
    elemental subroutine int32_as_c_int(f_val, c_val)
        integer(int32), intent(in) :: f_val
         !! the element containing the fortran variant of the number
        integer(c_int), intent(out) :: c_val
         !! the element that will hold the c_int representation of the int32

        c_val = int(f_val, kind=c_int)
    end subroutine int32_as_c_int

    !> Casts a c_int value into logical, elemental -> any shape
    elemental subroutine c_int_as_logical(c_val, f_val)
        integer(c_int), intent(in) :: c_val
         !! the element containing the c variant of the number
        logical, intent(out) :: f_val
         !! the element that will hold the logical representation of the c_int

        f_val = c_val /= 0
    end subroutine c_int_as_logical

    !> Casts a logical value into c_int, elemental -> any shape
    elemental subroutine logical_as_c_int(f_val, c_val)
        logical, intent(in) :: f_val
         !! the element containing the fortran variant of the number
        integer(c_int), intent(out) :: c_val
         !! the element that will hold the c_int representation of the logical

        if (f_val) then
           c_val = 1
        else
           c_val = 0
        end if
    end subroutine logical_as_c_int

    !> Casts a c_int value into character, elemental -> any shape
    elemental subroutine c_int_as_char(c_val, f_val)
        integer(c_int), intent(in) :: c_val
         !! the element containing the c variant of the character
        character(len=1), intent(out) :: f_val
         !! the element that will hold the char representation of the c_int

        if (c_val > 0_c_int .and. c_val < 128_c_int) then
            f_val = achar(c_val)
        else
            f_val = "?"
        end if
    end subroutine c_int_as_char

    !> Casts a character value into c_int, elemental -> any shape
    elemental subroutine char_as_c_int(f_val, c_val)
        character(len=1), intent(in) :: f_val
         !! the element containing the fortran variant of the character
        integer(c_int), intent(out) :: c_val
         !! the element that will hold the c_int (ASCII) representation of the character

        c_val = iachar(f_val)
    end subroutine char_as_c_int

    !> Casts a 1D c_int array into string
    pure subroutine c_int_1d_as_string(c_int_array, str_out)
        integer(c_int), dimension(:), intent(in) :: c_int_array
         !! c int array, representing characters
        character(len=:), allocatable, intent(out) :: str_out
         !! Fortran string, length determined by occuring null char in `c_int_array`

        integer(int32) :: i, str_len, array_len
        character(len=1) :: char

        array_len = size(c_int_array, 1)

        ! identify string length
        str_len = array_len
        do i = 1, array_len
            if (c_int_array(i) == iachar(c_null_char)) then
                str_len = i - 1
                exit
            end if
        end do

        ! create string
        allocate (character(len=str_len) :: str_out)
        do i = 1, str_len
            call c_int_as_char(c_int_array(i), char)
            str_out(i:i) = char
        end do
    end subroutine c_int_1d_as_string

    !> Casts a string into 1D c_int array
    pure subroutine string_as_c_int_1d(str, c_int_array)
        character(len=*), intent(in) :: str
         !! Fortran string to be converted
        integer(c_int), dimension(:), intent(out) :: c_int_array
         !! c int array, representing chars of `str`, will always end with null char. If array too small, it will hold fitting trimmed `str`.

        integer(int32) :: i_str, str_len

        ! determine string length to be converted, +1 because of the null char in c
        str_len = min(len_trim(str) + 1, size(c_int_array, 1)) - 1

        do i_str = 1, str_len
           c_int_array(i_str) = iachar(str(i_str:i_str))
        end do

        if (size(c_int_array, 1) > 0) then
           c_int_array(str_len + 1) = iachar(c_null_char)
        end if
    end subroutine string_as_c_int_1d

    !> Casts a 2D c_int array into 1D string array
    pure subroutine c_int_2d_as_string(c_int_array, str_out)
        integer(c_int), dimension(:, :), intent(in) :: c_int_array
         !! c int array, columns as ascii arrays
        character(len=:), dimension(:), allocatable, intent(out) :: str_out
         !! Fortran array of resulting strings

        integer(int32) :: i_str, n_rows, n_strings
        character(len=:), allocatable :: string

        n_rows = size(c_int_array, 1)
        n_strings = size(c_int_array, 2)

        allocate (character(len=n_rows) :: str_out(n_strings))
        ! create strings
        do i_str = 1, n_strings
            call c_int_1d_as_string(c_int_array(:, i_str), string)
            str_out(i_str) = string
        end do
    end subroutine c_int_2d_as_string

    !> Casts a 1D string array into 2D c_int array
    pure subroutine string_as_c_int_2d(strings, c_int_array)
        character(len=*), dimension(:), intent(in) :: strings
         !! Fortran array of strings
        integer(c_int), dimension(:, :), intent(out) :: c_int_array
         !! c int array, columns as ascii arrays

        integer(int32) :: i_str, n_strings

        do i_str = 1, size(strings, 1)
           call string_as_c_int_1d(strings(i_str), c_int_array(:, i_str))
        end do
    end subroutine string_as_c_int_2d

end module tox_conversions
