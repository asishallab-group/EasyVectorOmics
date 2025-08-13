module c_conversions
    use iso_fortran_env, only: int32, real64
    use iso_c_binding, only: c_ptr, c_int, c_double
    implicit none

contains

    !> Casts a c_double value into real64, elemental -> any shape
    elemental subroutine c_double_as_real64(c_val, f_real64)
        real(c_double), intent(in) :: c_val
         !! the element containing the c variant of the number
        real(real64), intent(out) :: f_real64
         !! the element that will hold the real64 representation of the c_double

        f_real64 = real(c_val, kind=real64)
    end subroutine c_double_as_real64

    !> Casts a c_int value into int32, elemental -> any shape
    elemental subroutine c_int_as_int32(c_val, f_int32)
        integer(c_int), intent(in) :: c_val
         !! the element containing the c variant of the number
        integer(int32), intent(out) :: f_int32
         !! the element that will hold the int32 representation of the c_int

        f_int32 = int(c_val, kind=int32)
    end subroutine c_int_as_int32

    !> Casts a c_int value into character, elemental -> any shape
    elemental subroutine c_int_as_char(c_val, f_char)
        integer(c_int), intent(in) :: c_val
         !! the element containing the c variant of the character
        character(len=1), intent(out) :: f_char
         !! the element that will hold the char representation of the c_int

        if (c_val > 0_c_int .and. c_val < 128_c_int) then
            f_char = achar(c_val)
        else
            f_char = "?"
        end if
    end subroutine c_int_as_char

    !> Treats a c_char array as ascii char array and translates it into Fortran character
    pure subroutine c_int_1d_as_string(c_int_array, array_len, str_out)
        integer(c_int), dimension(array_len), intent(in) :: c_int_array
         !! c int array, representing characters
        character(len=:), allocatable, intent(out) :: str_out
         !! Fortran string of input length
        integer(int32), intent(in) :: array_len
         !! string length

        integer(int32) :: i, str_len
        character(len=1) :: char

        ! identify string length
        str_len = array_len
        do i = 1, array_len
            if (c_int_array(i) == 0) then
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

    !> Treats a c_char array as ascii char array and translates it into Fortran character
    pure subroutine c_int_2d_as_string(c_int_array, n_rows, n_strings, str_out)
        integer(c_int), dimension(n_rows, n_strings), intent(in) :: c_int_array
         !! c int array, columns as ascii arrays
        character(len=:), dimension(:), allocatable, intent(out) :: str_out
         !! Fortran string of input length
        integer(int32), intent(in) :: n_rows
         !! number of rows in `c_int_array`, so maximum string length
        integer(int32), intent(in) :: n_strings
         !! number of columns/strings in `c_int_array`

        integer(int32) :: i_str
        character(len=:), allocatable :: string

        allocate (character(len=n_rows) :: str_out(n_strings))
        ! create strings
        do i_str = 1, n_strings
            call c_int_1d_as_string(c_int_array(:, i_str), n_rows, string)
            str_out(i_str) = string
        end do
    end subroutine c_int_2d_as_string

end module c_conversions
