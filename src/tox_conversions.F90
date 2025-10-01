module tox_conversions
    use iso_fortran_env, only: int32, real64
    use iso_c_binding, only: c_int, c_double, c_null_char, c_double_complex, c_char, c_ptr, c_associated, c_f_pointer
    use tox_errors, only: ERR_ALLOC_FAIL, is_err, set_ok, set_err, set_err_once, ERR_POINTER_NULL, ERR_INVALID_INPUT
    implicit none

    ! type guards to guarantee kind identity between fortran and c for correct interop in the c wrapper routines
    logical, parameter :: c_int_matches_int32 = 1 == 1 / merge(1, 0, c_int == int32)
    logical, parameter :: c_double_matches_real64 = 1 == 1 / merge(1, 0, c_double == real64)
    logical, parameter :: c_double_complex_matches_real64 = 1 == 1 / merge(1, 0, c_double_complex == real64)

contains

    subroutine check_fortran_pointer_inputs(c_pointer, dims, ierr)
        type(c_ptr), intent(in) :: c_pointer
            !! C pointer to be converted
        integer(int32), dimension(:), intent(in) :: dims
            !! dimensions as array
        integer(int32), intent(inout) :: ierr
            !! Error code

        if (.not. c_associated(c_pointer)) then
            call set_err_once(ierr, ERR_POINTER_NULL)
        else if (any(dims <= 0)) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
        end if
    end subroutine check_fortran_pointer_inputs

    subroutine fortran_pointer_int_1d(c_pointer, f_pointer, dims, ierr)
        type(c_ptr), intent(in) :: c_pointer
            !! C pointer to be converted
        integer(int32), dimension(:), pointer, intent(out) :: f_pointer
            !! Resulting fortran pointer
        integer(int32), dimension(1), intent(in) :: dims
            !! dimensions as array
        integer(int32), intent(inout) :: ierr
            !! Error code

        call check_fortran_pointer_inputs(c_pointer, dims, ierr)
        if (is_err(ierr)) then
            nullify(f_pointer)
        else
            call c_f_pointer(c_pointer, f_pointer, dims)
        end if
    end subroutine fortran_pointer_int_1d

    subroutine fortran_pointer_int_2d(c_pointer, f_pointer, dims, ierr)
        type(c_ptr), intent(in) :: c_pointer
            !! C pointer to be converted
        integer(int32), dimension(:, :), pointer, intent(out) :: f_pointer
            !! Resulting fortran pointer
        integer(int32), dimension(2), intent(in) :: dims
            !! dimensions as array
        integer(int32), intent(inout) :: ierr
            !! Error code

        call check_fortran_pointer_inputs(c_pointer, dims, ierr)
        if (is_err(ierr)) then
            nullify(f_pointer)
        else
            call c_f_pointer(c_pointer, f_pointer, dims)
        end if
    end subroutine fortran_pointer_int_2d

    subroutine fortran_pointer_real_1d(c_pointer, f_pointer, dims, ierr)
        type(c_ptr), intent(in) :: c_pointer
            !! C pointer to be converted
        real(real64), dimension(:), pointer, intent(out) :: f_pointer
            !! Resulting fortran pointer
        integer(int32), dimension(1), intent(in) :: dims
            !! dimensions as array
        integer(int32), intent(inout) :: ierr
            !! Error code

        call check_fortran_pointer_inputs(c_pointer, dims, ierr)
        if (is_err(ierr)) then
            nullify(f_pointer)
        else
            call c_f_pointer(c_pointer, f_pointer, dims)
        end if
    end subroutine fortran_pointer_real_1d

    subroutine fortran_pointer_real_2d(c_pointer, f_pointer, dims, ierr)
        type(c_ptr), intent(in) :: c_pointer
            !! C pointer to be converted
        real(real64), dimension(:, :), pointer, intent(out) :: f_pointer
            !! Resulting fortran pointer
        integer(int32), dimension(2), intent(in) :: dims
            !! dimensions as array
        integer(int32), intent(inout) :: ierr
            !! Error code

        call check_fortran_pointer_inputs(c_pointer, dims, ierr)
        if (is_err(ierr)) then
            nullify(f_pointer)
        else
            call c_f_pointer(c_pointer, f_pointer, dims)
        end if
    end subroutine fortran_pointer_real_2d

    subroutine fortran_pointer_complex_1d(c_pointer, f_pointer, dims, ierr)
        type(c_ptr), intent(in) :: c_pointer
            !! C pointer to be converted
        complex(real64), dimension(:), pointer, intent(out) :: f_pointer
            !! Resulting fortran pointer
        integer(int32), dimension(1), intent(in) :: dims
            !! dimensions as array
        integer(int32), intent(inout) :: ierr
            !! Error code

        call check_fortran_pointer_inputs(c_pointer, dims, ierr)
        if (is_err(ierr)) then
            nullify(f_pointer)
        else
            call c_f_pointer(c_pointer, f_pointer, dims)
        end if
    end subroutine fortran_pointer_complex_1d

    subroutine fortran_pointer_complex_2d(c_pointer, f_pointer, dims, ierr)
        type(c_ptr), intent(in) :: c_pointer
            !! C pointer to be converted
        complex(real64), dimension(:, :), pointer, intent(out) :: f_pointer
            !! Resulting fortran pointer
        integer(int32), dimension(2), intent(in) :: dims
            !! dimensions as array
        integer(int32), intent(inout) :: ierr
            !! Error code

        call check_fortran_pointer_inputs(c_pointer, dims, ierr)
        if (is_err(ierr)) then
            nullify(f_pointer)
        else
            call c_f_pointer(c_pointer, f_pointer, dims)
        end if
    end subroutine fortran_pointer_complex_2d

    !> Converts a c_int value to logical, elemental -> any shape
    elemental subroutine c_int_as_logical(c_val, f_val)
        integer(c_int), intent(in) :: c_val
         !! the element containing the c variant of the number
        logical, intent(out) :: f_val
         !! the element that will hold the logical representation of the c_int: `0` means `.false.`, else `.true.`

        f_val = c_val /= 0
    end subroutine c_int_as_logical

    !> Converts a logical value to c_int, elemental -> any shape
    elemental subroutine logical_as_c_int(f_val, c_val)
        logical, intent(in) :: f_val
         !! the element containing the fortran variant of the number
        integer(c_int), intent(out) :: c_val
         !! the element that will hold the c_int representation of the logical: `0` if `.false.`, else `1`

        if (f_val) then
           c_val = 1
        else
           c_val = 0
        end if
    end subroutine logical_as_c_int

    !> Converts a c_char value to fortran character, elemental -> any shape
    elemental subroutine c_char_as_char(c_val, f_val)
        character(kind=c_char, len=1), intent(in) :: c_val
         !! the element containing the c variant of the character
        character(len=1), intent(out) :: f_val
         !! the element that will hold the fortran char representation of the c_char

         f_val = c_val
    end subroutine c_char_as_char

    !> Converts a character value to c_char, elemental -> any shape
    elemental subroutine char_as_c_char(f_val, c_val)
        character(len=1), intent(in) :: f_val
         !! the element containing the fortran variant of the character
        character(kind=c_char, len=1), intent(out) :: c_val
         !! the element that will hold the c_char (ASCII) representation of the character

        c_val = f_val
    end subroutine char_as_c_char

    !> Converts a 1D c_char array to string
    pure subroutine c_char_1d_as_string(c_char_array, str_out, ierr)
        character(kind=c_char, len=1), dimension(:), intent(in) :: c_char_array
         !! c int array, representing characters
        character(len=:), allocatable, intent(out) :: str_out
         !! Fortran string, length determined by occuring null char in `c_char_array`
        integer(int32), intent(out) :: ierr
         !! Error code

        integer(int32) :: i, str_len, array_len

        call set_ok(ierr)

        array_len = size(c_char_array, 1)

        ! identify string length
        str_len = array_len
        do i = 1, array_len
            if (c_char_array(i) == c_null_char) then
                str_len = i - 1
                exit
            end if
        end do

        ! create string
        allocate(character(len=str_len) :: str_out, stat=ierr)
        if (is_err(ierr)) then
           call set_err(ierr, ERR_ALLOC_FAIL)
           return
        end if

        do i = 1, str_len
            call c_char_as_char(c_char_array(i), str_out(i:i))
        end do
    end subroutine c_char_1d_as_string

    !> Converts a string to 1D c_char array
    pure subroutine string_as_c_char_1d(str, c_char_array)
        character(len=*), intent(in) :: str
         !! Fortran string to be converted
        character(kind=c_char, len=1), dimension(:), intent(out) :: c_char_array
         !! c int array, representing chars of `str`, will always end with null char. If array too small, it will hold fitting trimmed `str`.

        integer(int32) :: i_str, str_len

        ! determine string length to be converted, +1 because of the null char in c
        str_len = min(len_trim(str) + 1, size(c_char_array, 1)) - 1

        do i_str = 1, str_len
           call c_char_as_char(str(i_str:i_str), c_char_array(i_str))
        end do

        if (size(c_char_array, 1) > 0) then
           c_char_array(str_len + 1) = c_null_char
        end if
    end subroutine string_as_c_char_1d

    !> Converts a 2D c_char array to 1D string array
    pure subroutine c_char_2d_as_string(c_char_array, str_out, ierr)
        character(kind=c_char, len=1), dimension(:, :), intent(in) :: c_char_array
         !! c int array, columns as ascii arrays
        character(len=:), dimension(:), allocatable, intent(out) :: str_out
         !! Fortran array of resulting strings
        integer(int32), intent(out) :: ierr
         !! Error code

        integer(int32) :: i_str, n_rows, n_strings
        character(len=:), allocatable :: string

        call set_ok(ierr)

        n_rows = size(c_char_array, 1)
        n_strings = size(c_char_array, 2)

        allocate(character(len=n_rows) :: str_out(n_strings), stat=ierr)
        if (is_err(ierr)) then
           call set_err(ierr, ERR_ALLOC_FAIL)
           return
        end if

        ! create strings
        do i_str = 1, n_strings
            call c_char_1d_as_string(c_char_array(:, i_str), string, ierr)
            if (is_err(ierr)) return
            str_out(i_str) = string
        end do
    end subroutine c_char_2d_as_string

    !> Converts a 1D string array to 2D c_char array
    pure subroutine string_as_c_char_2d(strings, c_char_array)
        character(len=*), dimension(:), intent(in) :: strings
         !! Fortran array of strings
        character(kind=c_char, len=1), dimension(:, :), intent(out) :: c_char_array
         !! c int array, columns as ascii arrays

        integer(int32) :: i_str, n_strings

        do i_str = 1, size(strings, 1)
           call string_as_c_char_1d(strings(i_str), c_char_array(:, i_str))
        end do
    end subroutine string_as_c_char_2d

end module tox_conversions
