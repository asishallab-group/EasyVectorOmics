!> Module for array utilities
module array_utils
    use, intrinsic :: iso_fortran_env, only: int32, real64
    implicit none

    PUBLIC :: get_array_dims, get_array_metadata_chars, ascii_to_string, read_file_header, write_file_header

    integer(int32), parameter :: ARRAY_FILE_MAGIC = int(z'46413230', int32) ! 'FA20' in hex

   contains

  subroutine write_file_header(filename, unit, type_code, ndim, dims, ierr, clen)
    use iso_fortran_env, only: int32
    implicit none
    character(len=*), intent(in) :: filename
    integer(int32), intent(in) :: type_code, ndim
    integer(int32), intent(in) :: dims(ndim)
    integer(int32), intent(in), optional :: clen
    integer(int32), intent(out) :: ierr
    integer(int32), INTENT(OUT) :: unit

    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace', iostat=ierr)
    if (ierr /= 0) then
      ierr = 400
      return
    end if
    write(unit, iostat=ierr) ARRAY_FILE_MAGIC
    if (ierr /= 0) then
      ierr = 401
      return
    end if
    write(unit, iostat=ierr) type_code
    if (ierr /= 0) then
      ierr = 402
      return
    end if
    write(unit, iostat=ierr) ndim
    if (ierr /= 0) then
      ierr = 403
      return
    end if
    write(unit, iostat=ierr) dims
    if (ierr /= 0) then
      ierr = 404
      return
    end if
    if (type_code == 3 .and. present(clen)) then
      write(unit, iostat=ierr) clen
      if (ierr /= 0) then
        ierr = 405
        return
      end if
    end if
  end subroutine write_file_header

  subroutine read_file_header(filename, unit, type_code, ndims, dims, clen, ierr)
    character(len=*), intent(in) :: filename
    integer(int32), intent(out) :: unit
    integer(int32), intent(out) :: type_code, ndims, clen
    integer(int32), intent(out) :: ierr
    integer(int32), allocatable :: dims(:)

    integer(int32) :: magic
    ierr = 0
    clen = 0

    
    if (ierr /= 0) then
      ierr = 101
      close(unit)
      return
    end if

    read(unit, iostat=ierr) magic
    if (ierr /= 0) then
      ierr = 102
      close(unit)
      return
    end if

    if (magic /= ARRAY_FILE_MAGIC) then
      ierr = 200
      close(unit)
      return
    end if

    read(unit, iostat=ierr) type_code
    if (ierr /= 0) then
      ierr = 103
      close(unit)
      return
    end if

    read(unit, iostat=ierr) ndims
    if (ierr /= 0) then
      ierr = 104
      close(unit)
      return
    end if

    allocate(dims(ndims))
    read(unit, iostat=ierr) dims
    if (ierr /= 0) then
      ierr = 105
      close(unit)
      return
    end if

    if(type_code==3) then
      read(unit, iostat=ierr) clen
      if (ierr /= 0) then
        ierr = 106
        close(unit)
        return
      end if
    else
      clen = 0 ! Not applicable for non-character types
    end if

  end subroutine

  !> Get the dimensions of an array file
  subroutine get_array_dims(filename, dims_out, ndims, ierr)
    implicit none

    character(len=*), intent(in) :: filename
    integer(int32), intent(out) :: dims_out(*)
    integer(int32), intent(out) :: ndims
    integer, intent(out) :: ierr

    integer :: unit, i
    integer(int32), allocatable :: dims(:)
    integer(int32) :: type_code, clen

    ! error handling
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='old', iostat=ierr)
    call read_file_header(filename, unit, type_code, ndims, dims, clen, ierr)
    close(unit)
    if (ierr /= 0) then
      return
    end if

    do i = 1, ndims
      dims_out(i) = dims(i)
    end do
  end subroutine

  !> subroutine to convert an ASCII array to a string
subroutine ascii_to_string(ascii_array, clen, str)
  implicit none

  integer(int32), intent(in) :: ascii_array(clen)
    !! Array of ASCII characters
  integer(int32), intent(in) :: clen
    !! Length of the ASCII array
  character(len=:), allocatable, INTENT(INOUT) :: str
  integer(int32) :: i

  allocate(character(len=clen) :: str)
  do i = 1, clen
    str(i:i) = char(ascii_array(i))
  end do

  ! Return the string (or use it as needed)
end subroutine ascii_to_string

!> Subroutine to get metadata of a char array file
  subroutine get_array_metadata_chars(filename, dims_out, ndims, type_code_out, clen_out, ierr)
    implicit none

    character(len=*), intent(in) :: filename
    !! Name of the file to read
    integer(int32), intent(out) :: dims_out(*)
    !! Output array for dimensions
    integer(int32), intent(out) :: ndims
    !! Output variable for the number of dimensions
    integer(int32), intent(out) :: type_code_out
    !! Output variable for the type code
    integer(int32), intent(out) :: clen_out
    !! Output variable for the character length
    integer(int32), intent(out) :: ierr
    !! Error code

    integer :: unit, i
    integer(int32), allocatable :: dims(:)

    ! error handling
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='old', iostat=ierr)
    call read_file_header(filename, unit, type_code_out, ndims, dims, clen_out, ierr)
    close(unit)
    if (ierr /= 0) then
      return
    end if

    do i = 1, ndims
      dims_out(i) = dims(i)
    end do
  end subroutine

end module array_utils

!> Subroutine to get the dimensions of an array file
subroutine get_array_dims_r(filename_ascii, fn_len, dims_out, ndims, ierr)
  use iso_fortran_env, only: int32
  use array_utils, only: ascii_to_string, get_array_dims
  implicit none

  ! Input
  integer(int32), intent(in) :: filename_ascii(fn_len)
    !! Array of ASCII characters representing the filename
  integer(int32), intent(in) :: fn_len
    !! Length of the filename array
  character(len=:), allocatable :: filename

  ! Output
  integer(int32), intent(out) :: dims_out(*)  ! R provides storage
    !! Output array for dimensions
  integer(int32), intent(out) :: ndims        ! Number of dimensions
    !! Output variable for the number of dimensions

  ! Local variables
  integer(int32) :: unit, magic, type_code, d, i
  integer(int32), allocatable :: dims(:)

  integer(int32), intent(out) :: ierr
  !! Error code

  ! ASCII → String
  call ascii_to_string(filename_ascii, fn_len, filename)

  call get_array_dims(filename, dims_out, ndims, ierr)

end subroutine get_array_dims_r

!> C binding for the subroutine to get the dimensions of an array file
subroutine get_array_dims_C(filename_ascii, fn_len, dims_out, ndims, ierr) bind(C, name="get_array_dims_C")
  use iso_c_binding, only: c_int
  use iso_fortran_env, only : int32
  use array_utils, only : ascii_to_string, get_array_dims
  implicit none

  ! Input
  integer(c_int), intent(in) :: filename_ascii(fn_len)
    !! Array of ASCII characters representing the filename
  integer(c_int), value :: fn_len
    !! Length of the filename array

  ! Output
  integer(c_int), intent(out) :: dims_out(*)
    !! Output array for dimensions
  integer(c_int), intent(out) :: ndims
    !! Output variable for the number of dimensions

  ! Local variables
  character(len=:), allocatable :: filename
  integer(c_int) :: unit, magic, type_code, d, i
  integer(c_int), allocatable :: dims(:)

  integer(c_int), intent(out) :: ierr
    !! Error code

  ! ASCII → String
  call ascii_to_string(filename_ascii, fn_len, filename)

  call get_array_dims(filename, dims_out, ndims, ierr)
end subroutine get_array_dims_C

!> Subroutine to get metadata of a char array file
subroutine get_array_metadata_chars_r(filename_ascii, fn_len, dims_out, ndims, type_code_out, clen_out, ierr)
  use iso_fortran_env, only: int32
  use array_utils, only : ascii_to_string, get_array_metadata_chars
  implicit none

  ! input
  integer(int32), intent(in) :: filename_ascii(fn_len)
    !! Array of ASCII characters representing the filename
  integer(int32), intent(in) :: fn_len
    !! Length of the filename array

  ! Output
  integer(int32), intent(out) :: dims_out(*)
    !! Output array for dimensions
  integer(int32), intent(out) :: ndims
    !! Output variable for the number of dimensions
  integer(int32), intent(out) :: type_code_out
    !! Output variable for the type code
  integer(int32), intent(out) :: clen_out
    !! Output variable for the character length
  integer(int32), intent(out) :: ierr
    !! Error code

  ! Internally
  character(len=:), allocatable :: filename
  integer(int32) :: unit, magic, type_code, d, i
  integer(int32), allocatable :: dims(:)

  ! ASCII → String
  call ascii_to_string(filename_ascii, fn_len, filename)

  call get_array_metadata_chars(filename, dims_out, ndims, type_code_out, clen_out, ierr)
end subroutine

!> C binding for the subroutine to get metadata of a char array file
subroutine get_array_metadata_chars_C(filename_ascii, fn_len, dims_out, dims_len, ndims, &
                          type_code_out, clen_out, ierr) bind(C, name="get_array_metadata_chars_C")
  use iso_c_binding, only: c_int
  use iso_fortran_env, only: int32
  use array_utils, only: ascii_to_string, get_array_metadata_chars
  implicit none

  ! Output
  integer(c_int), intent(in) :: filename_ascii(fn_len)
    !! Length of the filename array
  integer(c_int), value :: fn_len, dims_len

  ! Input
  integer(c_int), intent(out) :: dims_out(dims_len)
  !! Output array for dimensions
  integer(c_int), intent(out) :: ndims
  !! Output variable for the number of dimensions
  integer(c_int), intent(out) :: type_code_out
  !! Output variable for the type code
  integer(c_int), intent(out) :: clen_out
  !! Output variable for the character length
  integer(c_int), intent(out) :: ierr
  !! Error code

  ! Locally
  character(len=:), allocatable :: filename
  integer(int32) :: unit, magic, type_code, d, i
  integer(c_int), allocatable :: dims(:)

  ! ASCII → String
  call ascii_to_string(filename_ascii, fn_len, filename)

  call get_array_metadata_chars(filename, dims_out, ndims, type_code_out, clen_out, ierr)
end subroutine
