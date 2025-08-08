!> Module for array utilities
module array_utils
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use iso_c_binding
    implicit none

    public :: get_type_code

    integer(int32), parameter :: ARRAY_FILE_MAGIC = int(z'46413230', int32) ! 'FA20' in hex

   contains
  !> returns type code of the array file
  function get_type_code(filename) result(type_code)
    use iso_c_binding
    implicit none

    character(len=*), intent(in) :: filename
      !! name of the file to read
    integer :: type_code
    integer :: unit, magic

    open(newunit=unit, file=filename, form='unformatted', access='stream', status='old')
    read(unit) magic
    if (magic /= ARRAY_FILE_MAGIC) error stop "Invalid file format"
    read(unit) type_code
    close(unit)
  end function get_type_code

  !> Get the dimensions of an array file
  subroutine get_array_dims(filename, dims_out, ndims)
    use iso_c_binding
    implicit none

    character(len=*), intent(in) :: filename
      !! Name of the file to read
    integer(int32), intent(out) :: dims_out(*)
      !! Output array for dimensions
    integer(int32), intent(out) :: ndims
      !! Number of dimensions
    integer(int32), ALLOCATABLE:: dims(:)
      !! Output array for dimensions
    integer :: unit, magic, type_code, d, i

    open(newunit=unit, file=filename, form='unformatted', access='stream', status='old')
    read(unit) magic
    if (magic /= ARRAY_FILE_MAGIC) error stop "Invalid file format"
    read(unit) type_code
    read(unit) d
    allocate(dims(d))
    read(unit) dims
    close(unit)

    ndims = d
    do i = 1, d
      dims_out(i) = dims(i)
    end do
  end subroutine get_array_dims

  !> subroutine to convert an ASCII array to a string
subroutine ascii_to_string(ascii_array, clen, str)
  use iso_fortran_env, only: int32
  implicit none

  integer(int32), intent(in) :: ascii_array(clen)
    !! Array of ASCII characters
  integer(int32), intent(in) :: clen
    !! Length of the ASCII array
  character(len=:), allocatable, INTENT(INOUT) :: str
  integer :: i

  allocate(character(len=clen) :: str)
  do i = 1, clen
    str(i:i) = char(ascii_array(i))
  end do

  ! Return the string (or use it as needed)
end subroutine ascii_to_string

!> Subroutine to get metadata of a char array file
subroutine get_array_metadata_chars(filename, dims_out, ndims, type_code_out, clen_out)
  use iso_fortran_env, only: int32
  implicit none

  ! Input
  character(len=*), intent(in) :: filename
    !! Name of the file to read
  integer(int32), intent(out) :: dims_out(*)
    !! Output array for dimensions
  integer(int32), intent(out) :: ndims
    !! Number of dimensions
  integer(int32), intent(out) :: type_code_out
    !! Type code of the array
  integer(int32), intent(out) :: clen_out
    !! Character length if applicable

  ! Local variables
  integer :: unit, magic, type_code, d, i
  integer(int32), allocatable :: dims(:)

  ! Open file
  open(newunit=unit, file=filename, form='unformatted', access='stream', status='old')

  ! Read Header
  read(unit) magic
  if (magic /= ARRAY_FILE_MAGIC) error stop "Invalid file format"
  read(unit) type_code
  read(unit) d
  allocate(dims(d))
  read(unit) dims

  if (type_code == 3) then
    read(unit) clen_out
  else
    clen_out = 0
  end if

  close(unit)

  ! Return values
  type_code_out = type_code
  ndims = d

  do i = 1, d
    dims_out(i) = dims(i)
  end do

end subroutine get_array_metadata_chars

end module array_utils

!> Subroutine to get the dimensions of an array file
subroutine get_array_dims_r(filename_ascii, fn_len, dims_out, ndims)
  use iso_fortran_env, only: int32
  use array_utils
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
  integer :: unit, magic, type_code, d, i
  integer(int32), allocatable :: dims(:)

  ! ASCII → String
  call ascii_to_string(filename_ascii, fn_len, filename)

  call get_array_dims(filename, dims_out, ndims)

end subroutine get_array_dims_r

!> C binding for the subroutine to get the dimensions of an array file
subroutine get_array_dims_C(filename_ascii, fn_len, dims_out, ndims) bind(C, name="get_array_dims_C")
  use iso_c_binding
  use iso_fortran_env
  use array_utils
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


  ! ASCII → String
  call ascii_to_string(filename_ascii, fn_len, filename)

  call get_array_dims(filename, dims_out, ndims)
end subroutine get_array_dims_C

!> Subroutine to get metadata of a char array file
subroutine get_array_metadata_chars_r(filename_ascii, fn_len, dims_out, ndims, type_code_out, clen_out)
  use iso_fortran_env, only: int32
  use array_utils
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

  ! Internally
  character(len=:), allocatable :: filename
  integer :: unit, magic, type_code, d, i
  integer(int32), allocatable :: dims(:)

  ! ASCII → String
  call ascii_to_string(filename_ascii, fn_len, filename)

  call get_array_metadata_chars(filename, dims_out, ndims, type_code_out, clen_out)
end subroutine

!> C binding for the subroutine to get metadata of a char array file
subroutine get_array_metadata_chars_C(filename_ascii, fn_len, dims_out, dims_len, ndims, &
                          type_code_out, clen_out) bind(C, name="get_array_metadata_chars_C")
  use iso_c_binding
  use array_utils
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

  ! Locally
  character(len=:), allocatable :: filename
  integer :: unit, magic, type_code, d, i
  integer(c_int), allocatable :: dims(:)

  ! ASCII → String
  call ascii_to_string(filename_ascii, fn_len, filename)

  call get_array_metadata_chars(filename, dims_out, ndims, type_code_out, clen_out)
end subroutine
