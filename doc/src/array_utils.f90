!> Module for array utilities
module array_utils
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use iso_c_binding
    implicit none

    public :: get_type_code

    integer(int32), parameter :: ARRAY_FILE_MAGIC = int(z'46413230', int32) ! 'FA20' in hex

   contains
  !> @brief returns type code of the array file
  !> @param filename name of the file to read
  !> @return type code of the array file
  function get_type_code(filename) result(type_code)
    use iso_c_binding
    implicit none

    character(len=*), intent(in) :: filename
    integer :: type_code
    integer :: unit, magic

    open(newunit=unit, file=filename, form='unformatted', access='stream', status='old')
    read(unit) magic
    if (magic /= ARRAY_FILE_MAGIC) error stop "Invalid file format"
    read(unit) type_code
    close(unit)
  end function get_type_code

  !> @brief Get the dimensions of an array file
  !> @param filename Name of the file to read
  !> @return Allocatable array of dimensions
  function get_array_dims(filename) result(dims)
    use iso_c_binding
    implicit none

    character(len=*), intent(in) :: filename
    integer(int32), allocatable :: dims(:)
    integer :: unit, magic, type_code, d

    open(newunit=unit, file=filename, form='unformatted', access='stream', status='old')
    read(unit) magic
    if (magic /= ARRAY_FILE_MAGIC) error stop "Invalid file format"
    read(unit) type_code
    read(unit) d
    allocate(dims(d))
    read(unit) dims
    close(unit)
  end function get_array_dims
end module array_utils

!> @brief Subroutine to get the dimensions of an array file
!> @param filename_ascii Array of ASCII characters representing the filename
!> @param fn_len Length of the filename array
!> @param dims_out Output array for dimensions
!> @param ndims Output variable for the number of dimensions
subroutine get_array_dims_r(filename_ascii, fn_len, dims_out, ndims)
  use iso_fortran_env, only: int32
  implicit none

  ! Input
  integer(int32), intent(in) :: filename_ascii(fn_len)
  integer(int32), intent(in) :: fn_len

  ! Output
  integer(int32), intent(out) :: dims_out(*)  ! R gibt festen Speicher vor
  integer(int32), intent(out) :: ndims        ! Anzahl der Dimensionen

  ! Local variables
  character(len=:), allocatable :: filename
  integer :: unit, magic, type_code, d, i
  integer(int32), allocatable :: dims(:)

  ! ASCII → String
  allocate(character(len=fn_len) :: filename)
  do i = 1, fn_len
    filename(i:i) = char(filename_ascii(i))
  end do

  ! open file and read header
  open(newunit=unit, file=filename, form='unformatted', access='stream', status='old')
  read(unit) magic
  read(unit) type_code
  read(unit) d
  allocate(dims(d))
  read(unit) dims
  close(unit)

  ! return values
  ndims = d
  do i = 1, d
    dims_out(i) = dims(i)
  end do
end subroutine get_array_dims_r

!> @brief C binding for the subroutine to get the dimensions of an array file
!> @param filename_ascii Array of ASCII characters representing the filename
!> @param fn_len Length of the filename array
!> @param dims_out Output array for dimensions
!> @param ndims Output variable for the number of dimensions
subroutine get_array_dims_C(filename_ascii, fn_len, dims_out, ndims) bind(C, name="get_array_dims_C")
  use iso_c_binding
  use iso_fortran_env
  implicit none

  ! Input
  integer(c_int), intent(in) :: filename_ascii(fn_len)
  integer(c_int), value :: fn_len

  ! Output
  integer(c_int), intent(out) :: dims_out(*)
  integer(c_int), intent(out) :: ndims

  ! Local variables
  character(len=:), allocatable :: filename
  integer(c_int) :: unit, magic, type_code, d, i
  integer(c_int), allocatable :: dims(:)


  ! ASCII → String
  allocate(character(len=fn_len) :: filename)
  do i = 1, fn_len
    filename(i:i) = char(filename_ascii(i))
  end do

  ! open file and read header
  open(newunit=unit, file=filename, form='unformatted', access='stream', status='old')
  read(unit) magic
  read(unit) type_code
  read(unit) d
  allocate(dims(d))
  read(unit) dims
  close(unit)

  ! return values
  ndims = d
  do i = 1, d
    dims_out(i) = dims(i)
  end do
end subroutine get_array_dims_C

!> @brief Subroutine to get metadata of a char array file
!> @param filename_ascii Array of ASCII characters representing the filename
!> @param fn_len Length of the filename array
!> @param dims_out Output array for dimensions
!> @param ndims Output variable for the number of dimensions
!> @param type_code_out Output variable for the type code
!> @param clen_out Output variable for the character length
subroutine get_array_metadata_chars_r(filename_ascii, fn_len, dims_out, ndims, type_code_out, clen_out)
  use iso_fortran_env, only: int32
  implicit none

  integer(int32), parameter :: ARRAY_FILE_MAGIC = int(z'46413230', int32) ! 'FA20' in hex
  ! input
  integer(int32), intent(in) :: filename_ascii(fn_len)
  integer(int32), intent(in) :: fn_len

  ! Output
  integer(int32), intent(out) :: dims_out(*)
  integer(int32), intent(out) :: ndims
  integer(int32), intent(out) :: type_code_out
  integer(int32), intent(out) :: clen_out

  ! Internally
  character(len=:), allocatable :: filename
  integer :: unit, magic, type_code, d, i
  integer(int32), allocatable :: dims(:)

  ! ASCII → String
  allocate(character(len=fn_len) :: filename)
  do i = 1, fn_len
    filename(i:i) = char(filename_ascii(i))
  end do

  ! open file
  open(newunit=unit, file=filename, form='unformatted', access='stream', status='old')

  ! read Header
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

  ! return
  type_code_out = type_code
  ndims = d
  do i = 1, d
    dims_out(i) = dims(i)
  end do
end subroutine

!> @brief C binding for the subroutine to get metadata of a char array file
!> @param filename_ascii Array of ASCII characters representing the filename
!> @param fn_len Length of the filename array
!> @param dims_out Output array for dimensions
!> @param dims_len Length of the output dimensions array
!> @param ndims Output variable for the number of dimensions
!> @param type_code_out Output variable for the type code
!> @param clen_out Output variable for the character length
subroutine get_array_metadata_chars_C(filename_ascii, fn_len, dims_out, dims_len, ndims, &
                          type_code_out, clen_out) bind(C, name="get_array_metadata_chars_C")
  use iso_c_binding
  implicit none

  ! Output
  integer(c_int), intent(in) :: filename_ascii(fn_len)
  integer(c_int), value :: fn_len, dims_len

  ! Input
  integer(c_int), intent(out) :: dims_out(dims_len)
  integer(c_int), intent(out) :: ndims
  integer(c_int), intent(out) :: type_code_out
  integer(c_int), intent(out) :: clen_out

  ! Locally
  character(len=:), allocatable :: filename
  integer :: unit, magic, type_code, d, i
  integer(c_int), allocatable :: dims(:)
  integer(c_int), parameter :: ARRAY_FILE_MAGIC = int(z'46413230', c_int) ! 'FA20' in hex

  ! ASCII → String
  allocate(character(len=fn_len) :: filename)
  do i = 1, fn_len
    filename(i:i) = char(filename_ascii(i))
  end do

  ! open unit
  open(newunit=unit, file=filename, form='unformatted', access='stream', status='old')
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

  ! return
  type_code_out = type_code
  ndims = d 

  do i = 1, min(d, dims_len)
    dims_out(i) = dims(i)
  end do
end subroutine
