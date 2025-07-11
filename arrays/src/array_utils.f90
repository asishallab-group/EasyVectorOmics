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



end module array_utils

subroutine get_type_code_r(filename, type_code)
    use array_utils
    implicit none

    character(len=*), intent(in) :: filename
    integer, INTENT(OUT) :: type_code

    type_code = get_type_code(filename)
end subroutine get_type_code_r

subroutine get_array_dims(filename_ascii, fn_len, dims_out, ndims)
  use iso_fortran_env, only: int32
  implicit none

  ! Eingabeparameter
  integer(int32), intent(in) :: filename_ascii(fn_len)
  integer(int32), intent(in) :: fn_len

  ! Ausgabeparameter
  integer(int32), intent(out) :: dims_out(*)  ! R gibt festen Speicher vor
  integer(int32), intent(out) :: ndims        ! Anzahl der Dimensionen

  ! Interne Variablen
  character(len=:), allocatable :: filename
  integer :: unit, magic, type_code, d, i
  integer(int32), allocatable :: dims(:)

  ! ASCII → String
  allocate(character(len=fn_len) :: filename)
  do i = 1, fn_len
    filename(i:i) = char(filename_ascii(i))
  end do

  ! Datei öffnen und Header lesen
  open(newunit=unit, file=filename, form='unformatted', access='stream', status='old')
  read(unit) magic
  read(unit) type_code
  read(unit) d
  allocate(dims(d))
  read(unit) dims
  close(unit)

  ! Rückgabe
  ndims = d
  do i = 1, d
    dims_out(i) = dims(i)
  end do
end subroutine get_array_dims

subroutine get_array_metadata_chars(filename_ascii, fn_len, dims_out, ndims, type_code_out, clen_out)
  use iso_fortran_env, only: int32
  implicit none

  integer(int32), parameter :: ARRAY_FILE_MAGIC = int(z'46413230', int32) ! 'FA20' in hex
  ! Eingabe
  integer(int32), intent(in) :: filename_ascii(fn_len)
  integer(int32), intent(in) :: fn_len

  ! Ausgabe
  integer(int32), intent(out) :: dims_out(*)
  integer(int32), intent(out) :: ndims
  integer(int32), intent(out) :: type_code_out
  integer(int32), intent(out) :: clen_out

  ! Intern
  character(len=:), allocatable :: filename
  integer :: unit, magic, type_code, d, i
  integer(int32), allocatable :: dims(:)

  ! ASCII → String
  allocate(character(len=fn_len) :: filename)
  do i = 1, fn_len
    filename(i:i) = char(filename_ascii(i))
  end do

  ! Datei öffnen
  open(newunit=unit, file=filename, form='unformatted', access='stream', status='old')

  ! Header lesen
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

  ! Rückgabe
  type_code_out = type_code
  ndims = d
  do i = 1, d
    dims_out(i) = dims(i)
  end do
end subroutine
