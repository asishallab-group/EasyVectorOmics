!> Module for array utilities
module array_utils
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors
    implicit none

    PUBLIC :: get_array_metadata, ascii_to_string, read_file_header, write_file_header

    integer(int32), parameter :: ARRAY_FILE_MAGIC = int(z'46413230', int32) ! 'FA20' in hex
    !! Magic number for array files

   contains

  subroutine write_file_header(filename, unit, type_code, ndim, dims, ierr, clen)
    use iso_fortran_env, only: int32
    implicit none
    character(len=*), intent(in) :: filename
    !! filename to write to
    integer(int32), intent(in) :: type_code
    !! type code of the array (1=int, 2=real, 3=char)
    integer(int32), intent(in) :: ndim
    !! number of dimensions
    integer(int32), intent(in) :: dims(ndim)
    !! dimensions of the array
    integer(int32), intent(in), optional :: clen
    !! character length (only for character arrays)
    integer(int32), intent(inout) :: ierr
    !! error code
    integer(int32), INTENT(OUT) :: unit
    !! Fortran unit number for the file
    integer(int32) :: ioerror
    !! internal I/O error code
    
    call set_ok(ierr)
    call set_ok(ioerror)

    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace', iostat=ioerror)
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_FILE_OPEN)
      close(unit)
      return
    end if

    write(unit, iostat=ioerror) ARRAY_FILE_MAGIC
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_MAGIC)
      close(unit)
      return
    end if

    write(unit, iostat=ioerror) type_code
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_TYPE)
      close(unit)
      return
    end if

    write(unit, iostat=ierr) ndim
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_NDIMS)
      return
    end if

    write(unit, iostat=ierr) dims
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DIMS)
      close(unit)
      return
    end if

    if (type_code == 3 .and. present(clen)) then
      write(unit, iostat=ierr) clen
      if (.not. is_ok(ioerror)) then
        call set_err_once(ierr, ERR_WRITE_CHARLEN)
        return
      end if
    end if
  end subroutine write_file_header

  subroutine read_file_header(filename, unit, type_code, ndims, dims, clen, ierr)
    character(len=*), intent(in) :: filename
    !! filename to read from
    integer(int32), intent(out) :: unit
    !! Fortran unit number for the file
    integer(int32), intent(out) :: type_code
    !! type code of the array (1=int, 2=real, 3=char)
    integer(int32), INTENT(out) :: ndims
    !! number of dimensions
    integer(int32), intent(out) :: clen
    !! character length (only for character arrays)
    integer(int32), intent(out) :: ierr
    !! error code
    integer(int32), allocatable :: dims(:)
    !! dimensions of the array
    integer(int32) :: ioerror
    !! internal I/O error code

    integer(int32) :: magic
    !! magic number to verify file format

    clen = 0

    call set_ok(ioerror)
    call set_ok(ierr)

    open(newunit=unit, file=filename, form='unformatted', access='stream', status='old', iostat=ioerror)
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_FILE_OPEN)
      return
    end if

    read(unit, iostat=ioerror) magic
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_MAGIC)
      close(unit)
      return
    end if

    if (magic /= ARRAY_FILE_MAGIC) then
      call set_err_once(ierr, ERR_INVALID_FORMAT)
      close(unit)
      return
    end if

    read(unit, iostat=ioerror) type_code
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_TYPE)
      close(unit)
      return
    end if

    read(unit, iostat=ioerror) ndims
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_NDIMS)
      close(unit)
      return
    end if

    allocate(dims(ndims))
    read(unit, iostat=ioerror) dims
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DIMS)
      close(unit)
      return
    end if

    if(type_code==3) then
      read(unit, iostat=ioerror) clen
      if (.not. is_ok(ioerror)) then
        call set_err_once(ierr, ERR_READ_CHARLEN)
        close(unit)
        return
      end if
    else
      clen = 0 ! Not applicable for non-character types
    end if

  end subroutine

  !> Get the metadata of an array file
  subroutine get_array_metadata(filename, dims_out, ndims, ierr, clen)
    implicit none

    character(len=*), intent(in) :: filename
    !! Name of the file to read
    integer(int32), intent(out) :: dims_out(*)
    !! Output array for dimensions
    integer(int32), intent(out) :: ndims
    !! Output variable for the number of dimensions
    integer(int32), intent(out) :: ierr
    !! Error code

    integer(int32) :: unit
    !! Fortran unit number for the file
    integer(int32), allocatable :: dims(:)
    !! dimensions
    integer(int32) :: type_code
    !! type code of the array (1=int, 2=real, 3=char)
    integer(int32), INTENT(OUT), OPTIONAL :: clen
    !! character length (only for character arrays)
    integer(int32) :: i
    !! loop variable
    integer(int32) :: ioerror
    !! internal I/O error code

    call set_ok(ioerror)
    ! error handling

    call read_file_header(filename, unit, type_code, ndims, dims, clen, ierr)
    close(unit)
    if (.not. is_ok(ierr)) then
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
  !! Output string
  integer(int32) :: i
  !! loop variable

  allocate(character(len=clen) :: str)
  do i = 1, clen
    str(i:i) = char(ascii_array(i))
  end do
end subroutine ascii_to_string

end module array_utils

!> Subroutine to get the dimensions of an array file
subroutine get_array_metadata_r(filename_ascii, fn_len, dims_out, ndims, ierr, clen)
  use iso_fortran_env, only: int32
  use array_utils, only: ascii_to_string, get_array_metadata
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
  integer(int32), intent(out), optional :: clen
    !! Character length (only for character arrays)

  integer(int32), intent(out) :: ierr
  !! Error code

  call ascii_to_string(filename_ascii, fn_len, filename)

  call get_array_metadata(filename, dims_out, ndims, ierr, clen)

end subroutine get_array_metadata_r

!> C binding for the subroutine to get the dimensions of an array file
subroutine get_array_metadata_C(filename_ascii, fn_len, dims_out, ndims, ierr, clen) bind(C, name="get_array_metadata_C")
  use iso_c_binding, only: c_int
  use iso_fortran_env, only : int32
  use array_utils, only : ascii_to_string, get_array_metadata
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
   !! Filename as a string
  integer(c_int), intent(out) :: ierr
    !! Error code
  integer(c_int), intent(out) :: clen
    !! Character length (only for character arrays)
  
  clen = 0  ! Default value for non-character arrays

  call ascii_to_string(filename_ascii, fn_len, filename)

  call get_array_metadata(filename, dims_out, ndims, ierr, clen)
end subroutine get_array_metadata_C