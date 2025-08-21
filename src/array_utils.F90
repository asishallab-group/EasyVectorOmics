!> Module for array utilities
module array_utils
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors
    implicit none

    PUBLIC :: get_array_metadata, ascii_to_string, read_file_header, write_file_header, string_to_ascii_arr, check_okay_dims

    integer(int32), parameter :: ARRAY_FILE_MAGIC = int(z'46413230', int32) ! 'FA20' in hex
    !! Magic number for array files

   contains

  !> Check I/O error and set ierr accordingly
  subroutine check_okay_ioerror(ioerror, ierr, msg, unit)
    
    integer(int32), intent(in) :: ioerror
    !! IO error set by fortran
    integer(int32), intent(out) :: ierr
    !! Error code
    integer(int32), intent(in) :: msg
    !! Error code readable version used for setting
    integer(int32), intent(in) :: unit
    !! pass unit allowing it to be closed 

    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, msg)
      close(unit)
      return
    end if
    
  end subroutine

  subroutine check_okay_dims(dims, expected, ierr)
    integer(int32), intent(in) :: dims(:)
    integer(int32), intent(in) :: expected
    integer(int32), intent(out) :: ierr

    if (size(dims) /= expected) then
      call set_err_once(ierr, ERR_DIM_MISMATCH)
      RETURN
    end if

  end subroutine


  !> Opens unit and writes fileheader with all metadata to the given filename
  subroutine write_file_header(filename, unit, type_code, ndim, dims, ierr, clen)
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

    call check_okay_ioerror(ioerror, ierr, ERR_FILE_OPEN, unit)
    if (.not. is_ok(ierr)) return

    write(unit, iostat=ioerror) ARRAY_FILE_MAGIC
    call check_okay_ioerror(ioerror, ierr, ERR_WRITE_MAGIC, unit)
    if (.not. is_ok(ierr)) return

    write(unit, iostat=ioerror) type_code
    call check_okay_ioerror(ioerror, ierr, ERR_WRITE_TYPE, unit)
    if (.not. is_ok(ierr)) return

    write(unit, iostat=ioerror) ndim
    call check_okay_ioerror(ioerror, ierr, ERR_WRITE_NDIMS, unit)
    if (.not. is_ok(ierr)) return

    write(unit, iostat=ioerror) dims
    call check_okay_ioerror(ioerror, ierr, ERR_WRITE_DIMS, unit)
    if (.not. is_ok(ierr)) return

    if (type_code == 3 .and. present(clen)) then
      write(unit, iostat=ioerror) clen
      call check_okay_ioerror(ioerror, ierr, ERR_WRITE_CHARLEN, unit)
      if (.not. is_ok(ierr)) return
    end if
  end subroutine write_file_header

  !> Opens unit and reads file header with all metadata from given file
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


    call set_ok(ioerror)
    call set_ok(ierr)

    open(newunit=unit, file=filename, form='unformatted', access='stream', status='old', iostat=ioerror)
    call check_okay_ioerror(ioerror, ierr, ERR_FILE_OPEN, unit)
    if (.not. is_ok(ierr)) return

    read(unit, iostat=ioerror) magic
    call check_okay_ioerror(ioerror, ierr, ERR_READ_MAGIC, unit)
    if (.not. is_ok(ierr)) return

    if (magic /= ARRAY_FILE_MAGIC) then
      call set_err_once(ierr, ERR_INVALID_FORMAT)
      close(unit)
      return
    end if

    read(unit, iostat=ioerror) type_code
    call check_okay_ioerror(ioerror, ierr, ERR_READ_TYPE, unit)
    if (.not. is_ok(ierr)) return

    read(unit, iostat=ioerror) ndims
    call check_okay_ioerror(ioerror, ierr, ERR_READ_NDIMS, unit)
    if (.not. is_ok(ierr)) return

    allocate(dims(ndims))
    read(unit, iostat=ioerror) dims
    call check_okay_ioerror(ioerror, ierr, ERR_READ_DIMS, unit)
    if (.not. is_ok(ierr)) return

    if(type_code==3) then
      read(unit, iostat=ioerror) clen
      call check_okay_ioerror(ioerror, ierr, ERR_READ_CHARLEN, unit)
      if (.not. is_ok(ierr)) return
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
    if (.not. is_ok(ierr)) return
      
    do i = 1, ndims
      dims_out(i) = dims(i)
    end do
  end subroutine

  !> subroutine to convert an ASCII array to a string
  subroutine ascii_to_string(ascii_array, clen, str)
    implicit none
    integer(int32), intent(in) :: clen
      !! Length of the ASCII array
    integer(int32), intent(in) :: ascii_array(clen)
      !! Array of ASCII characters
    
    character(len=:), allocatable, INTENT(OUT) :: str
    !! Output string
    integer(int32) :: i
    !! loop variable

    allocate(character(len=clen) :: str)
    do i = 1, clen
      str(i:i) = char(ascii_array(i))
    end do
  end subroutine ascii_to_string

  !> converts a string array to an ascii array
  subroutine string_to_ascii_arr(flat, array_size, ascii_arr, clen)
    implicit none
    integer(int32), intent(in) :: array_size
    !! size of the input array
    integer(int32), intent(in) :: clen
    !! length of the longest string
    integer(int32), intent(out) :: ascii_arr(array_size*clen)
    !! ascii output array
    character(len=*), intent(in) :: flat(array_size)
    !! input array with characters
    
    integer(int32) :: i, j

    do i = 1, array_size
      do j = 1, clen
        if (j <= len_trim(flat(i))) then
          ascii_arr((i - 1) * clen + j) = iachar(flat(i)(j:j))
        else
          ascii_arr((i - 1) * clen + j) = 0
        end if
      end do
    end do
  end subroutine string_to_ascii_arr
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

  ! Output
  integer(int32), intent(out) :: dims_out(*)  ! R provides storage
    !! Output array for dimensions
  integer(int32), intent(out) :: ndims        ! Number of dimensions
    !! Output variable for the number of dimensions
  integer(int32), intent(out), optional :: clen
    !! Character length (only for character arrays)

  integer(int32), intent(out) :: ierr
  !! Error code
  character(len=:), allocatable :: filename

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
  integer(c_int), intent(out) :: ierr
    !! Error code
  integer(c_int), intent(out) :: clen
    !! Character length (only for character arrays)

  ! Local variables
  character(len=:), allocatable :: filename
   !! Filename as a string
  

  call ascii_to_string(filename_ascii, fn_len, filename)

  call get_array_metadata(filename, dims_out, ndims, ierr, clen)
end subroutine get_array_metadata_C