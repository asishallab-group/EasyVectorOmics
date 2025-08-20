!> Module for deserializing integer arrays from files
module int_deserialize_mod
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding, only : c_loc, c_f_pointer
  use array_utils, only: ascii_to_string, read_file_header
  use tox_errors
  implicit none

  private
  public :: deserialize_int_1d, deserialize_int_2d, &
           deserialize_int_3d, deserialize_int_4d, deserialize_int_5d, deserialize_int_flat

contains
  !> Deserialize a flat integer array from a file
  subroutine deserialize_int_flat(flat, dims, filename, ierr)
    integer(int32), pointer, intent(out) :: flat(:)
    !! Output flat array
    integer(int32), allocatable, intent(out) :: dims(:)
    !! Output dimensions array
    character(len=*), intent(in) :: filename
    !! Name of the file to read
    INTEGER(int32), INTENT(OUT) :: ierr
    !! Error code

    integer(int32) :: unit
    !! fortran representation of unit number
    integer(int32) :: ioerror
    !! internal error code
    integer(int32) :: magic, type_code, ndims, clen

    call set_ok(ierr)
    call set_ok(ioerror)
    ! Read file
    call read_file_header(filename, unit, type_code, ndims, dims, clen, ierr)
    if (.not. is_ok(ierr)) then
      return
    end if
    allocate(flat(product(dims)))
    read(unit, iostat=ioerror) flat
    close(unit)
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      deallocate(flat)
      return
    end if
  end subroutine deserialize_int_flat

  !> Deserialize a 1D integer array from a file
  !> @note The array is allocated by the deserialize flat routine, this is just for consistency 
  !> since a 1D array can not be deserialized to 1D
  subroutine deserialize_int_1d(arr, filename, ierr)
    integer(int32), pointer, intent(out) :: arr(:)
    !! Output array
    character(len=*), intent(in) :: filename
    !! Name of the file to read
    integer(int32), pointer :: flat(:)
    !! Output flat array
    integer(int32), allocatable :: dims(:)
    !! Output dimensions array
    integer(int32), intent(out) :: ierr
    !! Error code
    call set_ok(ierr)
    call deserialize_int_flat(flat, dims, filename, ierr)
    if (size(dims) /= 1) then
      call set_err_once(ierr, ERR_DIM_MISMATCH)
      RETURN
    end if
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1)])
  end subroutine deserialize_int_1d

  !> Deserialize a 2D integer array from a file
  !> @note The array is allocated by the deserialize flat routine
  subroutine deserialize_int_2d(arr, filename, ierr)
    integer(int32), pointer, intent(out) :: arr(:,:)
    !! Output array
    character(len=*), intent(in) :: filename
    !! Name of the file to read
    integer(int32), pointer :: flat(:)
    !! Output flat array
    integer(int32), allocatable :: dims(:)
    !! Output dimensions array
    integer(int32), intent(out) :: ierr
    !! Error code
    call set_ok(ierr)
    call deserialize_int_flat(flat, dims, filename, ierr)
    if (size(dims) /= 2) then
      call set_err_once(ierr, ERR_DIM_MISMATCH)
      RETURN
    end if
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1), dims(2)])
  end subroutine deserialize_int_2d

  !> Deserialize a 3D integer array from a file
  !> @note The array is allocated by the deserialize flat routine
  subroutine deserialize_int_3d(arr, filename, ierr)
    integer(int32), pointer, intent(out) :: arr(:,:,:)
    !! Output array
    character(len=*), intent(in) :: filename
    !! Name of the file to read
    integer(int32), pointer :: flat(:)
    !! Output flat array
    integer(int32), allocatable :: dims(:)
    !! Output dimensions array
    integer(int32), intent(out) :: ierr
    !! Error code
    call set_ok(ierr)
    call deserialize_int_flat(flat, dims, filename, ierr)
    if (size(dims) /= 3) then
      call set_err_once(ierr, ERR_DIM_MISMATCH)
      RETURN
    end if
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1), dims(2), dims(3)])
  end subroutine deserialize_int_3d

  !> Deserialize a 4D integer array from a file
  !> @note The array is allocated by the deserialize flat routine
  subroutine deserialize_int_4d(arr, filename, ierr)
    integer(int32), pointer, intent(out) :: arr(:,:,:,:)
    !! Output array
    character(len=*), intent(in) :: filename
    !! Name of the file to read
    integer(int32), pointer :: flat(:)
    !! Output flat array
    integer(int32), allocatable :: dims(:)
    !! Output dimensions array
    integer(int32), intent(out) :: ierr
    !! Error code

    call set_ok(ierr)
    call deserialize_int_flat(flat, dims, filename, ierr)
    if (size(dims) /= 4) then
      call set_err_once(ierr, ERR_DIM_MISMATCH)
      RETURN
    end if
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1), dims(2), dims(3), dims(4)])
  end subroutine deserialize_int_4d

  !> Deserialize a 5D integer array from a file
  !> @note The array is allocated by the deserialize flat routine
  subroutine deserialize_int_5d(arr, filename, ierr)
    integer(int32), pointer, intent(out) :: arr(:,:,:,:,:)
    !! Output array
    character(len=*), intent(in) :: filename
    !! Name of the file to read
    integer(int32), pointer :: flat(:)
    !! Output flat array
    integer(int32), allocatable :: dims(:)
    !! Output dimensions array
    integer(int32), intent(out) :: ierr
    !! Error code
    call set_ok(ierr)
    call deserialize_int_flat(flat, dims, filename, ierr)
    if (size(dims) /= 5) then
      call set_err_once(ierr, ERR_DIM_MISMATCH)
      RETURN
    end if
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1), dims(2), dims(3), dims(4), dims(5)])
  end subroutine deserialize_int_5d

end module int_deserialize_mod

!> R interface for deserializing an integer array from a file
!> @note The output array is handled and preallocated by R
subroutine deserialize_int_r(flat_arr, arr_size, filename_ascii, fn_len, ierr)
  use iso_fortran_env, only: int32
  use array_utils, only : ascii_to_string, read_file_header
  use tox_errors
  implicit none

  integer(int32), intent(out) :: flat_arr(arr_size)
  !! array passed by R
  integer(int32), intent(out) :: ierr
  !! error code
  integer(int32), intent(in)  :: filename_ascii(fn_len)
  !! filename to read from
  integer(int32), intent(in)  :: fn_len
  !! length of the filename
  integer(int32), intent(in)  :: arr_size
  !! size of the array
  integer(int32) :: ioerror

  character(len=:), allocatable :: filename
  !! filename in characters
  integer(int32), allocatable   :: dims(:)
  !! dimensions of the array
  integer(int32)                :: unit, type_code, ndims, clen

  call set_ok(ierr)
  call set_ok(ioerror)

  call ascii_to_string(filename_ascii, fn_len, filename)

  call read_file_header(filename, unit, type_code, ndims, dims, clen, ierr)
  if (.not. is_ok(ierr)) return

  if (product(dims) /= arr_size) then
    call set_err_once(ierr, ERR_SIZE_MISMATCH)
    return
  end if

  ! Read directly into R buffer
  read(unit, iostat=ioerror) flat_arr
  close(unit)
  if (.not. is_ok(ioerror)) then
    call set_err_once(ierr, ERR_READ_DATA)
    RETURN
  end if
end subroutine


!> C binding for the subroutine to deserialize an integer array from a file
!>@note It is assumed that the array is already allocated and passed together with its size
subroutine deserialize_int_C(arr, arr_size, filename_ascii, fn_len, ierr) bind(C, name="deserialize_int_C")
    use iso_c_binding, only: c_int
    use iso_fortran_env, only: int32
    use array_utils, only: ascii_to_string, read_file_header
    use tox_errors
    implicit none

    ! Inputs / Outputs
    integer(c_int), intent(out)   :: arr(arr_size)      ! Preallocated buffer from C/Python
    integer(c_int), value         :: arr_size           ! Buffer length
    integer(c_int), intent(in)    :: filename_ascii(fn_len)
    integer(c_int), value         :: fn_len
    integer(c_int), intent(out)   :: ierr

    ! Locals
    character(len=:), allocatable :: filename
    integer(int32), allocatable   :: dims(:)
    integer(int32)                :: unit
    integer(int32)                :: ioerror
    integer(int32)                :: type_code, ndims, clen

    call set_ok(ierr)
    call set_ok(ioerror)

    ! ASCII → String
    call ascii_to_string(filename_ascii, fn_len, filename)

    call read_file_header(filename, unit, type_code, ndims, dims, clen, ierr)
    if (.not. is_ok(ierr)) then
        return
    end if

    ! Safety check: ensure provided buffer matches size in file
    if (product(dims) /= arr_size) then
        call set_err_once(ierr, ERR_SIZE_MISMATCH)
        close(unit)
        return
    end if

    ! Read directly into C/Python-provided buffer → ZERO COPY
    read(unit, iostat=ioerror) arr
    close(unit)
    if (.not. is_ok(ioerror)) then
        call set_err_once(ierr, ERR_READ_DATA)
        return
    end if
end subroutine deserialize_int_C