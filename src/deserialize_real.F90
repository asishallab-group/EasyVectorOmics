!> Module for deserializing real (double precision) arrays from binary files
module real_deserialize_mod
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding, only : c_loc, c_f_pointer
  use array_utils, only: ascii_to_string, read_file_header
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use tox_errors
  implicit none

  private
  public :: deserialize_real_flat, deserialize_real_1d, deserialize_real_2d, &
           deserialize_real_3d, deserialize_real_4d, deserialize_real_5d

contains

  !> Deserialize a flat real array from a file
  subroutine deserialize_real_flat(flat, dims, filename, ierr)
    real(real64), pointer, intent(out) :: flat(:)
    !! Output flat array
    integer(int32), allocatable, intent(out), target :: dims(:)
    !! dimensions array
    character(len=*), intent(in) :: filename
    !! Name of the file to read
    integer(int32), intent(out) :: ierr
    !! Error code
    integer(int32) :: ioerror
    !! Fortran internal error code

    integer(int32) :: unit, magic, type_code, ndims, clen

    call set_ok(ierr)
    call set_ok(ioerror)
    call read_file_header(filename, unit, type_code, ndims, dims, clen, ierr)
    if (.not. is_ok(ierr)) return

    allocate(flat(product(dims)))
    read(unit, iostat=ioerror) flat
    close(unit)
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      deallocate(flat)
      return
    end if
  end subroutine deserialize_real_flat

  !> Deserialize a 1D real array from a file
  !> @note This file just moves a pointer, it exists for consistency
  subroutine deserialize_real_1d(arr, filename, ierr)
    real(real64), pointer, intent(out) :: arr(:)
    !! Output array
    character(len=*), intent(in) :: filename
    !! Name of the file to read

    real(real64), pointer :: flat(:)
    !! Output flat array
    integer(int32), allocatable :: dims(:)
    !! dimensions array
    integer(int32), intent(out) :: ierr
    !! Error code

    call set_ok(ierr)

    call deserialize_real_flat(flat, dims, filename, ierr)
    if(.not. is_ok(ierr)) return
    if (size(dims) /= 1) then
      call set_err_once(ierr, ERR_DIM_MISMATCH)
      RETURN
    end if
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1)])
  end subroutine deserialize_real_1d

  !> Deserialize a 2D real array from a file
  !> @note The array is allocated by the deserialize flat routine
  subroutine deserialize_real_2d(arr, filename, ierr)
    real(real64), pointer, intent(out) :: arr(:,:)
    !! Output array
    character(len=*), intent(in) :: filename
    !! Name of the file to read

    real(real64), pointer :: flat(:)
    !! Output flat array
    integer(int32), allocatable :: dims(:)
    !! dimensions array
    integer(int32), intent(out) :: ierr
    !! Error code

    call set_ok(ierr)

    call deserialize_real_flat(flat, dims, filename, ierr)
    if (.not. is_ok(ierr)) return
    if (size(dims) /= 2) then
      call set_err_once(ierr, ERR_DIM_MISMATCH)
      RETURN
    end if
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1), dims(2)])
  end subroutine deserialize_real_2d

  !> Deserialize a 3D real array from a file
  !> @note The array is allocated by the deserialize flat routine
  subroutine deserialize_real_3d(arr, filename, ierr)
    real(real64), pointer, intent(out) :: arr(:,:,:)
    !! Output array
    character(len=*), intent(in) :: filename
    !! Name of the file to read

    real(real64), pointer :: flat(:)
    !! Output flat array
    integer(int32), allocatable :: dims(:)
    !! dimensions array
    integer(int32), intent(out) :: ierr
    !! Error code

    call set_ok(ierr)

    call deserialize_real_flat(flat, dims, filename, ierr)
    if(.not. is_ok(ierr)) return
    if (size(dims) /= 3) then
      call set_err_once(ierr, ERR_DIM_MISMATCH)
      RETURN
    end if
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1), dims(2), dims(3)])
  end subroutine deserialize_real_3d

  !> Deserialize a 4D real array from a file
  !> @note The array is allocated by the deserialize flat routine
  subroutine deserialize_real_4d(arr, filename, ierr)
    real(real64), pointer, intent(out) :: arr(:,:,:,:)
    !! Output array
    character(len=*), intent(in) :: filename
    !! Name of the file to read

    real(real64), pointer :: flat(:)
    !! Output flat array
    integer(int32), allocatable :: dims(:)
    !! dimensions array
    integer(int32), intent(out) :: ierr
    !! Error code

    call set_ok(ierr)

    call deserialize_real_flat(flat, dims, filename, ierr)
    if(.not. is_ok(ierr)) return
    if (size(dims) /= 4) then
      call set_err_once(ierr, ERR_DIM_MISMATCH)
      RETURN
    end if
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1), dims(2), dims(3), dims(4)])
  end subroutine deserialize_real_4d

  !> Deserialize a 5D real array from a file
  !> @note The array is allocated by the deserialize flat routine
  subroutine deserialize_real_5d(arr, filename, ierr)
    real(real64), pointer, intent(out) :: arr(:,:,:,:,:)
    !! Output array
    character(len=*), intent(in) :: filename
    !! Name of the file to read#
    integer(int32), intent(out) :: ierr
    !! Error code

    real(real64), pointer :: flat(:)
    !! pointer to read into
    integer(int32), allocatable :: dims(:)
    !! dimensions

    call set_ok(ierr)

    call deserialize_real_flat(flat, dims, filename, ierr)
    if(.not. is_ok(ierr)) return
    if (size(dims) /= 5) then
      call set_err_once(ierr, ERR_DIM_MISMATCH)
      RETURN
    end if
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1), dims(2), dims(3), dims(4), dims(5)])
  end subroutine deserialize_real_5d

end module real_deserialize_mod

!> R binding for the subroutine to deserialize a flat real array from a file
!> @note It is assumed that the array is already allocated and passed together with its size
subroutine deserialize_real_flat_r(flat_arr, arr_size, filename_ascii, fn_len, ierr)
  use iso_fortran_env, only: real64, int32
  use array_utils, only : ascii_to_string, read_file_header
  use tox_errors
  implicit none

  real(real64), intent(out) :: flat_arr(arr_size)
  !! array provided by R
  integer(int32), intent(in) :: filename_ascii(fn_len)
  !! filename in ascii
  integer(int32), intent(in) :: fn_len
  !! length of the filename
  integer(int32), intent(in) :: arr_size
  !! size of the array
  integer(int32), intent(out) :: ierr
  !! error code
  integer(int32) :: ioerror
  !! internal fortran error

  character(len=:), allocatable :: filename
  !! filename
  integer(int32), allocatable :: dims(:)
  !! dimensions
  integer :: unit, type_code, ndims, clen

  call set_ok(ierr)
  call set_ok(ioerror)
  call ascii_to_string(filename_ascii, fn_len, filename)

  call read_file_header(filename, unit, type_code, ndims, dims, clen, ierr)
  if(.not. is_ok(ierr)) return

  if (product(dims) /= arr_size) then
    call set_err_once(ierr, ERR_DIM_MISMATCH)
    close(unit)
    return
  end if

  read(unit, iostat=ioerror) flat_arr
  close(unit)
  if (.not. is_ok(ioerror)) then
    call set_err_once(ierr, ERR_READ_DATA)
  end if

end subroutine


!> C binding for the subroutine to deserialize a real array from a file
!> @note It is assumed that the array is already allocated and passed together with its size
subroutine deserialize_real_C(arr, arr_size, filename_ascii, fn_len, ierr) bind(C, name="deserialize_real_C")
    use iso_c_binding, only : c_int, c_double
    use iso_fortran_env, only: int32, real64
    use array_utils, only: ascii_to_string, read_file_header
    use tox_errors
    implicit none

    ! Inputs / Outputs
    real(c_double), intent(out)   :: arr(arr_size)
    !! output array
    integer(c_int), value         :: arr_size  
    !! size of the output array
    integer(c_int), intent(in)    :: filename_ascii(fn_len)
    !! Filename in ascii
    integer(c_int), value         :: fn_len
    !! length of the filename
    integer(c_int), intent(out)   :: ierr
    !! error code
    integer(int32) :: ioerror
    !! internal fortran error

    ! Locals
    character(len=:), allocatable :: filename
    !! filename
    integer(int32), allocatable   :: dims(:)
    !! dimensions
    integer                       :: unit
    integer(int32)                :: type_code, ndims, clen

    ierr = 0

    ! ASCII to String
    call ascii_to_string(filename_ascii, fn_len, filename)

    call read_file_header(filename, unit, type_code, ndims, dims, clen, ierr)
    if(.not. is_ok(ierr)) return

    ! Safety check: ensure provided buffer matches size in file
    if (product(dims) /= arr_size) then
        call set_err_once(ierr, ERR_DIM_MISMATCH)
        close(unit)
        return
    end if

    ! Read directly into provided buffer
    read(unit, iostat=ioerror) arr
    close(unit)
    if (.not. is_ok(ioerror)) then
        call set_err_once(ierr, ERR_READ_DATA)
        return
    end if
end subroutine deserialize_real_C