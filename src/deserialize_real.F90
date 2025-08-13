!> Module for deserializing real (double precision) arrays from binary files
module real_deserialize_mod
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding, only : c_loc, c_f_pointer
  use array_utils, only: ascii_to_string, read_file_header
  use, intrinsic :: iso_fortran_env, only: real64, int32
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

    integer(int32) :: unit, magic, type_code, ndims, clen

    call read_file_header(filename, unit, type_code, ndims, dims, clen, ierr)
    if (ierr /= 0) then
      return
    end if
    allocate(flat(product(dims)))
    read(unit, iostat=ierr) flat
    close(unit)
    if (ierr /= 0) then
      ierr = 107
      return
    end if
  end subroutine deserialize_real_flat

  !> Deserialize a 1D real array from a file
  !> @note This file just moves a pointer, it exists for consistency
  subroutine deserialize_real_1d(arr, filename)
    real(real64), pointer, intent(out) :: arr(:)
    !! Output array
    character(len=*), intent(in) :: filename
    !! Name of the file to read

    real(real64), pointer :: flat(:)
    !! Output flat array
    integer(int32), allocatable :: dims(:)
    !! dimensions array
    integer(int32) :: ierr
    !! Error code

    call deserialize_real_flat(flat, dims, filename, ierr)
    if (size(dims) /= 1) error stop "Expected 1D array"
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1)])
  end subroutine deserialize_real_1d

  !> Deserialize a 2D real array from a file
  !> @note The array is allocated by the deserialize flat routine
  subroutine deserialize_real_2d(arr, filename)
    real(real64), pointer, intent(out) :: arr(:,:)
    !! Output array
    character(len=*), intent(in) :: filename
    !! Name of the file to read

    real(real64), pointer :: flat(:)
    !! Output flat array
    integer(int32), allocatable :: dims(:)
    !! dimensions array
    integer(int32) :: ierr
    !! Error code

    call deserialize_real_flat(flat, dims, filename, ierr)
    if (size(dims) /= 2) error stop "Expected 2D array"
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1), dims(2)])
  end subroutine deserialize_real_2d

  !> Deserialize a 3D real array from a file
  !> @note The array is allocated by the deserialize flat routine
  subroutine deserialize_real_3d(arr, filename)
    real(real64), pointer, intent(out) :: arr(:,:,:)
    !! Output array
    character(len=*), intent(in) :: filename
    !! Name of the file to read

    real(real64), pointer :: flat(:)
    !! Output flat array
    integer(int32), allocatable :: dims(:)
    !! dimensions array
    integer(int32) :: ierr
    !! Error code

    call deserialize_real_flat(flat, dims, filename, ierr)
    if (size(dims) /= 3) error stop "Expected 3D array"
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1), dims(2), dims(3)])
  end subroutine deserialize_real_3d

  !> Deserialize a 4D real array from a file
  !> @note The array is allocated by the deserialize flat routine
  subroutine deserialize_real_4d(arr, filename)
    real(real64), pointer, intent(out) :: arr(:,:,:,:)
    !! Output array
    character(len=*), intent(in) :: filename
    !! Name of the file to read

    real(real64), pointer :: flat(:)
    !! Output flat array
    integer(int32), allocatable :: dims(:)
    !! dimensions array
    integer(int32) :: ierr
    !! Error code

    call deserialize_real_flat(flat, dims, filename, ierr)
    if (size(dims) /= 4) error stop "Expected 4D array"
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1), dims(2), dims(3), dims(4)])
  end subroutine deserialize_real_4d

  !> Deserialize a 5D real array from a file
  !> @note The array is allocated by the deserialize flat routine
  subroutine deserialize_real_5d(arr, filename)
    real(real64), pointer, intent(out) :: arr(:,:,:,:,:)
    !! Output array
    character(len=*), intent(in) :: filename
    !! Name of the file to read#
    integer(int32) :: ierr
    !! Error code

    real(real64), pointer :: flat(:)
    integer(int32), allocatable :: dims(:)

    call deserialize_real_flat(flat, dims, filename, ierr)
    if (size(dims) /= 5) error stop "Expected 5D array"
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1), dims(2), dims(3), dims(4), dims(5)])
  end subroutine deserialize_real_5d

end module real_deserialize_mod

!> R binding for the subroutine to deserialize a flat real array from a file
!> @note It is assumed that the array is already allocated and passed together with its size
subroutine deserialize_real_flat_r(flat_arr, arr_size, filename_ascii, fn_len, ierr)
  use iso_fortran_env, only: real64, int32
  use array_utils, only : ascii_to_string, read_file_header
  implicit none

  real(real64), intent(out) :: flat_arr(arr_size)
  integer(int32), intent(in) :: filename_ascii(fn_len)
  integer(int32), intent(in) :: fn_len, arr_size
  integer(int32), intent(out) :: ierr

  character(len=:), allocatable :: filename
  integer(int32), allocatable :: dims(:)
  integer :: unit, type_code, ndims, clen

  call ascii_to_string(filename_ascii, fn_len, filename)

  call read_file_header(filename, unit, type_code, ndims, dims, clen, ierr)
  if (ierr /= 0) then
    close(unit)
    return
  end if

  if (product(dims) /= arr_size) then
    ierr = 201
    close(unit)
    return
  end if

  read(unit, iostat=ierr) flat_arr
  close(unit)
  if (ierr /= 0) ierr = 107

end subroutine


!> C binding for the subroutine to deserialize a real array from a file
!> @note It is assumed that the array is already allocated and passed together with its size
subroutine deserialize_real_C(arr, arr_size, filename_ascii, fn_len, ierr) bind(C, name="deserialize_real_C")
    use iso_c_binding, only : c_int, c_double
    use iso_fortran_env, only: int32, real64
    use array_utils, only: ascii_to_string, read_file_header
    implicit none

    ! Inputs / Outputs
    real(c_double), intent(inout) :: arr(arr_size)   ! Preallocated buffer from C/Python
    integer(c_int), value         :: arr_size        ! Buffer length
    integer(c_int), intent(in)    :: filename_ascii(fn_len)
    integer(c_int), value         :: fn_len
    integer(c_int), intent(out)   :: ierr

    ! Locals
    character(len=:), allocatable :: filename
    integer(int32), allocatable   :: dims(:)
    integer                       :: unit
    integer(int32)                :: type_code, ndims, clen

    ierr = 0

    ! ASCII → String
    call ascii_to_string(filename_ascii, fn_len, filename)

    call read_file_header(filename, unit, type_code, ndims, dims, clen, ierr)
    if (ierr /= 0) then
        close(unit)
        return
    end if

    ! Safety check: ensure provided buffer matches size in file
    if (product(dims) /= arr_size) then
        ierr = 302
        close(unit)
        return
    end if

    ! Read directly into provided buffer → ZERO COPY
    read(unit, iostat=ierr) arr
    close(unit)
    if (ierr /= 0) then
        ierr = 107
        return
    end if
end subroutine deserialize_real_C