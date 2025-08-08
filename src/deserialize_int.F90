!> Module for deserializing integer arrays from files
module int_deserialize_mod
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding, only: c_ptr, c_loc, c_char, c_null_char
  use array_utils
  implicit none

  private
  public :: deserialize_int_1d, deserialize_int_2d, &
           deserialize_int_3d, deserialize_int_4d, deserialize_int_5d, deserialize_int_flat

contains
  !> Deserialize a flat integer array from a file
  subroutine deserialize_int_flat(flat, dims, filename)
    use iso_c_binding
    integer(int32), pointer, intent(out) :: flat(:)
    !! Output flat array
    integer(int32), allocatable, intent(out), target :: dims(:)
    !! Output dimensions array
    character(len=*), intent(in) :: filename
    !! Name of the file to read

    integer(int32) :: unit, magic, type_code, ndims, ierr, clen

    ! Read file
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='old', iostat=ierr)
    call check_file_header(filename, unit, type_code, ndims, dims, clen, ierr)
    ! Allocate array of proper size
    allocate(flat(product(dims)))
    read(unit) flat
    close(unit)
  end subroutine deserialize_int_flat

  !> Deserialize a 1D integer array from a file
  !> @note The array must be allocated before calling this subroutine, this is just for consistency since 
  !> a 1D array can not be deserialized to 1D
  subroutine deserialize_int_1d(arr, filename)
    use iso_c_binding
    integer(int32), pointer, intent(out) :: arr(:)
    !! Output array
    character(len=*), intent(in) :: filename
    !! Name of the file to read
    integer(int32), pointer :: flat(:)
    !! Output flat array
    integer(int32), allocatable :: dims(:)
    !! Output dimensions array
    call deserialize_int_flat(flat, dims, filename)
    if (size(dims) /= 1) error stop "Expected 1D array"
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1)])
  end subroutine deserialize_int_1d

  !> Deserialize a 2D integer array from a file
  !> @note The array is allocated by the deserialize flat routine
  subroutine deserialize_int_2d(arr, filename)
    use iso_c_binding
    integer(int32), pointer, intent(out) :: arr(:,:)
    !! Output array
    character(len=*), intent(in) :: filename
    !! Name of the file to read
    integer(int32), pointer :: flat(:)
    !! Output flat array
    integer(int32), allocatable :: dims(:)
    !! Output dimensions array
    call deserialize_int_flat(flat, dims, filename)
    if (size(dims) /= 2) error stop "Expected 2D array"
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1), dims(2)])
  end subroutine deserialize_int_2d

  !> Deserialize a 3D integer array from a file
  !> @note The array is allocated by the deserialize flat routine
  subroutine deserialize_int_3d(arr, filename)
    use iso_c_binding
    integer(int32), pointer, intent(out) :: arr(:,:,:)
    !! Output array
    character(len=*), intent(in) :: filename
    !! Name of the file to read
    integer(int32), pointer :: flat(:)
    !! Output flat array
    integer(int32), allocatable :: dims(:)
    !! Output dimensions array
    call deserialize_int_flat(flat, dims, filename)
    if (size(dims) /= 3) error stop "Expected 3D array"
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1), dims(2), dims(3)])
  end subroutine deserialize_int_3d

  !> Deserialize a 4D integer array from a file
  !> @note The array is allocated by the deserialize flat routine
  subroutine deserialize_int_4d(arr, filename)
    use iso_c_binding
    integer(int32), pointer, intent(out) :: arr(:,:,:,:)
    !! Output array
    character(len=*), intent(in) :: filename
    !! Name of the file to read
    integer(int32), pointer :: flat(:)
    !! Output flat array
    integer(int32), allocatable :: dims(:)
    !! Output dimensions array
    call deserialize_int_flat(flat, dims, filename)
    if (size(dims) /= 4) error stop "Expected 4D array"
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1), dims(2), dims(3), dims(4)])
  end subroutine deserialize_int_4d

  !> Deserialize a 5D integer array from a file
  !> @note The array is allocated by the deserialize flat routine
  subroutine deserialize_int_5d(arr, filename)
    use iso_c_binding
    integer(int32), pointer, intent(out) :: arr(:,:,:,:,:)
    !! Output array
    character(len=*), intent(in) :: filename
    !! Name of the file to read
    integer(int32), pointer :: flat(:)
    !! Output flat array
    integer(int32), allocatable :: dims(:)
    !! Output dimensions array
    call deserialize_int_flat(flat, dims, filename)
    if (size(dims) /= 5) error stop "Expected 5D array"
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1), dims(2), dims(3), dims(4), dims(5)])
  end subroutine deserialize_int_5d

end module int_deserialize_mod

!> R interface for deserializing an integer array from a file
!> @note The output array is handled and preallocated by R
subroutine deserialize_int_r(flat_arr, arr_size, filename_ascii, fn_len)
  use iso_fortran_env, only: int32
  use int_deserialize_mod
  use array_utils
  implicit none

  ! Outputs
  integer(int32), intent(out) :: flat_arr(arr_size)
  !! Output flat array

  integer(int32), intent(in) :: filename_ascii(fn_len)
  !! ASCII representation of the filename
  integer(int32), intent(in) :: fn_len, arr_size
  !! Length of the filename array and size of the output array

  ! Local variables
  character(len=:), allocatable :: filename
  !! Filename as a string
  integer :: i, k
  integer(int32), pointer :: flat(:)
  integer(int32), allocatable, target :: dims(:)

  call ascii_to_string(filename_ascii, fn_len, filename)

  ! Read file
  call deserialize_int_flat(flat, dims, filename)

  do i = 1, product(dims)
    flat_arr(i) = flat(i)
  end do

  ! flat is no longer needed but manually created so it needs to be deallocated
  if (associated(flat)) deallocate(flat)
end subroutine

!> C binding for the subroutine to deserialize an integer array from a file
!>@note It is assumed that the array is already allocated and passed together with its size
subroutine deserialize_int_C(arr, arr_size, filename_ascii, fn_len) bind(C, name="deserialize_int_C")
  use iso_c_binding
  use int_deserialize_mod, only: deserialize_int_flat
  use iso_fortran_env, only: int32
  use array_utils
  implicit none

  integer(c_int), intent(inout) :: arr(arr_size)
  !! Output array, must be preallocated with the correct size
  integer(c_int), value :: arr_size
  !! Size of the output array
  integer(c_int), intent(in) :: filename_ascii(fn_len)
  !! ASCII representation of the filename
  integer(c_int), value :: fn_len
  !! Length of the filename array

  character(len=:), allocatable :: filename
  integer(int32) :: i, ierr

  integer(int32), pointer :: arr_f(:)
  integer(int32), allocatable :: dims(:)

  call ascii_to_string(filename_ascii, fn_len, filename)

  call deserialize_int_flat(arr_f, dims, filename)

  ! safety
  ierr = 0
  if (.not. associated(array_ptr)) then
    print *, "Error: arr_f not allocated"
    ierr = 301
  end if

  if (size(array_ptr) /= arr_size) then
    print *, "Error: Size does not match ", size(array_ptr), arr_size
    ierr = 302
  end if

  if (ierr /= 0) then
    print *, "Error in array pointer check ", ierr
    stop
  end if

  ! Move data in C buffer
  arr(:) = arr_f(:)
end subroutine