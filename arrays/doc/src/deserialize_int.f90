!> Module for deserializing integer arrays from files
module int_deserialize_mod
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding, only: c_ptr, c_loc, c_char, c_null_char
  implicit none

  private
  public :: deserialize_int_1d, deserialize_int_2d, &
           deserialize_int_3d, deserialize_int_4d, deserialize_int_5d, deserialize_int_flat

  integer(int32), parameter :: ARRAY_FILE_MAGIC = int(z'46413230', int32) ! 'FA20' in hex

contains
  !> @brief Deserialize a flat integer array from a file
  !> @param flat Pointer to the output flat array
  !> @param dims Output array for dimensions
  !> @param filename Name of the file to read
  subroutine deserialize_int_flat(flat, dims, filename)
    use iso_c_binding
    integer(int32), pointer, intent(out) :: flat(:)
    integer(int32), allocatable, intent(out), target :: dims(:)
    character(len=*), intent(in) :: filename

    integer :: unit, magic, type_code, d

    ! Read file
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='old')
    read(unit) magic
    if (magic /= ARRAY_FILE_MAGIC) error stop "Invalid file format"
    read(unit) type_code
    if (type_code /= 1) error stop "Expected int32 data"
    read(unit) d
    allocate(dims(d))
    read(unit) dims
    ! Allocate array of proper size
    allocate(flat(product(dims)))
    read(unit) flat
    close(unit)
  end subroutine deserialize_int_flat

  !> @brief Deserialize a 1D integer array from a file
  !> @param arr Pointer to the output array
  !> @param filename Name of the file to read
  !> @note The array must be allocated before calling this subroutine, this is just for consistency since 
  !> a 1D array can not be deserialized to 1D
  subroutine deserialize_int_1d(arr, filename)
    use iso_c_binding
    integer(int32), pointer, intent(out) :: arr(:)
    character(len=*), intent(in) :: filename
    integer(int32), pointer :: flat(:)
    integer(int32), allocatable :: dims(:)
    call deserialize_int_flat(flat, dims, filename)
    if (size(dims) /= 1) error stop "Expected 1D array"
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1)])
  end subroutine deserialize_int_1d

  !> @brief Deserialize a 2D integer array from a file
  !> @param arr Pointer to the output array
  !> @param filename Name of the file to read
  !> @note The array is allocated by the deserialize flat routine
  subroutine deserialize_int_2d(arr, filename)
    use iso_c_binding
    integer(int32), pointer, intent(out) :: arr(:,:)
    character(len=*), intent(in) :: filename
    integer(int32), pointer :: flat(:)
    integer(int32), allocatable :: dims(:)
    call deserialize_int_flat(flat, dims, filename)
    if (size(dims) /= 2) error stop "Expected 2D array"
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1), dims(2)])
  end subroutine deserialize_int_2d

  !> @brief Deserialize a 3D integer array from a file
  !> @param arr Pointer to the output array
  !> @param filename Name of the file to read
  !> @note The array is allocated by the deserialize flat routine
  subroutine deserialize_int_3d(arr, filename)
    use iso_c_binding
    integer(int32), pointer, intent(out) :: arr(:,:,:)
    character(len=*), intent(in) :: filename
    integer(int32), pointer :: flat(:)
    integer(int32), allocatable :: dims(:)
    call deserialize_int_flat(flat, dims, filename)
    if (size(dims) /= 3) error stop "Expected 3D array"
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1), dims(2), dims(3)])
  end subroutine deserialize_int_3d

  !> @brief Deserialize a 4D integer array from a file
  !> @param arr Pointer to the output array
  !> @param filename Name of the file to read
  !> @note The array is allocated by the deserialize flat routine
  subroutine deserialize_int_4d(arr, filename)
    use iso_c_binding
    integer(int32), pointer, intent(out) :: arr(:,:,:,:)
    character(len=*), intent(in) :: filename
    integer(int32), pointer :: flat(:)
    integer(int32), allocatable :: dims(:)
    call deserialize_int_flat(flat, dims, filename)
    if (size(dims) /= 4) error stop "Expected 4D array"
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1), dims(2), dims(3), dims(4)])
  end subroutine deserialize_int_4d

  !> @brief Deserialize a 5D integer array from a file
  !> @param arr Pointer to the output array
  !> @param filename Name of the file to read
  !> @note The array is allocated by the deserialize flat routine
  subroutine deserialize_int_5d(arr, filename)
    use iso_c_binding
    integer(int32), pointer, intent(out) :: arr(:,:,:,:,:)
    character(len=*), intent(in) :: filename
    integer(int32), pointer :: flat(:)
    integer(int32), allocatable :: dims(:)
    call deserialize_int_flat(flat, dims, filename)
    if (size(dims) /= 5) error stop "Expected 5D array"
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1), dims(2), dims(3), dims(4), dims(5)])
  end subroutine deserialize_int_5d

end module int_deserialize_mod

!> @brief R interface for deserializing an integer array from a file
!> @param flat_arr Output flat array
!> @param dims_out Output dimensions array
!> @param ndim_out Output number of dimensions
!> @param filename_ascii ASCII representation of the filename
!> @param fn_len Length of the filename array
!> @note The output array is handled and preallocated by R
subroutine deserialize_int_r(flat_arr, dims_out, ndim_out, filename_ascii, fn_len)
  use iso_fortran_env, only: int32
  use int_deserialize_mod
  implicit none

  ! Outputs
  integer(int32), intent(out) :: flat_arr(*)
  integer(int32), intent(out) :: dims_out(*)
  integer, intent(out) :: ndim_out

  ! Filename
  integer(int32), intent(in) :: filename_ascii(fn_len)
  integer, intent(in) :: fn_len

  ! Local variables
  character(len=:), allocatable :: filename
  integer :: i, k
  integer(int32), pointer :: flat(:)
  integer(int32), allocatable, target :: dims(:)

  allocate(character(len=fn_len) :: filename)
  do i = 1, fn_len
    filename(i:i) = char(filename_ascii(i))
  end do

  ! Read file
  call deserialize_int_flat(flat, dims, filename)

  ndim_out = size(dims)
  do i = 1, ndim_out
    dims_out(i) = dims(i)
  end do

  do i = 1, product(dims)
    flat_arr(i) = flat(i)
  end do

  ! flat is no longer needed but manually created so it needs to be deallocated
  if (associated(flat)) deallocate(flat)
end subroutine

!>@brief C binding for the subroutine to deserialize an integer array from a file
!>@param arr Output array
!>@param arr_size Size of the output array
!>@param filename_ascii ASCII representation of the filename
!>@param fn_len Length of the filename array
!>@note It is assumed that the array is already allocated and passed together with its size
subroutine deserialize_int_C(arr, arr_size, filename_ascii, fn_len) bind(C, name="deserialize_int_C")
  use iso_c_binding
  use int_deserialize_mod, only: deserialize_int_flat
  use iso_fortran_env, only: int32
  implicit none

  integer(c_int), intent(inout) :: arr(arr_size)
  integer(c_int), value :: arr_size
  integer(c_int), intent(in) :: filename_ascii(fn_len)
  integer(c_int), value :: fn_len

  character(len=:), allocatable :: filename
  integer :: i

  ! arr_f ist das Ergebnis von deserialize_int_flat – muss kopiert werden
  integer(int32), pointer :: arr_f(:)
  integer(int32), allocatable :: dims(:)

  ! ASCII → String
  allocate(character(len=fn_len) :: filename)
  do i = 1, fn_len
    filename(i:i) = char(filename_ascii(i))
  end do

  ! Daten einlesen – arr_f wird intern allokiert
  call deserialize_int_flat(arr_f, dims, filename)

  ! Sicherheitschecks
  if (.not. associated(arr_f)) then
      print *, "Error: arr_f not allocated"
      stop 1
  end if

  if (size(arr_f) /= arr_size) then
      print *, "Error: Size does not match ", size(arr_f), arr_size
      stop 2
  end if

  ! Move data in C buffer
  arr(:) = arr_f(:)
end subroutine