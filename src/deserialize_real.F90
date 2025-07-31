!> Module for deserializing real (double precision) arrays from binary files
module real_deserialize_mod
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding
  use, intrinsic :: iso_fortran_env, only: real64, int32
  implicit none

  private
  public :: deserialize_real_flat, deserialize_real_1d, deserialize_real_2d, &
           deserialize_real_3d, deserialize_real_4d, deserialize_real_5d

  integer(int32), parameter :: ARRAY_FILE_MAGIC = int(z'46413230', int32) ! 'FA20' in hex

contains

  !> @brief Deserialize a flat real array from a file
  !> @param flat Pointer to the output flat array
  !> @param dims Output array for dimensions
  !> @param filename Name of the file to read
  subroutine deserialize_real_flat(flat, dims, filename)
    use iso_c_binding
    real(real64), pointer, intent(out) :: flat(:)
    integer(int32), allocatable, intent(out), target :: dims(:)
    character(len=*), intent(in) :: filename

    integer :: unit, magic, type_code, d

    open(newunit=unit, file=filename, form='unformatted', access='stream', status='old')
    read(unit) magic
    if (magic /= ARRAY_FILE_MAGIC) error stop "Invalid file format"
    read(unit) type_code
    if (type_code /= 2) error stop "Expected real64 data"
    read(unit) d
    allocate(dims(d))
    read(unit) dims
    allocate(flat(product(dims)))
    read(unit) flat
    close(unit)
  end subroutine deserialize_real_flat

  !> @brief Deserialize a 1D real array from a file
  !> @param arr Pointer to the output array
  !> @param filename Name of the file to read
  !> @note This file just moves a pointer, it exists for consistency
  subroutine deserialize_real_1d(arr, filename)
    use iso_c_binding
    real(real64), pointer, intent(out) :: arr(:)
    character(len=*), intent(in) :: filename

    real(real64), pointer :: flat(:)
    integer(int32), allocatable :: dims(:)

    call deserialize_real_flat(flat, dims, filename)
    if (size(dims) /= 1) error stop "Expected 1D array"
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1)])
  end subroutine deserialize_real_1d

  !> @brief Deserialize a 2D real array from a file
  !> @param arr Pointer to the output array
  !> @param filename Name of the file to read
  !> @note The array is allocated by the deserialize flat routine
  subroutine deserialize_real_2d(arr, filename)
    use iso_c_binding
    real(real64), pointer, intent(out) :: arr(:,:)
    character(len=*), intent(in) :: filename

    real(real64), pointer :: flat(:)
    integer(int32), allocatable :: dims(:)

    call deserialize_real_flat(flat, dims, filename)
    if (size(dims) /= 2) error stop "Expected 2D array"
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1), dims(2)])
  end subroutine deserialize_real_2d

  !> @brief Deserialize a 3D real array from a file
  !> @param arr Pointer to the output array
  !> @param filename Name of the file to read
  !> @note The array is allocated by the deserialize flat routine
  subroutine deserialize_real_3d(arr, filename)
    use iso_c_binding
    real(real64), pointer, intent(out) :: arr(:,:,:)
    character(len=*), intent(in) :: filename

    real(real64), pointer :: flat(:)
    integer(int32), allocatable :: dims(:)

    call deserialize_real_flat(flat, dims, filename)
    if (size(dims) /= 3) error stop "Expected 3D array"
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1), dims(2), dims(3)])
  end subroutine deserialize_real_3d

  !> @brief Deserialize a 4D real array from a file
  !> @param arr Pointer to the output array
  !> @param filename Name of the file to read
  !> @note The array is allocated by the deserialize flat routine
  subroutine deserialize_real_4d(arr, filename)
    use iso_c_binding
    real(real64), pointer, intent(out) :: arr(:,:,:,:)
    character(len=*), intent(in) :: filename

    real(real64), pointer :: flat(:)
    integer(int32), allocatable :: dims(:)

    call deserialize_real_flat(flat, dims, filename)
    if (size(dims) /= 4) error stop "Expected 4D array"
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1), dims(2), dims(3), dims(4)])
  end subroutine deserialize_real_4d

  !> @brief Deserialize a 5D real array from a file
  !> @param arr Pointer to the output array
  !> @param filename Name of the file to read
  !> @note The array is allocated by the deserialize flat routine
  subroutine deserialize_real_5d(arr, filename)
    use iso_c_binding
    real(real64), pointer, intent(out) :: arr(:,:,:,:,:)
    character(len=*), intent(in) :: filename

    real(real64), pointer :: flat(:)
    integer(int32), allocatable :: dims(:)

    call deserialize_real_flat(flat, dims, filename)
    if (size(dims) /= 5) error stop "Expected 5D array"
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1), dims(2), dims(3), dims(4), dims(5)])
  end subroutine deserialize_real_5d

end module real_deserialize_mod

!> @brief R binding for the subroutine to deserialize a flat real array from a file
!> @param flat_arr Output flat array
!> @param dims_out Output dimensions array
!> @param ndim_out Number of dimensions
!> @param filename_ascii Array of ASCII characters representing the filename
!> @param fn_len Length of the filename array
subroutine deserialize_real_flat_r(flat_arr, arr_size, dims_out, ndim_out, filename_ascii, fn_len, ndim_actual)
  use iso_fortran_env
  use real_deserialize_mod
  implicit none

  ! This needs fixed size and pass the size as parameter
  real(real64), intent(out) :: flat_arr(arr_size)
  integer(int32), intent(out) :: dims_out(ndim_actual)
  integer(int32), intent(out) :: ndim_out

  ! Filename
  integer(int32), intent(in) :: filename_ascii(fn_len)
  integer(int32), intent(in) :: fn_len, arr_size, ndim_actual

  ! Local
  character(len=:), allocatable :: filename
  integer(int32) :: i, k
  real(real64), pointer :: flat(:)
  integer(int32), allocatable, target :: dims(:)

  allocate(character(len=fn_len) :: filename)
  do i = 1, fn_len
    filename(i:i) = char(filename_ascii(i))
  end do

  ! Read file
  call deserialize_real_flat(flat, dims, filename)

  ndim_out = size(dims)
  do i = 1, ndim_out
    dims_out(i) = dims(i)
  end do

  do i = 1, product(dims)
    flat_arr(i) = flat(i)
  end do

  !Flat is no longer needed but created manually so it needs to be deallocated
  if (associated(flat)) deallocate(flat)
end subroutine

!> @brief C binding for the subroutine to deserialize a real array from a file
!> @param arr Output array
!> @param arr_size Size of the output array
!> @param filename_ascii Array of ASCII characters representing the filename
!> @param fn_len Length of the filename array
!> @note It is assumed that the array is already allocated and passed together with its size
subroutine deserialize_real_C(arr, arr_size, filename_ascii, fn_len) bind(C, name="deserialize_real_C")
  use iso_c_binding
  use real_deserialize_mod, only: deserialize_real_flat
  use iso_fortran_env
  implicit none

  real(c_double), intent(inout) :: arr(arr_size)
  integer(c_int), value :: arr_size
  integer(c_int), intent(in) :: filename_ascii(fn_len)
  integer(c_int), value :: fn_len

  character(len=:), allocatable :: filename
  integer :: i

  ! arr_f is a pointer to the Fortran array
  real(real64), pointer :: arr_f(:)
  integer(int32), allocatable :: dims(:)

  ! ASCII → String
  allocate(character(len=fn_len) :: filename)
  do i = 1, fn_len
    filename(i:i) = char(filename_ascii(i))
  end do

  ! Read data
  call deserialize_real_flat(arr_f, dims, filename)

  ! Checks
  if (.not. associated(arr_f)) then
      print *, "Error: arr_f not allocated"
      stop 1
  end if

  if (size(arr_f) /= arr_size) then
      print *, "Error: Size does not match ", size(arr_f), arr_size
      stop 2
  end if

  ! Move to buffer
  arr(:) = arr_f(:)
end subroutine