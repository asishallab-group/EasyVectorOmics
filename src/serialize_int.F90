!> Module for serializing integer arrays to binary files.
module serialize_int
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding, only: c_loc
  implicit none

  public:: serialize_int_1d, serialize_int_2d, serialize_int_3d, &
           serialize_int_4d, serialize_int_5d, serialize_int_nd

  integer(int32), parameter :: ARRAY_FILE_MAGIC = int(z'46413230', int32) ! 'FA20' in hex

contains

  subroutine write_int_array(unit, type_code, ndim, dims)
    use iso_fortran_env, only: int32
    implicit none
    integer, intent(in) :: unit
    integer, intent(in) :: type_code, ndim
    integer(int32), intent(in) :: dims(ndim)

    write(unit) ARRAY_FILE_MAGIC
    write(unit) type_code
    write(unit) ndim
    write(unit) dims
  end subroutine write_int_array
  !> Serialize a 1D integer(int32) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  subroutine serialize_int_1d(arr, filename)
    integer(int32), intent(in) :: arr(:)
    !! array to save
    character(len=*), intent(in) :: filename
    !! output filename
    integer :: unit
    integer(int32) :: dims(1)
    dims = shape(arr)
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    call write_int_array(unit, 1, 1, dims)
    write(unit) arr
    close(unit)
  end subroutine

  !> Serialize a 2D integer(int32) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  subroutine serialize_int_2d(arr, filename)
    integer(int32), intent(in) :: arr(:,:)
    !! array to save
    character(len=*), intent(in) :: filename
    !! output filename
    integer :: unit
    integer(int32) :: dims(2)
    dims = shape(arr)
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    call write_int_array(unit, 1, 2, dims)
    write(unit) arr
    close(unit)
  end subroutine

  !> Serialize a 3D integer(int32) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  subroutine serialize_int_3d(arr, filename)
    integer(int32), intent(in) :: arr(:,:,:)
    !! array to save
    character(len=*), intent(in) :: filename
    !! output filename
    integer :: unit
    integer(int32) :: dims(3)
    dims = shape(arr)
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    call write_int_array(unit, 1, 3, dims)
    write(unit) arr
    close(unit)
  end subroutine

  !> Serialize a 4D integer(int32) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  subroutine serialize_int_4d(arr, filename)
    integer(int32), intent(in) :: arr(:,:,:,:)
    !! array to save
    character(len=*), intent(in) :: filename
    !! output filename
    integer :: unit
    integer(int32) :: dims(4)
    dims = shape(arr)
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    call write_int_array(unit, 1, 4, dims)
    write(unit) arr
    close(unit)
  end subroutine

  !> Serialize a 5D integer(int32) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  subroutine serialize_int_5d(arr, filename)
    integer(int32), intent(in) :: arr(:,:,:,:,:)
    !! array to save
    character(len=*), intent(in) :: filename
    !! output filename
    integer :: unit
    integer(int32) :: dims(5)
    dims = shape(arr)
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    call write_int_array(unit, 1, 5, dims)
    write(unit) arr
    close(unit)
  end subroutine

  !> Serialize a flat integer array with specified dimensions and number of dimensions to a binary file.
  !> @note this is called by R
  subroutine serialize_int_nd(arr, dims, ndim, filename)
    integer(int32), intent(in) :: arr(:)
    !! Flat integer array to serialize
    integer(int32), intent(in) :: dims(:)
    !! Dimensions of the array
    integer(int32), intent(in) :: ndim
    !! Number of dimensions
    character(len=*), intent(in) :: filename
    !! filename
    integer :: unit

    if (size(dims) /= ndim) then
      error stop "Dimension mismatch in serialize_int_nd"
    end if

    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    call write_int_array(unit, 1, ndim, dims)
    write(unit) arr
    close(unit)
  end subroutine
end module serialize_int


!> Serialize a flat integer array with specified dimensions and number of dimensions to a binary file.
!! R can not pass a multidimensional array directly, so we use a flat array and dimensions. Therefore, exposing serialize_int_*d to R is not needed.
subroutine serialize_int_flat_r(arr, array_size, dims, ndim, filename_ascii, fn_len)
  use iso_fortran_env, only: int32
  use array_utils
  use serialize_int, only: serialize_int_nd
  implicit none

  integer(int32), intent(in) :: arr(array_size)
  !! Flat integer array to serialize
  integer(int32), intent(in) :: dims(ndim)
  !! Dimensions of the array
  integer(int32), intent(in) :: ndim 
  !! Number of dimensions
  integer(int32), intent(in) :: array_size
  !! Size of the flat array
  integer(int32), intent(in) :: filename_ascii(fn_len)
  !! Array of ASCII characters representing the filename
  integer(int32), intent(in) :: fn_len
  !! Length of the filename array

  character(len=:), allocatable :: filename
  integer :: i, total_len

  call ascii_to_string(filename_ascii, fn_len, filename)

  ! Gesamtgröße berechnen (z. B. für Sicherheit oder Logging)
  total_len = 1
  do i = 1, ndim
    total_len = total_len * dims(i)
  end do

  ! optional check
  ! print *, "Serializing array with ", total_len, " elements and ", ndim, " dimensions."

  call serialize_int_nd(arr(1:total_len), dims(1:ndim), ndim, filename)
end subroutine

!> C binding for the subroutine to serialize a flat integer array to a binary file.
subroutine serialize_int_nd_C(arr, dims, ndim, filename_ascii, fn_len) bind(C, name="serialize_int_nd_C")
  use iso_c_binding
  use array_utils
  use serialize_int
  implicit none

  ! input
  type(c_ptr), value :: arr
    !! Pointer to the flat integer array
  integer(c_int), intent(in) :: dims(ndim)
    !! Dimensions of the array
  integer(c_int), value :: ndim
    !! Number of dimensions
  integer(c_int), intent(in) :: filename_ascii(fn_len)
    !! Array of ASCII characters representing the filename
  integer(c_int), value :: fn_len
    !! Length of the filename array

  ! Local
  character(len=:), allocatable :: filename
  integer(c_int), pointer :: arr_f(:)
  integer :: i

  call ascii_to_string(filename_ascii, fn_len, filename)
  ! 1D-Array to Fortran pointer
  call c_f_pointer(arr, arr_f, [product(dims(1:ndim))])

  ! save
  call serialize_int_nd(arr_f, dims, ndim, filename)
end subroutine
