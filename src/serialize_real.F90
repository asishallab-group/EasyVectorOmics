module serialize_real
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding, only: c_loc
  implicit none

  public:: serialize_real_1d, serialize_real_2d, serialize_real_3d, &
           serialize_real_4d, serialize_real_5d, serialize_real_nd

  integer(int32), parameter :: ARRAY_FILE_MAGIC = int(z'46413230', int32) ! 'FA20' in hex

contains

  !> Writes the header for a real array file.
  subroutine write_real_array_header(unit, type_code, ndim, dims)
    use iso_fortran_env, only: int32
    implicit none
    integer, intent(in) :: unit
      !! Type code for the array (e.g., 2 for real64)
    integer, intent(in) :: type_code, ndim
      !! Number of dimensions
    integer(int32), intent(in) :: dims(ndim)
      !! Dimensions of the array

    write(unit) ARRAY_FILE_MAGIC
    write(unit) type_code
    write(unit) ndim
    write(unit) dims
  end subroutine write_real_array_header

  !> Serialize a 1D real(real64) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  subroutine serialize_real_1d(arr, filename)
    real(real64), intent(in) :: arr(:)
      !! array to save
    character(len=*), intent(in) :: filename
      !! output filename
    integer :: unit
    integer(int32) :: dims(1)
    dims = shape(arr)
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    call write_real_array_header(unit, 2, 1, dims)
    write(unit) arr
    close(unit)
  end subroutine

  !> Serialize a 2D real(real64) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  subroutine serialize_real_2d(arr, filename)
    real(real64), intent(in) :: arr(:,:)
      !! array to save
    character(len=*), intent(in) :: filename
      !! filename
    integer :: unit
    integer(int32) :: dims(2)
    dims = shape(arr)
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    call write_real_array_header(unit, 2, 2, dims)
    write(unit) arr
    close(unit)
  end subroutine

  !> Serialize a 3D real(real64) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  !! @param arr The input real array to serialize.
  !! @param filename The output filename.
  subroutine serialize_real_3d(arr, filename)
    real(real64), intent(in) :: arr(:,:,:)
      !! array to save
    character(len=*), intent(in) :: filename
      !! filename
    integer :: unit
    integer(int32) :: dims(3)
    dims = shape(arr)
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    call write_real_array_header(unit, 2, 3, dims)
    write(unit) arr
    close(unit)
  end subroutine

  !> Serialize a 4D real(real64) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  subroutine serialize_real_4d(arr, filename)
    real(real64), intent(in) :: arr(:,:,:,:)
      !! array to save
    character(len=*), intent(in) :: filename
      !! filename
    integer :: unit
    integer(int32) :: dims(4)
    dims = shape(arr)
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    call write_real_array_header(unit, 2, 4, dims)
    write(unit) arr
    close(unit)
  end subroutine

  !> Serialize a 5D real(real64) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  subroutine serialize_real_5d(arr, filename)
    real(real64), intent(in) :: arr(:,:,:,:,:)
      !! array to save
    character(len=*), intent(in) :: filename
      !! filename
    integer :: unit
    integer(int32) :: dims(5)
    dims = shape(arr)
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    call write_real_array_header(unit, 2, 5, dims)
    write(unit) arr
    close(unit)
  end subroutine

  !> @brief Writes serialized real array from R to file with metdata.
  subroutine serialize_real_nd(arr, dims, ndim, filename)
    real(real64), intent(in) :: arr(:)
      !! array to save
    integer(int32), intent(in) :: dims(:)
      !! Dimensions of the array
    integer(int32), intent(in) :: ndim
      !! Number of dimensions
    character(len=*), intent(in) :: filename
      !! filename
    integer :: unit

    if (size(dims) /= ndim) then
      error stop "Dimension mismatch in serialize_real_nd"
    end if

    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    call write_real_array_header(unit, 2, ndim, dims)
    write(unit) arr
    close(unit)
  end subroutine
 
end module serialize_real

!> Serialize a flat integer array with specified dimensions and number of dimensions to a binary file.
!! R can not pass a multidimensional array directly, so we use a flat array and dimensions. Therefore, exposing serialize_int_*d to R is not needed.
subroutine serialize_real_flat_r(arr, array_size, dims, ndim, filename_ascii, fn_len)
  use iso_fortran_env
  use array_utils
  use serialize_real, only: serialize_real_nd
  implicit none
  real(real64), intent(in) :: arr(array_size) 
    !! Flat real array to serialize
  integer(int32), intent(in) :: dims(ndim)
    !! Dimensions of the array
  integer(int32), intent(in) :: ndim
    !! Number of dimensions
  integer(int32), intent(in) :: array_size
    !! Size of the flat array
  integer(int32), intent(in) :: filename_ascii(fn_len)
    !! Array of ASCII characters representing the filename
  integer(int32), intent(in) :: fn_len
    !! length of the filename array
  character(len=:), allocatable :: filename

  integer :: i, total_len

  call ascii_to_string(filename_ascii, fn_len, filename)

  ! calculate total size
  total_len = 1
  do i = 1, ndim
    total_len = total_len * dims(i)
  end do

  call serialize_real_nd(arr(1:total_len), dims(1:ndim), ndim, filename)
end subroutine

subroutine serialize_real_nd_C(arr, dims, ndim, filename_ascii, fn_len) bind(C, name="serialize_real_nd_C")
  use iso_c_binding
  use array_utils
  use serialize_real
  implicit none

  ! Input parameters
  type(c_ptr), value :: arr
    !! Pointer to the flat real array
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
  real(c_double), pointer :: arr_f(:)
  integer :: i

  call ascii_to_string(filename_ascii, fn_len, filename)

  ! 1D-Array from C pointer
  call c_f_pointer(arr, arr_f, [product(dims(1:ndim))])

  ! save
  call serialize_real_nd(arr_f, dims, ndim, filename)
end subroutine