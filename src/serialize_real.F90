module serialize_real
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding, only: c_loc
  use array_utils, only: write_file_header
  implicit none

  public:: serialize_real_1d, serialize_real_2d, serialize_real_3d, &
           serialize_real_4d, serialize_real_5d, serialize_real_nd

  integer(int32), parameter :: ARRAY_FILE_MAGIC = int(z'46413230', int32) ! 'FA20' in hex

contains

  !> Serialize a 1D real(real64) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  subroutine serialize_real_1d(arr, filename)
    real(real64), intent(in) :: arr(:)
      !! array to save
    character(len=*), intent(in) :: filename
      !! output filename
    integer(int32) :: unit
    integer(int32) :: dims(1)
    integer(int32) :: ierr
    dims = shape(arr)
    ierr = 0
    call write_file_header(filename, unit, 2, 1, dims, ierr)
    if (ierr /= 0) then
      return
    end if
    write(unit, iostat=ierr) arr
    if (ierr /= 0) then
      ierr = 405
    end if
    close(unit)
  end subroutine

  !> Serialize a 2D real(real64) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  subroutine serialize_real_2d(arr, filename)
    real(real64), intent(in) :: arr(:,:)
      !! array to save
    character(len=*), intent(in) :: filename
      !! filename
    integer(int32) :: unit
    integer(int32) :: dims(2)
    integer(int32) :: ierr
    ierr = 0
    dims = shape(arr)
    call write_file_header(filename, unit, 2, 2, dims, ierr)
    write(unit, iostat=ierr) arr
    if (ierr /= 0) then
      ierr = 405
    end if
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
    integer(int32) :: ierr
    ierr = 0
    dims = shape(arr)
    call write_file_header(filename, unit, 2, 3, dims, ierr)
    write(unit, iostat=ierr) arr
    if (ierr /= 0) then
      ierr = 405
    end if
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
    integer(int32) :: ierr
    ierr = 0
    dims = shape(arr)
    call write_file_header(filename, unit, 2, 4, dims, ierr)
    write(unit, iostat=ierr) arr
    if (ierr /= 0) then
      ierr = 405
    end if
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
    integer(int32) :: ierr
    ierr = 0
    dims = shape(arr)
    call write_file_header(filename, unit, 2, 5, dims, ierr)
    write(unit, iostat=ierr) arr
    if (ierr /= 0) then
      ierr = 405
    end if
    close(unit)
  end subroutine

  !> @brief Writes serialized real array from R to file with metdata.
  subroutine serialize_real_nd(arr, dims, ndim, filename, ierr)
    real(real64), intent(in) :: arr(:)
      !! array to save
    integer(int32), intent(in) :: dims(:)
      !! Dimensions of the array
    integer(int32), intent(in) :: ndim
      !! Number of dimensions
    character(len=*), intent(in) :: filename
      !! filename
    integer(int32) :: unit
    integer(int32), INTENT(OUT) :: ierr
    ierr = 0

    if (size(dims) /= ndim) then
      error stop "Dimension mismatch in serialize_real_nd"
    end if

    call write_file_header(filename, unit, 2, ndim, dims, ierr)
    write(unit, iostat=ierr) arr
    if (ierr /= 0) then
      ierr = 405
    end if
    close(unit)
  end subroutine
 
end module serialize_real

!> Serialize a flat integer array with specified dimensions and number of dimensions to a binary file.
!! R can not pass a multidimensional array directly, so we use a flat array and dimensions. Therefore, exposing serialize_int_*d to R is not needed.
subroutine serialize_real_flat_r(arr, array_size, dims, ndim, filename_ascii, fn_len, ierr)
  use iso_fortran_env, only: int32, real64
  use array_utils, only: ascii_to_string
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
  integer(int32), intent(out) :: ierr
    !! Error code

  integer :: i, total_len

  call ascii_to_string(filename_ascii, fn_len, filename)

  ! calculate total size
  total_len = 1
  do i = 1, ndim
    total_len = total_len * dims(i)
  end do

  call serialize_real_nd(arr(1:total_len), dims(1:ndim), ndim, filename, ierr)
end subroutine

subroutine serialize_real_nd_C(arr, dims, ndim, filename_ascii, fn_len, ierr) bind(C, name="serialize_real_nd_C")
  use iso_c_binding, only: c_ptr, c_int, c_f_pointer, c_double
  use array_utils, only: ascii_to_string
  use serialize_real, only: serialize_real_nd
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
  integer(c_int), intent(out) :: ierr

  ! Local
  character(len=:), allocatable :: filename
  real(c_double), pointer :: arr_f(:)
  integer :: i

  call ascii_to_string(filename_ascii, fn_len, filename)

  ! 1D-Array from C pointer
  call c_f_pointer(arr, arr_f, [product(dims(1:ndim))])

  ! save
  call serialize_real_nd(arr_f, dims, ndim, filename, ierr)
end subroutine