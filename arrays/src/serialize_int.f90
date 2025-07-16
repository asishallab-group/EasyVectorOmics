module serialize_int
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding, only: c_loc
  implicit none

  public:: serialize_int_1d, serialize_int_2d, serialize_int_3d, &
           serialize_int_4d, serialize_int_5d, serialize_int_nd

  integer(int32), parameter :: ARRAY_FILE_MAGIC = int(z'46413230', int32) ! 'FA20' in hex

contains

  !> Serialize a 1D integer(int32) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  !! @param arr The input integer array to serialize.
  !! @param filename The output filename.
  subroutine serialize_int_1d(arr, filename)
    integer(int32), intent(in) :: arr(:)
    character(len=*), intent(in) :: filename
    integer :: unit
    integer(int32) :: dims(1)
    dims = shape(arr)
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    write(unit) ARRAY_FILE_MAGIC
    write(unit) 1
    write(unit) 1
    write(unit) dims
    write(unit) arr
    close(unit)
  end subroutine

  !> Serialize a 2D integer(int32) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  !! @param arr The input integer array to serialize.
  !! @param filename The output filename.
  subroutine serialize_int_2d(arr, filename)
    integer(int32), intent(in) :: arr(:,:)
    character(len=*), intent(in) :: filename
    integer :: unit
    integer(int32) :: dims(2)
    dims = shape(arr)
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    write(unit) ARRAY_FILE_MAGIC
    write(unit) 1
    write(unit) 2
    write(unit) dims
    write(unit) arr
    close(unit)
  end subroutine

  !> Serialize a 3D integer(int32) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  !! @param arr The input integer array to serialize.
  !! @param filename The output filename.
  subroutine serialize_int_3d(arr, filename)
    integer(int32), intent(in) :: arr(:,:,:)
    character(len=*), intent(in) :: filename
    integer :: unit
    integer(int32) :: dims(3)
    dims = shape(arr)
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    write(unit) ARRAY_FILE_MAGIC
    write(unit) 1
    write(unit) 3
    write(unit) dims
    write(unit) arr
    close(unit)
  end subroutine

  !> Serialize a 4D integer(int32) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  !! @param arr The input integer array to serialize.
  !! @param filename The output filename.
  subroutine serialize_int_4d(arr, filename)
    integer(int32), intent(in) :: arr(:,:,:,:)
    character(len=*), intent(in) :: filename
    integer :: unit
    integer(int32) :: dims(4)
    dims = shape(arr)
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    write(unit) ARRAY_FILE_MAGIC
    write(unit) 1
    write(unit) 4
    write(unit) dims
    write(unit) arr
    close(unit)
  end subroutine

  !> Serialize a 5D integer(int32) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  !! @param arr The input integer array to serialize.
  !! @param filename The output filename.
  subroutine serialize_int_5d(arr, filename)
    integer(int32), intent(in) :: arr(:,:,:,:,:)
    character(len=*), intent(in) :: filename
    integer :: unit
    integer(int32) :: dims(5)
    dims = shape(arr)
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    write(unit) ARRAY_FILE_MAGIC
    write(unit) 1
    write(unit) 5
    write(unit) dims
    write(unit) arr
    close(unit)
  end subroutine

  subroutine serialize_int_nd(arr, dims, ndim, filename)
    integer(int32), intent(in) :: arr(:)
    integer(int32), intent(in) :: dims(:)
    integer(int32), intent(in) :: ndim
    character(len=*), intent(in) :: filename
    integer :: unit

    if (size(dims) /= ndim) then
      error stop "Dimension mismatch in serialize_int_nd"
    end if

    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    write(unit) ARRAY_FILE_MAGIC
    write(unit) 1
    write(unit) ndim
    write(unit) dims
    write(unit) arr
    close(unit)
  end subroutine
end module serialize_int


!> Serialize a flat integer array with specified dimensions and number of dimensions to a binary file.
!! R can not pass a multidimensional array directly, so we use a flat array and dimensions. Therefore, exposing serialize_int_*d to R is not needed.
!! @param arr The input integer array to serialize.
!! @param dims The dimensions of the array.
!! @param ndim The number of dimensions.
!! @param filename_ascii The output filename as an ASCII character array.
!! @param fn_len The length of the filename ASCII character array.
subroutine serialize_int_flat_r(arr, dims, ndim, filename_ascii, fn_len)
  use iso_fortran_env, only: int32
  use serialize_int, only: serialize_int_nd
  implicit none
  integer(int32), intent(in) :: arr(*)         ! assumed-size array
  integer(int32), intent(in) :: dims(*)
  integer(int32), intent(in) :: ndim
  integer(int32), intent(in) :: filename_ascii(fn_len)
  integer(int32), intent(in) :: fn_len

  character(len=:), allocatable :: filename
  integer :: i, total_len

  allocate(character(len=fn_len) :: filename)
  do i = 1, fn_len
    filename(i:i) = char(filename_ascii(i))
  end do

  ! Gesamtgröße berechnen (z. B. für Sicherheit oder Logging)
  total_len = 1
  do i = 1, ndim
    total_len = total_len * dims(i)
  end do

  ! optional check
  ! print *, "Serializing array with ", total_len, " elements and ", ndim, " dimensions."

  call serialize_int_nd(arr(1:total_len), dims(1:ndim), ndim, filename)
end subroutine


! --- C-Bindings für serialize_int_* ---

subroutine serialize_int_nd_C(arr, dims, ndim, filename_ascii, fn_len) bind(C, name="serialize_int_nd_C")
  use iso_c_binding
  use serialize_int
  implicit none

  ! Eingabeparameter
  type(c_ptr), value :: arr
  integer(c_int), intent(in) :: dims(ndim)
  integer(c_int), value :: ndim
  integer(c_int), intent(in) :: filename_ascii(fn_len)
  integer(c_int), value :: fn_len

  ! Lokale Variablen
  character(len=:), allocatable :: filename
  integer(c_int), pointer :: arr_f(:)
  integer :: i

  ! Filename umwandeln
  allocate(character(len=fn_len) :: filename)
  do i = 1, fn_len
    filename(i:i) = char(filename_ascii(i))
  end do

  ! 1D-Array aus dem C-Pointer
  call c_f_pointer(arr, arr_f, [product(dims(1:ndim))])

  ! Speichern
  call serialize_int_nd(arr_f, dims, ndim, filename)
end subroutine
