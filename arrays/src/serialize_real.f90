module serialize_real
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding, only: c_loc
  implicit none

  public:: serialize_real_1d, serialize_real_2d, serialize_real_3d, &
           serialize_real_4d, serialize_real_5d, serialize_real_nd

  integer(int32), parameter :: ARRAY_FILE_MAGIC = int(z'46413230', int32) ! 'FA20' in hex

contains

  !> Serialize a 1D real(real64) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  !! @param arr The input real array to serialize.
  !! @param filename The output filename.
  subroutine serialize_real_1d(arr, filename)
    real(real64), intent(in) :: arr(:)
    character(len=*), intent(in) :: filename
    integer :: unit
    integer(int32) :: dims(1)
    dims = shape(arr)
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    write(unit) ARRAY_FILE_MAGIC
    write(unit) 2
    write(unit) 1
    write(unit) dims
    write(unit) arr
    close(unit)
  end subroutine

  !> Serialize a 2D real(real64) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  !! @param arr The input real array to serialize.
  !! @param filename The output filename.
  subroutine serialize_real_2d(arr, filename)
    real(real64), intent(in) :: arr(:,:)
    character(len=*), intent(in) :: filename
    integer :: unit
    integer(int32) :: dims(2)
    dims = shape(arr)
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    write(unit) ARRAY_FILE_MAGIC
    write(unit) 2
    write(unit) 2
    write(unit) dims
    write(unit) arr
    close(unit)
  end subroutine

  !> Serialize a 3D real(real64) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  !! @param arr The input real array to serialize.
  !! @param filename The output filename.
  subroutine serialize_real_3d(arr, filename)
    real(real64), intent(in) :: arr(:,:,:)
    character(len=*), intent(in) :: filename
    integer :: unit
    integer(int32) :: dims(3)
    dims = shape(arr)
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    write(unit) ARRAY_FILE_MAGIC
    write(unit) 2
    write(unit) 3
    write(unit) dims
    write(unit) arr
    close(unit)
  end subroutine

  !> Serialize a 4D real(real64) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  !! @param arr The input real array to serialize.
  !! @param filename The output filename.
  subroutine serialize_real_4d(arr, filename)
    real(real64), intent(in) :: arr(:,:,:,:)
    character(len=*), intent(in) :: filename
    integer :: unit
    integer(int32) :: dims(4)
    dims = shape(arr)
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    write(unit) ARRAY_FILE_MAGIC
    write(unit) 2
    write(unit) 4
    write(unit) dims
    write(unit) arr
    close(unit)
  end subroutine

  !> Serialize a 5D real(real64) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  !! @param arr The input real array to serialize.
  !! @param filename The output filename.
  subroutine serialize_real_5d(arr, filename)
    real(real64), intent(in) :: arr(:,:,:,:,:)
    character(len=*), intent(in) :: filename
    integer :: unit
    integer(int32) :: dims(5)
    dims = shape(arr)
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    write(unit) ARRAY_FILE_MAGIC
    write(unit) 2
    write(unit) 5
    write(unit) dims
    write(unit) arr
    close(unit)
  end subroutine

    !> Writes serialized real array from R to file with metdata.
  subroutine serialize_real_nd(arr, dims, ndim, filename)
    real(real64), intent(in) :: arr(:)
    integer(int32), intent(in) :: dims(:)
    integer(int32), intent(in) :: ndim
    character(len=*), intent(in) :: filename
    integer :: unit

    if (size(dims) /= ndim) then
      error stop "Dimension mismatch in serialize_real_nd"
    end if

    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    write(unit) ARRAY_FILE_MAGIC
    write(unit) 2
    write(unit) ndim
    write(unit) dims
    write(unit) arr
    close(unit)
  end subroutine
 
end module serialize_real

!> Serialize a flat integer array with specified dimensions and number of dimensions to a binary file.
!! R can not pass a multidimensional array directly, so we use a flat array and dimensions. Therefore, exposing serialize_int_*d to R is not needed.
!! @param arr The input integer array to serialize.
!! @param dims The dimensions of the array.
!! @param ndim The number of dimensions.
!! @param filename_ascii The output filename as an ASCII character array.
!! @param fn_len The length of the filename ASCII character array.
subroutine serialize_real_flat_r(arr, dims, ndim, filename_ascii, fn_len)
  use iso_fortran_env
  use serialize_real, only: serialize_real_nd
  implicit none
  real(real64), intent(in) :: arr(*)         ! assumed-size array
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

  call serialize_real_nd(arr(1:total_len), dims(1:ndim), ndim, filename)
end subroutine

  ! --- C-Bindings für serialize_real_* ---

  subroutine serialize_real_1d_C(arr, n1, filename) bind(C, name="serialize_real_1d_C")
    use iso_c_binding
    use serialize_real
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    type(c_ptr), value :: arr
    integer, value :: n1
    character(kind=c_char), intent(in) :: filename(*)
    real(real64), pointer :: arr_f(:)
    character(len=:), allocatable :: fname
    integer :: i

    call c_f_pointer(arr, arr_f, [n1])

    i = 1
    do while (filename(i) /= c_null_char)
      i = i + 1
    end do
    fname = transfer(filename(1:i-1), fname)
    call serialize_real_1d(arr_f, fname)
  end subroutine

  subroutine serialize_real_2d_C(arr, n1, n2, filename) bind(C, name="serialize_real_2d_C")
    use iso_c_binding
    use serialize_real
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    type(c_ptr), value :: arr
    integer, value :: n1, n2
    character(kind=c_char), intent(in) :: filename(*)
    real(real64), pointer :: arr_f(:,:)
    character(len=:), allocatable :: fname
    integer :: i

    call c_f_pointer(arr, arr_f, [n1, n2])

    i = 1
    do while (filename(i) /= c_null_char)
      i = i + 1
    end do
    fname = transfer(filename(1:i-1), fname)
    call serialize_real_2d(arr_f, fname)
  end subroutine

  subroutine serialize_real_3d_C(arr, n1, n2, n3, filename) bind(C, name="serialize_real_3d_C")
    use iso_c_binding
    use serialize_real
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    type(c_ptr), value :: arr
    integer, value :: n1, n2, n3
    character(kind=c_char), intent(in) :: filename(*)
    real(real64), pointer :: arr_f(:,:,:)
    character(len=:), allocatable :: fname
    integer :: i

    call c_f_pointer(arr, arr_f, [n1, n2, n3])

    i = 1
    do while (filename(i) /= c_null_char)
      i = i + 1
    end do
    fname = transfer(filename(1:i-1), fname)
    call serialize_real_3d(arr_f, fname)
  end subroutine

  subroutine serialize_real_4d_C(arr, n1, n2, n3, n4, filename) bind(C, name="serialize_real_4d_C")
    use iso_c_binding
    use serialize_real
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    type(c_ptr), value :: arr
    integer, value :: n1, n2, n3, n4
    character(kind=c_char), intent(in) :: filename(*)
    real(real64), pointer :: arr_f(:,:,:,:)
    character(len=:), allocatable :: fname
    integer :: i

    call c_f_pointer(arr, arr_f, [n1, n2, n3, n4])

    i = 1
    do while (filename(i) /= c_null_char)
      i = i + 1
    end do
    fname = transfer(filename(1:i-1), fname)
    call serialize_real_4d(arr_f, fname)
  end subroutine

  subroutine serialize_real_5d_C(arr, n1, n2, n3, n4, n5, filename) bind(C, name="serialize_real_5d_C")
    use iso_c_binding
    use serialize_real
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    type(c_ptr), value :: arr
    integer, value :: n1, n2, n3, n4, n5
    character(kind=c_char), intent(in) :: filename(*)
    real(real64), pointer :: arr_f(:,:,:,:,:)
    character(len=:), allocatable :: fname
    integer :: i

    call c_f_pointer(arr, arr_f, [n1, n2, n3, n4, n5])

    i = 1
    do while (filename(i) /= c_null_char)
        i = i + 1
    end do
    fname = transfer(filename(1:i-1), fname)
    call serialize_real_5d(arr_f, fname)
  end subroutine