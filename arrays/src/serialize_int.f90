module serialize_int
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding, only: c_loc
  implicit none

  public:: serialize_int_1d, serialize_int_2d, serialize_int_3d, &
           serialize_int_4d, serialize_int_5d

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
end module serialize_int

!> R-Interface: 1D Integer-Array serialisieren
subroutine serialize_int_1d_r(arr, n1, filename_ascii, fn_len)
  use serialize_int
  implicit none
  integer(int32), intent(in) :: arr(n1)
  integer(int32), intent(in) :: n1
  integer(int32), intent(in) :: filename_ascii(fn_len)
  integer(int32), intent(in) :: fn_len

  character(len=:), allocatable :: filename
  integer :: i

  allocate(character(len=fn_len) :: filename)

  do i = 1, fn_len
    filename(i:i) = char(filename_ascii(i))
  end do

  call serialize_int_1d(arr, filename)
end subroutine


!> R-Interface: 2D Integer-Array serialisieren
subroutine serialize_int_2d_r(arr, n1, n2, filename)
  use iso_c_binding
  use serialize_int
  implicit none
  integer(int32), intent(in) :: arr(*)
  integer, intent(in) :: n1, n2
  character(len=*), intent(in) :: filename
  ! Siehe oben: kein c_f_pointer!
  call serialize_int_2d(reshape(arr(1:n1*n2), [n1, n2]), filename)
end subroutine

!> R-Interface: 3D Integer-Array serialisieren
subroutine serialize_int_3d_r(arr, n1, n2, n3, filename)
  use iso_c_binding
  use serialize_int
  implicit none
  integer(int32), intent(in) :: arr(*)
  integer, intent(in) :: n1, n2, n3
  character(len=*), intent(in) :: filename
  call serialize_int_3d(reshape(arr(1:n1*n2*n3), [n1, n2, n3]), filename)
end subroutine

!> R-Interface: 4D Integer-Array serialisieren
subroutine serialize_int_4d_r(arr, n1, n2, n3, n4, filename)
  use iso_c_binding
  use serialize_int
  implicit none
  integer(int32), intent(in) :: arr(*)
  integer, intent(in) :: n1, n2, n3, n4
  character(len=*), intent(in) :: filename
  call serialize_int_4d(reshape(arr(1:n1*n2*n3*n4), [n1, n2, n3, n4]), filename)
end subroutine

!> R-Interface: 5D Integer-Array serialisieren
subroutine serialize_int_5d_r(arr, n1, n2, n3, n4, n5, filename)
  use iso_c_binding
  use serialize_int
  implicit none
  integer(int32), intent(in) :: arr(*)
  integer, intent(in) :: n1, n2, n3, n4, n5
  character(len=*), intent(in) :: filename
  call serialize_int_5d(reshape(arr(1:n1*n2*n3*n4*n5), [n1, n2, n3, n4, n5]), filename)
end subroutine

! --- C-Bindings für serialize_int_* ---

subroutine serialize_int_1d_C(arr, n1, filename) bind(C, name="serialize_int_1d_C")
  use iso_c_binding
  use serialize_int
  use, intrinsic :: iso_fortran_env, only: int32
  implicit none
  type(c_ptr), value :: arr
  integer, value :: n1
  character(kind=c_char), intent(in) :: filename(*)
  integer(int32), pointer :: arr_f(:)
  character(len=:), allocatable :: fname
  integer :: i

  call c_f_pointer(arr, arr_f, [n1])

  i = 1
  do while (filename(i) /= c_null_char)
    i = i + 1
  end do
  fname = transfer(filename(1:i-1), fname)
  call serialize_int_1d(arr_f, fname)
end subroutine

subroutine serialize_int_2d_C(arr, n1, n2, filename) bind(C, name="serialize_int_2d_C")
  use iso_c_binding
  use serialize_int
  use, intrinsic :: iso_fortran_env, only: int32
  implicit none
  type(c_ptr), value :: arr
  integer, value :: n1, n2
  character(kind=c_char), intent(in) :: filename(*)
  integer(int32), pointer :: arr_f(:,:)
  character(len=:), allocatable :: fname
  integer :: i

  call c_f_pointer(arr, arr_f, [n1, n2])

  i = 1
  do while (filename(i) /= c_null_char)
    i = i + 1
  end do
  fname = transfer(filename(1:i-1), fname)
  call serialize_int_2d(arr_f, fname)
end subroutine

subroutine serialize_int_3d_C(arr, n1, n2, n3, filename) bind(C, name="serialize_int_3d_C")
  use iso_c_binding
  use serialize_int
  use, intrinsic :: iso_fortran_env, only: int32
  implicit none
  type(c_ptr), value :: arr
  integer, value :: n1, n2, n3
  character(kind=c_char), intent(in) :: filename(*)
  integer(int32), pointer :: arr_f(:,:,:)
  character(len=:), allocatable :: fname
  integer :: i

  call c_f_pointer(arr, arr_f, [n1, n2, n3])

  i = 1
  do while (filename(i) /= c_null_char)
    i = i + 1
  end do
  fname = transfer(filename(1:i-1), fname)
  call serialize_int_3d(arr_f, fname)
end subroutine

subroutine serialize_int_4d_C(arr, n1, n2, n3, n4, filename) bind(C, name="serialize_int_4d_C")
  use iso_c_binding
  use serialize_int
  use, intrinsic :: iso_fortran_env, only: int32
  implicit none
  type(c_ptr), value :: arr
  integer, value :: n1, n2, n3, n4
  character(kind=c_char), intent(in) :: filename(*)
  integer(int32), pointer :: arr_f(:,:,:,:)
  character(len=:), allocatable :: fname
  integer :: i

  call c_f_pointer(arr, arr_f, [n1, n2, n3, n4])

  i = 1
  do while (filename(i) /= c_null_char)
    i = i + 1
  end do
  fname = transfer(filename(1:i-1), fname)
  call serialize_int_4d(arr_f, fname)
end subroutine

subroutine serialize_int_5d_C(arr, n1, n2, n3, n4, n5, filename) bind(C, name="serialize_int_5d_C")
  use iso_c_binding
  use serialize_int
  use, intrinsic :: iso_fortran_env, only: int32
  implicit none
  type(c_ptr), value :: arr
  integer, value :: n1, n2, n3, n4, n5
  character(kind=c_char), intent(in) :: filename(*)
  integer(int32), pointer :: arr_f(:,:,:,:,:)
  character(len=:), allocatable :: fname
  integer :: i

  call c_f_pointer(arr, arr_f, [n1, n2, n3, n4, n5])

  i = 1
  do while (filename(i) /= c_null_char)
    i = i + 1
  end do
  fname = transfer(filename(1:i-1), fname)
  call serialize_int_5d(arr_f, fname)
end subroutine