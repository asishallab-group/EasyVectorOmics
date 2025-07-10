module serialize_real
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding, only: c_loc
  implicit none

  public:: serialize_real_1d, serialize_real_2d, serialize_real_3d, &
           serialize_real_4d, serialize_real_5d

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

 
end module serialize_real

  !> R-Interface: 1D Real-Array serialisieren
  subroutine serialize_real_1d_r(arr, n1, filename)
    use iso_c_binding
    use serialize_real
    implicit none
    real(real64), intent(in), target :: arr(*)
    integer, intent(in) :: n1
    character(len=*), intent(in) :: filename
    real(real64), pointer :: arr_f(:)
    call c_f_pointer(c_loc(arr(1)), arr_f, [n1])
    call serialize_real_1d(arr_f, filename)
  end subroutine

  !> R-Interface: 2D Real-Array serialisieren
  subroutine serialize_real_2d_r(arr, n1, n2, filename)
    use iso_c_binding
    use serialize_real
    implicit none
    real(real64), intent(in), target :: arr(*)
    integer, intent(in) :: n1, n2
    character(len=*), intent(in) :: filename
    real(real64), pointer :: arr_f(:,:)
    call c_f_pointer(c_loc(arr(1)), arr_f, [n1, n2])
    call serialize_real_2d(arr_f, filename)
  end subroutine

  !> R-Interface: 3D Real-Array serialisieren
  subroutine serialize_real_3d_r(arr, n1, n2, n3, filename)
    use iso_c_binding
    use serialize_real
    implicit none
    real(real64), intent(in), target :: arr(*)
    integer, intent(in) :: n1, n2, n3
    character(len=*), intent(in) :: filename
    real(real64), pointer :: arr_f(:,:,:)
    call c_f_pointer(c_loc(arr(1)), arr_f, [n1, n2, n3])
    call serialize_real_3d(arr_f, filename)
  end subroutine

  !> R-Interface: 4D Real-Array serialisieren
  subroutine serialize_real_4d_r(arr, n1, n2, n3, n4, filename)
    use iso_c_binding
    use serialize_real
    implicit none
    real(real64), intent(in), target :: arr(*)
    integer, intent(in) :: n1, n2, n3, n4
    character(len=*), intent(in) :: filename
    real(real64), pointer :: arr_f(:,:,:,:)
    call c_f_pointer(c_loc(arr(1)), arr_f, [n1, n2, n3, n4])
    call serialize_real_4d(arr_f, filename)
  end subroutine

  !> R-Interface: 5D Real-Array serialisieren
  subroutine serialize_real_5d_r(arr, n1, n2, n3, n4, n5, filename)
    use iso_c_binding
    use serialize_real
    implicit none
    real(real64), intent(in), target :: arr(*)
    integer, intent(in) :: n1, n2, n3, n4, n5
    character(len=*), intent(in) :: filename
    real(real64), pointer :: arr_f(:,:,:,:,:)
    call c_f_pointer(c_loc(arr(1)), arr_f, [n1, n2, n3, n4, n5])
    call serialize_real_5d(arr_f, filename)
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