!> Module providing serialization and deserialization routines for character arrays
!! of up to 5 dimensions, arrays are serialized to a custom binary format with a magic number and type/dimension metadata.

module serialize_char
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding, only: c_loc
  implicit none

  public:: serialize_char_1d, serialize_char_2d, serialize_char_3d, &
           serialize_char_4d, serialize_char_5d

  integer(int32), parameter :: ARRAY_FILE_MAGIC = int(z'46413230', int32) ! 'FA20' in hex

contains

  !> Serialize a 1D character array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, character length, and the array data.
  !! @param arr The input character array to serialize.
  !! @param filename The output filename.
  subroutine serialize_char_1d(arr, filename)
    character(len=*), intent(in) :: arr(:)
    character(len=*), intent(in) :: filename
    integer :: unit, clen, i, str_len
    integer(int32) :: dims(1)
    dims = shape(arr)
    clen = len(arr)
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    write(unit) ARRAY_FILE_MAGIC
    write(unit) 3
    write(unit) 1
    write(unit) dims
    write(unit) clen
    do i = 1, dims(1)
      str_len = len_trim(arr(i))
      write(unit) str_len
      write(unit) arr(i)(1:str_len)
    end do
    close(unit)
  end subroutine

  !> Serialize a 2D character array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, character length, and the array data.
  !! @param arr The input character array to serialize.
  !! @param filename The output filename.
  subroutine serialize_char_2d(arr, filename)
    character(len=*), intent(in) :: arr(:,:)
    character(len=*), intent(in) :: filename
    integer :: unit, clen, i, j, str_len
    integer(int32) :: dims(2)
    dims = shape(arr)
    clen = len(arr)
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    write(unit) ARRAY_FILE_MAGIC
    write(unit) 3  ! Typcode für CHAR
    write(unit) 2  ! Dimension
    write(unit) dims
    write(unit) clen
    do j = 1, dims(2)
      do i = 1, dims(1)
        str_len = len_trim(arr(i,j))
        write(unit) str_len
        write(unit) arr(i,j)(1:str_len)
      end do
    end do
    close(unit)
  end subroutine

  !> Serialize a 3D character array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, character length, and the array data.
  !! @param arr The input character array to serialize.
  !! @param filename The output filename.
  subroutine serialize_char_3d(arr, filename)
    character(len=*), intent(in) :: arr(:,:,:)
    character(len=*), intent(in) :: filename
    integer :: unit, clen, i, j, k, str_len
    integer(int32) :: dims(3)
    dims = shape(arr)
    clen = len(arr)
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    write(unit) ARRAY_FILE_MAGIC
    write(unit) 3
    write(unit) 3
    write(unit) dims
    write(unit) clen
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          str_len = len_trim(arr(i,j,k))
          write(unit) str_len
          write(unit) arr(i,j,k)(1:str_len)
        end do
      end do
    end do
    close(unit)
  end subroutine

  !> Serialize a 4D character array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, character length, and the array data.
  !! @param arr The input character array to serialize.
  !! @param filename The output filename.
  subroutine serialize_char_4d(arr, filename)
    character(len=*), intent(in) :: arr(:,:,:,:)
    character(len=*), intent(in) :: filename
    integer :: unit, clen, i, j, k, l, str_len
    integer(int32) :: dims(4)
    dims = shape(arr)
    clen = len(arr)
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    write(unit) ARRAY_FILE_MAGIC
    write(unit) 3
    write(unit) 4
    write(unit) dims
    write(unit) clen
    do l = 1, dims(4)
      do k = 1, dims(3)
        do j = 1, dims(2)
          do i = 1, dims(1)
            str_len = len_trim(arr(i,j,k,l))
            write(unit) str_len
            write(unit) arr(i,j,k,l)(1:str_len)
          end do
        end do
      end do
    end do
    close(unit)
  end subroutine

  !> Serialize a 5D character array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, character length, and the array data.
  !! @param arr The input character array to serialize.
  !! @param filename The output filename.
  subroutine serialize_char_5d(arr, filename)
    character(len=*), intent(in) :: arr(:,:,:,:,:)
    character(len=*), intent(in) :: filename
    integer :: unit, clen, i, j, k, l, m, str_len
    integer(int32) :: dims(5)
    dims = shape(arr)
    clen = len(arr)
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    write(unit) ARRAY_FILE_MAGIC
    write(unit) 3                    ! Typecode für CHARACTER
    write(unit) 5                    ! Dimension
    write(unit) dims                 ! Shape
    write(unit) clen                 ! Maximale Zeichenlänge
    do m = 1, dims(5)
      do l = 1, dims(4)
        do k = 1, dims(3)
          do j = 1, dims(2)
            do i = 1, dims(1)
              str_len = len_trim(arr(i,j,k,l,m))
              write(unit) str_len
              write(unit) arr(i,j,k,l,m)(1:str_len)
            end do
          end do
        end do
      end do
    end do
    close(unit)
  end subroutine


end module serialize_char

!> R-Interface: 1D Character-Array serialisieren
subroutine serialize_char_1d_r(arr, n1, filename)
  use iso_c_binding
  use serialize_char
  implicit none
  character(len=*), intent(in), target :: arr(*)
  integer, intent(in) :: n1
  character(len=*), intent(in) :: filename
  character(len=:), pointer :: arr_f(:)
  integer :: clen
  clen = len(arr)
  call c_f_pointer(c_loc(arr(1)), arr_f, [n1])
  call serialize_char_1d(arr_f, filename)
end subroutine

!> R-Interface: 2D Character-Array serialisieren
subroutine serialize_char_2d_r(arr, n1, n2, filename)
  use iso_c_binding
  use serialize_char
  implicit none
  character(len=*), intent(in), target :: arr(*)
  integer, intent(in) :: n1, n2
  character(len=*), intent(in) :: filename
  character(len=:), pointer :: arr_f(:,:)
  integer :: clen
  clen = len(arr)
  call c_f_pointer(c_loc(arr(1)), arr_f, [n1, n2])
  call serialize_char_2d(arr_f, filename)
end subroutine

!> R-Interface: 3D Character-Array serialisieren
subroutine serialize_char_3d_r(arr, n1, n2, n3, filename)
  use iso_c_binding
  use serialize_char
  implicit none
  character(len=*), intent(in), target :: arr(*)
  integer, intent(in) :: n1, n2, n3
  character(len=*), intent(in) :: filename
  character(len=:), pointer :: arr_f(:,:,:)
  integer :: clen
  clen = len(arr)
  call c_f_pointer(c_loc(arr(1)), arr_f, [n1, n2, n3])
  call serialize_char_3d(arr_f, filename)
end subroutine

!> R-Interface: 4D Character-Array serialisieren
subroutine serialize_char_4d_r(arr, n1, n2, n3, n4, filename)
  use iso_c_binding
  use serialize_char
  implicit none
  character(len=*), intent(in), target :: arr(*)
  integer, intent(in) :: n1, n2, n3, n4
  character(len=*), intent(in) :: filename
  character(len=:), pointer :: arr_f(:,:,:,:)
  integer :: clen
  clen = len(arr)
  call c_f_pointer(c_loc(arr(1)), arr_f, [n1, n2, n3, n4])
  call serialize_char_4d(arr_f, filename)
end subroutine

!> R-Interface: 5D Character-Array serialisieren
subroutine serialize_char_5d_r(arr, n1, n2, n3, n4, n5, filename)
  use iso_c_binding
  use serialize_char
  implicit none
  character(len=*), intent(in), target :: arr(*)
  integer, intent(in) :: n1, n2, n3, n4, n5
  character(len=*), intent(in) :: filename
  character(len=:), pointer :: arr_f(:,:,:,:,:)
  integer :: clen
  clen = len(arr)
  call c_f_pointer(c_loc(arr(1)), arr_f, [n1, n2, n3, n4, n5])
  call serialize_char_5d(arr_f, filename)
end subroutine

! --- C-Bindings für serialize_char_* ---

subroutine serialize_char_1d_C(arr, n1, filename) bind(C, name="serialize_char_1d_C")
  use iso_c_binding
  use serialize_char
  implicit none
  type(c_ptr), value :: arr
  integer, value :: n1
  character(kind=c_char), intent(in) :: filename(*)
  character(len=:), pointer :: arr_f(:)
  character(len=:), allocatable :: fname
  integer :: i

  call c_f_pointer(arr, arr_f, [n1])

  i = 1
  do while (filename(i) /= c_null_char)
    i = i + 1
  end do
  fname = transfer(filename(1:i-1), fname)
  call serialize_char_1d(arr_f, fname)
end subroutine

subroutine serialize_char_2d_C(arr, n1, n2, filename) bind(C, name="serialize_char_2d_C")
  use iso_c_binding
  use serialize_char
  implicit none
  type(c_ptr), value :: arr
  integer, value :: n1, n2
  character(kind=c_char), intent(in) :: filename(*)
  character(len=:), pointer :: arr_f(:,:)
  character(len=:), allocatable :: fname
  integer :: i

  call c_f_pointer(arr, arr_f, [n1, n2])

  i = 1
  do while (filename(i) /= c_null_char)
    i = i + 1
  end do
  fname = transfer(filename(1:i-1), fname)
  call serialize_char_2d(arr_f, fname)
end subroutine

subroutine serialize_char_3d_C(arr, n1, n2, n3, filename) bind(C, name="serialize_char_3d_C")
  use iso_c_binding
  use serialize_char
  implicit none
  type(c_ptr), value :: arr
  integer, value :: n1, n2, n3
  character(kind=c_char), intent(in) :: filename(*)
  character(len=:), pointer :: arr_f(:,:,:)
  character(len=:), allocatable :: fname
  integer :: i

  call c_f_pointer(arr, arr_f, [n1, n2, n3])

  i = 1
  do while (filename(i) /= c_null_char)
    i = i + 1
  end do
  fname = transfer(filename(1:i-1), fname)
  call serialize_char_3d(arr_f, fname)
end subroutine

subroutine serialize_char_4d_C(arr, n1, n2, n3, n4, filename) bind(C, name="serialize_char_4d_C")
  use iso_c_binding
  use serialize_char
  implicit none
  type(c_ptr), value :: arr
  integer, value :: n1, n2, n3, n4
  character(kind=c_char), intent(in) :: filename(*)
  character(len=:), pointer :: arr_f(:,:,:,:)
  character(len=:), allocatable :: fname
  integer :: i

  call c_f_pointer(arr, arr_f, [n1, n2, n3, n4])

  i = 1
  do while (filename(i) /= c_null_char)
    i = i + 1
  end do
  fname = transfer(filename(1:i-1), fname)
  call serialize_char_4d(arr_f, fname)
end subroutine

subroutine serialize_char_5d_C(arr, n1, n2, n3, n4, n5, filename) bind(C, name="serialize_char_5d_C")
  use iso_c_binding
  use serialize_char
  implicit none
  type(c_ptr), value :: arr
  integer, value :: n1, n2, n3, n4, n5
  character(kind=c_char), intent(in) :: filename(*)
  character(len=:), pointer :: arr_f(:,:,:,:,:)
  character(len=:), allocatable :: fname
  integer :: i

  call c_f_pointer(arr, arr_f, [n1, n2, n3, n4, n5])

  i = 1
  do while (filename(i) /= c_null_char)
    i = i + 1
  end do
  fname = transfer(filename(1:i-1), fname)
  call serialize_char_5d(arr_f, fname)
end subroutine
