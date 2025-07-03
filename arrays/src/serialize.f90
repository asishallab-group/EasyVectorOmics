!> Module providing serialization and deserialization routines for integer, real, and character arrays
!! of up to 5 dimensions, as well as utility routines for extracting rows, columns, and cells from 2D arrays.
!! Arrays are serialized to a custom binary format with a magic number and type/dimension metadata.
module serialize_mod
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding, only: c_loc
  implicit none

  public::serialize

  interface serialize
    module procedure serialize_int_1d, serialize_int_2d, serialize_int_3d, serialize_int_4d, serialize_int_5d
    module procedure serialize_real_1d, serialize_real_2d, serialize_real_3d, serialize_real_4d, serialize_real_5d
    module procedure serialize_char_1d, serialize_char_2d, serialize_char_3d, serialize_char_4d, serialize_char_5d
  end interface serialize

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

end module serialize_mod