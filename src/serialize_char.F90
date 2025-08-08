!> Module providing serialization and deserialization routines for character arrays
!! of up to 5 dimensions, arrays are serialized to a custom binary format with a magic number and type/dimension metadata.

module serialize_char
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding, only: c_loc
  implicit none

  public:: serialize_char_1d, serialize_char_2d, serialize_char_3d, &
           serialize_char_4d, serialize_char_5d, serialize_char_nd

  integer(int32), parameter :: ARRAY_FILE_MAGIC = int(z'46413230', int32) ! 'FA20' in hex

contains

  subroutine write_char_array_header(unit, type_code, ndim, dims, clen)
    use iso_fortran_env, only: int32
    implicit none
    integer, intent(in) :: unit
    integer, intent(in) :: type_code, ndim, clen
    integer(int32), intent(in) :: dims(ndim)

    write(unit) ARRAY_FILE_MAGIC
    write(unit) type_code
    write(unit) ndim
    write(unit) dims
    write(unit) clen
  end subroutine write_char_array_header

  !> Serialize a 1D character array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, character length, and the array data.
  subroutine serialize_char_1d(arr, filename)
    character(len=*), intent(in) :: arr(:)
    !! array to save
    character(len=*), intent(in) :: filename
    !! output filename
    integer :: unit, clen, i, str_len
    integer(int32) :: dims(1)
    dims = shape(arr)
    clen = len(arr)
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    call write_char_array_header(unit, 3, 1, dims, clen)
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
    !! array to save
    character(len=*), intent(in) :: filename
    !! output filename
    integer :: unit, clen, i, j, str_len
    integer(int32) :: dims(2)
    dims = shape(arr)
    clen = len(arr)
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    call write_char_array_header(unit, 3, 2, dims, clen)
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
  subroutine serialize_char_3d(arr, filename)
    character(len=*), intent(in) :: arr(:,:,:)
    !! array to save
    character(len=*), intent(in) :: filename
    !! output filename
    integer :: unit, clen, i, j, k, str_len
    integer(int32) :: dims(3)
    dims = shape(arr)
    clen = len(arr)
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    call write_char_array_header(unit, 3, 3, dims, clen)
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
  subroutine serialize_char_4d(arr, filename)
    character(len=*), intent(in) :: arr(:,:,:,:)
    !! array to save
    character(len=*), intent(in) :: filename
    !! output filename
    integer :: unit, clen, i, j, k, l, str_len
    integer(int32) :: dims(4)
    dims = shape(arr)
    clen = len(arr)
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    call write_char_array_header(unit, 3, 4, dims, clen)
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
  subroutine serialize_char_5d(arr, filename)
    character(len=*), intent(in) :: arr(:,:,:,:,:)
    !! array to save
    character(len=*), intent(in) :: filename
    !! output filename
    integer(int32) :: unit, clen, i, j, k, l, m, str_len
    integer(int32) :: dims(5)
    dims = shape(arr)
    clen = len(arr)
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    call write_char_array_header(unit, 3, 5, dims, clen)
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

  !> Serialize a character array of arbitrary dimensions to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, character length, and the array data.
  !! @note This routine is onyl called by R and serializes only flat character arrays to the memory
  subroutine serialize_char_nd(flat, dims, ndim, clen, filename)
    use iso_c_binding
    implicit none
    character(len=*), intent(in) :: flat(:)
    !! flat array to save
    integer(int32), intent(in) :: dims(:)
    !! dimensions of the array
    integer(int32), intent(in) :: ndim
    integer(int32), intent(in) :: clen
    character(len=*), intent(in) :: filename
    !! output filename

    integer :: unit, i, str_len

    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    call write_char_array_header(unit, 3, ndim, dims, clen)

    do i = 1, size(flat)
      str_len = len_trim(flat(i))
      write(unit) str_len
      if (str_len > 0) then
        write(unit) flat(i)(1:str_len)
      end if
    end do

    close(unit)
  end subroutine serialize_char_nd

end module serialize_char

!> serializes a flat character array to a binary file.
subroutine serialize_char_flat_r(ascii_arr, array_size, dims, ndim, clen, filename_ascii, fn_len)
  use iso_fortran_env
  use serialize_char
  use array_utils
  implicit none

  ! change to fixed size
  integer(int32), intent(in) :: ascii_arr(clen, array_size)
  !! Flat character array in ASCII format
  integer(int32), intent(in) :: dims(ndim)
  !! Dimensions of the array
  integer(int32), intent(in) :: ndim, array_size
  !! Number of dimensions
  integer(int32), intent(in) :: clen
  integer(int32), intent(in) :: filename_ascii(fn_len)
  integer(int32), intent(in) :: fn_len

  character(len=:), allocatable :: filename
  character(len=clen), allocatable :: flat(:)
  integer :: i, j, total

  total = product(dims)
  allocate(flat(total))

  ! ASCII to character conversion
  do i = 1, total
    flat(i) = ""
    do j = 1, clen
      if (ascii_arr(j, i) > 0) then
        flat(i)(j:j) = char(ascii_arr(j, i))
      else
        exit
      end if
    end do
  end do

  call ascii_to_string(filename_ascii, fn_len, filename)

  call serialize_char_nd(flat, dims, ndim, clen, filename)

end subroutine serialize_char_flat_r

!> C binding for the subroutine to serialize a flat character array to a binary file.
subroutine serialize_char_flat_C(ascii_ptr, dims, ndim, clen, filename_ascii, fn_len) bind(C, name="serialize_char_flat_C")
  use iso_c_binding
  use serialize_char
  use array_utils
  implicit none

  type(c_ptr), value :: ascii_ptr
  integer(c_int), intent(in) :: dims(ndim)
  !! Dimensions of the array
  integer(c_int), value :: ndim
  !! Number of dimensions
  integer(c_int), value :: clen
  !! Character length
  integer(c_int), intent(in) :: filename_ascii(fn_len)
  !! Array of ASCII characters representing the filename
  integer(c_int), value :: fn_len
  !! Length of the filename array

  integer(c_int), pointer :: ascii_arr(:)
  character(len=:), allocatable :: filename
  character(len=clen), allocatable :: flat(:)
  integer :: i, j, total

  total = product(dims)
  call c_f_pointer(ascii_ptr, ascii_arr, [clen * total])
  allocate(flat(total))

  ! ASCII to Fortran character(len=clen)
  do i = 1, total
    flat(i) = ""
    do j = 1, clen
      if (ascii_arr((i - 1) * clen + j) > 0) then
        flat(i)(j:j) = char(ascii_arr((i - 1) * clen + j))
      else
        exit
      end if
    end do
  end do

  call ascii_to_string(filename_ascii, fn_len, filename)

  call serialize_char_nd(flat, dims, ndim, clen, filename)
end subroutine serialize_char_flat_C