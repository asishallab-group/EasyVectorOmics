module array_ops
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

  ! ============================================
  ! Serialization Functions (INT)
  ! ============================================
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

  ! ============================================
  ! Serialization Functions (REAL)
  ! ============================================
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

  ! ============================================
  ! Serialization Functions (CHARACTER)
  ! ============================================
  subroutine serialize_char_1d(arr, filename)
    character(len=*), intent(in) :: arr(:)
    character(len=*), intent(in) :: filename
    integer :: unit, clen, i
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
      write(unit) arr(i)(1:clen)
    end do
    close(unit)
  end subroutine

  subroutine serialize_char_2d(arr, filename)
    character(len=*), intent(in) :: arr(:,:)
    character(len=*), intent(in) :: filename
    integer :: unit, clen, i, j
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
        write(unit) arr(i,j)(1:clen)
      end do
    end do

    close(unit)
  end subroutine

  subroutine serialize_char_3d(arr, filename)
    character(len=*), intent(in) :: arr(:,:,:)
    character(len=*), intent(in) :: filename
    integer :: unit, clen, i, j, k
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
          write(unit) arr(i,j,k)(1:clen)
        end do
      end do
    end do
    close(unit)
  end subroutine

  subroutine serialize_char_4d(arr, filename)
    character(len=*), intent(in) :: arr(:,:,:,:)
    character(len=*), intent(in) :: filename
    integer :: unit, clen, i, j, k, l
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
            write(unit) arr(i,j,k,l)(1:clen)
          end do
        end do
      end do
    end do
    close(unit)
  end subroutine

  subroutine serialize_char_5d(arr, filename)
    character(len=*), intent(in) :: arr(:,:,:,:,:)
    character(len=*), intent(in) :: filename
    integer :: unit, clen, i, j, k, l, m
    integer(int32) :: dims(5)
    dims = shape(arr)
    clen = len(arr)
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    write(unit) ARRAY_FILE_MAGIC
    write(unit) 3
    write(unit) 5
    write(unit) dims
    write(unit) clen
    do m = 1, dims(5)
      do l = 1, dims(4)
        do k = 1, dims(3)
          do j = 1, dims(2)
            do i = 1, dims(1)
              write(unit) arr(i,j,k,l,m)(1:clen)
            end do
          end do
        end do
      end do
    end do
    close(unit)
  end subroutine

  ! ============================================
  ! Flat Deserialization Functions
  ! ============================================

  subroutine deserialize_real_flat(flat, dims, filename)
    use iso_c_binding
    real(real64), pointer, intent(out) :: flat(:)
    integer(int32), allocatable, intent(out) :: dims(:)
    character(len=*), intent(in) :: filename

    integer :: unit, magic, type_code, d

    open(newunit=unit, file=filename, form='unformatted', access='stream', status='old')
    read(unit) magic
    if (magic /= ARRAY_FILE_MAGIC) error stop "Invalid file format"
    read(unit) type_code
    if (type_code /= 2) error stop "Expected real64 data"
    read(unit) d
    allocate(dims(d))
    read(unit) dims
    allocate(flat(product(dims)))
    read(unit) flat
    close(unit)
  end subroutine deserialize_real_flat

  subroutine deserialize_int_flat(flat, dims, filename)
    use iso_c_binding
    integer(int32), pointer, intent(out) :: flat(:)
    integer(int32), allocatable, intent(out) :: dims(:)
    character(len=*), intent(in) :: filename

    integer :: unit, magic, type_code, d

    open(newunit=unit, file=filename, form='unformatted', access='stream', status='old')
    read(unit) magic
    if (magic /= ARRAY_FILE_MAGIC) error stop "Invalid file format"
    read(unit) type_code
    if (type_code /= 1) error stop "Expected int32 data"
    read(unit) d
    allocate(dims(d))
    read(unit) dims
    allocate(flat(product(dims)))
    read(unit) flat
    close(unit)
  end subroutine deserialize_int_flat

  subroutine deserialize_char_flat(flat, dims, clen, filename)
    use iso_c_binding
    character(len=:), pointer, intent(out) :: flat(:)
    integer(int32), allocatable, intent(out) :: dims(:)
    integer, intent(out) :: clen
    character(len=*), intent(in) :: filename

    integer :: unit, magic, type_code, d, i

    open(newunit=unit, file=filename, form='unformatted', access='stream', status='old')
    read(unit) magic
    if (magic /= ARRAY_FILE_MAGIC) error stop "Invalid file format"
    read(unit) type_code
    if (type_code /= 3) error stop "Expected character data"
    read(unit) d
    allocate(dims(d))
    read(unit) dims
    read(unit) clen
    allocate(character(len=clen) :: flat(product(dims)))
    do i = 1, product(dims)
      read(unit) flat(i)
    end do
    close(unit)
  end subroutine deserialize_char_flat
end module array_ops

module reshape_utils
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding, only: c_loc, c_f_pointer
  implicit none

  public :: reshape_real_to_1D, reshape_real_to_2D, reshape_real_to_3D, reshape_real_to_4D, reshape_real_to_5D
  public :: reshape_int_to_1D, reshape_int_to_2D, reshape_int_to_3D, reshape_int_to_4D, reshape_int_to_5D
  public :: reshape_char_to_1D, reshape_char_to_2D, reshape_char_to_3D, reshape_char_to_4D, reshape_char_to_5D
contains

  subroutine reshape_real_to_1D(flat, arr, dims)
    real(real64), pointer, intent(in) :: flat(:)
    real(real64), pointer, intent(out) :: arr(:)
    integer(int32), intent(in) :: dims(1)
    call c_f_pointer(c_loc(flat(1)), arr, dims)
  end subroutine reshape_real_to_1D

  subroutine reshape_real_to_2D(flat, arr, dims)
    real(real64), pointer, intent(in) :: flat(:)
    real(real64), pointer, intent(out) :: arr(:,:)
    integer(int32), intent(in) :: dims(2)
    call c_f_pointer(c_loc(flat(1)), arr, dims)
  end subroutine reshape_real_to_2D

  subroutine reshape_real_to_3D(flat, arr, dims)
    real(real64), pointer, intent(in) :: flat(:)
    real(real64), pointer, intent(out) :: arr(:,:,:)
    integer(int32), intent(in) :: dims(3)
    call c_f_pointer(c_loc(flat(1)), arr, dims)
  end subroutine reshape_real_to_3D

  subroutine reshape_real_to_4D(flat, arr, dims)
    real(real64), pointer, intent(in) :: flat(:)
    real(real64), pointer, intent(out) :: arr(:,:,:,:)
    integer(int32), intent(in) :: dims(4)
    call c_f_pointer(c_loc(flat(1)), arr, dims)
  end subroutine reshape_real_to_4D

  subroutine reshape_real_to_5D(flat, arr, dims)
    real(real64), pointer, intent(in) :: flat(:)
    real(real64), pointer, intent(out) :: arr(:,:,:,:,:)
    integer(int32), intent(in) :: dims(5)
    call c_f_pointer(c_loc(flat(1)), arr, dims)
  end subroutine reshape_real_to_5D

  subroutine reshape_int_to_1D(flat, arr, dims)
    integer(int32), pointer, intent(in) :: flat(:)
    integer(int32), pointer, intent(out) :: arr(:)
    integer(int32), intent(in) :: dims(1)
    call c_f_pointer(c_loc(flat(1)), arr, dims)
  end subroutine reshape_int_to_1D

  subroutine reshape_int_to_2D(flat, arr, dims)
    integer(int32), pointer, intent(in) :: flat(:)
    integer(int32), pointer, intent(out) :: arr(:,:)
    integer(int32), intent(in) :: dims(2)
    call c_f_pointer(c_loc(flat(1)), arr, dims)
  end subroutine reshape_int_to_2D

  subroutine reshape_int_to_3D(flat, arr, dims)
    integer(int32), pointer, intent(in) :: flat(:)
    integer(int32), pointer, intent(out) :: arr(:,:,:)
    integer(int32), intent(in) :: dims(3)
    call c_f_pointer(c_loc(flat(1)), arr, dims)
  end subroutine reshape_int_to_3D

  subroutine reshape_int_to_4D(flat, arr, dims)
    integer(int32), pointer, intent(in) :: flat(:)
    integer(int32), pointer, intent(out) :: arr(:,:,:,:)
    integer(int32), intent(in) :: dims(4)
    call c_f_pointer(c_loc(flat(1)), arr, dims)
  end subroutine reshape_int_to_4D

  subroutine reshape_int_to_5D(flat, arr, dims)
    integer(int32), pointer, intent(in) :: flat(:)
    integer(int32), pointer, intent(out) :: arr(:,:,:,:,:)
    integer(int32), intent(in) :: dims(5)
    call c_f_pointer(c_loc(flat(1)), arr, dims)
  end subroutine reshape_int_to_5D

  ! ============================================
  ! Reshape-Routinen für CHARACTER-Arrays
  ! ============================================
  subroutine reshape_char_to_1D(flat, arr, dims, clen)
    character(len=*), pointer, intent(in) :: flat(:)
    character(len=*), pointer, intent(out) :: arr(:)
    integer(int32), intent(in) :: dims(1)
    integer, intent(in) :: clen
    integer :: i
    allocate(character(len=clen) :: arr(dims(1)))
    do i = 1, dims(1)
      arr(i) = flat(i)
    end do
  end subroutine reshape_char_to_1D

  subroutine reshape_char_to_2D(flat, arr, dims, clen)
    character(len=*), pointer, intent(in) :: flat(:)
    character(len=*), pointer, intent(out) :: arr(:,:)
    integer(int32), intent(in) :: dims(2)
    integer, intent(in) :: clen
    integer :: i, j, idx
    allocate(character(len=clen) :: arr(dims(1), dims(2)))
    idx = 1
    do j = 1, dims(2)
      do i = 1, dims(1)
        arr(i,j) = flat(idx)
        idx = idx + 1
      end do
    end do
  end subroutine reshape_char_to_2D

  subroutine reshape_char_to_3D(flat, arr, dims, clen)
    character(len=*), pointer, intent(in) :: flat(:)
    character(len=*), pointer, intent(out) :: arr(:,:,:)
    integer(int32), intent(in) :: dims(3)
    integer, intent(in) :: clen
    integer :: i, j, k, idx
    allocate(character(len=clen) :: arr(dims(1), dims(2), dims(3)))
    idx = 1
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          arr(i,j,k) = flat(idx)
          idx = idx + 1
        end do
      end do
    end do
  end subroutine reshape_char_to_3D

  subroutine reshape_char_to_4D(flat, arr, dims, clen)
    character(len=*), pointer, intent(in) :: flat(:)
    character(len=*), pointer, intent(out) :: arr(:,:,:,:)
    integer(int32), intent(in) :: dims(4)
    integer, intent(in) :: clen
    integer :: i, j, k, l, idx
    allocate(character(len=clen) :: arr(dims(1), dims(2), dims(3), dims(4)))
    idx = 1
    do l = 1, dims(4)
      do k = 1, dims(3)
        do j = 1, dims(2)
          do i = 1, dims(1)
            arr(i,j,k,l) = flat(idx)
            idx = idx + 1
          end do
        end do
      end do
    end do
  end subroutine reshape_char_to_4D

  subroutine reshape_char_to_5D(flat, arr, dims, clen)
    character(len=*), pointer, intent(in) :: flat(:)
    character(len=*), pointer, intent(out) :: arr(:,:,:,:,:)
    integer(int32), intent(in) :: dims(5)
    integer, intent(in) :: clen
    integer :: i, j, k, l, m, idx
    allocate(character(len=clen) :: arr(dims(1), dims(2), dims(3), dims(4), dims(5)))
    idx = 1
    do m = 1, dims(5)
      do l = 1, dims(4)
        do k = 1, dims(3)
          do j = 1, dims(2)
            do i = 1, dims(1)
              arr(i,j,k,l,m) = flat(idx)
              idx = idx + 1
            end do
          end do
        end do
      end do
    end do
  end subroutine reshape_char_to_5D

end module reshape_utils
