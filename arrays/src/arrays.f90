module array_ops
  use, intrinsic :: iso_fortran_env, only: int32, real64
  implicit none

  public::serialize, deserialize

  interface serialize
    module procedure serialize_int_array
    module procedure serialize_real_array
    module procedure serialize_char_array
  end interface serialize

  interface deserialize
    module procedure deserialize_int_array
    module procedure deserialize_real_array
    module procedure deserialize_char_array
  end interface deserialize

  integer(int32), parameter :: ARRAY_FILE_MAGIC = int(z'46413230', int32) ! 'FA20' in hex

contains

  !============================================
  ! Serialization Functions
  !============================================

  subroutine serialize_int_array(arr, filename)
    integer(int32), intent(in) :: arr(..)
    character(len=*), intent(in) :: filename
    integer :: unit
    integer :: d
    integer(int32), allocatable :: dims(:)
    integer(int32), allocatable :: flat(:)

    d = rank(arr)
    allocate(dims(d))
    dims = shape(arr)

    select rank(arr)
    rank(1)
      allocate(flat(size(arr)))
      flat = arr
    class default
      allocate(flat(size(arr)))
      flat = reshape(arr, [size(arr)])
    end select

    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    write(unit) ARRAY_FILE_MAGIC
    write(unit) 1  ! Type code for integer
    write(unit) d
    write(unit) dims
    write(unit) flat
    close(unit)
  end subroutine serialize_int_array

  subroutine serialize_real_array(arr, filename)
    real(real64), intent(in) :: arr(..)
    character(len=*), intent(in) :: filename
    integer :: unit, d
    integer(int32), allocatable :: dims(:)
    real(real64), allocatable :: flat(:)

    d = rank(arr)
    allocate(dims(d))
    dims = shape(arr)

    select rank(arr)
    rank(1)
      allocate(flat(size(arr)))
      flat = arr
    class default
      allocate(flat(size(arr)))
      flat = reshape(arr, [size(arr)])
    end select

    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    write(unit) ARRAY_FILE_MAGIC
    write(unit) 2  ! Type code for real
    write(unit) d
    write(unit) dims
    write(unit) flat
    close(unit)
  end subroutine serialize_real_array

  subroutine serialize_char_array(arr, filename)
    character(len=*), intent(in) :: arr(..)
    character(len=*), intent(in) :: filename
    integer :: unit, d, clen, i
    integer(int32), allocatable :: dims(:)
    character(len=:), allocatable :: flat(:)

    d = rank(arr)
    allocate(dims(d))
    dims = shape(arr)
    clen = len(arr)

    select rank(arr)
    rank(1)
      allocate(character(len=clen) :: flat(size(arr)))
      flat = arr
    class default
      allocate(character(len=clen) :: flat(size(arr)))
      flat = reshape(arr, [size(arr)])
    end select

    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace')
    write(unit) ARRAY_FILE_MAGIC
    write(unit) 3  ! Type code for character
    write(unit) d
    write(unit) dims
    write(unit) clen
    do i = 1, product(dims)
      write(unit) flat(i)
    end do
    close(unit)
  end subroutine serialize_char_array

  !============================================
  ! Deserialization Functions
  !============================================

  subroutine deserialize_int_array(arr, filename)
    integer(int32), allocatable, intent(out) :: arr(..)
    character(len=*), intent(in) :: filename
    integer :: unit, magic, type_code, d, i
    integer(int32), allocatable :: dims(:)
    integer(int32), allocatable :: flat(:)

    open(newunit=unit, file=filename, form='unformatted', access='stream', status='old')
    read(unit) magic
    if (magic /= ARRAY_FILE_MAGIC) error stop "Invalid file format"
    read(unit) type_code
    if (type_code /= 1) error stop "File contains non-integer data"
    read(unit) d
    allocate(dims(d))
    read(unit) dims
    allocate(flat(product(dims)))
    read(unit) flat
    select case(d)
    case(1)
      allocate(arr(dims(1)))
      arr = reshape(flat, dims)
    case(2)
      allocate(arr(dims(1), dims(2)))
      arr = reshape(flat, dims)
    case default
      allocate(arr(dims))
      arr = reshape(flat, dims)
    end select
    close(unit)
  end subroutine deserialize_int_array

  subroutine deserialize_real_array(arr, filename)
    real(real64), allocatable, intent(out) :: arr(..)
    character(len=*), intent(in) :: filename
    integer :: unit, magic, type_code, d, i
    integer(int32), allocatable :: dims(:)
    real(real64), allocatable :: flat(:)

    open(newunit=unit, file=filename, form='unformatted', access='stream', status='old')
    read(unit) magic
    if (magic /= ARRAY_FILE_MAGIC) error stop "Invalid file format"
    read(unit) type_code
    if (type_code /= 2) error stop "File contains non-real data"
    read(unit) d
    allocate(dims(d))
    read(unit) dims
    allocate(flat(product(dims)))
    read(unit) flat
    select case(d)
    case(1)
      allocate(arr(dims(1)))
      arr = reshape(flat, dims)
    case(2)
      allocate(arr(dims(1), dims(2)))
      arr = reshape(flat, dims)
    case default
      allocate(arr(dims))
      arr = reshape(flat, dims)
    end select
    close(unit)
  end subroutine deserialize_real_array

  subroutine deserialize_char_array(arr, filename)
    character(len=:), allocatable, intent(out) :: arr(..)
    character(len=*), intent(in) :: filename
    integer :: unit, magic, type_code, d, i, clen
    integer(int32), allocatable :: dims(:)
    character(len=:), allocatable :: flat(:)

    open(newunit=unit, file=filename, form='unformatted', access='stream', status='old')
    read(unit) magic
    if (magic /= ARRAY_FILE_MAGIC) error stop "Invalid file format"
    read(unit) type_code
    if (type_code /= 3) error stop "File contains non-character data"
    read(unit) d
    allocate(dims(d))
    read(unit) dims
    read(unit) clen
    allocate(character(len=clen) :: flat(product(dims)))
    do i = 1, product(dims)
      read(unit) flat(i)
    end do
    select case(d)
    case(1)
      allocate(character(len=clen) :: arr(dims(1)))
      arr = reshape(flat, dims)
    case(2)
      allocate(character(len=clen) :: arr(dims(1), dims(2)))
      arr = reshape(flat, dims)
    case default
      allocate(character(len=clen) :: arr(dims))
      arr = reshape(flat, dims)
    end select
    close(unit)
  end subroutine deserialize_char_array

  !============================================
  ! Accessor Functions
  !============================================

  ! Get a full row (returns a 1D array)
  function get_row(arr, row_index) result(row)
    integer(int32), intent(in) :: arr(:,:)
    integer, intent(in) :: row_index
    integer(int32) :: row(size(arr, 2))

    row = arr(row_index, :)
  end function get_row

  ! Get a full column (returns a 1D array)
  function get_col(arr, col_index) result(col)
    integer(int32), intent(in) :: arr(:,:)
    integer, intent(in) :: col_index
    integer(int32) :: col(size(arr, 1))

    col = arr(:, col_index)
  end function get_col

  ! Get a specific cell value
  function get_cell(arr, row_index, col_index) result(val)
    integer(int32), intent(in) :: arr(:,:)
    integer, intent(in) :: row_index, col_index
    integer(int32) :: val

    val = arr(row_index, col_index)
  end function get_cell

  ! Get column by name (assuming column_names array exists)
  function get_col_for_name(arr, column_names, name) result(col)
    integer(int32), intent(in) :: arr(:,:)
    character(len=*), intent(in) :: column_names(:), name
    integer(int32) :: col(size(arr, 1))
    integer :: i

    do i = 1, size(column_names)
      if (trim(column_names(i)) == trim(name)) then
        col = arr(:, i)
        return
      end if
    end do
    error stop "Column name not found"
  end function get_col_for_name

  ! Get cell by column name
  function get_cell_for_column_name(arr, column_names, col_name, row_index) result(val)
    integer(int32), intent(in) :: arr(:,:)
    character(len=*), intent(in) :: column_names(:), col_name
    integer, intent(in) :: row_index
    integer(int32) :: val
    integer :: i

    do i = 1, size(column_names)
      if (trim(column_names(i)) == trim(col_name)) then
        val = arr(row_index, i)
        return
      end if
    end do
    error stop "Column name not found"
  end function get_cell_for_column_name

end module array_ops

subroutine serialize_int_array_r(arr, filename)
  use array_ops
  implicit none
  integer(int32), intent(in) :: arr(..)
  character(len=*), intent(in) :: filename
    
  call serialize_int_array(arr, filename)
end subroutine serialize_int_array_r
