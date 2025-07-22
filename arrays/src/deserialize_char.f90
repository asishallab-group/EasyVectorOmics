!> Module for deserializing character arrays from files
module char_deserialize_mod
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding
  implicit none

  private
  public :: deserialize_char, deserialize_char_flat

  integer(int32), parameter :: ARRAY_FILE_MAGIC = int(z'46413230', int32) ! 'FA20' in hex

  ! Generic interface for deserialization of character arrays
  interface deserialize_char
    module procedure deserialize_char_1d
    module procedure deserialize_char_2d
    module procedure deserialize_char_3d
    module procedure deserialize_char_4d
    module procedure deserialize_char_5d
  end interface

contains
  !> @brief Subroutine to deserialize a flat character array from a file
  !> @param flat Output flat character array
  !> @param dims Output dimensions of the array
  !> @param clen Output maximum length of character string
  !> @param filename Name of the file to read
  subroutine deserialize_char_flat(flat, dims, clen, filename)
    use iso_c_binding
    character(len=:), pointer, intent(out) :: flat(:)
    integer(int32), allocatable, intent(out) :: dims(:)
    integer, intent(out) :: clen
    character(len=*), intent(in) :: filename

    integer :: unit, magic, type_code, d, i, str_len
    character(len=:), allocatable :: temp_str

    ! open file and read header
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='old')
    read(unit) magic
    if (magic /= ARRAY_FILE_MAGIC) error stop "Invalid file format"
    read(unit) type_code
    if (type_code /= 3) error stop "Expected character data"
    read(unit) d
    allocate(dims(d))
    read(unit) dims
    read(unit) clen
    !allocate proper length for flat array
    allocate(character(len=clen) :: flat(product(dims)))
    do i = 1, product(dims)
      read(unit) str_len
      if (str_len > 0) then
        allocate(character(len=str_len) :: temp_str)
        read(unit) temp_str
        flat(i) = temp_str
        deallocate(temp_str)
      else
        flat(i) = ''
      end if
    end do
    close(unit)
  end subroutine deserialize_char_flat

  !> @brief Subroutine to deserialize a 1D character array from a file
  !> @param arr Output character array
  !> @param filename Name of the file to read
  !> @note The array is read is flat and therefore 1D, this is just for compatibility and the generic interface
  !> @note All it does is move a pointer
  subroutine deserialize_char_1d(arr, filename)
    use iso_c_binding
    character(len=:), pointer, intent(out) :: arr(:)
    character(len=*), intent(in) :: filename

    character(len=:), pointer :: flat(:)
    integer(int32), allocatable :: dims(:)
    integer :: clen

    if (associated(arr)) nullify(arr)
    call deserialize_char_flat(flat, dims, clen, filename)
    if (size(dims) /= 1) error stop "Expected 1D array"
    arr => flat
  end subroutine

  !> @brief Subroutine to deserialize a 2D character array from a file
  !> @param arr Output character array
  !> @param filename Name of the file to read
  !!!> The array is read as flat and then reshaped to 2D
  subroutine deserialize_char_2d(arr, filename)
    use iso_c_binding
    character(len=:), pointer, intent(out) :: arr(:,:)
    character(len=*), intent(in) :: filename

    character(len=:), pointer :: flat(:)
    integer(int32), allocatable :: dims(:)
    integer :: clen, i, j, idx

    call deserialize_char_flat(flat, dims, clen, filename)
    if (size(dims) /= 2) error stop "Expected 2D array"
    allocate(character(len=clen) :: arr(dims(1), dims(2)))
    ! Reshape flat array to 2D, c_f_pointer is not used here since it is note very stable for character arrays
    idx = 1
    do j = 1, dims(2)
      do i = 1, dims(1)
        arr(i,j) = flat(idx)
        idx = idx + 1
      end do
    end do
    ! Deallocate flat array to avoid memory leaks
    deallocate(flat)
  end subroutine

  !> @brief Subroutine to deserialize a 3D character array from a file
  !> @param arr Output character array
  !> @param filename Name of the file to read
  !!!> The array is read as flat and then reshaped to 3D
  subroutine deserialize_char_3d(arr, filename)
    use iso_c_binding
    character(len=:), pointer, intent(out) :: arr(:,:,:)
    character(len=*), intent(in) :: filename

    character(len=:), pointer :: flat(:)
    integer(int32), allocatable :: dims(:)
    integer :: clen, i, j, k, idx

    !Read file
    call deserialize_char_flat(flat, dims, clen, filename)
    if (size(dims) /= 3) error stop "Expected 3D array"
    !Allocate 3D array
    allocate(character(len=clen) :: arr(dims(1), dims(2), dims(3)))
    !Reshape flat array to 3D
    idx = 1
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          arr(i,j,k) = flat(idx)
          idx = idx + 1
        end do
      end do
    end do
    !Deallocate flat array to avoid memory leaks
    deallocate(flat)
  end subroutine

  !> @brief Subroutine to deserialize a 4D character array from a file
  !> @param arr Output character array
  !> @param filename Name of the file to read
  !!!> The array is read as flat and then reshaped to 4D
  subroutine deserialize_char_4d(arr, filename)
    use iso_c_binding
    character(len=:), pointer, intent(out) :: arr(:,:,:,:)
    character(len=*), intent(in) :: filename

    character(len=:), pointer :: flat(:)
    integer(int32), allocatable :: dims(:)
    integer :: clen, i, j, k, l, idx
    !Read file
    call deserialize_char_flat(flat, dims, clen, filename)
    if (size(dims) /= 4) error stop "Expected 4D array"
    !Allocate 4D array
    allocate(character(len=clen) :: arr(dims(1), dims(2), dims(3), dims(4)))
    !Reshape flat array to 4D
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
    !Deallocate flat array to avoid memory leaks
    deallocate(flat)
  end subroutine

  !> @brief Subroutine to deserialize a 5D character array from a file
  !> @param arr Output character array
  !> @param filename Name of the file to read
  !!!> The array is read as flat and then reshaped to 5D
  subroutine deserialize_char_5d(arr, filename)
      character(len=:), pointer, intent(out) :: arr(:,:,:,:,:)
      character(len=*), intent(in) :: filename

      character(len=:), pointer :: flat(:)
      integer(int32), allocatable :: dims(:)
      integer :: max_len, i, j, k, l, m, idx

      !Avoid memory leaks
      if (associated(arr)) nullify(arr)
      !Read file
      call deserialize_char_flat(flat, dims, max_len, filename)
      if (size(dims) /= 5) error stop "Expected 5D array"
      !Allocate 5D array
      allocate(character(len=max_len) :: arr(dims(1), dims(2), dims(3), dims(4), dims(5)))
      !RESHAPE
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

      !Deallocate flat array to avoid memory leaks
      deallocate(flat) 
  end subroutine
end module char_deserialize_mod

!> @brief Subroutine to deserialize a flat character array from a file and return it as an ASCII array callable by R
!> @param ascii_arr Output array of ASCII characters, preallocated by R
!> @param arr_size Size of the output array allocated by R
!> @param dims_out Output array for dimensions
!> @param ndim_out Output variable for the number of dimensions
!> @param clen_out Output variable for the character length
!> @param filename_ascii Array of ASCII characters representing the filename
!> @param fn_len Length of the filename array
!> @param ndim_actual Actual number of dimensions expected in the intput
!> @note The array is returned flat and needs to be reshaped in R
subroutine deserialize_char_flat_r(ascii_arr, arr_size, dims_out, ndim_out, clen_out, filename_ascii, fn_len, ndim_actual)
  use iso_fortran_env
  use char_deserialize_mod
  implicit none

  ! Arrays are allocated by R
  integer(int32), intent(out) :: ascii_arr(arr_size)
  integer(int32), intent(out) :: dims_out(ndim_actual)
  integer, intent(out) :: ndim_out, clen_out
  integer(int32), intent(in) :: filename_ascii(fn_len)
  integer(int32), intent(in) :: fn_len
  integer(int32), intent(in) :: arr_size, ndim_actual

  character(len=:), allocatable :: filename
  character(len=:), pointer :: flat(:)
  integer(int32), allocatable :: dims(:)
  integer :: i, j, clen, total

  ! ASCII → String
  allocate(character(len=fn_len) :: filename)
  do i = 1, fn_len
    filename(i:i) = char(filename_ascii(i))
  end do

  ! Deserialize flat character array
  call deserialize_char_flat(flat, dims, clen, filename)
  total = product(dims)
  clen_out = clen
  ndim_out = size(dims)

  do i = 1, ndim_out
    dims_out(i) = dims(i)
  end do

  ! Write data to ASCII array
  do i = 1, total
    do j = 1, clen
      if (j <= len_trim(flat(i))) then
        ascii_arr((i - 1) * clen + j) = iachar(flat(i)(j:j))
      else
        ascii_arr((i - 1) * clen + j) = 0
      end if
    end do
  end do
end subroutine deserialize_char_flat_r

!> @brief C binding for the subroutine to deserialize a flat character array from a file
!> @param ascii_arr Output array of ASCII characters, preallocated by R
!> @param clen maximum length of character string
!> @param total Total number of strings in the array
!> @param dims_out Output array for dimensions
!> @param ndim_out Output variable for the number of dimensions
!> @param clen_out Output variable for the character length
!> @param filename_ascii Array of ASCII characters representing the filename
!> @param fn_len Length of the filename array
!> @param ndim_actual Actual number of dimensions expected in the input
!> @note The array is returned flat and needs to be reshaped in C/python
subroutine deserialize_char_flat_C(ascii_arr, clen, total, dims_out, ndim_out, clen_out, &
                                   filename_ascii, fn_len, ndim_actual) bind(C, name="deserialize_char_flat_C")
  use iso_c_binding
  use char_deserialize_mod
  implicit none

  ! Arguments
  integer(c_int), intent(out) :: ascii_arr(clen, total)
  integer(c_int), value       :: clen, total
  integer(c_int), intent(out) :: dims_out(ndim_actual)
  integer(c_int), intent(out) :: ndim_out, clen_out
  integer(c_int), intent(in)  :: filename_ascii(fn_len)
  integer(c_int), value       :: fn_len, ndim_actual

  character(len=:), allocatable :: filename
  character(len=:), pointer     :: flat(:)
  integer(c_int), allocatable   :: dims(:)
  integer :: i, j, actual_clen

  ! Convert Filename
  allocate(character(len=fn_len) :: filename)
  do i = 1, fn_len
    filename(i:i) = char(filename_ascii(i))
  end do

  ! Deserialize
  call deserialize_char_flat(flat, dims, actual_clen, filename)
  clen_out = actual_clen
  ndim_out = size(dims)

  do i = 1, min(ndim_out, ndim_actual)
    dims_out(i) = dims(i)
  end do

  ! Convert to ASCII (Null-padding)
  do i = 1, total
    do j = 1, clen
      if (j <= len_trim(flat(i))) then
        ascii_arr(j, i) = iachar(flat(i)(j:j))
      else
        ascii_arr(j, i) = 0
      end if
    end do
  end do
end subroutine
