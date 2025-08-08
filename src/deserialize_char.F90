!> Module for deserializing character arrays from files
module char_deserialize_mod
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use array_utils
  use iso_c_binding
  implicit none

  private
  public :: deserialize_char_1d, deserialize_char_2d, deserialize_char_3d, deserialize_char_flat, &
          deserialize_char_4d, deserialize_char_5d

contains
  !> Subroutine to deserialize a flat character array from a file
  subroutine deserialize_char_flat(flat, dims, clen, filename)
    use iso_c_binding
    character(len=:), pointer, intent(out) :: flat(:)
      !! Output flat character array
    integer(int32), allocatable, intent(out) :: dims(:)
      !! Output dimensions of the array
    integer, intent(out) :: clen
      !! Maximum length of character string
    character(len=*), intent(in) :: filename
      !! Name of the file to read

    integer(int32) :: unit, magic, type_code, ndim, i, str_len, ierr
    character(len=:), allocatable :: temp_str

    ! open file and read header
    open(newunit=unit, file=filename, form='unformatted', access='stream', status='old', iostat=ierr)
    call check_file_header(filename, unit, type_code, ndim, dims, clen, ierr)
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

  !> Subroutine to deserialize a 1D character array from a file
  !> only moves a pointer for completeness
  subroutine deserialize_char_1d(arr, filename)
    use iso_c_binding
    character(len=:), pointer, intent(out) :: arr(:)
      !! Output character array
    character(len=*), intent(in) :: filename
      !! Name of the file to read

    character(len=:), pointer :: flat(:)
      !! Flat character array
    integer(int32), allocatable :: dims(:)
      !! Output dimensions of the array
    integer :: clen
      !! Maximum length of character string

    if (associated(arr)) nullify(arr)
    call deserialize_char_flat(flat, dims, clen, filename)
    if (size(dims) /= 1) error stop "Expected 1D array"
    arr => flat
  end subroutine

  !> Subroutine to deserialize a 2D character array from a file
  !!!> The array is read as flat and then reshaped to 2D
  subroutine deserialize_char_2d(arr, filename)
    use iso_c_binding
    character(len=:), pointer, intent(out) :: arr(:,:)
    !! Output character array
    character(len=*), intent(in) :: filename
    !! Name of the file to read

    character(len=:), pointer :: flat(:)
    !! Flat character array
    integer(int32), allocatable :: dims(:)
    !! Output dimensions of the array
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

  !> Subroutine to deserialize a 3D character array from a file
  subroutine deserialize_char_3d(arr, filename)
    use iso_c_binding
    character(len=:), pointer, intent(out) :: arr(:,:,:)
    !! Output character array
    character(len=*), intent(in) :: filename
    !! Name of the file to read

    character(len=:), pointer :: flat(:)
    !! Flat character array
    integer(int32), allocatable :: dims(:)
    !! Output dimensions of the array
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

  !> Subroutine to deserialize a 4D character array from a file
  subroutine deserialize_char_4d(arr, filename)
    use iso_c_binding
    character(len=:), pointer, intent(out) :: arr(:,:,:,:)
    !! Output character array
    character(len=*), intent(in) :: filename
    !! Name of the file to read

    character(len=:), pointer :: flat(:)
    !! Flat character array
    integer(int32), allocatable :: dims(:)
    !! Output dimensions of the array
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

  !> Subroutine to deserialize a 5D character array from a file
  subroutine deserialize_char_5d(arr, filename)
      character(len=:), pointer, intent(out) :: arr(:,:,:,:,:)
      !! Output character array
      character(len=*), intent(in) :: filename
      !! Name of the file to read

      character(len=:), pointer :: flat(:)
      !! Flat character array
      integer(int32), allocatable :: dims(:)
      !! Output dimensions of the array
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

!> Subroutine to deserialize a flat character array from a file and return it as an ASCII array callable by R
!> @note The array is returned flat and needs to be reshaped in R
subroutine deserialize_char_flat_r(ascii_arr, arr_size, filename_ascii, fn_len)
  use iso_fortran_env
  use char_deserialize_mod
  use array_utils
  implicit none

  ! Arrays are allocated by R
  integer(int32), intent(out) :: ascii_arr(arr_size)
  !! Output array of ASCII characters, preallocated by R
  integer(int32), intent(in) :: filename_ascii(fn_len)
  !! Array of ASCII characters representing the filename
  integer(int32), intent(in) :: fn_len
  !! Length of the filename array
  integer(int32), intent(in) :: arr_size
  !! Size of the ASCII array

  character(len=:), allocatable :: filename
  character(len=:), pointer :: flat(:)
  integer(int32), allocatable :: dims(:)
  integer :: i, j, clen, total_array_size

  call ascii_to_string(filename_ascii, fn_len, filename)

  ! Deserialize flat character array
  call deserialize_char_flat(flat, dims, clen, filename)
  total_array_size = product(dims)

  ! Write data to ASCII array
  do i = 1, total_array_size
    do j = 1, clen
      if (j <= len_trim(flat(i))) then
        ascii_arr((i - 1) * clen + j) = iachar(flat(i)(j:j))
      else
        ascii_arr((i - 1) * clen + j) = 0
      end if
    end do
  end do
end subroutine deserialize_char_flat_r

!> C binding for the subroutine to deserialize a flat character array from a file
!> @note The array is returned flat and needs to be reshaped in C/python
subroutine deserialize_char_flat_C(ascii_arr, clen, total_array_size, &
                                   filename_ascii, fn_len) bind(C, name="deserialize_char_flat_C")
  use iso_c_binding
  use char_deserialize_mod
  use array_utils
  implicit none

  ! Arguments
  integer(c_int), intent(out) :: ascii_arr(clen, total_array_size)
  !! Output array of ASCII characters, preallocated by C/Python
  integer(c_int), value       :: clen
  !! Length of each character string
  integer(c_int), value :: total_array_size
  !! Total size of the ASCII array
  integer(c_int), intent(in)  :: filename_ascii(fn_len)
  !! Array of ASCII characters representing the filename
  integer(c_int), value       :: fn_len
  !! Length of the filename array

  character(len=:), allocatable :: filename
  character(len=:), pointer     :: flat(:)
  integer(c_int), allocatable   :: dims(:)
  integer :: i, j, actual_clen

  call ascii_to_string(filename_ascii, fn_len, filename)

  ! Deserialize
  call deserialize_char_flat(flat, dims, actual_clen, filename)

  ! Convert to ASCII (Null-padding)
  do i = 1, total_array_size
    do j = 1, clen
      if (j <= len_trim(flat(i))) then
        ascii_arr(j, i) = iachar(flat(i)(j:j))
      else
        ascii_arr(j, i) = 0
      end if
    end do
  end do
end subroutine
