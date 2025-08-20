!> Module for deserializing character arrays from files
module char_deserialize_mod
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use array_utils, only : ascii_to_string, read_file_header
  use tox_errors
  implicit none

  private
  public :: deserialize_char_1d, deserialize_char_2d, deserialize_char_3d, deserialize_char_flat, &
          deserialize_char_4d, deserialize_char_5d

contains
  !> Subroutine to deserialize a flat character array from a file
  subroutine deserialize_char_flat(flat, dims, clen, filename, ierr)
    character(len=:), pointer, intent(out) :: flat(:)
      !! Output flat character array
    integer(int32), allocatable, intent(out) :: dims(:)
      !! Output dimensions of the array
    integer(int32), intent(out) :: clen
      !! Maximum length of character string
    character(len=*), intent(in) :: filename
      !! Name of the file to read
    integer(int32), intent(out) :: ierr
      !! Error code
    integer(int32) :: ioerror
      !! Internal I/O error code

    integer(int32) :: unit, magic, type_code, ndim, i, str_len
    character(len=:), allocatable :: temp_str

    call set_ok(ierr)
    call set_ok(ioerror)
    ! open file and read header
    call read_file_header(filename, unit, type_code, ndim, dims, clen, ierr)
    if (.not. is_ok(ierr)) then
      return
    end if
    !allocate proper length for flat array
    allocate(character(len=clen) :: flat(product(dims)))
    do i = 1, product(dims)
      read(unit, iostat=ioerror) str_len
      if (.not. is_ok(ioerror)) then
        call set_err_once(ierr, ERR_READ_CHARLEN)
        close(unit)
        return
      end if
      if (str_len > 0) then
        allocate(character(len=str_len) :: temp_str)
        read(unit, iostat=ioerror) temp_str
        if (.not. is_ok(ioerror)) then
          call set_err_once(ierr, ERR_READ_DATA)
          close(unit)
          return
        end if
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
  subroutine deserialize_char_1d(arr, filename, ierr)
    character(len=:), pointer, intent(out) :: arr(:)
      !! Output character array
    character(len=*), intent(in) :: filename
      !! Name of the file to read

    character(len=:), pointer :: flat(:)
      !! Flat character array
    integer(int32), allocatable :: dims(:)
      !! Output dimensions of the array
    integer(int32) :: clen
      !! Maximum length of character string
    integer(int32), intent(out) :: ierr
      !! Error code
    call set_ok(ierr)

    if (associated(arr)) nullify(arr)
    call deserialize_char_flat(flat, dims, clen, filename, ierr)
    if (size(dims) /= 1) then
      call set_err_once(ierr, ERR_DIM_MISMATCH)
      if (associated(flat)) nullify(flat)
      return
    end if
    arr => flat
  end subroutine

  !> Subroutine to deserialize a 2D character array from a file
  !!!> The array is read as flat and then reshaped to 2D
  subroutine deserialize_char_2d(arr, filename, ierr)
    character(len=:), pointer, intent(out) :: arr(:,:)
    !! Output character array
    character(len=*), intent(in) :: filename
    !! Name of the file to read

    character(len=:), pointer :: flat(:)
    !! Flat character array
    integer(int32), allocatable :: dims(:)
    !! Output dimensions of the array
    integer(int32) :: clen
    !! Maximum length of character string
    integer(int32) :: i, j, idx
    !! Loop indices
    integer(int32), intent(out) :: ierr
    !! Error code

    call set_ok(ierr)

    if (associated(arr)) nullify(arr)

    call deserialize_char_flat(flat, dims, clen, filename, ierr)

    if (size(dims) /= 2) then
      call set_err_once(ierr, ERR_DIM_MISMATCH)
      if (associated(flat)) nullify(flat)
      return
    end if

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
  subroutine deserialize_char_3d(arr, filename, ierr)
    character(len=:), pointer, intent(out) :: arr(:,:,:)
    !! Output character array
    character(len=*), intent(in) :: filename
    !! Name of the file to read

    character(len=:), pointer :: flat(:)
    !! Flat character array
    integer(int32), allocatable :: dims(:)
    !! Output dimensions of the array
    integer(int32) :: clen
    !! Maximum length of character string
    integer(int32) :: i, j, k, idx
    !! Loop indices
    integer(int32), intent(out) :: ierr
    !! Error code

    call set_ok(ierr)
    !Read file
    call deserialize_char_flat(flat, dims, clen, filename, ierr)
    if (size(dims) /= 3) then
      call set_err_once(ierr, ERR_DIM_MISMATCH)
      if (associated(flat)) nullify(flat)
      return
    end if
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
  subroutine deserialize_char_4d(arr, filename, ierr)
    character(len=:), pointer, intent(out) :: arr(:,:,:,:)
    !! Output character array
    character(len=*), intent(in) :: filename
    !! Name of the file to read

    character(len=:), pointer :: flat(:)
    !! Flat character array
    integer(int32), allocatable :: dims(:)
    !! Output dimensions of the array
    integer(int32) :: clen
    !! Maximum length of character string
    integer(int32) :: i, j, k, l, idx
    !! Loop indices
    integer(int32), intent(out) :: ierr
    !! Error code

    call set_ok(ierr)

    !Read file
    call deserialize_char_flat(flat, dims, clen, filename, ierr)
    if (size(dims) /= 4) then
      call set_err_once(ierr, ERR_DIM_MISMATCH)
      if (associated(flat)) nullify(flat)
      return
    end if
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
  subroutine deserialize_char_5d(arr, filename, ierr)
      character(len=:), pointer, intent(out) :: arr(:,:,:,:,:)
      !! Output character array
      character(len=*), intent(in) :: filename
      !! Name of the file to read

      character(len=:), pointer :: flat(:)
      !! Flat character array
      integer(int32), allocatable :: dims(:)
      !! Output dimensions of the array
      integer(int32) :: clen 
      !! Maximum length of character string
      integer(int32) :: i, j, k, l, m, idx
      !! Loop indices
      integer(int32), intent(out) :: ierr
      !! Error code

      call set_ok(ierr)
      !Avoid memory leaks
      if (associated(arr)) nullify(arr)
      !Read file
      call deserialize_char_flat(flat, dims, clen, filename, ierr)
      if (size(dims) /= 5) then
          call set_err_once(ierr, ERR_DIM_MISMATCH)
          if (associated(flat)) nullify(flat)
          return
      end if
      !Allocate 5D array
      allocate(character(len=clen) :: arr(dims(1), dims(2), dims(3), dims(4), dims(5)))
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
subroutine deserialize_char_flat_r(ascii_arr, arr_size, filename_ascii, fn_len, ierr)
  use iso_fortran_env, only: int32
  use char_deserialize_mod, only: deserialize_char_flat
  use array_utils, only: ascii_to_string, string_to_ascii_arr
  use tox_errors, only : set_ok, is_ok
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
  integer(int32), intent(out) :: ierr
  !! Error code

  character(len=:), allocatable :: filename
  !! Filename as a string
  character(len=:), pointer :: flat(:)
  !! Flat character array
  integer(int32), allocatable :: dims(:)
  !! Output dimensions of the array
  integer(int32) :: clen
  !! Maximum length of character string
  integer(int32) :: total_array_size
  !! Total size of the ASCII array
  integer(int32) :: i, j
  !! Loop indices

  call set_ok(ierr)

  call ascii_to_string(filename_ascii, fn_len, filename)

  ! Deserialize flat character array
  call deserialize_char_flat(flat, dims, clen, filename, ierr)
  if(.not. is_ok(ierr)) then
    DEALLOCATE(flat)
    return
  end if
  total_array_size = product(dims)

  ! Write data to ASCII array
  call string_to_ascii_arr(flat, total_array_size, ascii_arr, clen)

  deallocate(flat)
end subroutine deserialize_char_flat_r

!> C binding for the subroutine to deserialize a flat character array from a file
!> @note The array is returned flat and needs to be reshaped in C/python
subroutine deserialize_char_flat_C(ascii_arr, clen, total_array_size, &
                                   filename_ascii, fn_len, ierr) bind(C, name="deserialize_char_flat_C")
  use iso_c_binding, only: c_int
  use iso_fortran_env, only: int32
  use char_deserialize_mod, only: deserialize_char_flat
  use array_utils, only : ascii_to_string, string_to_ascii_arr
  use tox_errors, only : is_ok, set_ok
  implicit none

  ! Arguments
  integer(c_int), intent(out) :: ascii_arr(clen*total_array_size)
  !! Output array of ASCII characters, preallocated by C/Python (flat)
  integer(c_int), value       :: clen
  !! Length of each character string
  integer(c_int), value       :: total_array_size
  !! Total size of the ASCII array
  integer(c_int), intent(in)  :: filename_ascii(fn_len)
  !! Array of ASCII characters representing the filename
  integer(c_int), value       :: fn_len
  !! Length of the filename array
  integer(c_int), intent(out) :: ierr
  !! error code

  character(len=:), allocatable :: filename
  character(len=:), pointer     :: flat(:)
  integer(c_int), allocatable   :: dims(:)
  integer(int32) :: i, j, actual_clen

  call set_ok(ierr)
  call ascii_to_string(filename_ascii, fn_len, filename)

  ! Deserialize
  call deserialize_char_flat(flat, dims, actual_clen, filename, ierr)
  if (.not. is_ok(ierr)) then
    deallocate(flat)
    return
  end if

  ! Convert to ASCII (null-padding), flatten manually
  call string_to_ascii_arr(flat, total_array_size, ascii_arr, clen)

  deallocate(flat)
end subroutine