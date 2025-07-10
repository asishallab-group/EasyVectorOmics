module int_deserialize_mod
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding, only: c_ptr, c_loc, c_char, c_null_char
  implicit none

  private
  public :: deserialize_int, deserialize_int_1d, deserialize_int_2d, &
           deserialize_int_3d, deserialize_int_4d, deserialize_int_5d, deserialize_int_flat

  integer(int32), parameter :: ARRAY_FILE_MAGIC = int(z'46413230', int32) ! 'FA20' in hex

  interface deserialize_int
    module procedure deserialize_int_1d
    module procedure deserialize_int_2d
    module procedure deserialize_int_3d
    module procedure deserialize_int_4d
    module procedure deserialize_int_5d
  end interface

contains

  ! Hilfsroutine: Flaches Array + Dimensionen lesen
  subroutine deserialize_int_flat(flat, dims, filename)
    use iso_c_binding
    integer(int32), pointer, intent(out) :: flat(:)
    integer(int32), allocatable, intent(out), target :: dims(:)
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

  ! Überladene Routinen: 1D
  subroutine deserialize_int_1d(arr, filename)
    use iso_c_binding
    integer(int32), pointer, intent(out) :: arr(:)
    character(len=*), intent(in) :: filename
    integer(int32), pointer :: flat(:)
    integer(int32), allocatable :: dims(:)
    call deserialize_int_flat(flat, dims, filename)
    if (size(dims) /= 1) error stop "Expected 1D array"
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1)])
  end subroutine deserialize_int_1d

  ! 2D
  subroutine deserialize_int_2d(arr, filename)
    use iso_c_binding
    integer(int32), pointer, intent(out) :: arr(:,:)
    character(len=*), intent(in) :: filename
    integer(int32), pointer :: flat(:)
    integer(int32), allocatable :: dims(:)
    call deserialize_int_flat(flat, dims, filename)
    if (size(dims) /= 2) error stop "Expected 2D array"
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1), dims(2)])
  end subroutine deserialize_int_2d

  ! 3D
  subroutine deserialize_int_3d(arr, filename)
    use iso_c_binding
    integer(int32), pointer, intent(out) :: arr(:,:,:)
    character(len=*), intent(in) :: filename
    integer(int32), pointer :: flat(:)
    integer(int32), allocatable :: dims(:)
    call deserialize_int_flat(flat, dims, filename)
    if (size(dims) /= 3) error stop "Expected 3D array"
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1), dims(2), dims(3)])
  end subroutine deserialize_int_3d

  ! 4D
  subroutine deserialize_int_4d(arr, filename)
    use iso_c_binding
    integer(int32), pointer, intent(out) :: arr(:,:,:,:)
    character(len=*), intent(in) :: filename
    integer(int32), pointer :: flat(:)
    integer(int32), allocatable :: dims(:)
    call deserialize_int_flat(flat, dims, filename)
    if (size(dims) /= 4) error stop "Expected 4D array"
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1), dims(2), dims(3), dims(4)])
  end subroutine deserialize_int_4d

  ! 5D
  subroutine deserialize_int_5d(arr, filename)
    use iso_c_binding
    integer(int32), pointer, intent(out) :: arr(:,:,:,:,:)
    character(len=*), intent(in) :: filename
    integer(int32), pointer :: flat(:)
    integer(int32), allocatable :: dims(:)
    call deserialize_int_flat(flat, dims, filename)
    if (size(dims) /= 5) error stop "Expected 5D array"
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1), dims(2), dims(3), dims(4), dims(5)])
  end subroutine deserialize_int_5d

end module int_deserialize_mod

subroutine deserialize_int_1d_r(arr, size, filename_ascii, fn_len)
  use int_deserialize_mod, only: deserialize_int_1d
  use iso_fortran_env, only: int32
  implicit none
  integer , INTENT(IN) :: size
  integer(int32), intent(out) :: arr(size)
  integer(int32), intent(in) :: filename_ascii(fn_len)
  integer(int32), pointer :: arr_f(:)
  integer(int32), intent(in) :: fn_len

  character(len=:), allocatable :: filename
  integer :: i

  allocate(character(len=fn_len) :: filename)

  do i = 1, fn_len
    filename(i:i) = char(filename_ascii(i))
  end do
  
  call deserialize_int_1d(arr_f, filename)

  arr = arr_f
end subroutine deserialize_int_1d_r

subroutine deserialize_int_flat_r(flat_arr, dims_out, ndim_out, filename_ascii, fn_len)
  use iso_fortran_env, only: int32
  use int_deserialize_mod
  implicit none

  ! Rückgabe an R
  integer(int32), intent(out) :: flat_arr(*)
  integer(int32), intent(out) :: dims_out(*)
  integer, intent(out) :: ndim_out

  ! Dateiname
  integer(int32), intent(in) :: filename_ascii(fn_len)
  integer, intent(in) :: fn_len

  ! Intern
  character(len=:), allocatable :: filename
  integer :: i, k
  integer(int32), pointer :: flat(:)
  integer(int32), allocatable, target :: dims(:)

  allocate(character(len=fn_len) :: filename)
  do i = 1, fn_len
    filename(i:i) = char(filename_ascii(i))
  end do

  ! Lese Daten
  call deserialize_int_flat(flat, dims, filename)

  ndim_out = size(dims)
  do i = 1, ndim_out
    dims_out(i) = dims(i)
  end do

  do i = 1, product(dims)
    flat_arr(i) = flat(i)
  end do
end subroutine


subroutine deserialize_int_1d_C(arr, filename) bind(C, name="deserialize_int_1d_C")
  use iso_c_binding, only: c_ptr, c_loc, c_char, c_null_char
  use int_deserialize_mod, only: deserialize_int_1d
  use iso_fortran_env, only: int32
  implicit none
  type(c_ptr), intent(out) :: arr
  character(kind=c_char), intent(in) :: filename(*)
  integer(int32), pointer :: arr_f(:)
  character(len=:), allocatable :: fname
  integer :: i
  arr_f => null()
  i = 1
  do while (filename(i) /= c_null_char)
    i = i + 1
  end do
  fname = transfer(filename(1:i-1), fname)
  call deserialize_int_1d(arr_f, fname)
  arr = c_loc(arr_f)
end subroutine deserialize_int_1d_C

subroutine deserialize_int_2d_C(arr, filename) bind(C, name="deserialize_int_2d_C")
  use iso_c_binding, only: c_ptr, c_loc, c_char, c_null_char
  use int_deserialize_mod, only: deserialize_int_2d
  use iso_fortran_env, only: int32
  implicit none
  type(c_ptr), intent(out) :: arr
  character(kind=c_char), intent(in) :: filename(*)
  integer(int32), pointer :: arr_f(:,:)
  character(len=:), allocatable :: fname
  integer :: i
  arr_f => null()
  i = 1
  do while (filename(i) /= c_null_char)
    i = i + 1
  end do
  fname = transfer(filename(1:i-1), fname)
  call deserialize_int_2d(arr_f, fname)
  arr = c_loc(arr_f)
end subroutine deserialize_int_2d_C

subroutine deserialize_int_3d_C(arr, filename) bind(C, name="deserialize_int_3d_C")
  use iso_c_binding, only: c_ptr, c_loc, c_char, c_null_char
  use int_deserialize_mod, only: deserialize_int_3d
  use iso_fortran_env, only: int32
  implicit none
  type(c_ptr), intent(out) :: arr
  character(kind=c_char), intent(in) :: filename(*)
  integer(int32), pointer :: arr_f(:,:,:)
  character(len=:), allocatable :: fname
  integer :: i
  arr_f => null()
  i = 1
  do while (filename(i) /= c_null_char)
    i = i + 1
  end do
  fname = transfer(filename(1:i-1), fname)
  call deserialize_int_3d(arr_f, fname)
  arr = c_loc(arr_f)
end subroutine deserialize_int_3d_C

subroutine deserialize_int_4d_C(arr, filename) bind(C, name="deserialize_int_4d_C")
  use iso_c_binding, only: c_ptr, c_loc, c_char, c_null_char
  use int_deserialize_mod, only: deserialize_int_4d
  use iso_fortran_env, only: int32
  implicit none
  type(c_ptr), intent(out) :: arr
  character(kind=c_char), intent(in) :: filename(*)
  integer(int32), pointer :: arr_f(:,:,:,:)
  character(len=:), allocatable :: fname
  integer :: i
  arr_f => null()
  i = 1
  do while (filename(i) /= c_null_char)
    i = i + 1
  end do
  fname = transfer(filename(1:i-1), fname)
  call deserialize_int_4d(arr_f, fname)
  arr = c_loc(arr_f)
end subroutine deserialize_int_4d_C

subroutine deserialize_int_5d_C(arr, filename) bind(C, name="deserialize_int_5d_C")
  use iso_c_binding, only: c_ptr, c_loc, c_char, c_null_char
  use int_deserialize_mod, only: deserialize_int_5d
  use iso_fortran_env, only: int32
  implicit none
  type(c_ptr), intent(out) :: arr
  character(kind=c_char), intent(in) :: filename(*)
  integer(int32), pointer :: arr_f(:,:,:,:,:)
  character(len=:), allocatable :: fname
  integer :: i
  arr_f => null()
  i = 1
  do while (filename(i) /= c_null_char)
    i = i + 1
  end do
  fname = transfer(filename(1:i-1), fname)
  call deserialize_int_5d(arr_f, fname)
  arr = c_loc(arr_f)
end subroutine deserialize_int_5d_C