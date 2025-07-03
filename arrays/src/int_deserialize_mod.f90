module int_deserialize_mod
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding
  implicit none

  private
  public :: deserialize_int

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
