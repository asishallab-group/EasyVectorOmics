module real_deserialize_mod
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding
  use, intrinsic :: iso_fortran_env, only: real64, int32
  implicit none

  private
  public :: deserialize_real

  integer(int32), parameter :: ARRAY_FILE_MAGIC = int(z'46413230', int32) ! 'FA20' in hex

  interface deserialize_real
    module procedure deserialize_real_1d
    module procedure deserialize_real_2d
    module procedure deserialize_real_3d
    module procedure deserialize_real_4d
    module procedure deserialize_real_5d
  end interface

contains

  ! Hilfsroutine: Flaches Array + Dimensionen lesen
  subroutine deserialize_real_flat(flat, dims, filename)
    use iso_c_binding
    real(real64), pointer, intent(out) :: flat(:)
    integer(int32), allocatable, intent(out), target :: dims(:)
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

  ! 1D
  subroutine deserialize_real_1d(arr, filename)
    use iso_c_binding
    real(real64), pointer, intent(out) :: arr(:)
    character(len=*), intent(in) :: filename

    real(real64), pointer :: flat(:)
    integer(int32), allocatable :: dims(:)

    call deserialize_real_flat(flat, dims, filename)
    if (size(dims) /= 1) error stop "Expected 1D array"
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1)])
  end subroutine deserialize_real_1d

  ! 2D
  subroutine deserialize_real_2d(arr, filename)
    use iso_c_binding
    real(real64), pointer, intent(out) :: arr(:,:)
    character(len=*), intent(in) :: filename

    real(real64), pointer :: flat(:)
    integer(int32), allocatable :: dims(:)

    call deserialize_real_flat(flat, dims, filename)
    if (size(dims) /= 2) error stop "Expected 2D array"
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1), dims(2)])
  end subroutine deserialize_real_2d

  ! 3D
  subroutine deserialize_real_3d(arr, filename)
    use iso_c_binding
    real(real64), pointer, intent(out) :: arr(:,:,:)
    character(len=*), intent(in) :: filename

    real(real64), pointer :: flat(:)
    integer(int32), allocatable :: dims(:)

    call deserialize_real_flat(flat, dims, filename)
    if (size(dims) /= 3) error stop "Expected 3D array"
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1), dims(2), dims(3)])
  end subroutine deserialize_real_3d

  ! 4D
  subroutine deserialize_real_4d(arr, filename)
    use iso_c_binding
    real(real64), pointer, intent(out) :: arr(:,:,:,:)
    character(len=*), intent(in) :: filename

    real(real64), pointer :: flat(:)
    integer(int32), allocatable :: dims(:)

    call deserialize_real_flat(flat, dims, filename)
    if (size(dims) /= 4) error stop "Expected 4D array"
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1), dims(2), dims(3), dims(4)])
  end subroutine deserialize_real_4d

  ! 5D
  subroutine deserialize_real_5d(arr, filename)
    use iso_c_binding
    real(real64), pointer, intent(out) :: arr(:,:,:,:,:)
    character(len=*), intent(in) :: filename

    real(real64), pointer :: flat(:)
    integer(int32), allocatable :: dims(:)

    call deserialize_real_flat(flat, dims, filename)
    if (size(dims) /= 5) error stop "Expected 5D array"
    call c_f_pointer(c_loc(flat(1)), arr, shape=[dims(1), dims(2), dims(3), dims(4), dims(5)])
  end subroutine deserialize_real_5d

end module real_deserialize_mod
