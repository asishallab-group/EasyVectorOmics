module real_deserialize_mod
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding
  use, intrinsic :: iso_fortran_env, only: real64, int32
  implicit none

  private
  public :: deserialize_real, deserialize_real_flat

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

subroutine deserialize_real_flat_r(flat_arr, dims_out, ndim_out, filename_ascii, fn_len)
  use iso_fortran_env
  use real_deserialize_mod
  implicit none

  ! Rückgabe an R
  real(real64), intent(out) :: flat_arr(*)
  integer(int32), intent(out) :: dims_out(*)
  integer, intent(out) :: ndim_out

  ! Dateiname
  integer(int32), intent(in) :: filename_ascii(fn_len)
  integer, intent(in) :: fn_len

  ! Intern
  character(len=:), allocatable :: filename
  integer :: i, k
  real(real64), pointer :: flat(:)
  integer(int32), allocatable, target :: dims(:)

  allocate(character(len=fn_len) :: filename)
  do i = 1, fn_len
    filename(i:i) = char(filename_ascii(i))
  end do

  ! Lese Daten
  call deserialize_real_flat(flat, dims, filename)

  ndim_out = size(dims)
  do i = 1, ndim_out
    dims_out(i) = dims(i)
  end do

  do i = 1, product(dims)
    flat_arr(i) = flat(i)
  end do

  !Flat is no longer needed but created manually so it needs to be deallocated
  if (associated(flat)) deallocate(flat)
end subroutine

subroutine deserialize_real_C(arr, arr_size, filename_ascii, fn_len) bind(C, name="deserialize_real_C")
  use iso_c_binding
  use real_deserialize_mod, only: deserialize_real_flat
  use iso_fortran_env
  implicit none

  real(c_double), intent(inout) :: arr(arr_size)
  integer(c_int), value :: arr_size
  integer(c_int), intent(in) :: filename_ascii(fn_len)
  integer(c_int), value :: fn_len

  character(len=:), allocatable :: filename
  integer :: i

  ! arr_f ist das Ergebnis von deserialize_int_flat – muss kopiert werden
  real(real64), pointer :: arr_f(:)
  integer(int32), allocatable :: dims(:)

  ! ASCII → String
  allocate(character(len=fn_len) :: filename)
  do i = 1, fn_len
    filename(i:i) = char(filename_ascii(i))
  end do

  ! Daten einlesen – arr_f wird intern allokiert
  call deserialize_real_flat(arr_f, dims, filename)

  ! Sicherheitschecks
  if (.not. associated(arr_f)) then
      print *, "Fehler: arr_f nicht allokiert"
      stop 1
  end if

  if (size(arr_f) /= arr_size) then
      print *, "Fehler: Größe passt nicht: ", size(arr_f), arr_size
      stop 2
  end if

  ! Jetzt in den Python-Puffer kopieren
  arr(:) = arr_f(:)
end subroutine