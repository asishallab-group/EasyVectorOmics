module char_deserialize_mod
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding
  implicit none

  private
  public :: deserialize_char, deserialize_char_flat

  integer(int32), parameter :: ARRAY_FILE_MAGIC = int(z'46413230', int32) ! 'FA20' in hex

  interface deserialize_char
    module procedure deserialize_char_1d
    module procedure deserialize_char_2d
    module procedure deserialize_char_3d
    module procedure deserialize_char_4d
    module procedure deserialize_char_5d
  end interface

contains

  ! Hilfsroutine: Flaches Array + Dimensionen lesen (angepasst für variable Stringlängen)
  subroutine deserialize_char_flat(flat, dims, clen, filename)
    use iso_c_binding
    character(len=:), pointer, intent(out) :: flat(:)
    integer(int32), allocatable, intent(out) :: dims(:)
    integer, intent(out) :: clen
    character(len=*), intent(in) :: filename

    integer :: unit, magic, type_code, d, i, str_len
    character(len=:), allocatable :: temp_str

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

  ! 1D
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

  ! 2D
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
    idx = 1
    do j = 1, dims(2)
      do i = 1, dims(1)
        arr(i,j) = flat(idx)
        idx = idx + 1
      end do
    end do
    deallocate(flat)
  end subroutine

  ! 3D
  subroutine deserialize_char_3d(arr, filename)
    use iso_c_binding
    character(len=:), pointer, intent(out) :: arr(:,:,:)
    character(len=*), intent(in) :: filename

    character(len=:), pointer :: flat(:)
    integer(int32), allocatable :: dims(:)
    integer :: clen, i, j, k, idx

    call deserialize_char_flat(flat, dims, clen, filename)
    if (size(dims) /= 3) error stop "Expected 3D array"
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
    deallocate(flat)
  end subroutine

  ! 4D
  subroutine deserialize_char_4d(arr, filename)
    use iso_c_binding
    character(len=:), pointer, intent(out) :: arr(:,:,:,:)
    character(len=*), intent(in) :: filename

    character(len=:), pointer :: flat(:)
    integer(int32), allocatable :: dims(:)
    integer :: clen, i, j, k, l, idx

    call deserialize_char_flat(flat, dims, clen, filename)
    if (size(dims) /= 4) error stop "Expected 4D array"
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
    deallocate(flat)
  end subroutine

    ! 5D Deserialisierung (konsistent mit 1D-4D Routinen)
    subroutine deserialize_char_5d(arr, filename)
        character(len=:), pointer, intent(out) :: arr(:,:,:,:,:)
        character(len=*), intent(in) :: filename

        character(len=:), pointer :: flat(:)
        integer(int32), allocatable :: dims(:)
        integer :: max_len, i, j, k, l, m, idx

        if (associated(arr)) nullify(arr)

        call deserialize_char_flat(flat, dims, max_len, filename)
        if (size(dims) /= 5) error stop "Expected 5D array"

        allocate(character(len=max_len) :: arr(dims(1), dims(2), dims(3), dims(4), dims(5)))

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

        deallocate(flat)  ! Wichtig: Freigabe des temporären Arrays
    end subroutine

end module char_deserialize_mod

subroutine deserialize_char_flat_r(ascii_arr, dims_out, ndim_out, clen_out, filename_ascii, fn_len)
  use iso_fortran_env
  use char_deserialize_mod
  implicit none

  integer(int32), intent(out) :: ascii_arr(*)
  integer(int32), intent(out) :: dims_out(*)
  integer, intent(out) :: ndim_out, clen_out
  integer(int32), intent(in) :: filename_ascii(fn_len)
  integer(int32), intent(in) :: fn_len

  character(len=:), allocatable :: filename
  character(len=:), pointer :: flat(:)
  integer(int32), allocatable :: dims(:)
  integer :: i, j, clen, total

  allocate(character(len=fn_len) :: filename)
  do i = 1, fn_len
    filename(i:i) = char(filename_ascii(i))
  end do

  call deserialize_char_flat(flat, dims, clen, filename)
  total = product(dims)
  clen_out = clen
  ndim_out = size(dims)

  do i = 1, ndim_out
    dims_out(i) = dims(i)
  end do

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

