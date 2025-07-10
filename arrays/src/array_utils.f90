module array_utils
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use iso_c_binding
    implicit none

    public :: get_type_code

    integer(int32), parameter :: ARRAY_FILE_MAGIC = int(z'46413230', int32) ! 'FA20' in hex

   contains
  !> @brief returns type code of the array file
  !> @param filename name of the file to read
  !> @return type code of the array file
  function get_type_code(filename) result(type_code)
    use iso_c_binding
    implicit none

    character(len=*), intent(in) :: filename
    integer :: type_code
    integer :: unit, magic

    open(newunit=unit, file=filename, form='unformatted', access='stream', status='old')
    read(unit) magic
    if (magic /= ARRAY_FILE_MAGIC) error stop "Invalid file format"
    read(unit) type_code
    close(unit)
  end function get_type_code
end module array_utils

subroutine get_type_code_r(filename, type_code)
    use array_utils
    implicit none

    character(len=*), intent(in) :: filename
    integer, INTENT(OUT) :: type_code

    type_code = get_type_code(filename)
end subroutine get_type_code_r