module array_utils
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use iso_c_binding
    implicit none

    public :: get_row, get_col, get_cell, get_type_code
    integer(int32), parameter :: ARRAY_FILE_MAGIC = int(z'46413230', int32) ! 'FA20' in hex

   contains
  !> Extract a row from a 2D integer array.
  !! @param arr Input 2D integer array.
  !! @param idx Row index to extract.
  !! @param row Output row as a 1D integer array.
  subroutine get_row(arr, idx, row)
    integer, intent(in) :: arr(:,:), idx
    integer, intent(out) :: row(size(arr,2))
    row = arr(idx,:)
  end subroutine get_row

  !> Extract a column from a 2D integer array.
  !! @param arr Input 2D integer array.
  !! @param idx Column index to extract.
  !! @param col Output column as a 1D integer array.
  subroutine get_col(arr, idx, col)
    integer, intent(in) :: arr(:,:), idx
    integer, intent(out) :: col(size(arr,1))
    col = arr(:,idx)
  end subroutine get_col

  !> Extract a single cell from a 2D integer array.
  !! @param arr Input 2D integer array.
  !! @param row_idx Row index.
  !! @param col_idx Column index.
  !! @param cell Output cell value.
  subroutine get_cell(arr, row_idx, col_idx, cell)
    integer, intent(in) :: arr(:,:)
    integer, intent(in) :: row_idx, col_idx
    integer, intent(out) :: cell
    cell = arr(row_idx, col_idx)
  end subroutine get_cell

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