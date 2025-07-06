module array_utils
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use iso_c_binding
    implicit none

    public :: get_type_code, get_cell, get_row, get_col
    interface get_cell
        module procedure get_cell_int, get_cell_real, get_cell_char
    end interface get_cell

    interface get_row
        module procedure get_row_int, get_row_real, get_row_char
    end interface get_row

    interface get_col
        module procedure get_col_int, get_col_real, get_col_char
    end interface get_col
    integer(int32), parameter :: ARRAY_FILE_MAGIC = int(z'46413230', int32) ! 'FA20' in hex

   contains
  !> Extract a row from a 2D integer array.
  !! @param arr Input 2D integer array.
  !! @param idx Row index to extract.
  !! @param row Output row as a 1D integer array.
  subroutine get_row_int(arr, idx, row)
    integer, intent(in) :: arr(:,:), idx
    integer, intent(out) :: row(size(arr,2))
    row = arr(idx,:)
  end subroutine get_row_int

  subroutine get_row_real(arr, idx, row)
    real, intent(in) :: arr(:,:)
    integer, intent(in) :: idx
    real, intent(out) :: row(size(arr,2))
    row = arr(idx,:)
  end subroutine get_row_real

  subroutine get_row_char(arr, idx, row, ncol, strlen)
    integer, intent(in) :: idx, ncol, strlen
    character(len=*), intent(in) :: arr(:, :)
    character(len=*), intent(out) :: row(:)
    integer :: j
    do j = 1, ncol
        row(j) = arr(idx, j)(1:strlen)
    end do
  end subroutine get_row_char

  !> Extract a column from a 2D integer array.
  !! @param arr Input 2D integer array.
  !! @param idx Column index to extract.
  !! @param col Output column as a 1D integer array.
  subroutine get_col_int(arr, idx, col)
    integer(int32), intent(in) :: arr(:,:), idx
    integer(int32), intent(out) :: col(size(arr,1))
    col = arr(:,idx)
  end subroutine get_col_int

  subroutine get_col_real(arr, idx, col)
    real(real64), intent(in) :: arr(:,:), idx
    real(real64), intent(out) :: col(size(arr,1))
    col = arr(:,idx)
  end subroutine get_col_real

  subroutine get_col_char(arr, idx, col, nrow, strlen)
    integer(int32), intent(in) :: idx, nrow, strlen
    character(len=*), intent(in) :: arr(:, :)
    character(len=*), intent(out) :: col(:)
    integer :: i
    do i = 1, nrow
        col(i) = arr(i, idx)(1:strlen)
    end do
  end subroutine get_col_char

  !> Extract a single cell from a 2D integer array.
  !! @param arr Input 2D integer array.
  !! @param row_idx Row index.
  !! @param col_idx Column index.
  !! @param cell Output cell value.
  subroutine get_cell_int(arr, row_idx, col_idx, cell)
    integer, intent(in) :: arr(:,:)
    integer, intent(in) :: row_idx, col_idx
    integer, intent(out) :: cell
    cell = arr(row_idx, col_idx)
  end subroutine get_cell_int

  subroutine get_cell_real(arr, row_idx, col_idx, cell)
    real, intent(in) :: arr(:,:)
    integer, intent(in) :: row_idx, col_idx
    integer, intent(out) :: cell
    cell = arr(row_idx, col_idx)
  end subroutine get_cell_real

  subroutine get_cell_char(arr, row_idx, col_idx, cell, strlen)
    character(len=*), intent(in) :: arr(:,:)
    integer, intent(in) :: row_idx, col_idx, strlen
    character(len=*), intent(out) :: cell
    cell = arr(row_idx, col_idx)(1:strlen)
  end subroutine get_cell_char

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

subroutine get_row_int_r(arr, n_rows, n_cols, idx, row)
    use array_utils
    implicit none
    integer(int32), intent(in) :: n_rows, n_cols
    integer(int32), intent(in) :: arr(n_rows, n_cols)
    integer(int32), intent(in) :: idx
    integer(int32), intent(out) :: row(n_cols)
    call get_row_int(arr, idx, row)
end subroutine get_row_int_r

subroutine get_row_real_r(arr, n_rows, n_cols, idx, row)
    use array_utils
    implicit none
    integer(int32), intent(in) :: n_rows, n_cols
    real(real64), intent(in) :: arr(n_rows, n_cols)
    integer(int32), intent(in) :: idx
    real(real64), intent(out) :: row(n_cols)
    call get_row_real(arr, idx, row)
end subroutine get_row_real_r

subroutine get_row_char_r(arr, nrow, ncol, idx, row, strlen)
    !String needs fixed length strlen
    use array_utils
    implicit none
    integer(int32), intent(in) :: nrow, ncol, idx, strlen
    character(len=*), intent(in) :: arr(nrow, ncol)
    character(len=*), intent(out) :: row(ncol)
    call get_row_char(arr, idx, row, ncol, strlen)
end subroutine get_row_char_r

subroutine get_col_int_r(arr, n_rows, n_cols, idx, col)
    use array_utils
    implicit none
    integer(int32), intent(in) :: n_rows, n_cols
    integer(int32), intent(in) :: arr(n_rows, n_cols)
    integer(int32), intent(in) :: idx
    integer(int32), intent(out) :: col(n_rows)
    call get_row_int(arr, idx, col)
end subroutine get_col_int_r

subroutine get_col_real_r(arr, n_rows, n_cols, idx, col)
    use array_utils
    implicit none
    integer(int32), intent(in) :: n_rows, n_cols
    real(real64), intent(in) :: arr(n_rows, n_cols)
    integer(int32), intent(in) :: idx
    real(real64), intent(out) :: col(n_cols)
    call get_row_real(arr, idx, col)
end subroutine get_col_real_r

subroutine get_col_char_r(arr, nrow, ncol, idx, col, strlen)
  !String needs fixed length strlen
  use array_utils
  implicit none
  integer(int32), intent(in) :: nrow, ncol, idx, strlen
  character(len=*), intent(in) :: arr(nrow, ncol)
  character(len=*), intent(out) :: col(ncol)
  call get_row_char(arr, idx, col, ncol, strlen)
end subroutine get_col_char_r

subroutine get_cell_int_r(arr, n_rows, n_cols, row_idx, col_idx, cell)
    use array_utils
    implicit none
    integer(int32), intent(in) :: n_rows, n_cols
    integer(int32), intent(in) :: arr(n_rows, n_cols)
    integer(int32), intent(in) :: row_idx, col_idx
    integer(int32), intent(out) :: cell
    call get_cell_int(arr, row_idx, col_idx, cell)
end subroutine get_cell_int_r

subroutine get_cell_real_r(arr, n_rows, n_cols, row_idx, col_idx, cell)
    use array_utils
    implicit none
    integer(int32), intent(in) :: n_rows, n_cols
    real(real64), intent(in) :: arr(n_rows, n_cols)
    integer(int32), intent(in) :: row_idx, col_idx
    real(real64), intent(out) :: cell
    call get_cell_real(arr, row_idx, col_idx, cell)
end subroutine get_cell_real_r

subroutine get_cell_char_r(arr, nrow, ncol, row_idx, col_idx, cell, strlen)
    !String needs fixed length strlen
    use array_utils
    implicit none
    integer(int32), intent(in) :: nrow, ncol, row_idx, col_idx, strlen
    character(len=*), intent(in) :: arr(nrow, ncol)
    character(len=*), intent(out) :: cell
    call get_cell_char(arr, row_idx, col_idx, cell, strlen)
end subroutine get_cell_char_r

