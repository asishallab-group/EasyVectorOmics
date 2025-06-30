program test_arrays
  use array_ops
  use reshape_utils
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding
  implicit none

  integer(int32), allocatable :: row(:), col(:)
  integer(int32), allocatable :: iarr(:,:), iarr1d(:), iarr3d(:,:,:), iarr4d(:,:,:,:), iarr5d(:,:,:,:,:)
  integer(int32), pointer :: iarr2(:,:), iarr1d2(:), iarr3d2(:,:,:), iarr4d2(:,:,:,:)
  integer(int32), pointer :: iarr_flat(:)
  integer(int32), allocatable :: idims(:)

  real(real64), allocatable :: rarr(:,:), rarr1d(:), rarr3d(:,:,:), rarr4d(:,:,:,:), rarr5d(:,:,:,:,:)
  real(real64), pointer :: rarr2(:,:), rarr1d2(:), rarr3d2(:,:,:), rarr4d2(:,:,:,:)

  real(real64), pointer :: rarr_flat(:)
  integer(int32), allocatable :: rdims(:)
  character(len=:), pointer :: carr2(:,:) ! Korrigiert: deferred length
  character(len=:), pointer :: carr1d2(:)
  character(len=:), pointer :: carr3d2(:,:,:)
  character(len=:), pointer :: carr4d2(:,:,:,:)
  character(len=:), pointer :: carr5d(:,:,:,:,:)
  character(len=:), allocatable :: carr(:,:), carr1d(:), carr3d(:,:,:), carr4d(:,:,:,:)
  character(len=:), pointer :: carr5d_ref(:,:,:,:,:)
  character(len=:), pointer :: carr_flat(:)
  integer(int32), allocatable :: cdims(:)
  integer :: clen
  character(len=5), allocatable :: crow(:), ccol(:)
  character(len=5), allocatable :: cnames(:)
  integer(int32) :: val, i, j
  character(len=100) :: fname
  logical :: ok

  ! Funktionsprototypen für Accessor-Tests
  interface
    function get_row(arr, idx) result(row)
      import :: int32
      integer(int32), intent(in) :: arr(:,:), idx
      integer(int32) :: row(size(arr,2))
    end function
    function get_col(arr, idx) result(col)
      import :: int32
      integer(int32), intent(in) :: arr(:,:), idx
      integer(int32) :: col(size(arr,1))
    end function
    function get_cell(arr, i, j) result(val)
      import :: int32
      integer(int32), intent(in) :: arr(:,:), i, j
      integer(int32) :: val
    end function
    function get_col_for_name(arr, cnames, name) result(col)
      import :: int32
      integer(int32), intent(in) :: arr(:,:)
      character(len=*), intent(in) :: cnames(:), name
      integer(int32) :: col(size(arr,1))
    end function
    function get_cell_for_column_name(arr, cnames, name, rowidx) result(val)
      import :: int32
      integer(int32), intent(in) :: arr(:,:)
      character(len=*), intent(in) :: cnames(:), name
      integer(int32), intent(in) :: rowidx
      integer(int32) :: val
    end function
  end interface

  print *, "==== Test: serialize/deserialize integer array (reshape) ===="
  allocate(iarr(2,3))
  iarr = reshape([1,2,3,4,5,6], [2,3])
  fname = "test_iarr.bin"
  call serialize(iarr, fname)
  call deserialize_int_flat(iarr_flat, idims, fname)
  allocate(iarr2(idims(1), idims(2)))
  call reshape_int_to_2D(iarr_flat, iarr2, idims)
  if (all(iarr == iarr2)) then
    print *, "PASS: Integer array serialization/reshape"
  else
    print *, "FAIL: Integer array serialization/reshape"
  end if

  print *, "==== Test: serialize/deserialize real array (reshape) ===="
  allocate(rarr(2,2))
  rarr = reshape([1.5_real64, 2.5_real64, 3.5_real64, 4.5_real64], [2,2])
  fname = "test_rarr.bin"
  call serialize(rarr, fname)
  call deserialize_real_flat(rarr_flat, rdims, fname)
  allocate(rarr2(rdims(1), rdims(2)))
  call reshape_real_to_2D(rarr_flat, rarr2, rdims)
  if (all(abs(rarr - rarr2) < 1e-12_real64)) then
    print *, "PASS: Real array serialization/reshape"
  else
    print *, "FAIL: Real array serialization/reshape"
  end if

  print *, "==== Test: serialize/deserialize character array (reshape) ===="
  allocate(character(len=clen):: carr(2,2))
  carr = reshape(['foo  ','bar  ','baz  ','qux  '], [2,2])
  fname = "test_carr.bin"
  call serialize(carr, fname)
  call deserialize_char_flat(carr_flat, cdims, clen, fname)
  allocate(character(len=clen) :: carr2(cdims(1), cdims(2)))
  call reshape_char_to_2D(carr_flat, carr2, cdims, clen)
  if (all(carr == carr2)) then
    print *, "PASS: Character array serialization/reshape"
  else
    print *, "FAIL: Character array serialization/reshape"
  end if

  print *, "==== Edge Case: 1x1 array serialization ===="
  if (allocated(iarr)) deallocate(iarr)
  allocate(iarr(1,1))
  iarr = 42
  fname = "test_iarr_1x1.bin"
  call serialize(iarr, fname)
  call deserialize_int_flat(iarr_flat, idims, fname)
  allocate(iarr2(idims(1), idims(2)))
  call reshape_int_to_2D(iarr_flat, iarr2, idims)
  if (iarr2(1,1) == 42) then
    print *, "PASS: 1x1 array serialization"
  else
    print *, "FAIL: 1x1 array serialization"
  end if

  print *, "==== Edge Case: Empty array (0 rows) ===="
  if (allocated(iarr)) deallocate(iarr)
  allocate(iarr(0,3))
  fname = "test_iarr_empty.bin"
  call serialize(iarr, fname)
  call deserialize_int_flat(iarr_flat, idims, fname)
  allocate(iarr2(0,3))
  call reshape_int_to_2D(iarr_flat, iarr2, idims)
  if (size(iarr2,1) == 0 .and. size(iarr2,2) == 3) then
    print *, "PASS: Empty array serialization"
  else
    print *, "FAIL: Empty array serialization"
  end if

  print *, "==== Test: serialize/deserialize 1D integer array ===="
  if (allocated(iarr1d)) deallocate(iarr1d)
  allocate(iarr1d(5))
  iarr1d = [10,20,30,40,50]
  fname = "test_iarr1d.bin"
  call serialize(iarr1d, fname)
  call deserialize_int_flat(iarr_flat, idims, fname)
  allocate(iarr1d2(idims(1)))
  call reshape_int_to_1D(iarr_flat, iarr1d2, idims)
  if (all(iarr1d == iarr1d2)) then
    print *, "PASS: 1D integer array serialization"
  else
    print *, "FAIL: 1D integer array serialization"
  end if

  print *, "==== Test: serialize/deserialize 3D integer array ===="
  if (allocated(iarr3d)) deallocate(iarr3d)
  allocate(iarr3d(2,2,2))
  iarr3d = reshape([1,2,3,4,5,6,7,8], [2,2,2])
  fname = "test_iarr3d.bin"
  call serialize(iarr3d, fname)
  call deserialize_int_flat(iarr_flat, idims, fname)
  allocate(iarr3d2(idims(1), idims(2), idims(3)))
  call reshape_int_to_3D(iarr_flat, iarr3d2, idims)
  if (all(iarr3d == iarr3d2)) then
    print *, "PASS: 3D integer array serialization"
  else
    print *, "FAIL: 3D integer array serialization"
  end if

  print *, "==== Test: serialize/deserialize 4D integer array ===="
  if (allocated(iarr4d)) deallocate(iarr4d)
  allocate(iarr4d(2,2,1,2))
  iarr4d = reshape([(i, i=1,8)], [2,2,1,2])
  fname = "test_iarr4d.bin"
  call serialize(iarr4d, fname)
  call deserialize_int_flat(iarr_flat, idims, fname)
  allocate(iarr4d2(idims(1), idims(2), idims(3), idims(4)))
  call reshape_int_to_4D(iarr_flat, iarr4d2, idims)
  if (all(iarr4d == iarr4d2)) then
    print *, "PASS: 4D integer array serialization"
  else
    print *, "FAIL: 4D integer array serialization"
  end if

  print *, "==== Test: serialize/deserialize 1D real array ===="
  if (allocated(rarr1d)) deallocate(rarr1d)
  allocate(rarr1d(4))
  rarr1d = [1.1_real64, 2.2_real64, 3.3_real64, 4.4_real64]
  fname = "test_rarr1d.bin"
  call serialize(rarr1d, fname)
  call deserialize_real_flat(rarr_flat, rdims, fname)
  allocate(rarr1d2(rdims(1)))
  call reshape_real_to_1D(rarr_flat, rarr1d2, rdims)
  if (all(abs(rarr1d - rarr1d2) < 1e-12_real64)) then
    print *, "PASS: 1D real array serialization"
  else
    print *, "FAIL: 1D real array serialization"
  end if

  print *, "==== Test: serialize/deserialize 3D real array ===="
  if (allocated(rarr3d)) deallocate(rarr3d)
  allocate(rarr3d(2,2,2))
  rarr3d = reshape([1.0_real64,2.0_real64,3.0_real64,4.0_real64,5.0_real64,6.0_real64,7.0_real64,8.0_real64], [2,2,2])
  fname = "test_rarr3d.bin"
  call serialize(rarr3d, fname)
  call deserialize_real_flat(rarr_flat, rdims, fname)
  allocate(rarr3d2(rdims(1), rdims(2), rdims(3)))
  call reshape_real_to_3D(rarr_flat, rarr3d2, rdims)
  if (all(abs(rarr3d - rarr3d2) < 1e-12_real64)) then
    print *, "PASS: 3D real array serialization"
  else
    print *, "FAIL: 3D real array serialization"
  end if

  print *, "==== Test: serialize/deserialize 4D real array ===="
  if (allocated(rarr4d)) deallocate(rarr4d)
  allocate(rarr4d(2,2,1,2))
  rarr4d = reshape([(real(i,real64), i=1,8)], [2,2,1,2])
  fname = "test_rarr4d.bin"
  call serialize(rarr4d, fname)
  call deserialize_real_flat(rarr_flat, rdims, fname)
  allocate(rarr4d2(rdims(1), rdims(2), rdims(3), rdims(4)))
  call reshape_real_to_4D(rarr_flat, rarr4d2, rdims)
  if (all(abs(rarr4d - rarr4d2) < 1e-12_real64)) then
    print *, "PASS: 4D real array serialization"
  else
    print *, "FAIL: 4D real array serialization"
  end if

  print *, "==== Test: serialize/deserialize 1D character array ===="
  if (allocated(carr1d)) deallocate(carr1d)
  allocate(character(len=clen):: carr1d(3))
  carr1d = ['foo  ','bar  ','baz  ']
  fname = "test_carr1d.bin"
  call serialize(carr1d, fname)
  call deserialize_char_flat(carr_flat, cdims, clen, fname)
  allocate(character(len=clen) :: carr1d2(cdims(1)))
  call reshape_char_to_1D(carr_flat, carr1d2, cdims, clen)
  if (all(carr1d == carr1d2)) then
    print *, "PASS: 1D character array serialization"
  else
    print *, "FAIL: 1D character array serialization"
  end if

  print *, "==== Test: serialize/deserialize 3D character array ===="
  if (allocated(carr3d)) deallocate(carr3d)
  allocate(character(len=clen):: carr3d(2,2,1))
  carr3d = reshape(['foo  ','bar  ','baz  ','qux  '], [2,2,1])
  fname = "test_carr3d.bin"
  call serialize(carr3d, fname)
  call deserialize_char_flat(carr_flat, cdims, clen, fname)
  allocate(character(len=clen) :: carr3d2(cdims(1), cdims(2), cdims(3)))
  call reshape_char_to_3D(carr_flat, carr3d2, cdims, clen)
  if (all(carr3d == carr3d2)) then
    print *, "PASS: 3D character array serialization"
  else
    print *, "FAIL: 3D character array serialization"
  end if

  print *, "==== Test: serialize/deserialize 4D character array ===="
  if (allocated(carr4d)) deallocate(carr4d)
  allocate(character(len=clen):: carr4d(2,1,1,2))
  carr4d = reshape(['foo  ','bar  ','baz  ','qux  '], [2,1,1,2])
  fname = "test_carr4d.bin"
  call serialize(carr4d, fname)
  call deserialize_char_flat(carr_flat, cdims, clen, fname)
  allocate(character(len=clen) :: carr4d2(cdims(1), cdims(2), cdims(3), cdims(4)))
  call reshape_char_to_4D(carr_flat, carr4d2, cdims, clen)
  if (all(carr4d == carr4d2)) then
    print *, "PASS: 4D character array serialization"
  else
    print *, "FAIL: 4D character array serialization"
  end if

  ! Arrays für Accessor-Tests wiederherstellen
  if (allocated(iarr)) deallocate(iarr)
  if (allocated(row)) deallocate(row)
  if (allocated(col)) deallocate(col)
  allocate(iarr(2,3))
  iarr = reshape([1,2,3,4,5,6], [2,3])

  print *, "==== Test: get_row ===="
  if (allocated(row)) deallocate(row)
  allocate(row(3))
  row = get_row(iarr, 2)
  print *, "row=", row
  if (all(row == [2_int32,4_int32,6_int32])) then
    print *, "PASS: get_row"
  else
    print *, "FAIL: get_row"
  end if

  print *, "==== Test: get_col ===="
  if (allocated(col)) deallocate(col)
  allocate(col(2))
  col = get_col(iarr, 1)
  print *, "col=", col
  if (all(col == [1_int32,2_int32])) then
    print *, "PASS: get_col"
  else
    print *, "FAIL: get_col"
  end if

  print *, "==== Test: get_cell ===="
  val = get_cell(iarr, 2, 3)
  if (val == 6) then
    print *, "PASS: get_cell"
  else
    print *, "FAIL: get_cell"
  end if

  print *, "==== Test: get_col_for_name ===="
  if (allocated(cnames)) deallocate(cnames)
  allocate(cnames(3))
  cnames = ['A','B','C']
  col = get_col_for_name(iarr, cnames, 'B')
  print *, "col_for_name=", col
  if (all(col == [3_int32,4_int32])) then
    print *, "PASS: get_col_for_name"
  else
    print *, "FAIL: get_col_for_name"
  end if

  print *, "==== Test: get_cell_for_column_name ===="
  val = get_cell_for_column_name(iarr, cnames, 'C', 2)
  print *, "cell_for_column_name=", val
  if (val == 6_int32) then
    print *, "PASS: get_cell_for_column_name"
  else
    print *, "FAIL: get_cell_for_column_name"
  end if

  print *, "==== Edge Case: 1x1 array serialization ===="
  if (allocated(iarr)) deallocate(iarr)
  allocate(iarr(1,1))
  iarr = 42
  fname = "test_iarr_1x1.bin"
  call serialize(iarr, fname)
  call deserialize(iarr2, fname)
  if (iarr2(1,1) == 42) then
    print *, "PASS: 1x1 array serialization"
  else
    print *, "FAIL: 1x1 array serialization"
  end if

  print *, "==== Edge Case: Empty array (0 rows) ===="
  if (allocated(iarr)) deallocate(iarr)
  allocate(iarr(0,3))
  fname = "test_iarr_empty.bin"
  call serialize(iarr, fname)
  call deserialize(iarr2, fname)
  if (size(iarr2,1) == 0 .and. size(iarr2,2) == 3) then
    print *, "PASS: Empty array serialization"
  else
    print *, "FAIL: Empty array serialization"
  end if

  print *, "==== All tests done. ===="

end program test_arrays
