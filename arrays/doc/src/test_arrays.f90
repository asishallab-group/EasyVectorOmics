!>@brief Test program for deserializing integer, real, and character arrays
!> @note This program tests the serialization/deserialization routines for integer, real, and character arrays
program test_arrays
  use array_utils
  use int_deserialize_mod
  use real_deserialize_mod
  use char_deserialize_mod
  use serialize_mod
  use serialize_char
  use serialize_int
  use serialize_real
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding
  implicit none

  integer(int32), allocatable :: row(:), col(:)
  integer(int32), allocatable :: iarr(:,:), iarr1d(:), iarr3d(:,:,:), iarr4d(:,:,:,:), iarr5d(:,:,:,:,:)
  integer(int32), pointer :: iarr2(:,:), iarr1d2(:), iarr3d2(:,:,:), iarr4d2(:,:,:,:), iarr5d2(:,:,:,:,:)
  integer(int32), pointer :: iarr_flat(:)
  integer(int32), allocatable :: idims(:)

  real(real64), allocatable :: rarr(:,:), rarr1d(:), rarr3d(:,:,:), rarr4d(:,:,:,:), rarr5d(:,:,:,:,:)
  real(real64), pointer :: rarr2(:,:), rarr1d2(:), rarr3d2(:,:,:), rarr4d2(:,:,:,:), rarr5d2(:,:,:,:,:)
  real(real64), pointer :: rarr_flat(:)
  integer(int32), allocatable :: rdims(:)

  character(len=:), pointer :: carr2(:,:), carr1d2(:), carr3d2(:,:,:), carr4d2(:,:,:,:), carr5d2(:,:,:,:,:)
  character(len=:), allocatable :: carr(:,:), carr1d(:), carr3d(:,:,:), carr4d(:,:,:,:), carr5d(:,:,:,:,:)
  character(len=:), pointer :: carr_flat(:)
  integer(int32), allocatable :: cdims(:)
  integer :: clen, cell
  character(len=5), allocatable :: crow(:), ccol(:)
  character(len=5), allocatable :: cnames(:)
  integer(int32) :: val, i, j
  character(len=100) :: fname
  logical :: ok

  character(len=:), pointer :: protein_data(:,:,:,:,:) => null()
  character(len=:), pointer :: protein_data_loaded(:,:,:,:,:) => null()
  character(len=*), parameter :: test_file = "test_proteins.bin"
  integer :: k, l, m
  logical :: test_passed = .true.

  integer(int32), ALLOCATABLE :: int_arr_6d(:,:,:,:,:,:)
  integer(int32), pointer :: int_arr_6d2(:,:,:,:,:,:)

  print *, "==== Test: deserialize_int 1D/2D/3D/4D/5D ===="
  ! 1D
  allocate(iarr1d(5)); iarr1d = [10,20,30,40,50]
  fname = "test_iarr1d.bin"
  call serialize_int_1d(iarr1d, fname)
  call deserialize_int_1d(iarr1d2, fname)
  if (all(iarr1d == iarr1d2)) then
    print *, "PASS: deserialize_int_1d"
  else
    print *, "FAIL: deserialize_int_1d"
  end if

  ! 2D
  allocate(iarr(2,3)); iarr = reshape([1,2,3,4,5,6], [2,3])
  fname = "test_iarr2d.bin"
  call serialize_int_2d(iarr, fname)
  call deserialize_int_2d(iarr2, fname)
  if (all(iarr == iarr2)) then
    print *, "PASS: deserialize_int_2d"
  else
    print *, "FAIL: deserialize_int_2d"
  end if

  ! 3D
  allocate(iarr3d(2,2,2)); iarr3d = reshape([1,2,3,4,5,6,7,8], [2,2,2])
  fname = "test_iarr3d.bin"
  call serialize_int_3d(iarr3d, fname)
  call deserialize_int_3d(iarr3d2, fname)
  if (all(iarr3d == iarr3d2)) then
    print *, "PASS: deserialize_int_3d"
  else
    print *, "FAIL: deserialize_int_3d"
  end if

  ! 4D
  allocate(iarr4d(2,2,1,2)); iarr4d = reshape([(i, i=1,8)], [2,2,1,2])
  fname = "test_iarr4d.bin"
  call serialize_int_4d(iarr4d, fname)
  call deserialize_int_4d(iarr4d2, fname)
  if (all(iarr4d == iarr4d2)) then
    print *, "PASS: deserialize_int_4d"
  else
    print *, "FAIL: deserialize_int_4d"
  end if

  ! 5D
  allocate(iarr5d(2,1,2,1,2)); iarr5d = reshape([(i, i=1,8)], [2,1,2,1,2])
  fname = "test_iarr5d.bin"
  call serialize_int_5d(iarr5d, fname)
  call deserialize_int_5d(iarr5d2, fname)
  if (all(iarr5d == iarr5d2)) then
    print *, "PASS: deserialize_int_5d"
  else
    print *, "FAIL: deserialize_int_5d"
  end if

  print *, "==== Test: deserialize_real 1D/2D/3D/4D/5D ===="
  ! 1D
  allocate(rarr1d(4)); rarr1d = [1.1_real64, 2.2_real64, 3.3_real64, 4.4_real64]
  fname = "test_rarr1d.bin"
  call serialize_real_1d(rarr1d, fname)
  call deserialize_real_1d(rarr1d2, fname)
  if (all(abs(rarr1d - rarr1d2) < 1e-12_real64)) then
    print *, "PASS: deserialize_real_1d"
  else
    print *, "FAIL: deserialize_real_1d"
  end if

  ! 2D
  allocate(rarr(2,2)); rarr = reshape([1.5_real64, 2.5_real64, 3.5_real64, 4.5_real64], [2,2])
  fname = "test_rarr2d.bin"
  call serialize_real_2d(rarr, fname)
  call deserialize_real_2d(rarr2, fname)
  if (all(abs(rarr - rarr2) < 1e-12_real64)) then
    print *, "PASS: deserialize_real_2d"
  else
    print *, "FAIL: deserialize_real_2d"
  end if

  ! 3D
  allocate(rarr3d(2,2,2)); rarr3d = reshape([1.0_real64,2.0_real64,3.0_real64,4.0_real64, &
    5.0_real64,6.0_real64,7.0_real64,8.0_real64], [2,2,2])
  fname = "test_rarr3d.bin"
  call serialize_real_3d(rarr3d, fname)
  call deserialize_real_3d(rarr3d2, fname)
  if (all(abs(rarr3d - rarr3d2) < 1e-12_real64)) then
    print *, "PASS: deserialize_real_3d"
  else
    print *, "FAIL: deserialize_real_3d"
  end if

  ! 4D
  allocate(rarr4d(2,2,1,2)); rarr4d = reshape([(real(i,real64), i=1,8)], [2,2,1,2])
  fname = "test_rarr4d.bin"
  call serialize_real_4d(rarr4d, fname)
  call deserialize_real_4d(rarr4d2, fname)
  if (all(abs(rarr4d - rarr4d2) < 1e-12_real64)) then
    print *, "PASS: deserialize_real_4d"
  else
    print *, "FAIL: deserialize_real_4d"
  end if

  ! 5D
  allocate(rarr5d(2,1,2,1,2)); rarr5d = reshape([(real(i,real64), i=1,8)], [2,1,2,1,2])
  fname = "test_rarr5d.bin"
  call serialize_real_5d(rarr5d, fname)
  call deserialize_real_5d(rarr5d2, fname)
  if (all(abs(rarr5d - rarr5d2) < 1e-12_real64)) then
    print *, "PASS: deserialize_real_5d"
  else
    print *, "FAIL: deserialize_real_5d"
  end if

  print *, "==== Test: deserialize_char 1D/2D/3D/4D/5D ===="
  ! 1D
  clen = 5
  allocate(character(len=clen):: carr1d(3))
  carr1d = ['foo  ','bar  ','baz  ']
  fname = "test_carr1d.bin"
  call serialize_char_1d(carr1d, fname)
  call deserialize_char_1d(carr1d2, fname)
  if (all(carr1d == carr1d2)) then
    print *, "PASS: deserialize_char_1d"
  else
    print *, "FAIL: deserialize_char_1d"
  end if

  ! 2D
  clen = 5
  allocate(character(len=clen):: carr(2,2))
  carr = reshape(['foo  ','bar  ','baz  ','qux  '], [2,2])
  fname = "test_carr2d.bin"
  call serialize_char_2d(carr, fname)
  call deserialize_char_2d(carr2, fname)
  if (all(carr == carr2)) then
    print *, "PASS: deserialize_char_2d"
  else
    print *, "FAIL: deserialize_char_2d"
  end if

  ! 3D
  clen = 5
  allocate(character(len=clen):: carr3d(2,2,1))
  carr3d = reshape(['foo  ','bar  ','baz  ','qux  '], [2,2,1])
  fname = "test_carr3d.bin"
  call serialize_char_3d(carr3d, fname)
  call deserialize_char_3d(carr3d2, fname)
  if (all(carr3d == carr3d2)) then
    print *, "PASS: deserialize_char_3d"
  else
    print *, "FAIL: deserialize_char_3d"
  end if

  ! 4D
  clen = 5
  allocate(character(len=clen):: carr4d(2,1,1,2))
  carr4d = reshape(['foo  ','bar  ','baz  ','qux  '], [2,1,1,2])
  fname = "test_carr4d.bin"
  call serialize_char_4d(carr4d, fname)
  call deserialize_char_4d(carr4d2, fname)
  if (all(carr4d == carr4d2)) then
    print *, "PASS: deserialize_char_4d"
  else
    print *, "FAIL: deserialize_char_4d"
  end if

  ! 5D
  clen = 5
  allocate(character(len=clen):: carr5d(2,1,2,1,2))
  carr5d = reshape(['foo  ','bar  ','baz  ','qux  ','aaa  ','bbb  ','ccc  ','ddd  '], [2,1,2,1,2])
  fname = "test_carr5d.bin"
  call serialize_char_5d(carr5d, fname)
  call deserialize_char_5d(carr5d2, fname)
  if (all(carr5d == carr5d2)) then
    print *, "PASS: deserialize_char_5d"
  else
    print *, "FAIL: deserialize_char_5d"
  end if

  ! Arrays für Accessor-Tests wiederherstellen
  if (allocated(iarr)) deallocate(iarr)
  if (allocated(row)) deallocate(row)
  if (allocated(col)) deallocate(col)
  allocate(iarr(2,3))
  iarr = reshape([1,2,3,4,5,6], [2,3])


  print *, "==== Edge Case: 1x1 array serialization ===="
  if (allocated(iarr)) deallocate(iarr)
  allocate(iarr(1,1))
  iarr = 42
  fname = "test_iarr_1x1.bin"
  call serialize_int_2d(iarr, fname)
  call deserialize_int_2d(iarr2, fname)
  print *, "iarr2=", iarr2
  if (iarr2(1,1) == 42) then
    print *, "PASS: 1x1 array serialization"
  else
    print *, "FAIL: 1x1 array serialization"
  end if

  print *, "==== Edge Case: Empty array (0 rows) ===="
  if (allocated(iarr)) deallocate(iarr)
  allocate(iarr(0,3))
  fname = "test_iarr_empty.bin"
  call serialize_int_2d(iarr, fname)
  call deserialize_int_2d(iarr2, fname)
  if (size(iarr2,1) == 0 .and. size(iarr2,2) == 3) then
    print *, "PASS: Empty array serialization"
  else
    print *, "FAIL: Empty array serialization"
  end if

  ! Test variable string length arrays

  print *, "Starting variable length string test for 5D array..."

  ! 1. prep test data (2x1x2x1x2 Array)
  allocate(character(len=12) :: protein_data(2,1,2,1,2))  ! Maximallänge 12 Zeichen

  protein_data(1,1,1,1,1) = "P12345"         ! 6 Zeichen
  protein_data(2,1,1,1,1) = "Q9Y6K8-2"       ! 8 Zeichen
  protein_data(1,1,2,1,1) = "A0A024R3R5"     ! 10 Zeichen
  protein_data(2,1,2,1,1) = "O14773"         ! 6 Zeichen 
  protein_data(1,1,1,1,2) = "P0DTD1-PROT"    ! 11 Zeichen
  protein_data(2,1,1,1,2) = ""               ! Leerer String
  protein_data(1,1,2,1,2) = "X6R8Y4"         ! 6 Zeichen
  protein_data(2,1,2,1,2) = "A0A0B4J2F5-1"   ! 12 Zeichen

  ! 2. serialize and deserialize
  print *, "Serializing data..."
  call serialize_char_5d(protein_data, test_file)

  print *, "Deserializing data..."
  call deserialize_char_5d(protein_data_loaded, test_file)

  ! 3. verify results
  print *, "Verifying data..."
  do m = 1, 2
    do l = 1, 1
      do k = 1, 2
        do j = 1, 1
          do i = 1, 2
            if (protein_data(i,j,k,l,m) /= protein_data_loaded(i,j,k,l,m)) then
              print *, "Mismatch at position (", i, ",", j, ",", k, ",", l, ",", m, ")"
              print *, "  Original: '", protein_data(i,j,k,l,m), "' (Length:", len_trim(protein_data(i,j,k,l,m)), ")"
              print *, "  Loaded:   '", protein_data_loaded(i,j,k,l,m), "' (Length:", len_trim(protein_data_loaded(i,j,k,l,m)), ")"
              test_passed = .false.
            end if
          end do
        end do
      end do
    end do
  end do

  ! 4. print result
  if (test_passed) then
    print *, "PASS: 5D variable length string test successful"
    print *, "All", size(protein_data), "elements match perfectly!"
  else
    print *, "FAIL: 5D variable length string test failed"
  end if

  ! 5. free memory
  if (associated(protein_data)) deallocate(protein_data)
  if (associated(protein_data_loaded)) deallocate(protein_data_loaded)
  print *, "Test completed. Cleaned up memory."


  print *, "==== All tests done. ===="

  ! Test what happens when interface gets 6D array
  ! allocate(int_arr_6d(2,2,1,2,1,1)); int_arr_6d = reshape([(i, i=1,8)], [2,2,1,2,1,1])
  ! call serialize(int_arr_6d, "error_test.bin")

end program test_arrays
