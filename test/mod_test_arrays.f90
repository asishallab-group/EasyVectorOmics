module mod_test_arrays
    use asserts
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
    PUBLIC

    abstract interface
        subroutine test_interface()
        end subroutine test_interface
    end interface

type :: test_case
    character(len=64) :: name
    procedure(test_interface), pointer, nopass :: test_proc => null()
  end type test_case

contains

  !> @brief Get array of all available tests.
  function get_all_tests() result(all_tests)
    type(test_case) :: all_tests(18)  ! Number of tests
    
    all_tests(1) = test_case("test_integer_array_1d", test_integer_array_1d)
    all_tests(2) = test_case("test_integer_array_2d", test_integer_array_2d)
    all_tests(3) = test_case("test_integer_array_3d", test_integer_array_3d)
    all_tests(4) = test_case("test_integer_array_4d", test_integer_array_4d)
    all_tests(5) = test_case("test_integer_array_5d", test_integer_array_5d)
    all_tests(6) = test_case("test_real_array_1d", test_real_array_1d)
    all_tests(7) = test_case("test_real_array_2d", test_real_array_2d)
    all_tests(8) = test_case("test_real_array_3d", test_real_array_3d)
    all_tests(9) = test_case("test_real_array_4d", test_real_array_4d)
    all_tests(10) = test_case("test_real_array_5d", test_real_array_5d)
    all_tests(11) = test_case("test_char_array_1d", test_char_array_1d)
    all_tests(12) = test_case("test_char_array_2d", test_char_array_2d)
    all_tests(13) = test_case("test_char_array_3d", test_char_array_3d)
    all_tests(14) = test_case("test_char_array_4d", test_char_array_4d)
    all_tests(15) = test_case("test_char_array_5d", test_char_array_5d)
    all_tests(16) = test_case("test_integer_array_1x1", test_integer_array_1x1)
    all_tests(17) = test_case("test_integer_array_empty", test_integer_array_empty)
    all_tests(18) = test_case("test_char_array_protein", test_char_array_protein)
  end function get_all_tests

  !> @brief Run all tests in this module.
  subroutine run_all_tests_array()
    type(test_case) :: all_tests(18)
    integer :: i
    
    all_tests = get_all_tests()
    
    do i = 1, size(all_tests)
      call all_tests(i)%test_proc()
      print *, trim(all_tests(i)%name), " passed."
    end do
    print *, "All array tests passed successfully."
  end subroutine run_all_tests_array

  !> @brief Run specific tests by name.
  subroutine run_named_tests_array(test_names)
    character(len=*), intent(in) :: test_names(:)
    type(test_case) :: all_tests(18)
    integer :: i, j
    logical :: found
    
    all_tests = get_all_tests()
    
    do i = 1, size(test_names)
      found = .false.
      do j = 1, size(all_tests)
        if (trim(test_names(i)) == trim(all_tests(j)%name)) then
          call all_tests(j)%test_proc()
          print *, trim(test_names(i)), " passed."
          found = .true.
          exit
        end if
      end do
      if (.not. found) then
        print *, "Unknown test: ", trim(test_names(i))
      end if
    end do
  end subroutine run_named_tests_array

  !> @brief Test for 1D integer array
  subroutine test_integer_array_1d()
    integer(int32), ALLOCATABLE :: iarr1d(:)
    integer(int32), pointer :: iarr1d2(:)
    character(len=100) :: fname

    allocate(iarr1d(5)) 
    iarr1d = [10,20,30,40,50]
    fname = "test_iarr1d.bin"

    call serialize_int_1d(iarr1d, fname)
    call deserialize_int_1d(iarr1d2, fname)

    call assert_equal_array_int(iarr1d, iarr1d2, size(iarr1d), "Original array does not match deserialized array")
  end subroutine test_integer_array_1d

  !> @brief Test for 2D integer array
  subroutine test_integer_array_2d()
    INTEGER(int32), ALLOCATABLE :: iarr(:,:)
    INTEGER(int32), POINTER :: iarr2(:,:)
    CHARACTER(len=100) :: fname

    allocate(iarr(2,3))
    iarr = reshape([1,2,3,4,5,6], [2,3])
    fname = "test_iarr2d.bin"
    call serialize_int_2d(iarr, fname)
    call deserialize_int_2d(iarr2, fname)
    call assert_equal_array_int(iarr, iarr2, size(iarr), "Original array does not match deserialized array")
  end subroutine test_integer_array_2d

  !> @brief Test for 3D integer array
  subroutine test_integer_array_3d()
    integer(int32), allocatable :: iarr3d(:,:,:)
    integer(int32), pointer :: iarr3d2(:,:,:)
    character(len=100) :: fname

    allocate(iarr3d(2,2,2))
    iarr3d = reshape([1,2,3,4,5,6,7,8], [2,2,2])
    fname = "test_iarr3d.bin"

    call serialize_int_3d(iarr3d, fname)
    call deserialize_int_3d(iarr3d2, fname)

    call assert_equal_array_int(iarr3d, iarr3d2, size(iarr3d), "Original array does not match deserialized array")
  end subroutine test_integer_array_3d

  !> @brief Test for 4D integer array
  subroutine test_integer_array_4d()
    integer(int32), allocatable :: iarr4d(:,:,:,:)
    integer(int32), pointer :: iarr4d2(:,:,:,:)
    character(len=100) :: fname
    integer :: i

    allocate(iarr4d(2,2,1,2))
    iarr4d = reshape([(i, i=1,8)], [2,2,1,2])
    fname = "test_iarr4d.bin"

    call serialize_int_4d(iarr4d, fname)
    call deserialize_int_4d(iarr4d2, fname)

    call assert_equal_array_int(iarr4d, iarr4d2, size(iarr4d), "Original array does not match deserialized array")
  end subroutine test_integer_array_4d

  !> @brief Test for 5D integer array
  subroutine test_integer_array_5d()
    integer(int32), allocatable :: iarr5d(:,:,:,:,:)
    integer(int32), pointer :: iarr5d2(:,:,:,:,:)
    character(len=100) :: fname
    integer :: i

    allocate(iarr5d(2,1,2,1,2))
    iarr5d = reshape([(i, i=1,8)], [2,1,2,1,2])
    fname = "test_iarr5d.bin"

    call serialize_int_5d(iarr5d, fname)
    call deserialize_int_5d(iarr5d2, fname)

    call assert_equal_array_int(iarr5d, iarr5d2, size(iarr5d), "Original array does not match deserialized array")
  end subroutine test_integer_array_5d

  !> @brief Test for 1D real array
  subroutine test_real_array_1d()
    real(real64), ALLOCATABLE :: rarr1d(:)
    real(real64), POINTER :: rarr1d2(:)
    CHARACTER(len=100) :: fname

    allocate(rarr1d(4)) 
    rarr1d = [1.1_real64, 2.2_real64, 3.3_real64, 4.4_real64]
    fname = "test_rarr1d.bin"

    call serialize_real_1d(rarr1d, fname)
    call deserialize_real_1d(rarr1d2, fname)

    call assert_equal_array_real(rarr1d, rarr1d2, size(rarr1d), 1d-12, "Original array does not match deserialized array")
  end subroutine test_real_array_1d

  !> @brief Test for 2D real array
  subroutine test_real_array_2d()
    real(real64), ALLOCATABLE :: rarr2d(:,:)
    real(real64), POINTER :: rarr2d2(:,:)
    CHARACTER(len=100) :: fname

    allocate(rarr2d(2,2))
    rarr2d = reshape([1.5_real64, 2.5_real64, 3.5_real64, 4.5_real64], [2,2])
    fname = "test_rarr2d.bin"

    call serialize_real_2d(rarr2d, fname)
    call deserialize_real_2d(rarr2d2, fname)

    call assert_equal_array_real(rarr2d, rarr2d2, size(rarr2d), 1d-12, "Original array does not match deserialized array")
  end subroutine test_real_array_2d

  !> @brief Test for 3D real array
  subroutine test_real_array_3d()
    real(real64), ALLOCATABLE :: rarr3d(:,:,:)
    real(real64), POINTER :: rarr3d2(:,:,:)
    CHARACTER(len=100) :: fname

    allocate(rarr3d(2,2,2))
    rarr3d = reshape([1.0_real64,2.0_real64,3.0_real64,4.0_real64,5.0_real64,6.0_real64,7.0_real64,8.0_real64], [2,2,2])
    fname = "test_rarr3d.bin"
    
    call serialize_real_3d(rarr3d, fname)
    call deserialize_real_3d(rarr3d2, fname)

    call assert_equal_array_real(rarr3d, rarr3d2, size(rarr3d), 1d-12, "Original array does not match deserialized array")
  end subroutine test_real_array_3d

  !> @brief Test for 4D real array
  subroutine test_real_array_4d()
    real(real64), ALLOCATABLE :: rarr4d(:,:,:,:)
    real(real64), POINTER :: rarr4d2(:,:,:,:)
    CHARACTER(len=100) :: fname
    integer :: i

    allocate(rarr4d(2,2,1,2)) 
    rarr4d = reshape([(real(i,real64), i=1,8)], [2,2,1,2])
    fname = "test_rarr4d.bin"
    call serialize_real_4d(rarr4d, fname)
    call deserialize_real_4d(rarr4d2, fname)

    call assert_equal_array_real(rarr4d, rarr4d2, size(rarr4d), 1d-12, "Original array does not match deserialized array")
  end subroutine test_real_array_4d

  !> @brief Test for 5D real array
  subroutine test_real_array_5d()
    real(real64), ALLOCATABLE :: rarr5d(:,:,:,:,:)
    real(real64), POINTER :: rarr5d2(:,:,:,:,:)
    CHARACTER(len=100) :: fname
    integer :: i

    allocate(rarr5d(2,1,2,1,2)) 
    rarr5d = reshape([(real(i,real64), i=1,8)], [2,1,2,1,2])
    fname = "test_rarr5d.bin"
    call serialize_real_5d(rarr5d, fname)
    call deserialize_real_5d(rarr5d2, fname)

    call assert_equal_array_real(rarr5d, rarr5d2, size(rarr5d), 1d-12, "Original array does not match deserialized array")
  end subroutine test_real_array_5d

  !> @brief Test for 1D character array
  subroutine test_char_array_1d()
    character(len=:), ALLOCATABLE :: carr1d(:)
    character(len=:), pointer :: carr1d2(:)
    character(len=100) :: fname
    integer :: clen

    clen = 5
    allocate(character(len=clen):: carr1d(3))
    carr1d = ['foo  ','bar  ','baz  ']
    fname = "test_carr1d.bin"

    call serialize_char_1d(carr1d, fname)
    call deserialize_char_1d(carr1d2, fname)

    call assert_equal_array_char(carr1d, carr1d2, size(carr1d), "Original array does not match deserialized array")
  end subroutine test_char_array_1d

  !> @brief Test for 2D character array
  subroutine test_char_array_2d()
    character(len=:), ALLOCATABLE :: carr2d(:,:)
    character(len=:), pointer :: carr2d2(:,:)
    character(len=100) :: fname
    integer :: clen

    clen = 5
    allocate(character(len=clen):: carr2d(2,2))
    carr2d = reshape(['foo  ','bar  ','baz  ','qux  '], [2,2])
    fname = "test_carr2d.bin"

    call serialize_char_2d(carr2d, fname)
    call deserialize_char_2d(carr2d2, fname)

    call assert_equal_array_char(carr2d, carr2d2, size(carr2d), "Original array does not match deserialized array")
  end subroutine test_char_array_2d

  !> @brief Test for 3D character array
  subroutine test_char_array_3d()
    character(len=:), ALLOCATABLE :: carr3d(:,:,:)
    character(len=:), pointer :: carr3d2(:,:,:)
    character(len=100) :: fname
    integer :: clen

    clen = 5
    allocate(character(len=clen):: carr3d(2,2,1))
    carr3d = reshape(['foo  ','bar  ','baz  ','qux  '], [2,2,1])
    fname = "test_carr3d.bin"
    
    call serialize_char_3d(carr3d, fname)
    call deserialize_char_3d(carr3d2, fname)

    call assert_equal_array_char(carr3d, carr3d2, size(carr3d), "Original array does not match deserialized array")
  end subroutine test_char_array_3d

  !> @brief Test for 4D character array
  subroutine test_char_array_4d()
    character(len=:), ALLOCATABLE :: carr4d(:,:,:,:)
    character(len=:), pointer :: carr4d2(:,:,:,:)
    character(len=100) :: fname
    integer :: clen

    clen = 5
    allocate(character(len=clen):: carr4d(2,1,1,2))
    carr4d = reshape(['foo  ','bar  ','baz  ','qux  '], [2,1,1,2])
    fname = "test_carr4d.bin"

    call serialize_char_4d(carr4d, fname)
    call deserialize_char_4d(carr4d2, fname)

    call assert_equal_array_char(carr4d, carr4d2, size(carr4d), "Original array does not match deserialized array")
  end subroutine test_char_array_4d

  !> @brief Test for 5D character array
  subroutine test_char_array_5d()
    character(len=:), ALLOCATABLE :: carr5d(:,:,:,:,:)
    character(len=:), pointer :: carr5d2(:,:,:,:,:)
    character(len=100) :: fname
    integer :: clen

    clen = 5
    allocate(character(len=clen):: carr5d(2,1,2,1,2))
    carr5d = reshape(['foo  ','bar  ','baz  ','qux  ','aaa  ','bbb  ','ccc  ','ddd  '], [2,1,2,1,2])
    fname = "test_carr5d.bin"

    call serialize_char_5d(carr5d, fname)
    call deserialize_char_5d(carr5d2, fname)

    call assert_equal_array_char(carr5d, carr5d2, size(carr5d), "Original array does not match deserialized array")
  end subroutine test_char_array_5d

  ! Edge case tests
  !> @brief Test for 1x1 matrix
  subroutine test_integer_array_1x1()
    integer(int32), allocatable :: iarr(:,:)
    integer(int32), pointer :: iarr2(:,:)
    character(len=100) :: fname

    if (allocated(iarr)) deallocate(iarr)
    allocate(iarr(1,1))
    iarr = 42
    fname = "test_iarr_1x1.bin"

    call serialize_int_2d(iarr, fname)
    call deserialize_int_2d(iarr2, fname)

    call assert_equal_array_int(iarr, iarr2, size(iarr), "Original array does not match deserialized array")
  end subroutine test_integer_array_1x1

  !> @brief Test for empty array
  subroutine test_integer_array_empty()
    integer(int32), allocatable :: iarr(:,:)
    integer(int32), pointer :: iarr2(:,:)
    character(len=100) :: fname

    if (allocated(iarr)) deallocate(iarr)
    allocate(iarr(0,3))
    fname = "test_iarr_empty.bin"

    call serialize_int_2d(iarr, fname)
    call deserialize_int_2d(iarr2, fname)

    call assert_equal_array_int(iarr, iarr2, size(iarr), "Original array does not match deserialized array")
  end subroutine test_integer_array_empty

  !> @brief Test for protein data
  subroutine test_char_array_protein()
    character(len=:), pointer :: protein_data(:,:,:,:,:) => null()
    character(len=:), pointer :: protein_data_loaded(:,:,:,:,:) => null()
    character(len=*), parameter :: test_file = "test_proteins.bin"
    allocate(character(len=12) :: protein_data(2,1,2,1,2))  ! max 12 Symbols

    protein_data(1,1,1,1,1) = "P12345"         ! 6 Symbols
    protein_data(2,1,1,1,1) = "Q9Y6K8-2"       ! 8 Symbols
    protein_data(1,1,2,1,1) = "A0A024R3R5"     ! 10 Symbols
    protein_data(2,1,2,1,1) = "O14773"         ! 6 Symbols 
    protein_data(1,1,1,1,2) = "P0DTD1-PROT"    ! 11 Symbols
    protein_data(2,1,1,1,2) = ""               ! empty String
    protein_data(1,1,2,1,2) = "X6R8Y4"         ! 6 Symbols
    protein_data(2,1,2,1,2) = "A0A0B4J2F5-1"   ! 12 Symbols

    call serialize_char_5d(protein_data, test_file)
    call deserialize_char_5d(protein_data_loaded, test_file)
  end subroutine test_char_array_protein
end module mod_test_arrays
