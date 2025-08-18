! filepath: test/mod_test_shift_vectors.f90
!> Unit test suite for compute shift vectors routines.
module mod_test_shift_vectors
   use asserts
   use tox_shift_vectors
   use tox_errors, only: ERR_OK, ERR_EMPTY_INPUT, ERR_DIM_MISMATCH
   use, intrinsic :: iso_fortran_env, only: real64, int32
   implicit none
   public

   ! Abstract interface for all test procedures
   abstract interface
      subroutine test_interface()
      end subroutine test_interface
   end interface

   ! Type to hold test name and procedure pointer
   type :: test_case
      character(len=64) :: name
      procedure(test_interface), pointer, nopass :: test_proc => null()
   end type test_case

contains

   !> Get array of all available tests.
   function get_all_tests() result(all_tests)
      type(test_case) :: all_tests(6)
      all_tests(1) = test_case("test_correct_family_mapping", test_correct_family_mapping)
      all_tests(2) = test_case("test_invalid_family_mapping", test_invalid_family_mapping)
      all_tests(3) = test_case("test_zero_distance", test_zero_distance)
      all_tests(4) = test_case("test_multiple_genes_per_family", test_multiple_genes_per_family)
      all_tests(5) = test_case("test_single_gene_per_family", test_single_gene_per_family)
      all_tests(6) = test_case("test_dimension_edge_cases", test_dimension_edge_cases)
   end function get_all_tests

   !> Run all shift vector tests.
   subroutine run_all_tests_shift_vectors()
      type(test_case) :: all_tests(6)
      integer(int32) :: i
      all_tests = get_all_tests()
      do i = 1, size(all_tests)
         call all_tests(i)%test_proc()
         print *, trim(all_tests(i)%name), " passed."
      end do
      print *, "All shift vector tests passed successfully."
   end subroutine run_all_tests_shift_vectors

   !> Run specific shift vector tests by name.
   subroutine run_named_tests_shift_vectors(test_names)
      character(len=*), intent(in) :: test_names(:)
      type(test_case) :: all_tests(6)
      integer(int32) :: i, j
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
   end subroutine run_named_tests_shift_vectors

   !> Test correct mapping between families and genes. //TODO change to simple test without douplicates on families
   subroutine test_correct_family_mapping()
      real(real64) :: expression_vectors(3, 5), family_centroids(3, 3), shift_vectors(6, 5)
      integer(int32) :: gene_to_family(5), family_ids(3), ierr, i
      expression_vectors = reshape([(real(i, real64), i=1, 15)], [3, 5])
      family_centroids = reshape([(real(i, real64), i=5, -3, -1)], [3, 3])
      gene_to_family = [3, 2, 1, 3, 1]
      family_ids = [2, 3, 1]

     call compute_shift_vector_field(3, 5, 3, expression_vectors, family_centroids, gene_to_family, family_ids, shift_vectors, ierr)

      ! Check array size stays correct
      call assert_true(size(shift_vectors, 1) == 6 .and. size(shift_vectors, 2) == 5, "Shift_vectors shape is incorrect")

      !//TODO change to assert real array
      ! Check values in shift_vectors
      call assert_equal_real(shift_vectors(1, 1), 2.0_real64, 1e-12_real64, "Shift_vectors(1, 1) should be 2")
      call assert_equal_real(shift_vectors(2, 1), 1.0_real64, 1e-12_real64, "Shift_vectors(2, 1) should be 1")
      call assert_equal_real(shift_vectors(3, 1), 0.0_real64, 1e-12_real64, "Shift_vectors(3, 1) should be 0")
      call assert_equal_real(shift_vectors(4, 1), -1.0_real64, 1e-12_real64, "Shift_vectors(4, 1) should be -1")
      call assert_equal_real(shift_vectors(5, 1), 1.0_real64, 1e-12_real64, "Shift_vectors(5, 1) should be 1")
      call assert_equal_real(shift_vectors(6, 1), 3.0_real64, 1e-12_real64, "Shift_vectors(6, 1) should be 3")

      call assert_equal_real(shift_vectors(1, 2), 5.0_real64, 1e-12_real64, "Shift_vectors(1, 2) should be 5")
      call assert_equal_real(shift_vectors(2, 2), 4.0_real64, 1e-12_real64, "Shift_vectors(2, 2) should be 4")
      call assert_equal_real(shift_vectors(3, 2), 3.0_real64, 1e-12_real64, "Shift_vectors(3, 2) should be 3")
      call assert_equal_real(shift_vectors(4, 2), -1.0_real64, 1e-12_real64, "Shift_vectors(4, 2) should be -1")
      call assert_equal_real(shift_vectors(5, 2), 1.0_real64, 1e-12_real64, "Shift_vectors(5, 2) should be 1")
      call assert_equal_real(shift_vectors(6, 2), 3.0_real64, 1e-12_real64, "Shift_vectors(6, 2) should be 3")

      call assert_equal_real(shift_vectors(1, 3), -1.0_real64, 1e-12_real64, "Shift_vectors(1, 3) should be -1")
      call assert_equal_real(shift_vectors(2, 3), -2.0_real64, 1e-12_real64, "Shift_vectors(2, 3) should be -2")
      call assert_equal_real(shift_vectors(3, 3), -3.0_real64, 1e-12_real64, "Shift_vectors(3, 3) should be -3")
      call assert_equal_real(shift_vectors(4, 3), 8.0_real64, 1e-12_real64, "Shift_vectors(4, 3) should be 8")
      call assert_equal_real(shift_vectors(5, 3), 10.0_real64, 1e-12_real64, "Shift_vectors(5, 3) should be 10")
      call assert_equal_real(shift_vectors(6, 3), 12.0_real64, 1e-12_real64, "Shift_vectors(6, 3) should be 12")

      call assert_equal_real(shift_vectors(1, 4), 2.0_real64, 1e-12_real64, "Shift_vectors(1, 4) should be 2")
      call assert_equal_real(shift_vectors(2, 4), 1.0_real64, 1e-12_real64, "Shift_vectors(2, 4) should be 1")
      call assert_equal_real(shift_vectors(3, 4), 0.0_real64, 1e-12_real64, "Shift_vectors(3, 4) should be 0")
      call assert_equal_real(shift_vectors(4, 4), 8.0_real64, 1e-12_real64, "Shift_vectors(4, 4) should be 8")
      call assert_equal_real(shift_vectors(5, 4), 10.0_real64, 1e-12_real64, "Shift_vectors(5, 4) should be 10")
      call assert_equal_real(shift_vectors(6, 4), 12.0_real64, 1e-12_real64, "Shift_vectors(6, 4) should be 12")

      call assert_equal_real(shift_vectors(1, 5), -1.0_real64, 1e-12_real64, "Shift_vectors(1, 5) should be -1")
      call assert_equal_real(shift_vectors(2, 5), -2.0_real64, 1e-12_real64, "Shift_vectors(2, 5) should be -2")
      call assert_equal_real(shift_vectors(3, 5), -3.0_real64, 1e-12_real64, "Shift_vectors(3, 5) should be -3")
      call assert_equal_real(shift_vectors(4, 5), 14.0_real64, 1e-12_real64, "Shift_vectors(4, 5) should be 14")
      call assert_equal_real(shift_vectors(5, 5), 16.0_real64, 1e-12_real64, "Shift_vectors(5, 5) should be 16")
      call assert_equal_real(shift_vectors(6, 5), 18.0_real64, 1e-12_real64, "Shift_vectors(6, 5) should be 18")
   end subroutine test_correct_family_mapping

   !> Test for invalid family id mapping raising error
   subroutine test_invalid_family_mapping()
      real(real64) :: expression_vectors(3, 2), family_centroids(3, 3), shift_vectors(6, 2)
      integer(int32) :: gene_to_family(2), family_ids(3), ierr, i
      expression_vectors = reshape([(real(i, real64), i=1, 6)], [3, 2])
      family_centroids = reshape([(real(i, real64), i=5, -3, -1)], [3, 3])
      ! Invalid gene_to_family id (4) should raise an error, as there are only 3 family_ids.
      gene_to_family = [3, 4]
      family_ids = [2, 3, 1]

      ! Call the function with invalid mapping
     call compute_shift_vector_field(3, 2, 3, expression_vectors, family_centroids, gene_to_family, family_ids, shift_vectors, ierr)

      ! Check for expected ERR_INVALID_INPUT error
      call assert_equal_int(ierr, ERR_INVALID_INPUT, "Invalid family id mapping should return ERR_INVALID_INPUT")
   end subroutine test_invalid_family_mapping

   !> Test for zero distance between paralog and centroid
   subroutine test_zero_distance()
      real(real64) :: expression_vectors(3, 2), family_centroids(3, 2), shift_vectors(6, 2), expected_shift_vectors(6, 2)
      integer(int32) :: gene_to_family(2), family_ids(2), ierr

      ! Expression vectors and family_centroids have the same values and get mapped 1 to 1 and 2 to 2
      expression_vectors = reshape([1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64, 6.0_real64], [3, 2])
      family_centroids = reshape([1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64, 6.0_real64], [3, 2])
      gene_to_family = [1, 2]
      family_ids = [1, 2]

      ! Call the function with zero distance values
     call compute_shift_vector_field(3, 2, 2, expression_vectors, family_centroids, gene_to_family, family_ids, shift_vectors, ierr)

      ! Check for exprected 0 values in shift_vectors
      expected_shift_vectors = reshape([1.0_real64,2.0_real64,3.0_real64,0.0_real64,0.0_real64,0.0_real64,4.0_real64,5.0_real64,6.0_real64,0.0_real64,0.0_real64,0.0_real64], [6,2])
      call assert_equal_array_real(shift_vectors, expected_shift_vectors, 6, 1e-12_real64, "Shift vectors should be zero distance from centroids")
   end subroutine

   !> Test for multiple genes per family
   subroutine test_multiple_genes_per_family()
      real(real64) :: expression_vectors(2, 4), family_centroids(2, 2), shift_vectors(4, 4), expected_shift_vectors(4, 4)
      integer(int32) :: gene_to_family(4), family_ids(2), ierr, i

      ! Expression vectors belong to 2 different families
      expression_vectors = reshape([(real(i, real64), i=1, 8)], [2, 4])
      family_centroids = reshape([10.0_real64, 20.0_real64, 30.0_real64, 40.0_real64], [2, 2])
      gene_to_family = [1, 2, 1, 2]
      family_ids = [1, 2]

      expected_shift_vectors = reshape([10.0_real64, 20.0_real64, -9.0_real64, -18.0_real64, 30.0_real64, 40.0_real64, -27.0_real64, -36.0_real64, 50.0_real64, 60.0_real64, -5.0_real64, -14.0_real64, 70.0_real64, 80.0_real64, -23.0_real64, -23.0_real64], [4, 4])

      ! Call the function with multiple genes per family
     call compute_shift_vector_field(2, 4, 2, expression_vectors, family_centroids, gene_to_family, family_ids, shift_vectors, ierr)

  call assert_equal_array_real(shift_vectors, expected_shift_vectors, 4, 1e-12_real64, "Shift vectors should match expected values")

   end subroutine test_multiple_genes_per_family

   !> Test for single gene per family
   subroutine test_single_gene_per_family()
      real(real64) :: expression_vectors(2, 4), family_centroids(2, 4), shift_vectors(4, 4), expected_shift_vectors(4, 4)
      integer(int32) :: gene_to_family(4), family_ids(4), ierr, i

      ! Each expression vector belongs to one single family
      expression_vectors = reshape([(real(i, real64), i=1, 8)], [2, 4])
      family_centroids = reshape([(real(i, real64), i=10, 80, 10)], [2, 4])
      gene_to_family = [1, 2, 3, 4]
      family_ids = [1, 2, 3, 4]

      expected_shift_vectors = reshape([10.0_real64, 20.0_real64, -9.0_real64, -18.0_real64, 30.0_real64, 40.0_real64, -27.0_real64, -36.0_real64, 50.0_real64, 60.0_real64, -45.0_real64, -54.0_real64, 70.0_real64, 80.0_real64, -63.0_real64, -72.0_real64], [4, 4])

      ! Call the function with single genes per family
     call compute_shift_vector_field(2, 4, 4, expression_vectors, family_centroids, gene_to_family, family_ids, shift_vectors, ierr)

  call assert_equal_array_real(shift_vectors, expected_shift_vectors, 4, 1e-12_real64, "Shift vectors should match expected values")

   end subroutine test_single_gene_per_family

   !> Test for dimension edge cases (0 genes with dimension 1 and 1 family)
   subroutine test_dimension_edge_cases()
      real(real64) :: expression_vectors(1, 0), family_centroids(1, 1), shift_vectors(2, 0)
      integer(int32) :: gene_to_family(0), family_ids(1), ierr

      ! Call the function with edge case arrays
     call compute_shift_vector_field(1, 0, 1, expression_vectors, family_centroids, gene_to_family, family_ids, shift_vectors, ierr)

      !//TODO Should this raise a empty input error?

      ! Check for expected 0 length of shift_vectors array
      call assert_equal_int(size(shift_vectors), 0, "Shift vectors array length should be 0")

   end subroutine test_dimension_edge_cases

end module mod_test_shift_vectors

!//TODO check back in code guidelines for correct testing
!//TODO check auto formatting