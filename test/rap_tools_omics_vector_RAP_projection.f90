module mod_test_rap_tools_omics_vector_RAP_projection
   use asserts
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

   integer, parameter :: TEST_COUNT = 9

contains

   !> Get array of all available tests.
   function get_all_tests() result(all_tests)
      type(test_case) :: all_tests(TEST_COUNT)

      all_tests(1) = test_case("test_omics_vector_RAP_projection_all_selected", test_omics_vector_RAP_projection_all_selected)
      all_tests(2) = test_case("test_omics_vector_RAP_projection_one_axis_selected", test_omics_vector_RAP_projection_one_axis_selected)
      all_tests(3) = test_case("test_omics_vector_RAP_projection_one_vector_selected", test_omics_vector_RAP_projection_one_vector_selected)
      all_tests(4) = test_case("test_omics_vector_RAP_projection_constant_vector", test_omics_vector_RAP_projection_constant_vector)
      all_tests(5) = test_case("test_omics_vector_RAP_projection_orthogonal_vector", test_omics_vector_RAP_projection_orthogonal_vector)
      all_tests(6) = test_case("test_omics_vector_RAP_projection_no_axes", test_omics_vector_RAP_projection_no_axes)
      all_tests(7) = test_case("test_omics_vector_RAP_projection_no_vectors", test_omics_vector_RAP_projection_no_vectors)
      all_tests(8) = test_case("test_omics_vector_RAP_projection_mixed_selection", test_omics_vector_RAP_projection_mixed_selection)
      all_tests(9) = test_case("test_omics_vector_RAP_projection_non_square_vecs", test_omics_vector_RAP_projection_non_square_vecs)
   end function get_all_tests

   subroutine call_omics_vector_RAP_projection(test_name, vecs, axes_mask, vecs_mask)
      implicit none

      character(len=*), intent(in) :: test_name
      real(real64), dimension(:,:), intent(in) :: vecs
      logical, dimension(:), intent(in) :: axes_mask, vecs_mask

      integer(int32) :: n_axes, n_vecs, n_selected_axes, n_selected_vecs
      real(real64), allocatable :: projections(:,:)
      integer(int32) :: i_vec

      n_axes = size(axes_mask)
      n_vecs = size(vecs_mask)
      n_selected_axes = count(axes_mask)
      n_selected_vecs = count(vecs_mask)

      allocate(projections(n_selected_axes, n_selected_vecs))
      projections = 1

      call omics_vector_RAP_projection_r(vecs, n_axes, n_vecs, vecs_mask, n_selected_vecs, axes_mask, n_selected_axes, projections)

      do i_vec = 1, n_selected_vecs
         call assert_equal_real(&
            sum(projections(:, i_vec)) / real(n_selected_axes, real64),&
            0.0_real64,&
            1d-12,&
            "test_omics_vector_RAP_projection: Case '" // trim(test_name) // "' failed"&
         )
      end do

      deallocate(projections)
   end subroutine call_omics_vector_RAP_projection

   subroutine test_omics_vector_RAP_projection_all_selected()
      implicit none

      real(real64), dimension(3,3) :: vecs
      logical :: axes_mask(3), vecs_mask(3)

      vecs = reshape([1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0], [3,3])
      axes_mask = [.true., .true., .true.]
      vecs_mask = [.true., .true., .true.]
      call call_omics_vector_RAP_projection("All selected", vecs, axes_mask, vecs_mask)
   end subroutine test_omics_vector_RAP_projection_all_selected

   subroutine test_omics_vector_RAP_projection_one_axis_selected()
      implicit none

      real(real64), dimension(3,3) :: vecs
      logical :: axes_mask(3), vecs_mask(3)

      vecs = reshape([1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0], [3,3])
      axes_mask = [.false., .true., .false.]
      vecs_mask = [.true., .true., .true.]
      call call_omics_vector_RAP_projection("One axis selected", vecs, axes_mask, vecs_mask)
   end subroutine test_omics_vector_RAP_projection_one_axis_selected

   subroutine test_omics_vector_RAP_projection_one_vector_selected()
      implicit none

      real(real64), dimension(3,3) :: vecs
      logical :: axes_mask(3), vecs_mask(3)

      vecs = reshape([1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0], [3,3])
      axes_mask = [.true., .true., .true.]
      vecs_mask = [.false., .true., .false.]
      call call_omics_vector_RAP_projection("One vector selected", vecs, axes_mask, vecs_mask)
   end subroutine test_omics_vector_RAP_projection_one_vector_selected

   subroutine test_omics_vector_RAP_projection_constant_vector()
      implicit none

      real(real64), dimension(3,3) :: vecs
      logical :: axes_mask(3), vecs_mask(3)

      vecs = 0.0_real64
      vecs(:,1) = [5.0, 5.0, 5.0]
      axes_mask = [.true., .true., .true.]
      vecs_mask = [.true., .false., .false.]
      call call_omics_vector_RAP_projection("Constant vector", vecs, axes_mask, vecs_mask)
   end subroutine test_omics_vector_RAP_projection_constant_vector

   subroutine test_omics_vector_RAP_projection_orthogonal_vector()
      implicit none

      real(real64), dimension(3,3) :: vecs
      logical :: axes_mask(3), vecs_mask(3)

      vecs = 0.0_real64
      vecs(:,1) = [1.0, 0.0, -1.0]
      axes_mask = [.true., .true., .true.]
      vecs_mask = [.true., .false., .false.]
      call call_omics_vector_RAP_projection("Orthogonal vector", vecs, axes_mask, vecs_mask)
   end subroutine test_omics_vector_RAP_projection_orthogonal_vector

   subroutine test_omics_vector_RAP_projection_no_axes()
      implicit none

      real(real64), dimension(3,3) :: vecs
      logical :: axes_mask(3), vecs_mask(3)

      vecs = reshape([1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0], [3,3])
      axes_mask = [.false., .false., .false.]
      vecs_mask = [.true., .true., .true.]
      call call_omics_vector_RAP_projection("No axes selected", vecs, axes_mask, vecs_mask)
   end subroutine test_omics_vector_RAP_projection_no_axes

   subroutine test_omics_vector_RAP_projection_no_vectors()
      implicit none

      real(real64), dimension(3,3) :: vecs
      logical :: axes_mask(3), vecs_mask(3)

      vecs = reshape([1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0], [3,3])
      axes_mask = [.true., .true., .true.]
      vecs_mask = [.false., .false., .false.]
      call call_omics_vector_RAP_projection("No vectors selected", vecs, axes_mask, vecs_mask)
   end subroutine test_omics_vector_RAP_projection_no_vectors

   subroutine test_omics_vector_RAP_projection_mixed_selection()
      implicit none

      real(real64), dimension(3,3) :: vecs
      logical :: axes_mask(3), vecs_mask(3)

      vecs = reshape([1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0], [3,3])
      axes_mask = [.true., .false., .true.]
      vecs_mask = [.true., .false., .true.]
      call call_omics_vector_RAP_projection("Mixed selection", vecs, axes_mask, vecs_mask)
   end subroutine test_omics_vector_RAP_projection_mixed_selection

   subroutine test_omics_vector_RAP_projection_non_square_vecs()
      real(real64), dimension(4,2) :: vecs
      logical :: axes_mask(4), vecs_mask(2)

      ! Define two 4D vectors
      vecs(:,1) = [1.0, 2.0, 3.0, 4.0]
      vecs(:,2) = [4.0, 3.0, 2.0, 1.0]

      ! Select all axes and both vectors
      axes_mask = [.true., .true., .true., .true.]
      vecs_mask = [.true., .true.]

      call call_omics_vector_RAP_projection("Non-square vecs (4D x 2)", vecs, axes_mask, vecs_mask)
   end subroutine test_omics_vector_RAP_projection_non_square_vecs




   !> Run all omics_vector_RAP_projection test.
   subroutine run_all_tests_rap_tools_omics_vector_RAP_projection()
      type(test_case) :: all_tests(TEST_COUNT)
      integer :: i

      all_tests = get_all_tests()

      do i = 1, size(all_tests)
         call all_tests(i)%test_proc()
         print *, trim(all_tests(i)%name), " passed."
      end do
      print *, "All rap_tools_omics_vector_RAP_projection tests passed successfully."
   end subroutine run_all_tests_rap_tools_omics_vector_RAP_projection

   !> Run specific omics_vector_RAP_projection tests by name.
   subroutine run_named_tests_rap_tools_omics_vector_RAP_projection(test_names)
      character(len=*), intent(in) :: test_names(:)
      type(test_case) :: all_tests(TEST_COUNT)
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
   end subroutine run_named_tests_rap_tools_omics_vector_RAP_projection
end module mod_test_rap_tools_omics_vector_RAP_projection
