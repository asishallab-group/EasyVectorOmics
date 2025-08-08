!> Module for tools related to relative axis planes (RAPs), so planes in the higher dimensional gene expression space
module relative_axis_plane_tools
   use, intrinsic :: iso_fortran_env, only: real64, int32

   private
   public :: omics_vector_RAP_projection
contains

   !> Project selected vectors (e.g. expression vectors) onto the RAP constructed from a selected set of axes.
   pure subroutine omics_vector_RAP_projection(vecs, n_axes, n_vecs, vecs_selection_mask, n_selected_vecs, axes_selection_mask, n_selected_axes, projections)
      implicit none

      real(real64), dimension(n_axes, n_vecs), intent(in) :: vecs
         !! matrix with expression vectors
      integer(int32), intent(in) :: n_axes
         !! number of axes
      integer(int32), intent(in) :: n_vecs
         !! number of vectors per axis
      logical, dimension(n_vecs), intent(in) :: vecs_selection_mask
         !! `.true.` for vectors where projection is to be computed
      integer(int32), intent(in) :: n_selected_vecs
         !! count of `.true.` values in `vecs_selection_mask`
      logical, dimension(n_axes), intent(in) :: axes_selection_mask
         !! `.true.` for axes to be included in RAP
      integer(int32), intent(in) :: n_selected_axes
         !! count of `.true.` values in `axes_selection_mask`
      real(real64), dimension(n_selected_axes, n_selected_vecs), intent(out) :: projections
         !! projected vectors

      ! fill projections matrix with selected vectors
      integer(int32) :: i_vec, i_axis, i_vec_proj, i_axis_proj
      i_vec_proj = 1
      do i_vec = 1, n_vecs
         if (vecs_selection_mask(i_vec)) then
            i_axis_proj = 1
            do i_axis = 1, n_axes
               if (axes_selection_mask(i_axis)) then
                  projections(i_axis_proj, i_vec_proj) = vecs(i_axis, i_vec)

                  i_axis_proj = i_axis_proj + 1
               end if
            end do

            i_vec_proj = i_vec_proj + 1
         end if
      end do

      call project_selected_vecs_onto_rap(projections, n_selected_axes, n_selected_vecs)
   end subroutine omics_vector_RAP_projection

   !> Projects selected vectors onto its RAP
   pure subroutine project_selected_vecs_onto_rap(selected_vecs, n_selected_axes, n_selected_vecs)
      implicit none

      real(real64), dimension(n_selected_axes, n_selected_vecs), intent(inout) :: selected_vecs
         !! matrix with vectors for selected axes
      integer(int32), intent(in) :: n_selected_axes
         !! number of selected axes
      integer(int32), intent(in) :: n_selected_vecs
         !! number of selected vectors per axis

      ! project selected vectors onto RAP
      integer(int32) :: i_vec, i_axis
      real(real64) :: diagonal_component

      do i_vec = 1, n_selected_vecs

         ! calculate diagonal component to be subtracted from vectors for projection
         diagonal_component = 0.0_real64
         do i_axis = 1, n_selected_axes
            diagonal_component = diagonal_component + selected_vecs(i_axis, i_vec)
         end do
         diagonal_component = diagonal_component / n_selected_axes

         ! transform vector to its projection
         do i_axis = 1, n_selected_axes
            selected_vecs(i_axis, i_vec) = selected_vecs(i_axis, i_vec) - diagonal_component
         end do
      end do
   end subroutine project_selected_vecs_onto_rap
end module relative_axis_plane_tools

subroutine omics_vector_RAP_projection_r(vecs, n_axes, n_vecs, vecs_selection_mask, n_selected_vecs, axes_selection_mask, n_selected_axes, projections)
   use relative_axis_plane_tools
   use, intrinsic :: iso_fortran_env, only: real64, int32
   implicit none

   real(real64), dimension(n_axes, n_vecs), intent(in) :: vecs
      !! matrix with expression vectors
   integer(int32), intent(in) :: n_axes
      !! number of axes
   integer(int32), intent(in) :: n_vecs
      !! number of vectors per axis
   logical, dimension(n_vecs), intent(in) :: vecs_selection_mask
      !! `.true.` for vectors where projection is to be computed
   integer(int32), intent(in) :: n_selected_vecs
      !! count of `.true.` values in `vecs_selection_mask`
   logical, dimension(n_axes), intent(in) :: axes_selection_mask
      !! `.true.` for axes to be included in RAP
   integer(int32), intent(in) :: n_selected_axes
      !! count of `.true.` values in `axes_selection_mask`
   real(real64), dimension(n_selected_axes, n_selected_vecs), intent(out) :: projections
      !! projected vectors

   call omics_vector_RAP_projection(vecs, n_axes, n_vecs, vecs_selection_mask, n_selected_vecs, axes_selection_mask, n_selected_axes, projections)
end subroutine omics_vector_RAP_projection_r

subroutine omics_vector_RAP_projection_c(vecs, n_axes, n_vecs, vecs_selection_mask, n_selected_vecs, axes_selection_mask, n_selected_axes, projections) bind(C, name="omics_vector_RAP_projection_c")
   use iso_c_binding
   use relative_axis_plane_tools
   implicit none

   real(c_double), dimension(n_axes, n_vecs), intent(in) :: vecs
      !! matrix with expression vectors
   integer(c_int), intent(in), value :: n_axes
      !! number of axes
   integer(c_int), intent(in), value :: n_vecs
      !! number of vectors per axis
   integer(c_int), dimension(n_vecs), intent(in) :: vecs_selection_mask
      !! `.true.` for vectors where projection is to be computed
   integer(c_int), intent(in), value :: n_selected_vecs
      !! count of `.true.` values in `vecs_selection_mask`
   integer(c_int), dimension(n_axes), intent(in) :: axes_selection_mask
      !! `.true.` for axes to be included in RAP
   integer(c_int), intent(in), value :: n_selected_axes
      !! count of `.true.` values in `axes_selection_mask`
   real(c_double), dimension(n_selected_axes, n_selected_vecs), intent(out) :: projections
      !! projected vectors

   call omics_vector_RAP_projection(vecs, n_axes, n_vecs, vecs_selection_mask /= 0, n_selected_vecs, axes_selection_mask /= 0, n_selected_axes, projections)
end subroutine omics_vector_RAP_projection_c