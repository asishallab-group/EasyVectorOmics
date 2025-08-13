!> Module for tools related to relative axis planes (RAPs), so planes in the higher dimensional gene expression space
module relative_axis_plane_tools
   use, intrinsic :: iso_fortran_env, only: real64, int32

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

   !> Compute the signed clock hand angle between two RAP-projected and normalized vectors.
   !| Calculates the signed rotation angle between two normalized vectors in RAP space.
   !| For 2D/3D: automatic directionality calculation. For >3D: uses selected axes for directionality.
   pure subroutine clock_hand_angle_between_vectors(v1, v2, n_dims, signed_angle, selected_axes_for_signed)
      implicit none

      real(real64), dimension(n_dims), intent(in) :: v1
         !! First normalized vector in RAP space
      real(real64), dimension(n_dims), intent(in) :: v2
         !! Second normalized vector in RAP space  
      integer(int32), intent(in) :: n_dims
         !! Dimension of both vectors
      real(real64), intent(out) :: signed_angle
         !! Signed angle between vectors in radians [-π, π]
      integer(int32), dimension(3), intent(in) :: selected_axes_for_signed
         !! Indices of 3 axes to use for directionality calculation (ignored if n_dims <= 3)

      real(real64) :: dot_product, unsigned_angle, orientation_sign
      integer(int32) :: i

      ! Calculate dot product of normalized vectors
      dot_product = 0.0_real64
      do i = 1, n_dims
         dot_product = dot_product + v1(i) * v2(i)
      end do

      ! Clamp dot product to [-1, 1] to handle numerical precision issues
      dot_product = max(-1.0_real64, min(1.0_real64, dot_product))

      ! Calculate unsigned angle using arccos
      unsigned_angle = acos(dot_product)

      ! Calculate orientation sign for directionality
      if (n_dims == 2) then
         ! For 2D: use determinant directly
         orientation_sign = sign(1.0_real64, v1(1) * v2(2) - v1(2) * v2(1))
      else if (n_dims == 3) then
         ! For 3D: use z-component of cross product
         orientation_sign = sign(1.0_real64, v1(1) * v2(2) - v1(2) * v2(1))
      else
         ! For >3D: use selected axes to form 3D subspace
         orientation_sign = sign(1.0_real64, &
            v1(selected_axes_for_signed(1)) * v2(selected_axes_for_signed(2)) - &
            v1(selected_axes_for_signed(2)) * v2(selected_axes_for_signed(1)))
      end if

      ! Apply sign to unsigned angle
      signed_angle = orientation_sign * unsigned_angle
   end subroutine clock_hand_angle_between_vectors

   !> Compute signed rotation angles between RAP-projected and normalized vector pairs.
   !| Takes separate arrays of RAP-projected and normalized vectors (e.g. expression 
   !| centroids and paralogs) and computes the signed rotation angle between corresponding pairs.
   !| This measures both magnitude and directionality of angular separation in RAP space.
   pure subroutine clock_hand_angles_for_shift_vectors(origins, targets, n_dims, n_vecs, &
                                                      vecs_selection_mask, &
                                                      n_selected_vecs, selected_axes_for_signed, &
                                                      signed_angles)
      implicit none

      real(real64), dimension(n_dims, n_vecs), intent(in) :: origins
         !! First set of RAP-projected, normalized vectors (e.g. expression centroids)
      real(real64), dimension(n_dims, n_vecs), intent(in) :: targets
         !! Second set of RAP-projected, normalized vectors (e.g. paralogs)
      integer(int32), intent(in) :: n_dims
         !! Dimension of each vector in RAP space
      integer(int32), intent(in) :: n_vecs
         !! Number of vector pairs
      logical, dimension(n_vecs), intent(in) :: vecs_selection_mask
         !! .true. for vector pairs where angle should be computed
      integer(int32), intent(in) :: n_selected_vecs
         !! Count of .true. values in vecs_selection_mask
      integer(int32), dimension(3), intent(in) :: selected_axes_for_signed
         !! Indices of 3 axes to use for directionality calculation (ignored if n_dims <= 3)
      real(real64), dimension(n_selected_vecs), intent(out) :: signed_angles
         !! Signed rotation angles between vector pairs in radians [-π, π]

      integer(int32) :: i_vec, i_result

      i_result = 1
      do i_vec = 1, n_vecs
         if (vecs_selection_mask(i_vec)) then
            call clock_hand_angle_between_vectors(origins(:, i_vec), targets(:, i_vec), n_dims, &
                                                 signed_angles(i_result), selected_axes_for_signed)
            
            i_result = i_result + 1
         end if
      end do
   end subroutine clock_hand_angles_for_shift_vectors


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

subroutine clock_hand_angle_between_vectors_r(v1, v2, n_dims, signed_angle, selected_axes_for_signed)
   use relative_axis_plane_tools
   use, intrinsic :: iso_fortran_env, only: real64, int32
   implicit none

   real(real64), dimension(n_dims), intent(in) :: v1
      !! First normalized vector in RAP space
   real(real64), dimension(n_dims), intent(in) :: v2
      !! Second normalized vector in RAP space  
   integer(int32), intent(in) :: n_dims
      !! Dimension of both vectors
   real(real64), intent(out) :: signed_angle
      !! Signed angle between vectors in radians [-π, π]
   integer(int32), dimension(3), intent(in) :: selected_axes_for_signed
      !! Indices of 3 axes to use for directionality calculation (ignored if n_dims <= 3)

   call clock_hand_angle_between_vectors(v1, v2, n_dims, signed_angle, selected_axes_for_signed)
end subroutine clock_hand_angle_between_vectors_r

subroutine clock_hand_angle_between_vectors_c(v1, v2, n_dims, signed_angle, selected_axes_for_signed) bind(C, name="clock_hand_angle_between_vectors_c")
   use iso_c_binding
   use relative_axis_plane_tools
   implicit none

   real(c_double), dimension(n_dims), intent(in) :: v1
      !! First normalized vector in RAP space
   real(c_double), dimension(n_dims), intent(in) :: v2
      !! Second normalized vector in RAP space  
   integer(c_int), intent(in), value :: n_dims
      !! Dimension of both vectors
   real(c_double), intent(out) :: signed_angle
      !! Signed angle between vectors in radians [-π, π]
   integer(c_int), dimension(3), intent(in) :: selected_axes_for_signed
      !! Indices of 3 axes to use for directionality calculation (ignored if n_dims <= 3)

   call clock_hand_angle_between_vectors(v1, v2, n_dims, signed_angle, selected_axes_for_signed)
end subroutine clock_hand_angle_between_vectors_c

subroutine clock_hand_angles_for_shift_vectors_r(origins, targets, n_dims, n_vecs, vecs_selection_mask, n_selected_vecs, selected_axes_for_signed, signed_angles)
   use relative_axis_plane_tools
   use, intrinsic :: iso_fortran_env, only: real64, int32
   implicit none

   real(real64), dimension(n_dims, n_vecs), intent(in) :: origins
      !! First set of RAP-projected, normalized vectors (e.g. expression centroids)
   real(real64), dimension(n_dims, n_vecs), intent(in) :: targets
      !! Second set of RAP-projected, normalized vectors (e.g. paralogs)
   integer(int32), intent(in) :: n_dims
      !! Dimension of each vector in RAP space
   integer(int32), intent(in) :: n_vecs
      !! Number of vector pairs
   logical, dimension(n_vecs), intent(in) :: vecs_selection_mask
      !! .true. for vector pairs where angle should be computed
   integer(int32), intent(in) :: n_selected_vecs
      !! Count of .true. values in vecs_selection_mask
   integer(int32), dimension(3), intent(in) :: selected_axes_for_signed
      !! Indices of 3 axes to use for directionality calculation (ignored if n_dims <= 3)
   real(real64), dimension(n_selected_vecs), intent(out) :: signed_angles
      !! Signed rotation angles between vector pairs in radians [-π, π]

   call clock_hand_angles_for_shift_vectors(origins, targets, n_dims, n_vecs, vecs_selection_mask, n_selected_vecs, selected_axes_for_signed, signed_angles)
end subroutine clock_hand_angles_for_shift_vectors_r

subroutine clock_hand_angles_for_shift_vectors_c(origins, targets, n_dims, n_vecs, vecs_selection_mask, n_selected_vecs, selected_axes_for_signed, signed_angles) bind(C, name="clock_hand_angles_for_shift_vectors_c")
   use iso_c_binding
   use relative_axis_plane_tools
   implicit none

   real(c_double), dimension(n_dims, n_vecs), intent(in) :: origins
      !! First set of RAP-projected, normalized vectors (e.g. expression centroids)
   real(c_double), dimension(n_dims, n_vecs), intent(in) :: targets
      !! Second set of RAP-projected, normalized vectors (e.g. paralogs)
   integer(c_int), intent(in), value :: n_dims
      !! Dimension of each vector in RAP space
   integer(c_int), intent(in), value :: n_vecs
      !! Number of vector pairs
   integer(c_int), dimension(n_vecs), intent(in) :: vecs_selection_mask
      !! .true. for vector pairs where angle should be computed
   integer(c_int), intent(in), value :: n_selected_vecs
      !! Count of .true. values in vecs_selection_mask
   integer(c_int), dimension(3), intent(in) :: selected_axes_for_signed
      !! Indices of 3 axes to use for directionality calculation (ignored if n_dims <= 3)
   real(c_double), dimension(n_selected_vecs), intent(out) :: signed_angles
      !! Signed rotation angles between vector pairs in radians [-π, π]

   call clock_hand_angles_for_shift_vectors(origins, targets, n_dims, n_vecs, vecs_selection_mask /= 0, n_selected_vecs, selected_axes_for_signed, signed_angles)
end subroutine clock_hand_angles_for_shift_vectors_c
