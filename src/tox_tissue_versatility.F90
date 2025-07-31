!> Module for calculating normalized tissue (axis) versatility.
!! This module implements the angle-based metric for tissue versatility,
!! quantifying how uniformly a gene is expressed across selected axes (tissues).
module avmod
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
contains

  !> Computes normalized tissue versatility for selected expression vectors.
  !! The metric is based on the angle between each gene expression vector and the space diagonal.
  !! Versatility is normalized to [0, 1], where 0 means uniform expression and 1 means expression in only one axis.
  !!
  !! @param n_axes Number of axes (tissues/dimensions).
  !! @param n_vectors Number of expression vectors (genes).
  !! @param expression_vectors 2D array (n_axes, n_vectors), each column is a gene expression vector.
  !! @param exp_vecs_selection_index Logical array (n_vectors), .TRUE. for vectors to process.
  !! @param axes_selection Logical array (n_axes), .TRUE. for axes to include in calculation.
  !! @param tissue_versatilities Output: real array, length = count(exp_vecs_selection_index), stores the calculated tissue versatilities.
  pure subroutine compute_tissue_versatility(n_axes, n_vectors, expression_vectors, &
                                             exp_vecs_selection_index, axes_selection, &
                                             tissue_versatilities, tissue_angles_deg)
    integer, intent(in) :: n_axes, n_vectors
    real(real64), intent(in) :: expression_vectors(n_axes, n_vectors)
    logical, intent(in) :: exp_vecs_selection_index(n_vectors)
    logical, intent(in) :: axes_selection(n_axes)
    real(real64), intent(out) :: tissue_versatilities(count(exp_vecs_selection_index))
    real(real64), intent(out) :: tissue_angles_deg(count(exp_vecs_selection_index))

    ! Local variables
    integer :: i_vec, i_axis, n_selected_axes, out_idx
    real(real64) :: norm_diag, dot_prod, norm_v, cos_phi, angle_rad, norm_factor
    real(real64), parameter :: rad2deg = 180.0_real64 / acos(-1.0_real64)

    ! Compute the norm of the space diagonal (only active axes)
    n_selected_axes = count(axes_selection)
    if (n_selected_axes <= 0) then
      tissue_versatilities = -1.0_real64   ! Error: no axes selected
      return
    end if
    norm_diag = sqrt(real(n_selected_axes, real64))
    ! Precompute normalization factor for tissue versatility
    norm_factor = 1.0_real64 - 1.0_real64 / norm_diag

    ! Loop over selected expression vectors
    ! Note: If the expression vector is zero in all selected axes, tissue versatility (TV) is set to 1 (maximum specificity) and the angle is set to acos(0) = 90 degrees.
    out_idx = 0
    do i_vec = 1, n_vectors
      if (.not. exp_vecs_selection_index(i_vec)) cycle  ! Skip if not selected

      ! Compute dot product and norm for the vector in active axes
      dot_prod = 0.0_real64
      norm_v = 0.0_real64

      do i_axis = 1, n_axes
        if (axes_selection(i_axis)) then
          dot_prod = dot_prod + expression_vectors(i_axis, i_vec) 
          norm_v = norm_v + expression_vectors(i_axis, i_vec)**2
        end if
      end do

      out_idx = out_idx + 1
      ! If the vector is zero, set TV = 1 (maximum specificity)
      if (norm_v <= 0.0_real64) then
        tissue_versatilities(out_idx) = 1.0_real64
        tissue_angles_deg(out_idx) = 90.0_real64
        cycle
      else
        cos_phi = dot_prod / (sqrt(norm_v) * norm_diag)
      end if
      ! Clamp cos_phi for numerical safety
      cos_phi = max(-1.0_real64, min(1.0_real64, cos_phi))
      if (abs(cos_phi - 1.0_real64) < 1e-12_real64) cos_phi = 1.0_real64
      angle_rad = acos(cos_phi)
      tissue_versatilities(out_idx) = (1.0_real64 - cos_phi) / norm_factor
      tissue_angles_deg(out_idx) = angle_rad * rad2deg
      if (abs(tissue_angles_deg(out_idx)) < 1e-12_real64) tissue_angles_deg(out_idx) = 0.0_real64
    end do

  end subroutine compute_tissue_versatility

end