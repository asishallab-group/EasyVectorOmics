!> @file f42_loess_smoothing.F90
!> @brief LOESS smoothing module for omics data.
!> @details This module implements the LOESS smoothing algorithm, supporting optional masking for robust smoothing and cross-language compatibility.

module loess_module
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
contains

  !> Finds the indices of the true values in a logical mask.
  !! @param mask Logical array of size n.
  !! @param n Size of the mask.
  !! @param idx_out Integer array to store the indices of true values.
  !! @param m_max Maximum size of idx_out.
  !! @param m_out Actual size of idx_out (number of true values found).
  pure subroutine which(mask, n, idx_out, m_max, m_out)
    logical, intent(in) :: mask(:)
    integer, intent(in) :: n
    integer, intent(out) :: idx_out(:)
    integer, intent(in) :: m_max
    integer, intent(out) :: m_out
    integer :: i, count
    count = 0
    idx_out = 0  ! Initialize to avoid garbage values
    do i = 1, n
      if (mask(i)) then
        count = count + 1
        if (count <= m_max) then
          idx_out(count) = i
        end if
      end if
    end do
    m_out = count
  end subroutine which

  !> Performs LOESS smoothing on a set of data points.
  !! @param n_total Total number of reference points.
  !! @param n_target Number of target points to smooth.
  !! @param d Dimensionality of the data.
  !! @param x_ref Reference x-coordinates.
  !! @param y_ref Reference y-coordinates (d x n_total).
  !! @param indices_used Indices of reference points used for smoothing.
  !! @param x_query Target x-coordinates to smooth.
  !! @param kernel_sigma Bandwidth parameter for the kernel.
  !! @param kernel_cutoff Cutoff for the kernel.
  !! @param y_out Output smoothed values (d x n_target).
  !! @param workspace_weights Temporary array for weights.
  !! @param workspace_values Temporary array for values.
  !! @param mask_in Optional logical mask for explicit exclusion.
  pure subroutine loess_smooth(n_total, n_target, d, x_ref, y_ref, indices_used, x_query, &
                          kernel_sigma, kernel_cutoff, y_out, workspace_weights, &
                          workspace_values, mask_in)
    integer, intent(in) :: n_total, n_target, d
    real(real64), intent(in) :: x_ref(n_total), y_ref(d, n_total), x_query(n_target)
    integer, intent(in) :: indices_used(n_total)
    real(real64), intent(in) :: kernel_sigma, kernel_cutoff
    real(real64), intent(out) :: y_out(d, n_target)
    real(real64), intent(inout) :: workspace_weights(n_total), workspace_values(d, n_total)
    logical, intent(in), optional :: mask_in(n_total)

    integer :: q, i, k, idx, m_out, valid_indices(n_total)
    real(real64) :: query_x, ref_x, delta, sum_weights, weight
    logical :: mask(n_total)

    do q = 1, n_target
      query_x = x_query(q)
      sum_weights = 0.0
      y_out(:, q) = 0.0
      mask = .false.
      do i = 1, n_total
        idx = indices_used(i)
        ref_x = x_ref(idx)
        delta = abs(query_x - ref_x)
        if (present(mask_in)) then
          mask(i) = mask_in(idx) .and. (delta <= kernel_cutoff * kernel_sigma)
        else
          mask(i) = (delta <= kernel_cutoff * kernel_sigma)
        end if
      end do
      call which(mask, n_total, valid_indices, n_total, m_out)
      do i = 1, m_out
        idx = indices_used(valid_indices(i))
        ref_x = x_ref(idx)
        delta = abs(query_x - ref_x)
        weight = exp(-(delta / kernel_sigma)**2)
        sum_weights = sum_weights + weight
        do k = 1, d
          y_out(k, q) = y_out(k, q) + weight * y_ref(k, idx)
        end do
      end do
      if (sum_weights > 0.0) then
        y_out(:, q) = y_out(:, q) / sum_weights
      else
        y_out(:, q) = y_ref(:, indices_used(q))
      end if
    end do
  end subroutine loess_smooth

end module loess_module
