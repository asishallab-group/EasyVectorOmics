#include "macros.h"

module tox_trajectory_contribution_analysis
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, is_err, validate_dimension_size, validate_all_in_range_real, validate_all_in_range_int
    implicit none

contains

    pure subroutine compute_baselines_factor_dependent(factor, dependent, n_dims, mode, factor_baseline, dependent_baseline, ierr)
        integer(int32), intent(in) :: n_dims
        real(real64), dimension(n_dims), intent(in) :: factor
        real(real64), dimension(n_dims), intent(in) :: dependent
        integer(int32), intent(in) :: mode
        real(real64), intent(out) :: factor_baseline
        real(real64), intent(out) :: dependent_baseline
        integer(int32), intent(out) :: ierr

        call set_ok(ierr)
        factor_baseline = 0.0_real64
        dependent_baseline = 0.0_real64
    end subroutine compute_baselines_factor_dependent

    pure subroutine compute_contributions_helper(factor, dependent, n_dims, mode, local_contributions, total_contribution, ierr)
        integer(int32), intent(in) :: n_dims
        real(real64), dimension(n_dims), intent(in) :: factor
        real(real64), dimension(n_dims), intent(in) :: dependent
        integer(int32), intent(in) :: mode
        real(real64), dimension(n_dims), intent(out) :: local_contributions
        real(real64), intent(out) :: total_contribution
        integer(int32), intent(out) :: ierr

        integer(int32) :: i_dim
        real(real64) :: factor_baseline, dependent_baseline

        call set_ok(ierr)

        total_contribution = 0.0_real64
        do i_dim = 1, n_dims
            call compute_baselines_factor_dependent(factor, dependent, n_dims, mode, factor_baseline, dependent_baseline, ierr)
            if (is_err(ierr)) return

            local_contributions(i_dim) = (factor(i_dim) - factor_baseline) * (dependent(i_dim) - dependent_baseline)
            total_contribution = total_contribution + local_contributions(i_dim)
        end do
    end subroutine compute_contributions_helper

    pure subroutine compute_contributions(factor, dependent, n_dims, mode, local_contributions, total_contribution, ierr)
        integer(int32), intent(in) :: n_dims
        real(real64), dimension(n_dims), intent(in) :: factor
        real(real64), dimension(n_dims), intent(in) :: dependent
        integer(int32), intent(in) :: mode
        real(real64), dimension(n_dims), intent(out) :: local_contributions
        real(real64), intent(out) :: total_contribution
        integer(int32), intent(out) :: ierr

        call set_ok(ierr)

        call validate_dimension_size(n_dims, ierr)
        call validate_all_in_range_real(factor, n_dims, ierr)
        call validate_all_in_range_real(dependent, n_dims, ierr)

        if (is_err(ierr)) return

        call compute_contributions_helper(factor, dependent, n_dims, mode, local_contributions, total_contribution, ierr)
    end subroutine compute_contributions

    pure subroutine compute_all_contributions(trajectories, n_factors, n_samples, n_timepoints, factor_indices, n_selected_factors, dependent_indices, n_selected_dependents, mode, local_contributions, total_contributions, temp_factor, temp_dependent, ierr)
        integer(int32), intent(in) :: n_factors
        integer(int32), intent(in) :: n_samples
        integer(int32), intent(in) :: n_timepoints
        integer(int32), intent(in) :: n_selected_factors
        integer(int32), intent(in) :: n_selected_dependents
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
        integer(int32), dimension(n_selected_factors), intent(in) :: factor_indices
        integer(int32), dimension(n_selected_dependents), intent(in) :: dependent_indices
        integer(int32), intent(in) :: mode
        real(real64), dimension(n_samples, n_selected_dependents, n_timepoints), intent(out) :: local_contributions
        real(real64), dimension(n_selected_dependents, n_timepoints), intent(out) :: total_contributions
        real(real64), dimension(n_samples), intent(out) :: temp_factor
        real(real64), dimension(n_samples), intent(out) :: temp_dependent
        integer(int32), intent(out) :: ierr

        integer(int32) :: i_timepoint, i_dependent, i_factor, i_sel_factor, i_sel_dependent, i_sample

        call set_ok(ierr)

        call validate_dimension_size(n_factors, ierr)
        call validate_dimension_size(n_samples, ierr)
        call validate_dimension_size(n_timepoints, ierr)
        call validate_dimension_size(n_selected_factors, ierr)
        call validate_dimension_size(n_selected_dependents, ierr)
        call validate_all_in_range_real(trajectories, n_timepoints * n_samples * n_samples, ierr)
        call validate_all_in_range_int(factor_indices, n_selected_factors, ierr, min=1, max=n_factors)
        call validate_all_in_range_int(dependent_indices, n_selected_dependents, ierr, min=1, max=n_factors)

        if (is_err(ierr)) return

        do i_timepoint = 1, n_timepoints
            do i_sel_dependent = 1, n_selected_dependents
                i_dependent = dependent_indices(i_sel_dependent)
                do i_sample = 1, n_samples
                    temp_dependent = trajectories(i_dependent, i_sample, i_timepoint)
                end do
                do i_sel_factor = 1, n_selected_factors
                    i_factor = factor_indices(i_sel_factor)
                    do i_sample = 1, n_samples
                        temp_factor = trajectories(i_factor, i_sample, i_timepoint)
                    end do

                    call compute_contributions_helper(temp_factor, temp_dependent, n_samples, mode, local_contributions(:, i_dependent, i_timepoint), total_contributions(i_dependent, i_timepoint), ierr)
                    if (is_err(ierr)) return
                end do
            end do
        end do
    end subroutine compute_all_contributions

end module tox_trajectory_contribution_analysis