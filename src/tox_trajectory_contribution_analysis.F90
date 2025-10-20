#include "macros.h"

module tox_trajectory_contribution_analysis
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use f42_utils, only: is_close
    use tox_errors, only: set_ok, set_err, is_err, ERR_IDX_OUT_OF_BOUNDS, ERR_EMPTY_INPUT, ERR_DIVISION_BY_ZERO, ERR_INVALID_INPUT, ERR_NAN_INF
    implicit none

    integer(int32), parameter :: MODE_NORMAL = 1
    integer(int32), parameter :: MODE_RAP    = 2
contains

    !> Calculates trajectory contribution using cosine similarity
    !! $$ contribution = \frac{factor \cdot dependent}{\left \Vert factor \right \Vert \cdot \left \Vert dependent \right \Vert} $$
    !! There are two modes:
    !!
    !! 1. **Normal Mode**: meant for normal vectors, it returns the cosine similarity
    !! 2. **RAP Mode**: meant for RAP projected vectors, it focusses only on angular similarity and returns the `arccos` of the cosine similarity
    pure subroutine trajectory_contribution(factor, dependent, n_timepoints, mode, contribution, ierr)
        integer(int32), intent(in) :: n_timepoints
            !! Number of timepoints included in `factor` and `dependent`
        real(real64), dimension(n_timepoints), intent(in) :: factor
            !! Vector of independent variable that contributes to `dependent`
        real(real64), dimension(n_timepoints), intent(in) :: dependent
            !! Vector of dependent variable
        integer(int32), intent(in) :: mode
            !! | Mode   | Value |
            !! |--------|-------|
            !! | Normal | 1     |
            !! | RAP    | 2     |
        real(real64), intent(out) :: contribution
            !! Contribution value
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i_timepoint
        real(real64) :: norm_factor, norm_dependent, magnitude, dot_prod

        call set_ok(ierr)

        dot_prod = 0.0_real64
        norm_factor = 0.0_real64
        norm_dependent = 0.0_real64
        do i_timepoint = 1, n_timepoints
            norm_factor = norm_factor + factor(i_timepoint) ** 2
            norm_dependent = norm_dependent + dependent(i_timepoint) ** 2
            dot_prod = dot_prod + factor(i_timepoint) * dependent(i_timepoint)
        end do

        if (ieee_is_nan(magnitude) .or. ieee_is_nan(dot_prod)) then
            call set_err(ierr, ERR_NAN_INF)
            return
        end if

        magnitude = sqrt(norm_factor) * sqrt(norm_dependent)
        if (is_close(magnitude, 0.0_real64)) then
            call set_err(ierr, ERR_DIVISION_BY_ZERO)
            return
        end if

        select case (mode)
        case (MODE_NORMAL)
            ! -1 <= x <= 1 ; the closer to one, the more directionally similar
            contribution = dot_prod / magnitude
        case (MODE_RAP)
            ! the smaller the angle, the more aligned
            contribution = acos(dot_prod / magnitude)
        case default
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end select
    end subroutine trajectory_contribution

    !> Calculates spike contributions using cosine similarity for specific indices
    !! $$ contributions(i) =
    !!         \frac{factor(i) \cdot dependent(i)}
    !!              {\left \Vert factor \right \Vert \cdot \left \Vert dependent \right \Vert}
    !!         \; | \; 1 \le i \le n\_timepoints $$
    !! There are two modes:
    !!
    !! 1. **Normal Mode**: meant for normal vectors, it returns the cosine similarity
    !! 2. **RAP Mode**: meant for RAP projected vectors, it focusses only on angular similarity and returns the `arccos` of the cosine similarity
    pure subroutine spike_contribution(factor, dependent, n_timepoints, mode, contribution, ierr)
        integer(int32), intent(in) :: n_timepoints
            !! Number of timepoints included in `factor` and `dependent`
        real(real64), dimension(n_timepoints), intent(in) :: factor
            !! Vector of independent variable that contributes to `dependent`
        real(real64), dimension(n_timepoints), intent(in) :: dependent
            !! Vector of independent variable that contributes  `dependent`
        integer(int32), intent(in) :: mode
            !! | Mode   | Value |
            !! |--------|-------|
            !! | Normal | 1     |
            !! | RAP    | 2     |
        real(real64), dimension(n_timepoints), intent(out) :: contribution
            !! Contribution value
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i_timepoint
        real(real64) :: norm_factor, norm_dependent, magnitude

        call set_ok(ierr)

        norm_factor = 0.0_real64
        norm_dependent = 0.0_real64
        do i_timepoint = 1, n_timepoints
            norm_factor = norm_factor + factor(i_timepoint) ** 2
            norm_dependent = norm_dependent + dependent(i_timepoint) ** 2
        end do
        magnitude = sqrt(norm_factor) * sqrt(norm_dependent)

        if (ieee_is_nan(magnitude)) then
            call set_err(ierr, ERR_NAN_INF)
            return
        end if

        if (is_close(magnitude, 0.0_real64)) then
            call set_err(ierr, ERR_DIVISION_BY_ZERO)
            return
        end if

        select case (mode)
        case (MODE_NORMAL)
            do i_timepoint = 1, n_timepoints
                ! -1 <= x <= 1 ; the closer to one, the more directionally similar for specific timepoint
                contribution(i_timepoint) = (factor(i_timepoint) * dependent(i_timepoint)) / magnitude
            end do
        case (MODE_RAP)
            do i_timepoint = 1, n_timepoints
                ! the smaller the angle, the more aligned for specific timepoint
                contribution(i_timepoint) = acos((factor(i_timepoint) * dependent(i_timepoint)) / magnitude)
            end do
        case default
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end select
    end subroutine spike_contribution

    !> Accessor routine to extract sample data for a specific factor and timepoint
    pure subroutine get_vec_across_samples(trajectories, n_factors, n_samples, n_timepoints, i_factor, i_timepoint, result_vector, ierr)
        integer(int32), intent(in) :: n_factors
            !! Number of factors in `trajectories`
        integer(int32), intent(in) :: n_samples
            !! Number of samples in `trajectories`
        integer(int32), intent(in) :: n_timepoints
            !! Number of timepoints in `trajectories`
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
            !! Trajectories of expression data
        integer(int32), intent(in) :: i_factor
            !! Factor index the retrieved samples should be related to
        integer(int32), intent(in) :: i_timepoint
            !! Timepoint index the retrieved samples should be related to
        real(real64), dimension(n_samples), intent(out) :: result_vector
            !! Result vector, equivalent to `trajectories(i_factor, :, i_timepoint)`
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i_sample

        call set_ok(ierr)

        if (n_factors <= 0 .or. n_samples <= 0 .or. n_timepoints <= 0) then
            call set_err(ierr, ERR_EMPTY_INPUT)
            return
        end if

        if (i_factor < 1 .or. i_timepoint < 1 .or. i_factor > n_factors .or. i_timepoint > n_timepoints) then
            call set_err(ierr, ERR_IDX_OUT_OF_BOUNDS)
            return
        end if

        do i_sample = 1, n_samples
            result_vector(i_sample) = trajectories(i_factor, i_sample, i_timepoint)
        end do
    end subroutine get_vec_across_samples

    !> Accessor routine to extract timepoint data for a specific factor and sample
    pure subroutine get_vec_across_timepoints(trajectories, n_factors, n_samples, n_timepoints, i_factor, i_sample, result_vector, ierr)
        integer(int32), intent(in) :: n_factors
            !! Number of factors in `trajectories`
        integer(int32), intent(in) :: n_samples
            !! Number of samples in `trajectories`
        integer(int32), intent(in) :: n_timepoints
            !! Number of timepoints in `trajectories`
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
            !! Trajectories of expression data
        integer(int32), intent(in) :: i_factor
            !! Factor index the retrieved timepoints should be related to
        integer(int32), intent(in) :: i_sample
            !! Timepoint index the retrieved timepoints should be related to
        real(real64), dimension(n_timepoints), intent(out) :: result_vector
            !! Result vector, equivalent to `trajectories(i_factor, i_sample, :)`
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i_timepoint

        call set_ok(ierr)

        if (n_factors <= 0 .or. n_samples <= 0 .or. n_timepoints <= 0) then
            call set_err(ierr, ERR_EMPTY_INPUT)
            return
        end if

        if (i_factor < 1 .or. i_sample < 1 .or. i_factor > n_factors .or. i_sample > n_samples) then
            call set_err(ierr, ERR_IDX_OUT_OF_BOUNDS)
            return
        end if

        do i_timepoint = 1, n_timepoints
            result_vector(i_timepoint) = trajectories(i_factor, i_sample, i_timepoint)
        end do
    end subroutine get_vec_across_timepoints

end module tox_trajectory_contribution_analysis

pure subroutine trajectory_contribution_c(factor, dependent, n_timepoints, mode, contribution, ierr) bind(C, name="trajectory_contribution_c")
    use tox_trajectory_contribution_analysis, only: trajectory_contribution
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    M_USE_NULL_VALIDATION
    integer(c_int), intent(in), target :: n_timepoints
        !! Number of timepoints included in `factor` and `dependent`
    real(c_double), dimension(n_timepoints), intent(in), target :: factor
        !! Vector of independent variable that contributes to `dependent`
    real(c_double), dimension(n_timepoints), intent(in), target :: dependent
        !! Vector of dependent variable
    integer(c_int), intent(in), target :: mode
        !! | Mode   | Value |
        !! |--------|-------|
        !! | Normal | 1     |
        !! | RAP    | 2     |
    real(c_double), intent(out), target :: contribution
        !! Contribution value
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(factor)
    M_CHECK_NON_NULL(dependent)
    M_CHECK_NON_NULL(mode)
    M_CHECK_NON_NULL(contribution)

    call trajectory_contribution(factor, dependent, n_timepoints, mode, contribution, ierr)
end subroutine trajectory_contribution_c

!> Calculates spike contributions using cosine similarity for specific indices
!! $$ contributions(i) =
!!         \frac{factor(i) \cdot dependent(i)}
!!              {\left \Vert factor \right \Vert \cdot \left \Vert dependent \right \Vert}
!!         \; | \; 1 \le i \le n\_timepoints $$
!! There are two modes:
!!
!! 1. **Normal Mode**: meant for normal vectors, it returns the cosine similarity
!! 2. **RAP Mode**: meant for RAP projected vectors, it focusses only on angular similarity and returns the `arccos` of the cosine similarity
pure subroutine spike_contribution_c(factor, dependent, n_timepoints, mode, contribution, ierr) bind(C, name="spike_contribution_c")
    use tox_trajectory_contribution_analysis, only: spike_contribution
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    M_USE_NULL_VALIDATION
    integer(c_int), intent(in), target :: n_timepoints
        !! Number of timepoints included in `factor` and `dependent`
    real(c_double), dimension(n_timepoints), intent(in), target :: factor
        !! Vector of independent variable that contributes to `dependent`
    real(c_double), dimension(n_timepoints), intent(in), target :: dependent
        !! Vector of independent variable that contributes  `dependent`
    integer(c_int), intent(in), target :: mode
        !! | Mode   | Value |
        !! |--------|-------|
        !! | Normal | 1     |
        !! | RAP    | 2     |
    real(c_double), dimension(n_timepoints), intent(out), target :: contribution
        !! Contribution value
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(factor)
    M_CHECK_NON_NULL(dependent)
    M_CHECK_NON_NULL(mode)
    M_CHECK_NON_NULL(contribution)

    call spike_contribution(factor, dependent, n_timepoints, mode, contribution, ierr)
end subroutine spike_contribution_c