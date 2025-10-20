module tox_trajectory_contribution_analysis
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, set_err, is_err, ERR_INVALID_INPUT, ERR_EMPTY_INPUT
    implicit none
contains

    !> Accessor routine to extract sample data for a specific factor and timepoint
    pure subroutine get_vec_across_samples(trajectories, n_factors, n_samples, n_timepoints, i_factor, i_timepoint, result_vector, ierr)
        integer(int32), intent(in) :: n_factors
            !! Number of factors in trajectories
        integer(int32), intent(in) :: n_samples
            !! Number of samples in trajectories
        integer(int32), intent(in) :: n_timepoints
            !! Number of timepoints in trajectories
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
            !! Trajectories of expression data
        integer(int32), intent(in) :: i_factor
            !! Factor index the retrieved samples should be related to
        integer(int32), intent(in) :: i_timepoint
            !! Timepoint index the retrieved samples should be related to
        real(real64), dimension(n_samples), intent(out) :: result_vector
            !! Result vector, equivalent to trajectories(i_factor, :, i_timepoint)
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i_sample

        call set_ok(ierr)

        if (n_factors <= 0 .or. n_samples <= 0 .or. n_timepoints <= 0) then
            call set_err(ierr, ERR_EMPTY_INPUT)
            return
        end if

        if (i_factor < 1 .or. i_timepoint < 1 .or. i_factor > n_factors .or. i_timepoint > n_timepoints) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if

        do i_sample = 1, n_samples
            result_vector(i_sample) = trajectories(i_factor, i_sample, i_timepoint)
        end do
    end subroutine get_vec_across_samples

    !> Accessor routine to extract timepoint data for a specific factor and sample
    pure subroutine get_vec_across_timepoints(trajectories, n_factors, n_samples, n_timepoints, i_factor, i_sample, result_vector, ierr)
        integer(int32), intent(in) :: n_factors
            !! Number of factors in trajectories
        integer(int32), intent(in) :: n_samples
            !! Number of samples in trajectories
        integer(int32), intent(in) :: n_timepoints
            !! Number of timepoints in trajectories
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
            !! Trajectories of expression data
        integer(int32), intent(in) :: i_factor
            !! Factor index the retrieved timepoints should be related to
        integer(int32), intent(in) :: i_sample
            !! Timepoint index the retrieved timepoints should be related to
        real(real64), dimension(n_timepoints), intent(out) :: result_vector
            !! Result vector, equivalent to trajectories(i_factor, i_sample, :)
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i_timepoint

        call set_ok(ierr)

        if (n_factors <= 0 .or. n_samples <= 0 .or. n_timepoints <= 0) then
            call set_err(ierr, ERR_EMPTY_INPUT)
            return
        end if

        if (i_factor < 1 .or. i_sample < 1 .or. i_factor > n_factors .or. i_sample > n_samples) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if

        do i_timepoint = 1, n_timepoints
            result_vector(i_timepoint) = trajectories(i_factor, i_sample, i_timepoint)
        end do
    end subroutine get_vec_across_timepoints

end module tox_trajectory_contribution_analysis
