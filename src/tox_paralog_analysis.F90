module tox_paralog_analysis
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, set_err, is_err, ERR_INVALID_INPUT
    use f42_utils, only: nth_percentile, add_vector, subtract_vector, norm
    implicit none

    private
    public :: detect_patterns, angle_between, mask_check_active
    public :: DOSAGE_PATTERN, SUBFUNC_PATTERN

    integer(int32), parameter :: DOSAGE_PATTERN = 0
    integer(int32), parameter :: SUBFUNC_PATTERN = 1

contains

    pure subroutine detect_patterns(ancestor, paralogs, n_paralogs, n_dims, rdi_threshold, pattern, paralog_angles, sorted_paralog_angles_perm, work_arr_paralog_subsets, n_paralog_subsets, temp_paralog_vector, ierr)

        integer(int32), intent(in) :: n_dims
            !! size of `ancestor` vector and vectors in `paralogs`
        integer(int32), intent(in) :: n_paralogs
            !! number of vectors in `paralogs`
        real(real64), dimension(n_dims), intent(in) :: ancestor
            !! expression vector of ancestral ortholog
        real(real64), dimension(n_dims, n_paralogs), intent(in) :: paralogs
            !! expression vectors of paralogs
        real(real64), intent(in) :: rdi_threshold
            !! max allowed residual distance from `ancestor`
        integer(int32), intent(in) :: pattern
            !! used pattern for detection
            !!
            !! |       Pattern        | Value |
            !! |----------------------|-------|
            !! |    Dosage Effect     |   0   |
            !! | Subfunctionalization |   1   |
            !!
        real(real64), dimension(n_paralogs), intent(in) :: paralog_angles
            !! vector, holding the ascending sorted angles between ancestor and paralogs. Needed for filtering paralogs
        integer(int32), dimension(n_paralogs), intent(in) :: sorted_paralog_angles_perm
            !! vector, holding the indices for the `paralog_angles` for ascending sorted order
        integer(int32), intent(in) :: n_paralog_subsets
            !! number of paralog subsets that can be stored in `work_arr_paralog_subsets`
        integer(int32), dimension((n_paralogs + 31) / 32, n_paralog_subsets), intent(out) :: work_arr_paralog_subsets
            !! working array to hold bitmask encoded subsets for detection.
            !! @note
            !! Each bitmask is built of 32 bit chunks. `(n_paralogs + 31) / 32` is equivalent to `ceil(n_paralogs / 32)`
            !! @endnote
        real(real64), dimension(n_dims), intent(out) :: temp_paralog_vector
            !! vector used for pruning subsets
        integer(int32), intent(out) :: ierr
            !! error code

        ! Locals
        integer(int32), dimension((n_paralogs + 31) / 32) :: active_mask, global_mask, candidate_mask
        real(real64) :: percentile_angle, median_angle, min_norm, max_angle, residual_norm, residual_angle
        integer(int32) :: i_paralog, i_subset_size, i_dim, max_subset_size, n_active_masks, n_new_active_masks, n_results

        call set_ok(ierr)

        call detect_max_subset_size(n_paralog_subsets, n_paralogs, max_subset_size)

        if (max_subset_size < 2) return

        max_angle = maxval(paralog_angles)
        min_norm = 0_real64

        do i_paralog = 1, n_paralogs
            min_norm = min(min_norm, norm(paralogs(:, i_paralog), n_dims))
        end do

        n_active_masks = 0_int32
        work_arr_paralog_subsets = 0_int32
        global_mask = 0_int32

        select case (pattern)
        case (DOSAGE_PATTERN)
            percentile_angle = nth_percentile(5_int32, paralog_angles, sorted_paralog_angles_perm)
            median_angle = nth_percentile(50_int32, paralog_angles, sorted_paralog_angles_perm)

            ! only paralogs with angles below the gene-family median or lower five percentile are marked active
            do i_paralog = 1, n_paralogs
                if (paralog_angles(i_paralog) < median_angle .and. paralog_angles(i_paralog) > percentile_angle) then
                    n_active_masks = n_active_masks + 1
                    call mask_mark_active(work_arr_paralog_subsets(:, n_active_masks), i_paralog, ierr)
                    call mask_mark_active(global_mask, i_paralog, ierr)
                    if (is_err(ierr)) return
                end if
            end do
        case (SUBFUNC_PATTERN)
            median_angle = nth_percentile(50_int32, paralog_angles, sorted_paralog_angles_perm)

            ! only paralogs with angles greater than the gene-family median angle are marked active
            do i_paralog = 1, n_paralogs
                if (paralog_angles(i_paralog) > median_angle) then
                    n_active_masks = n_active_masks + 1
                    call mask_mark_active(work_arr_paralog_subsets(:, n_active_masks), i_paralog, ierr)
                    call mask_mark_active(global_mask, i_paralog, ierr)
                    if (is_err(ierr)) return
                end if
            end do
        case default
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end select

        n_results = 0_int32
        do i_subset_size = 2, max_subset_size
            n_new_active_masks = 0_int32
            do while (n_active_masks > 0)
                ! take handled mask, always the first one
                active_mask = work_arr_paralog_subsets(:, n_results + 1)

                ! move last active mask to first index
                work_arr_paralog_subsets(:, 1) = work_arr_paralog_subsets(:, n_results + n_active_masks)
                ! replace moved mask by last new active mask
                work_arr_paralog_subsets(:, n_results + n_active_masks) = work_arr_paralog_subsets(:, n_results + n_active_masks + n_new_active_masks)
                n_active_masks = n_active_masks - 1

                temp_paralog_vector = ancestor
                do i_paralog = 1, n_paralogs
                    if (mask_check_active(active_mask, i_paralog)) then
                        call subtract_vector(temp_paralog_vector, n_dims, paralogs(:, i_paralog))
                    end if
                    if (is_err(ierr)) return
                end do

                ! generate extended subsets by adding successing paralogs of the first active paralog if suitable.
                ! Predecessors already paired with the first active one.
                ! `mask_count_leftmost_zero_bits` returns the index of inactive paralog before first active one, so +2 gives the first potential non-active
                do i_paralog = mask_count_leftmost_zero_bits(active_mask) + 2, n_paralogs
                    if (mask_check_active(global_mask, i_paralog) .and. .not. mask_check_active(active_mask, i_paralog)) then
                        candidate_mask = active_mask
                        call mask_mark_active(candidate_mask, i_paralog, ierr)
                        if (is_err(ierr)) return

                        call subtract_vector(temp_paralog_vector, n_dims, paralogs(:, i_paralog))
                        residual_norm = norm(temp_paralog_vector, n_dims)
                        if (residual_norm < rdi_threshold) then
                            ! move first new active mask to end
                            work_arr_paralog_subsets(:, n_results + n_active_masks + n_new_active_masks + 1) = work_arr_paralog_subsets(:, n_results + n_active_masks + 1)
                            ! move first active mask before new active masks
                            work_arr_paralog_subsets(:, n_results + n_active_masks + 1) = work_arr_paralog_subsets(:, n_results + 1)
                            ! store result
                            n_results = n_results + 1
                            work_arr_paralog_subsets(:, n_results) = candidate_mask
                        else
                            call angle_between(temp_paralog_vector, ancestor, n_dims, residual_angle)
                            if (&
                              (pattern == DOSAGE_PATTERN .and. max_angle > residual_angle)&
                              .or.&
                              (pattern == SUBFUNC_PATTERN .and. min_norm < residual_norm)&
                            ) then
                                n_new_active_masks = n_new_active_masks + 1
                                work_arr_paralog_subsets(:, n_results + n_active_masks + n_new_active_masks) = candidate_mask
                            end if
                        end if

                        call add_vector(temp_paralog_vector, n_dims, paralogs(:, i_paralog))
                    end if
                    if (is_err(ierr)) return
                end do
            end do

            n_active_masks = n_new_active_masks
        end do
    end subroutine detect_patterns

    pure subroutine detect_max_subset_size(work_array_size, n_paralogs, max_subset_size)
        integer(int32), intent(in) :: work_array_size, n_paralogs
        integer(int32), intent(out) :: max_subset_size

        integer(int32), parameter :: max_int32 = huge(0_int32)
        integer(int32) :: total_size, n_to_choose, ierr, combination_count

        max_subset_size = 0
        total_size = 0
        n_to_choose = 0
        combination_count = 1
        do while (total_size < work_array_size)
            max_subset_size = n_to_choose
            n_to_choose = n_to_choose + 1
            ! Overflow check
            if (combination_count > max_int32 / (n_paralogs - n_to_choose + 1)) exit
            combination_count = combination_count * (n_paralogs - n_to_choose + 1)
            combination_count = combination_count / n_to_choose

            ! Overflow check
            if (total_size > max_int32 - combination_count) exit
            total_size = total_size + combination_count
        end do
    end subroutine detect_max_subset_size

    pure function mask_count_leftmost_zero_bits(bit_mask) result(n_zeros)
        integer(int32), dimension(:), intent(in) :: bit_mask
            !! chunked mask to mark active paralogs
        logical :: n_zeros
            !! number of leftmost zeros in the bitmask

        integer(int32) :: i_mask_chunk, i_bit

        n_zeros = 0
        i_mask_chunk = 1
        do while (i_mask_chunk < size(bit_mask))
            n_zeros = n_zeros + 32
            i_mask_chunk = i_mask_chunk + 1
            if (bit_mask(i_mask_chunk) > 0) exit
        end do

        do i_bit = 0, 31
            if (.not. btest(bit_mask(i_mask_chunk), i_bit)) exit
            n_zeros = n_zeros + 1
        end do
    end function mask_count_leftmost_zero_bits

    pure subroutine mask_mark_active(bit_mask, i_paralog, ierr)
        integer(int32), dimension(:), intent(out) :: bit_mask
            !! chunked mask to mark active paralogs
        integer(int32), intent(in) :: i_paralog
            !! index of paralog top be marked active
        integer(int32), intent(out) :: ierr
            !! error code

        integer(int32) :: i_mask_chunk

        call set_ok(ierr)

        i_mask_chunk = (i_paralog - 1) / 32 + 1
        if (i_mask_chunk > size(bit_mask)) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if

        bit_mask(i_mask_chunk) = ibset(bit_mask(i_mask_chunk), mod(i_paralog - 1, 32))
    end subroutine mask_mark_active

    pure function mask_check_active(bit_mask, i_paralog) result(is_active)
        integer(int32), dimension(:), intent(in) :: bit_mask
            !! chunked mask to mark active paralogs
        integer(int32), intent(in) :: i_paralog
            !! index of paralog top be marked active
        logical :: is_active
            !! check result

        integer(int32) :: i_mask_chunk

        i_mask_chunk = (i_paralog - 1) / 32 + 1
        if (i_mask_chunk > size(bit_mask)) then
            is_active = .false.
        else
            is_active = btest(bit_mask(i_mask_chunk), mod(i_paralog - 1, 32))
        end if

    end function mask_check_active

    pure subroutine angle_between(v1, v2, n_dims, angle)
        integer(int32), intent(in) :: n_dims
            !! number of elements in `v1` and `v2`
        real(real64), dimension(n_dims), intent(in) :: v1
            !! first vector for angle calculation
        real(real64), dimension(n_dims), intent(in) :: v2
            !! second vector for angle calculation
        real(real64), intent(out) :: angle
            !! will hold calculated angle

        integer(int32) :: i
        real(real64) :: theta, dot_product, norm1, norm2

        dot_product = 0
        norm1 = 0
        norm2 = 0
        do i = 1, size(v1)
            dot_product = dot_product + v1(i) * v2(i)
            norm1 = norm1 + v1(i) ** 2
            norm2 = norm2 + v2(i) ** 2
        end do
        angle = acos(dot_product / (sqrt(norm1) * sqrt(norm2)))
    end subroutine angle_between


end module tox_paralog_analysis