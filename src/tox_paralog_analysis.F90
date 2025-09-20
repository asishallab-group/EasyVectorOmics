module tox_paralog_analysis
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, set_err, is_err, ERR_INVALID_INPUT
    use f42_utils, only: nth_percentile, add_vector, subtract_vector, norm
    implicit none

    integer(int32), parameter :: DOSAGE_PATTERN = 0
    integer(int32), parameter :: SUBFUNC_PATTERN = 1

contains

    pure subroutine detect_patterns(ancestor, paralogs, n_paralogs, n_dims, rdi_threshold, pattern, paralog_angles, sorted_paralog_angles_perm, n_results, max_subset_size, work_arr_paralog_subsets, n_paralog_subsets, global_mask, candidate_mask, temp_paralog_vector, ierr)

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
        integer(int32), intent(out) :: n_results
            !! number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`
        integer(int32), intent(in) :: max_subset_size
            !! maximum subset size of checked paralog subsets. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
        integer(int32), intent(in) :: n_paralog_subsets
            !! number of paralog subsets that can be stored in `work_arr_paralog_subsets`. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
        integer(int32), dimension((n_paralogs + 31) / 32, n_paralog_subsets), intent(out) :: work_arr_paralog_subsets
            !! working array to hold bitmask encoded subsets for detection.
            !! @note
            !! Each bitmask is built of 32 bit chunks. `(n_paralogs + 31) / 32` is equivalent to `ceil(n_paralogs / 32)` and represents the number of chunks
            !! @endnote
        integer(int32), dimension((n_paralogs + 31) / 32), intent(out) :: global_mask
            !! bit mask that will have indices of paralogs kept by pattern set to 1, else 0
        integer(int32), dimension((n_paralogs + 31) / 32), intent(out) :: candidate_mask
            !! working array to hold a subset that is a potential result candidate
        real(real64), dimension(n_dims), intent(out) :: temp_paralog_vector
            !! vector used for pruning subsets
        integer(int32), intent(out) :: ierr
            !! error code

        ! Locals
        real(real64) :: min_norm, max_angle, residual_norm, residual_angle
        integer(int32) :: i_paralog, i_subset_size, i_dim, n_active_masks, n_new_active_masks

        call set_ok(ierr)

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
            block
                real(real64) :: percentile_angle, median_angle
                percentile_angle = nth_percentile(5_int32, paralog_angles, sorted_paralog_angles_perm)
                median_angle = nth_percentile(50_int32, paralog_angles, sorted_paralog_angles_perm)

                ! only paralogs with angles below the gene-family median or lower five percentile are marked active
                do i_paralog = 1, n_paralogs
                    if (paralog_angles(i_paralog) < median_angle .and. paralog_angles(i_paralog) > percentile_angle) then
                        n_active_masks = n_active_masks + 1
                        call mask_set_state(work_arr_paralog_subsets(:, n_active_masks), i_paralog, .true., ierr)
                        call mask_set_state(global_mask, i_paralog, .true., ierr)
                        if (is_err(ierr)) return
                    end if
                end do
            end block
        case (SUBFUNC_PATTERN)
            block
                real(real64) :: median_angle
                median_angle = nth_percentile(50_int32, paralog_angles, sorted_paralog_angles_perm)

                ! only paralogs with angles greater than the gene-family median angle are marked active
                do i_paralog = 1, n_paralogs
                    if (paralog_angles(i_paralog) > median_angle) then
                        n_active_masks = n_active_masks + 1
                        call mask_set_state(work_arr_paralog_subsets(:, n_active_masks), i_paralog, .true., ierr)
                        call mask_set_state(global_mask, i_paralog, .true., ierr)
                        if (is_err(ierr)) return
                    end if
                end do
            end block
        case default
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end select

        n_results = 0_int32
        do i_subset_size = 2, max_subset_size
            n_new_active_masks = 0_int32
            do while (n_active_masks > 0)
                ! take handled mask, always the first one
                candidate_mask = work_arr_paralog_subsets(:, n_results + 1)

                ! move last active mask to first index
                work_arr_paralog_subsets(:, 1) = work_arr_paralog_subsets(:, n_results + n_active_masks)
                ! replace moved mask by last new active mask
                work_arr_paralog_subsets(:, n_results + n_active_masks) = work_arr_paralog_subsets(:, n_results + n_active_masks + n_new_active_masks)
                n_active_masks = n_active_masks - 1

                temp_paralog_vector = ancestor
                do i_paralog = 1, n_paralogs
                    if (mask_check_state(candidate_mask, i_paralog)) then
                        call subtract_vector(temp_paralog_vector, n_dims, paralogs(:, i_paralog))
                    end if
                    if (is_err(ierr)) return
                end do

                ! generate extended subsets by adding successing paralogs of the first active paralog if suitable.
                ! `mask_count_leftmost_zero_bits` returns the index of inactive paralog before first active one, so +2 gives the first potential non-active
                do i_paralog = mask_count_leftmost_zero_bits(candidate_mask) + 2, n_paralogs
                    if (mask_check_state(global_mask, i_paralog) .and. .not. mask_check_state(candidate_mask, i_paralog)) then
                        call mask_set_state(candidate_mask, i_paralog, .true., ierr)
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
                        call mask_set_state(candidate_mask, i_paralog, .false., ierr)
                    end if
                    if (is_err(ierr)) return
                end do
            end do

            n_active_masks = n_new_active_masks
        end do
    end subroutine detect_patterns

    pure subroutine calc_work_arr_paralog_subsets_size(max_subset_size, n_paralogs, work_array_size)
        integer(int32), intent(in) :: max_subset_size, n_paralogs
        integer(int32), intent(out) :: work_array_size

        integer(int32), parameter :: max_int32 = huge(0_int32)
        integer(int32) :: subset_size, extensions_count

        work_array_size = 0
        extensions_count = 1
        do subset_size = 1, max_subset_size

            ! calculate the number of extensions of current subsets

            ! overflow check
            if (extensions_count > max_int32 / (n_paralogs - subset_size + 1)) exit
            extensions_count = extensions_count * (n_paralogs - subset_size + 1)
            extensions_count = extensions_count / subset_size


            ! The current subsets will be replaced by their extensions.
            ! In worst case all extended subsets won't be pruned.
            ! Thus, the extensions count will be the work array size.
            ! As every subset will be extended by one paralog of its successors, the last subset doesn't have extensions, but could be a result.
            ! As this applies to each iteration and the results are stored in the work array, the current subset size holds the number of those possible results.

            ! overflow check
            if (extensions_count > max_int32 - subset_size) exit
            ! if there are less extensions than before, the work array size won't grow anymore
            if (extensions_count + subset_size < work_array_size) exit
            work_array_size = extensions_count + subset_size
        end do
    end subroutine calc_work_arr_paralog_subsets_size

    pure function mask_count_leftmost_zero_bits(bit_mask) result(n_zeros)
        integer(int32), dimension(:), intent(in) :: bit_mask
            !! chunked mask to mark active paralogs
        integer(int32) :: n_zeros
            !! number of leftmost zeros in the bitmask

        integer(int32) :: i_mask_chunk, i_bit

        n_zeros = 0
        i_mask_chunk = size(bit_mask)
        do while (i_mask_chunk > 0)
            if (bit_mask(i_mask_chunk) /= 0) exit
            n_zeros = n_zeros + 32
            i_mask_chunk = i_mask_chunk - 1
        end do

        if (i_mask_chunk > 0) then
            do i_bit = 31, 0, -1
                if (btest(bit_mask(i_mask_chunk), i_bit)) exit
                n_zeros = n_zeros + 1
            end do
        end if
    end function mask_count_leftmost_zero_bits

    pure subroutine mask_set_state(bit_mask, i_paralog, state, ierr)
        integer(int32), dimension(:), intent(out) :: bit_mask
            !! chunked mask to mark active paralogs
        integer(int32), intent(in) :: i_paralog
            !! index of paralog top be marked active
        logical, intent(in) :: state
            !! state the bit should be set to
        integer(int32), intent(out) :: ierr
            !! error code

        integer(int32) :: i_mask_chunk

        call set_ok(ierr)

        i_mask_chunk = (i_paralog - 1) / 32 + 1
        if (i_mask_chunk > size(bit_mask)) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if

        if (state) then
            bit_mask(i_mask_chunk) = ibset(bit_mask(i_mask_chunk), mod(i_paralog - 1, 32))
        else
            bit_mask(i_mask_chunk) = ibclr(bit_mask(i_mask_chunk), mod(i_paralog - 1, 32))
        end if
    end subroutine mask_set_state

    pure function mask_check_state(bit_mask, i_paralog) result(state)
        integer(int32), dimension(:), intent(in) :: bit_mask
            !! chunked mask to mark active paralogs
        integer(int32), intent(in) :: i_paralog
            !! index of paralog top be marked active
        logical :: state
            !! check result

        integer(int32) :: i_mask_chunk

        i_mask_chunk = (i_paralog - 1) / 32 + 1
        if (i_mask_chunk > size(bit_mask)) then
            state = .false.
        else
            state = btest(bit_mask(i_mask_chunk), mod(i_paralog - 1, 32))
        end if

    end function mask_check_state

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