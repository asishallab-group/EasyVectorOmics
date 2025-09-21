module tox_paralog_analysis
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, set_err, is_err, ERR_INVALID_INPUT, ERR_SIZE_MISMATCH
    use f42_utils, only: calc_percentile, add_vector, subtract_vector, norm, PI
    implicit none

    integer(int32), parameter :: DOSAGE_PATTERN = 0
    integer(int32), parameter :: SUBFUNC_PATTERN = 1

    real(real64), parameter :: DOSAGE_MAX_ANGLE_DEG = PI
    real(real64), parameter :: DOSAGE_GAIN_GAMMA = 0.1_real64

contains

    pure subroutine detect_patterns(ancestor, paralogs, n_paralogs, n_dims, rdi_threshold, pattern, filtered_paralogs_mask, n_mask_chunks, paralog_angles, sorted_paralog_angles_perm, n_results, max_subset_size, work_arr_paralog_subsets, n_paralog_subsets, active_mask, temp_paralog_vector, ierr)
        integer(int32), intent(in) :: n_dims
            !! size of `ancestor` vector and vectors in `paralogs`
        integer(int32), intent(in) :: n_paralogs
            !! number of vectors in `paralogs`
        integer(int32), intent(in) :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_paralogs` paralogs. Use subroutine `mask_chunk_count` for calculation
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
        integer(int32), dimension(n_mask_chunks, n_paralog_subsets), intent(out) :: work_arr_paralog_subsets
            !! working array to hold bitmask encoded subsets for detection.
            !! @note
            !! Each bitmask is built of 32 bit chunks. `(n_paralogs + 31) / 32` is equivalent to `ceil(n_paralogs / 32)` and represents the number of chunks
            !! @endnote
        integer(int32), dimension(n_mask_chunks), intent(in) :: filtered_paralogs_mask
            !! bit mask with paralogs' indices kept by pattern set to 1, else 0. Use `filter_paralogs_by_pattern` for its calculation
        integer(int32), dimension(n_mask_chunks), intent(out) :: active_mask
            !! working array to hold the extended subsets
        real(real64), dimension(n_dims), intent(out) :: temp_paralog_vector
            !! vector used for pruning subsets
        integer(int32), intent(out) :: ierr
            !! error code

        ! Locals
        real(real64) :: residual_norm, residual_angle
        integer(int32) :: i_paralog, subset_size, n_active_masks, n_new_active_masks

        call set_ok(ierr)

        work_arr_paralog_subsets = 0_int32
        n_active_masks = 0_int32

        ! initialize first `n_paralogs` subsets of size 1 to be extended
        do i_paralog = 1, n_paralogs
            if (mask_check_state(filtered_paralogs_mask, i_paralog)) then
                n_active_masks = n_active_masks + 1
                call mask_set_state(work_arr_paralog_subsets(:, n_active_masks), i_paralog, .true., ierr)
                if (is_err(ierr)) return
            end if
        end do

        n_results = 0_int32
        do subset_size = 2, max_subset_size
            n_new_active_masks = 0_int32
            do while (n_active_masks > 0)
                call take_active_mask(work_arr_paralog_subsets, n_mask_chunks, n_paralog_subsets, n_results, n_active_masks, n_new_active_masks, active_mask, ierr)
                if (is_err(ierr)) return

                call generate_subsets(active_mask, filtered_paralogs_mask, n_mask_chunks, pattern, ancestor, paralogs, n_paralogs, n_dims, temp_paralog_vector, rdi_threshold, work_arr_paralog_subsets, n_paralog_subsets, n_results, n_active_masks, n_new_active_masks, ierr)
                if (is_err(ierr)) return
            end do

            n_active_masks = n_new_active_masks
        end do
    end subroutine detect_patterns

    pure subroutine generate_subsets(candidate_mask, filtered_paralogs_mask, n_mask_chunks, pattern, ancestor, paralogs, n_paralogs, n_dims, temp_paralog_vector, rdi_threshold, work_arr_paralog_subsets, n_paralog_subsets, n_results, n_active_masks, n_new_active_masks, ierr)
        integer(int32), intent(in) :: n_dims
            !! size of `ancestor` vector and vectors in `paralogs`
        integer(int32), intent(in) :: n_paralogs
            !! number of vectors in `paralogs`
        integer(int32), intent(in) :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_paralogs` paralogs
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
        integer(int32), intent(in) :: n_paralog_subsets
            !! number of paralog subsets that can be stored in `work_arr_paralog_subsets`. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
        integer(int32), dimension(n_mask_chunks, n_paralog_subsets), intent(out) :: work_arr_paralog_subsets
            !! working array to hold bitmask encoded subsets for detection.
        integer(int32), dimension(n_mask_chunks), intent(in) :: filtered_paralogs_mask
            !! bit mask that will have indices of paralogs kept by pattern set to 1, else 0
        integer(int32), dimension(n_mask_chunks), intent(inout) :: candidate_mask
            !! working array to hold a subset that is a potential result candidate
        real(real64), dimension(n_dims), intent(inout) :: temp_paralog_vector
            !! vector used for pruning subsets
        integer(int32), intent(inout) :: n_results
            !! number of results in `work_arr_paralog_subsets`
        integer(int32), intent(in) :: n_active_masks
            !! number of active subsets in `work_arr_paralog_subsets`
        integer(int32), intent(inout) :: n_new_active_masks
            !! number of new active subsets in `work_arr_paralog_subsets`
        integer(int32), intent(out) :: ierr
            !! error code

        integer(int32) :: i_paralog

        select case (pattern)
        case (DOSAGE_PATTERN)
            block
                real(real64) :: subset_angle

                temp_paralog_vector = 0
                do i_paralog = 1, n_paralogs
                    if (mask_check_state(candidate_mask, i_paralog)) then
                        call add_vector(temp_paralog_vector, n_dims, paralogs(:, i_paralog))
                    end if
                    if (is_err(ierr)) return
                end do

                ! generate extended subsets by adding successing paralogs of the first active paralog if suitable.
                do i_paralog = mask_get_first_successor_idx(candidate_mask, n_paralogs), n_paralogs
                    if (mask_check_state(filtered_paralogs_mask, i_paralog) .and. .not. mask_check_state(candidate_mask, i_paralog)) then
                        call mask_set_state(candidate_mask, i_paralog, .true., ierr)
                        if (is_err(ierr)) return

                        call add_vector(temp_paralog_vector, n_dims, paralogs(:, i_paralog))

                        call angle_between(temp_paralog_vector, ancestor, n_dims, subset_angle)
                        if (subset_angle <= DOSAGE_MAX_ANGLE_DEG) then
                            if (norm(temp_paralog_vector, n_dims) >= (1 + DOSAGE_GAIN_GAMMA) * norm(ancestor, n_dims)) then
                                call add_to_results(work_arr_paralog_subsets, n_mask_chunks, n_paralog_subsets, n_results, n_active_masks, n_new_active_masks, candidate_mask, ierr)
                            else
                                call add_new_active_mask(work_arr_paralog_subsets, n_mask_chunks, n_paralog_subsets, n_results, n_active_masks, n_new_active_masks, candidate_mask, ierr)
                            end if
                            if (is_err(ierr)) return
                        end if

                        call subtract_vector(temp_paralog_vector, n_dims, paralogs(:, i_paralog))
                        call mask_set_state(candidate_mask, i_paralog, .false., ierr)
                        if (is_err(ierr)) return
                    end if
                end do
            end block
        case (SUBFUNC_PATTERN)
            block
                real(real64) :: residual_norm

                temp_paralog_vector = ancestor
                do i_paralog = 1, n_paralogs
                    if (mask_check_state(candidate_mask, i_paralog)) then
                        call subtract_vector(temp_paralog_vector, n_dims, paralogs(:, i_paralog))
                    end if
                    if (is_err(ierr)) return
                end do

                ! generate extended subsets by adding successing paralogs of the first active paralog if suitable.
                do i_paralog = mask_get_first_successor_idx(candidate_mask, n_paralogs), n_paralogs
                    if (mask_check_state(filtered_paralogs_mask, i_paralog) .and. .not. mask_check_state(candidate_mask, i_paralog)) then
                        call mask_set_state(candidate_mask, i_paralog, .true., ierr)
                        if (is_err(ierr)) return

                        call subtract_vector(temp_paralog_vector, n_dims, paralogs(:, i_paralog))

                        residual_norm = norm(temp_paralog_vector, n_dims)
                        if (residual_norm < rdi_threshold) then
                            call add_to_results(work_arr_paralog_subsets, n_mask_chunks, n_paralog_subsets, n_results, n_active_masks, n_new_active_masks, candidate_mask, ierr)
                        else ! if (min_norm_among_remaining < residual_norm)
                            call add_new_active_mask(work_arr_paralog_subsets, n_mask_chunks, n_paralog_subsets, n_results, n_active_masks, n_new_active_masks, candidate_mask, ierr)
                        end if
                        if (is_err(ierr)) return

                        call add_vector(temp_paralog_vector, n_dims, paralogs(:, i_paralog))
                        call mask_set_state(candidate_mask, i_paralog, .false., ierr)
                        if (is_err(ierr)) return
                    end if
                end do
            end block
        case default
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end select

    end subroutine generate_subsets

    pure subroutine take_active_mask(subsets, n_mask_chunks, n_subsets, n_results, n_active_masks, n_new_active_masks, active_mask, ierr)
        integer(int32), intent(in) :: n_mask_chunks
            !! number of 32 bit chunks in a mask
        integer(int32), intent(in) :: n_subsets
            !! number of subsets that can be hold
        integer(int32), intent(in) :: n_results
            !! number of results in `subsets`
        integer(int32), intent(inout) :: n_active_masks
            !! number of active masks in `subsets`
        integer(int32), intent(in) :: n_new_active_masks
            !! number of new active masks in `subsets`
        integer(int32), dimension(n_mask_chunks, n_subsets), intent(inout) :: subsets
            !! working array to hold bitmask encoded subsets for detection.
        integer(int32), dimension(n_mask_chunks), intent(out) :: active_mask
            !! taken active mask
        integer(int32), intent(out) :: ierr
            !! error code

        if (n_active_masks <= 0 .or. n_results < 0 .or. n_new_active_masks < 0) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if

        if (n_subsets < n_results + n_active_masks + n_new_active_masks) then
            call set_err(ierr, ERR_SIZE_MISMATCH)
            return
        end if

        ! take handled mask, always the first one after results
        active_mask = subsets(:, n_results + 1)

        ! move last active mask to first index after results
        subsets(:, n_results + 1) = subsets(:, n_results + n_active_masks)
        ! replace moved mask by last new active mask
        subsets(:, n_results + n_active_masks) = subsets(:, n_results + n_active_masks + n_new_active_masks)
        n_active_masks = n_active_masks - 1
    end subroutine take_active_mask

    pure subroutine add_to_results(subsets, n_mask_chunks, n_subsets, n_results, n_active_masks, n_new_active_masks, result, ierr)
        integer(int32), intent(in) :: n_mask_chunks
            !! number of 32 bit chunks in a mask
        integer(int32), intent(in) :: n_subsets
            !! number of subsets that can be hold
        integer(int32), intent(inout) :: n_results
            !! number of results in `subsets`
        integer(int32), intent(in) :: n_active_masks
            !! number of active masks in `subsets`
        integer(int32), intent(in) :: n_new_active_masks
            !! number of new active masks in `subsets`
        integer(int32), dimension(n_mask_chunks, n_subsets), intent(inout) :: subsets
            !! working array to hold bitmask encoded subsets for detection.
        integer(int32), dimension(n_mask_chunks), intent(in) :: result
            !! result to add
        integer(int32), intent(out) :: ierr
            !! error code

        if (n_active_masks < 0 .or. n_results < 0 .or. n_new_active_masks < 0) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if

        if (n_subsets < n_results + n_active_masks + n_new_active_masks + 1) then
            call set_err(ierr, ERR_SIZE_MISMATCH)
            return
        end if

        ! move first new active mask to end
        subsets(:, n_results + n_active_masks + n_new_active_masks + 1) = subsets(:, n_results + n_active_masks + 1)
        ! move first active mask before new active masks
        subsets(:, n_results + n_active_masks + 1) = subsets(:, n_results + 1)
        ! store result
        n_results = n_results + 1
        subsets(:, n_results) = result
    end subroutine add_to_results

    pure subroutine add_new_active_mask(subsets, n_mask_chunks, n_subsets, n_results, n_active_masks, n_new_active_masks, new_active_mask, ierr)
        integer(int32), intent(in) :: n_mask_chunks
            !! number of 32 bit chunks in a mask
        integer(int32), intent(in) :: n_subsets
            !! number of subsets that can be hold
        integer(int32), intent(in) :: n_results
            !! number of results in `subsets`
        integer(int32), intent(in) :: n_active_masks
            !! number of active masks in `subsets`
        integer(int32), intent(inout) :: n_new_active_masks
            !! number of new active masks in `subsets`
        integer(int32), dimension(n_mask_chunks, n_subsets), intent(out) :: subsets
            !! working array to hold bitmask encoded subsets for detection.
        integer(int32), dimension(n_mask_chunks), intent(in) :: new_active_mask
            !! new active mask to add
        integer(int32), intent(out) :: ierr
            !! error code

        if (n_active_masks < 0 .or. n_results < 0 .or. n_new_active_masks < 0) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if

        if (n_subsets < n_results + n_active_masks + n_new_active_masks + 1) then
            call set_err(ierr, ERR_SIZE_MISMATCH)
            return
        end if

        n_new_active_masks = n_new_active_masks + 1
        subsets(:, n_results + n_active_masks + n_new_active_masks) = new_active_mask
    end subroutine add_new_active_mask

    pure subroutine mask_chunk_count(n_paralogs, count)
        integer(int32), intent(in) :: n_paralogs
            !! number of vectors in `paralogs`
        integer(int32), intent(out) :: count
            !! number of 32 bit chunks a mask needs to encode `n_paralogs` paralogs

        !! Each bitmask is built of 32 bit chunks. `(n_paralogs + 31) / 32` is equivalent to `ceil(n_paralogs / 32)` and represents the number of chunks
        count = (n_paralogs + 31) / 32
    end subroutine mask_chunk_count

    pure subroutine filter_paralogs_by_pattern(pattern, pattern_percentile, paralog_angles, sorted_paralog_angles_perm, n_paralogs, mask, n_mask_chunks, ierr)
        integer(int32), intent(in) :: n_paralogs
            !! number of vectors in `paralogs`
        integer(int32), intent(in) :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_paralogs` paralogs
        integer(int32), intent(in) :: pattern
            !! used pattern for detection
            !!
            !! |       Pattern        | Value |
            !! |----------------------|-------|
            !! |    Dosage Effect     |   0   |
            !! | Subfunctionalization |   1   |
            !!
        real(real64), intent(in) :: pattern_percentile
            !! value between 0.0-1.0 defining the percentile that is used for filtering
        real(real64), dimension(n_paralogs), intent(in) :: paralog_angles
            !! vector, holding the ascending sorted angles between ancestor and paralogs. Needed for filtering paralogs
        integer(int32), dimension(n_paralogs), intent(in) :: sorted_paralog_angles_perm
            !! vector, holding the indices for the `paralog_angles` for ascending sorted order
        integer(int32), dimension(n_mask_chunks), intent(out) :: mask
            !! bit mask that will have indices of paralogs kept by pattern set to 1, else 0
        integer(int32), intent(out) :: ierr
            !! error code

        real(real64) :: percentile_value
        integer(int32) :: i_paralog

        call set_ok(ierr)

        mask = 0_int32

        select case (pattern)
        case (DOSAGE_PATTERN)
            call calc_percentile(pattern_percentile, paralog_angles, sorted_paralog_angles_perm, percentile_value, ierr)
            if (is_err(ierr)) return

            ! only paralogs with angles below the gene-family median or lower five percentile are marked active
            do i_paralog = 1, n_paralogs
                if (paralog_angles(i_paralog) <= percentile_value) then
                    call mask_set_state(mask, i_paralog, .true., ierr)
                    if (is_err(ierr)) return
                end if
            end do
        case (SUBFUNC_PATTERN)
            call calc_percentile(pattern_percentile, paralog_angles, sorted_paralog_angles_perm, percentile_value, ierr)
            if (is_err(ierr)) return

            ! only paralogs with angles greater than the gene-family median angle are marked active
            do i_paralog = 1, n_paralogs
                if (paralog_angles(i_paralog) >= percentile_value) then
                    call mask_set_state(mask, i_paralog, .true., ierr)
                    if (is_err(ierr)) return
                end if
            end do
        case default
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end select
    end subroutine filter_paralogs_by_pattern

    pure subroutine calc_work_arr_paralog_subsets_size(max_subset_size, n_paralogs, work_array_size, filtered_paralogs_mask, n_mask_chunks, ierr)
        integer(int32), intent(in) :: n_paralogs
        integer(int32), intent(in) :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_paralogs` paralogs
        integer(int32), intent(inout) :: max_subset_size
        integer(int32), intent(out) :: work_array_size
        integer(int32), dimension(n_mask_chunks), intent(in) :: filtered_paralogs_mask
        integer(int32), intent(out) :: ierr

        integer(int32), parameter :: max_int32 = huge(0_int32)
        integer(int32) :: i_paralog, subset_size, extensions_count, results, previous_results, n_paralogs_filtered
        real(real64) :: percentile_rank

        n_paralogs_filtered = 0
        do i_paralog = 1, n_paralogs
            if (mask_check_state(filtered_paralogs_mask, i_paralog)) then
                n_paralogs_filtered = n_paralogs_filtered + 1
            end if
        end do

        work_array_size = 0
        extensions_count = 1
        previous_results = 0
        results = 0
        do subset_size = 1, max_subset_size
            ! results holds the number of subsets of current subset size that don't have any successing paralogs to be extended with -> result in worst case
            ! previous_results holds the number of results that come from previous subset sizes
            if (previous_results > max_int32 - results) then
                max_subset_size = subset_size - 1
                exit
            end if
            previous_results = previous_results + results
            results = extensions_count - results

            ! calculate the number of extensions of current subsets
            ! overflow check
            if (extensions_count > max_int32 / (n_paralogs - subset_size + 1)) then
                max_subset_size = subset_size - 1
                exit
            end if
            extensions_count = extensions_count * (n_paralogs - subset_size + 1)
            extensions_count = extensions_count / subset_size

            ! The current subsets will be replaced by their extensions.
            ! In worst case all extended subsets won't be pruned.
            ! Thus, the extensions count will be the work array size, plus the subsets from previous iterations that became results.

            ! if there are less extensions than before, the work array size won't grow anymore
            if (extensions_count > max_int32 - previous_results) then
                max_subset_size = subset_size - 1
                exit
            end if
            if (extensions_count + previous_results < work_array_size) exit
            work_array_size = extensions_count + previous_results
        end do
    end subroutine calc_work_arr_paralog_subsets_size

    pure function mask_get_first_successor_idx(bit_mask, n_paralogs) result(idx)
        integer(int32), dimension(:), intent(in) :: bit_mask
            !! chunked mask to mark active paralogs
        integer(int32), intent(in) :: n_paralogs
            !! number of paralogs `bit_mask` is made for
        integer(int32) :: idx
            !! index of last active paralog

        integer(int32) :: i_mask_chunk, i_bit

        idx = size(bit_mask) * 32
        do i_mask_chunk = size(bit_mask), 1, -1
            if (bit_mask(i_mask_chunk) == 0) then
                idx = idx - 32
            else                
                do i_bit = 31, 0, -1
                    if (btest(bit_mask(i_mask_chunk), i_bit)) exit
                    idx = idx - 1
                end do
            end if
        end do
        idx = idx + 1
    end function mask_get_first_successor_idx

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
        theta = dot_product / (sqrt(norm1) * sqrt(norm2))
        theta = max(-1.0_real64, min(1.0_real64, theta))
        angle = acos(theta)
    end subroutine angle_between
end module tox_paralog_analysis