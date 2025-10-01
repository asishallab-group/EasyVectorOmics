module tox_paralog_analysis
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, set_err, is_err, ERR_INVALID_INPUT, ERR_SIZE_MISMATCH
    use f42_utils, only: add_vector, subtract_vector, norm
    implicit none

    integer(int32), parameter :: DOSAGE_PATTERN = 0
    integer(int32), parameter :: SUBFUNC_PATTERN = 1

contains

    pure subroutine detect_dosage_effect(ancestor, paralogs, n_paralogs, n_dims, filtered_paralogs_mask, n_mask_chunks, n_results, max_subset_size, work_arr_paralog_subsets, n_paralog_subsets, active_mask, temp_paralog_vector, ierr, max_angle, gain_gamma)
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
        integer(int32), intent(out) :: n_results
            !! number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`
        integer(int32), intent(in) :: max_subset_size
            !! maximum subset size of checked paralog subsets. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
        integer(int32), intent(in) :: n_paralog_subsets
            !! number of paralog subsets that can be stored in `work_arr_paralog_subsets`. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
        integer(int32), dimension(n_mask_chunks, n_paralog_subsets), intent(out) :: work_arr_paralog_subsets
            !! working array to hold bitmask encoded subsets for detection.
        integer(int32), dimension(n_mask_chunks), intent(in) :: filtered_paralogs_mask
            !! bit mask with paralogs' indices kept by pattern set to 1, else 0. Use `filter_paralogs_by_pattern` for its calculation
        integer(int32), dimension(n_mask_chunks), intent(out) :: active_mask
            !! working array to hold the extended subsets
        real(real64), dimension(n_dims), intent(out) :: temp_paralog_vector
            !! vector used for pruning subsets
        integer(int32), intent(out) :: ierr
            !! error code
        real(real64), intent(in), optional :: gain_gamma
            !! in dosage mode required true positive magnitude gain for dosage, default 0.1
        real(real64), intent(in), optional :: max_angle
            !! in dosage mode maximum angle in radians that a subset candidate must not exceed, otherwise pruned, default is Pi

        call detect_patterns(ancestor, paralogs, n_paralogs, n_dims, DOSAGE_PATTERN, filtered_paralogs_mask, n_mask_chunks, n_results, max_subset_size, work_arr_paralog_subsets, n_paralog_subsets, active_mask, temp_paralog_vector, dosage_max_angle=max_angle, dosage_gain_gamma=gain_gamma, ierr=ierr)
    end subroutine detect_dosage_effect

    pure subroutine detect_subfunctionalization(ancestor, paralogs, n_paralogs, n_dims, rdi_threshold, filtered_paralogs_mask, n_mask_chunks, n_results, max_subset_size, work_arr_paralog_subsets, n_paralog_subsets, active_mask, temp_paralog_vector, paralog_norms, sorted_paralog_norms_perm, temp_work_array, ierr)
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
        integer(int32), intent(out) :: n_results
            !! number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`
        integer(int32), intent(in) :: max_subset_size
            !! maximum subset size of checked paralog subsets. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
        integer(int32), intent(in) :: n_paralog_subsets
            !! number of paralog subsets that can be stored in `work_arr_paralog_subsets`. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
        integer(int32), dimension(n_mask_chunks, n_paralog_subsets), intent(out) :: work_arr_paralog_subsets
            !! working array to hold bitmask encoded subsets for detection.
        integer(int32), dimension(n_mask_chunks), intent(in) :: filtered_paralogs_mask
            !! bit mask with paralogs' indices kept by pattern set to 1, else 0. Use `filter_paralogs_by_pattern` for its calculation
        integer(int32), dimension(n_mask_chunks), intent(out) :: active_mask
            !! working array to hold the extended subsets
        real(real64), dimension(n_dims), intent(out) :: temp_paralog_vector
            !! vector used for pruning subsets
        integer(int32), intent(out) :: ierr
            !! error code
        real(real64), dimension(n_paralogs), intent(in) :: paralog_norms
            !! in subfunctionalization mode needed for subset pruning, holds the euclidean norms of paralogs (you can use the `norm` function from `f42_utils` function for this)
        integer(int32), dimension(n_paralogs), intent(in) :: sorted_paralog_norms_perm
            !! in subfunctionalization mode needed for subset pruning, as the minimum norm of the paralogs that could extend a subset should not be lower than the subset angle to the ancestor
        real(real64), dimension(n_paralogs), intent(out) :: temp_work_array
            !! in subfunctionalization mode needed for efficient check of minimum value after a certain index

        call detect_patterns(ancestor, paralogs, n_paralogs, n_dims, SUBFUNC_PATTERN, filtered_paralogs_mask, n_mask_chunks, n_results, max_subset_size, work_arr_paralog_subsets, n_paralog_subsets, active_mask, temp_paralog_vector, subfunc_rdi_threshold=rdi_threshold, subfunc_paralog_norms=paralog_norms, subfunc_sorted_paralog_norms_perm=sorted_paralog_norms_perm, subfunc_temp_work_array=temp_work_array, ierr=ierr)
    end subroutine detect_subfunctionalization

    pure subroutine detect_patterns(ancestor, paralogs, n_paralogs, n_dims, pattern, filtered_paralogs_mask, n_mask_chunks, n_results, max_subset_size, work_arr_paralog_subsets, n_paralog_subsets, active_mask, temp_paralog_vector, dosage_max_angle, dosage_gain_gamma, subfunc_rdi_threshold, subfunc_paralog_norms, subfunc_sorted_paralog_norms_perm, subfunc_temp_work_array, ierr)
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
        integer(int32), intent(in) :: pattern
            !! used pattern for detection
            !!
            !! |       Pattern        | Value |
            !! |----------------------|-------|
            !! |    Dosage Effect     |   0   |
            !! | Subfunctionalization |   1   |
            !!
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
        real(real64), intent(in), optional :: dosage_gain_gamma
            !! in dosage mode required true positive magnitude gain for dosage, default 0.1
        real(real64), intent(in), optional :: dosage_max_angle
            !! in dosage mode maximum angle in radians that a subset candidate must not exceed, otherwise pruned, default is Pi
        real(real64), dimension(n_paralogs), intent(in), optional :: subfunc_paralog_norms
            !! in subfunctionalization mode needed for subset pruning, holds the euclidean norms of paralogs (you can use the `norm` from `f42_utils` function for this)
        integer(int32), dimension(n_paralogs), intent(in), optional :: subfunc_sorted_paralog_norms_perm
            !! in subfunctionalization mode needed for subset pruning, as the minimum norm of the paralogs that could extend a subset should not be lower than the subset angle to the ancestor
        real(real64), dimension(n_paralogs), intent(out), optional :: subfunc_temp_work_array
            !! in subfunctionalization mode needed for efficient check of minimum value after a certain index
        real(real64), intent(in), optional :: subfunc_rdi_threshold
            !! max allowed residual distance from `ancestor`

        ! Locals
        integer(int32) :: i_paralog, subset_size, n_active_masks, n_new_active_masks

        call set_ok(ierr)

        work_arr_paralog_subsets = 0_int32
        n_active_masks = 0_int32

        ! initialize first `n_paralogs - 1` subsets of size 1 to be extended
        ! -1 because the subset with last paralog set cannot be extended, as it doesn't have successors
        do i_paralog = 1, n_paralogs - 1
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

                call generate_subsets(active_mask, filtered_paralogs_mask, n_mask_chunks, pattern, ancestor, paralogs, n_paralogs, n_dims, temp_paralog_vector, work_arr_paralog_subsets, n_paralog_subsets, n_results, n_active_masks, n_new_active_masks, dosage_max_angle, dosage_gain_gamma, subfunc_rdi_threshold, subfunc_paralog_norms, subfunc_sorted_paralog_norms_perm, subfunc_temp_work_array, ierr)
                if (is_err(ierr)) return
            end do

            n_active_masks = n_new_active_masks
        end do
    end subroutine detect_patterns

    pure subroutine generate_subsets(candidate_mask, filtered_paralogs_mask, n_mask_chunks, pattern, ancestor, paralogs, n_paralogs, n_dims, temp_paralog_vector, work_arr_paralog_subsets, n_paralog_subsets, n_results, n_active_masks, n_new_active_masks, dosage_max_angle, dosage_gain_gamma, subfunc_rdi_threshold, subfunc_paralog_norms, subfunc_sorted_paralog_norms_perm, subfunc_temp_work_array, ierr)
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
        real(real64), intent(in), optional :: subfunc_rdi_threshold
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
        real(real64), intent(in), optional :: dosage_gain_gamma
            !! in dosage mode required true positive magnitude gain for dosage, default 0.1
        real(real64), intent(in), optional :: dosage_max_angle
            !! in dosage mode maximum angle in radians that a subset candidate must not exceed, otherwise pruned, default is Pi
        real(real64), dimension(n_paralogs), intent(in), optional :: subfunc_paralog_norms
            !! in subfunctionalization mode needed for subset pruning, holds the euclidean norms of paralogs (you can use the `norm` from `f42_utils` function for this)
        integer(int32), dimension(n_paralogs), intent(in), optional :: subfunc_sorted_paralog_norms_perm
            !! in subfunctionalization mode needed for subset pruning, as the minimum norm of the paralogs that could extend a subset should not be lower than the subset angle to the ancestor
        real(real64), dimension(n_paralogs), intent(out), optional :: subfunc_temp_work_array
            !! in subfunctionalization mode needed for efficient check of minimum value after a certain index

        integer(int32) :: i_paralog

        select case (pattern)
        case (DOSAGE_PATTERN)
            block
                use f42_utils, only: PI
                real(real64) :: subset_angle, gain, max_angle

                gain = merge(dosage_gain_gamma, 0.1_real64, present(dosage_gain_gamma))
                max_angle = merge(dosage_max_angle, PI, present(dosage_max_angle))

                !! prepare residual, so the extending paralog just needs to be included in one operation and excluded after calculation
                temp_paralog_vector = 0
                do i_paralog = 1, n_paralogs
                    if (mask_check_state(candidate_mask, i_paralog)) then
                        call add_vector(temp_paralog_vector, paralogs(:, i_paralog))
                    end if
                    if (is_err(ierr)) return
                end do

                ! generate extended subsets by adding successing paralogs of the first active paralog if suitable.
                do i_paralog = mask_get_first_successor_idx(candidate_mask), n_paralogs
                    if (mask_check_state(filtered_paralogs_mask, i_paralog) .and. .not. mask_check_state(candidate_mask, i_paralog)) then
                        call mask_set_state(candidate_mask, i_paralog, .true., ierr)
                        if (is_err(ierr)) return

                        call add_vector(temp_paralog_vector, paralogs(:, i_paralog))

                        call angle_between(temp_paralog_vector, ancestor, n_dims, subset_angle)
                        if (subset_angle <= max_angle) then
                            if (norm(temp_paralog_vector) >= (1 + gain) * norm(ancestor)) then
                                call add_to_results(work_arr_paralog_subsets, n_mask_chunks, n_paralog_subsets, n_results, n_active_masks, n_new_active_masks, candidate_mask, ierr)
                            else
                                call add_new_active_mask(work_arr_paralog_subsets, n_mask_chunks, n_paralog_subsets, n_results, n_active_masks, n_new_active_masks, candidate_mask, ierr)
                            end if
                            if (is_err(ierr)) return
                        end if

                        call subtract_vector(temp_paralog_vector, paralogs(:, i_paralog))
                        call mask_set_state(candidate_mask, i_paralog, .false., ierr)
                        if (is_err(ierr)) return
                    end if
                end do
            end block
        case (SUBFUNC_PATTERN)
            block
                real(real64) :: residual_norm

                if (.not. (present(subfunc_paralog_norms) .and. present(subfunc_sorted_paralog_norms_perm) .and. present(subfunc_temp_work_array) .and. present(subfunc_rdi_threshold))) then
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
                end if

                !! initialize work array with min values, so each index i holds the min value in subarray subfunc_paralog_norms(i:n_paralogs)
                call fill_array_with_minvals_for_each_idx(subfunc_temp_work_array, subfunc_paralog_norms, subfunc_sorted_paralog_norms_perm, n_paralogs, ierr)
                if (is_err(ierr)) return

                !! also, prepare residual, so the extending paralog just needs to be included in one operation and excluded after calculation
                temp_paralog_vector = ancestor
                do i_paralog = 1, n_paralogs
                    !! residual
                    if (mask_check_state(candidate_mask, i_paralog)) then
                        call subtract_vector(temp_paralog_vector, paralogs(:, i_paralog))
                    end if
                    if (is_err(ierr)) return

                end do

                ! generate extended subsets by adding successing paralogs of the first active paralog if suitable.
                do i_paralog = mask_get_first_successor_idx(candidate_mask), n_paralogs
                    if (mask_check_state(filtered_paralogs_mask, i_paralog) .and. .not. mask_check_state(candidate_mask, i_paralog)) then
                        call mask_set_state(candidate_mask, i_paralog, .true., ierr)
                        if (is_err(ierr)) return

                        call subtract_vector(temp_paralog_vector, paralogs(:, i_paralog))

                        residual_norm = norm(temp_paralog_vector)
                        if (residual_norm <= subfunc_rdi_threshold) then
                            call add_to_results(work_arr_paralog_subsets, n_mask_chunks, n_paralog_subsets, n_results, n_active_masks, n_new_active_masks, candidate_mask, ierr)
                        else if (i_paralog < n_paralogs) then
                            ! subfunc_temp_work_array(i_paralog+1) is min(nor) for i in i_paralog+1:n_paralogs )
                            ! so if the minimum norm of the remaining paralogs is not lower the residual, prune this subset branch
                            if (subfunc_temp_work_array(i_paralog + 1) <= residual_norm) then
                                call add_new_active_mask(work_arr_paralog_subsets, n_mask_chunks, n_paralog_subsets, n_results, n_active_masks, n_new_active_masks, candidate_mask, ierr)
                            end if
                        end if
                        if (is_err(ierr)) return

                        call add_vector(temp_paralog_vector, paralogs(:, i_paralog))
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

    !> Helper for subfunctionalization pruning. It initializes a working array with each index i holding the min value in subarray src_arr(i:src_arr_len)
    pure subroutine fill_array_with_minvals_for_each_idx(out_arr, src_arr, sorted_src_arr_perm, src_arr_len, ierr)
        integer(int32), intent(in) :: src_arr_len
            !! number elements in `src_arr`
        real(real64), dimension(src_arr_len), intent(in), optional :: src_arr
            !! array that holds the original values
        integer(int32), dimension(src_arr_len), intent(in), optional :: sorted_src_arr_perm
            !! sorted permutation vector, so each index holds the index of the value in `src_arr` as if `src_arr` was soted ascending
        real(real64), dimension(src_arr_len), intent(out), optional :: out_arr
            !! output array, e.g. source [50, 75, 0, 100, 25] would lead to `out_arr` [0, 0, 0, 25, 25]
        integer(int32), intent(out) :: ierr
            !! error code

        integer(int32) :: i_out, i_perm, last_min_index, current_min_idx

        call set_ok(ierr)

        last_min_index = 0
        do i_perm = 1, src_arr_len
            current_min_idx = sorted_src_arr_perm(i_perm)
            if (current_min_idx > src_arr_len) then
                call set_err(ierr, ERR_INVALID_INPUT)
                exit
            end if
            if (current_min_idx > last_min_index) then
                do i_out = last_min_index + 1, current_min_idx
                    out_arr(i_out) = src_arr(current_min_idx)
                end do

                if (current_min_idx == src_arr_len) exit
                last_min_index = current_min_idx
            end if
        end do
    end subroutine fill_array_with_minvals_for_each_idx

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

        call set_ok(ierr)

        if (n_active_masks <= 0 .or. n_results < 0 .or. n_new_active_masks < 0) then
            call set_err(ierr, ERR_INVALID_INPUT)
            return
        end if

        if (n_subsets < n_results + n_active_masks + n_new_active_masks) then
            call set_err(ierr, ERR_SIZE_MISMATCH)
            return
        end if

        ! take handled mask, always the first one after results
        active_mask = subsets(:, n_results + n_active_masks)

        ! replace taken mask by last new active mask
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

        call set_ok(ierr)

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

        call set_ok(ierr)

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

    pure subroutine filter_paralogs_by_pattern_subfunctionalization(paralog_angles, threshold, n_paralogs, mask, n_mask_chunks, ierr)
        integer(int32), intent(in) :: n_paralogs
            !! number of vectors in `paralogs`
        integer(int32), intent(in) :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_paralogs` paralogs
        real(real64), dimension(n_paralogs), intent(in) :: paralog_angles
            !! vector, holding the ascending sorted angles between ancestor and paralogs. Needed for filtering paralogs
        real(real64), intent(in) :: threshold
            !! filter threshold
        integer(int32), dimension(n_mask_chunks), intent(out) :: mask
            !! bit mask that will have indices of paralogs kept by pattern set to 1, else 0
        integer(int32), intent(out) :: ierr
            !! error code

        call filter_paralogs_by_pattern(SUBFUNC_PATTERN, paralog_angles, threshold, n_paralogs, mask, n_mask_chunks, ierr)
    end subroutine filter_paralogs_by_pattern_subfunctionalization

    pure subroutine filter_paralogs_by_pattern_dosage_effect(paralog_angles, threshold, n_paralogs, mask, n_mask_chunks, ierr)
        integer(int32), intent(in) :: n_paralogs
            !! number of vectors in `paralogs`
        integer(int32), intent(in) :: n_mask_chunks
            !! number of 32 bit chunks a mask needs to encode `n_paralogs` paralogs
        real(real64), dimension(n_paralogs), intent(in) :: paralog_angles
            !! vector, holding the ascending sorted angles between ancestor and paralogs. Needed for filtering paralogs
        real(real64), intent(in) :: threshold
            !! filter threshold
        integer(int32), dimension(n_mask_chunks), intent(out) :: mask
            !! bit mask that will have indices of paralogs kept by pattern set to 1, else 0
        integer(int32), intent(out) :: ierr
            !! error code

        call filter_paralogs_by_pattern(DOSAGE_PATTERN, paralog_angles, threshold, n_paralogs, mask, n_mask_chunks, ierr)
    end subroutine filter_paralogs_by_pattern_dosage_effect

    pure subroutine filter_paralogs_by_pattern(pattern, paralog_angles, threshold, n_paralogs, mask, n_mask_chunks, ierr)
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
        real(real64), dimension(n_paralogs), intent(in) :: paralog_angles
            !! vector, holding the ascending sorted angles between ancestor and paralogs. Needed for filtering paralogs
        real(real64), intent(in) :: threshold
            !! filter threshold
        integer(int32), dimension(n_mask_chunks), intent(out) :: mask
            !! bit mask that will have indices of paralogs kept by pattern set to 1, else 0
        integer(int32), intent(out) :: ierr
            !! error code

        integer(int32) :: i_paralog

        call set_ok(ierr)

        mask = 0_int32

        select case (pattern)
        case (DOSAGE_PATTERN)
            ! only paralogs with angles below the gene-family median or lower five percentile are marked active
            do i_paralog = 1, n_paralogs
                if (paralog_angles(i_paralog) <= threshold) then
                    call mask_set_state(mask, i_paralog, .true., ierr)
                    if (is_err(ierr)) return
                end if
            end do
        case (SUBFUNC_PATTERN)
            ! only paralogs with angles greater than the gene-family median angle are marked active
            do i_paralog = 1, n_paralogs
                if (paralog_angles(i_paralog) >= threshold) then
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
        integer(int32) :: i_paralog, subset_size, extensions_count, results, previous_results, n_paralogsiltered

        n_paralogsiltered = 0
        do i_paralog = 1, n_paralogs
            if (mask_check_state(filtered_paralogs_mask, i_paralog)) then
                n_paralogsiltered = n_paralogsiltered + 1
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

        ! all subsets with last paralog enabled are counted as a result.
        ! as the subset of size 1 with last paralog is not a valid subset, remove it (can not be extended, thus also not part of initialization)
        work_array_size = work_array_size - 1
    end subroutine calc_work_arr_paralog_subsets_size

    pure function mask_get_first_successor_idx(bit_mask) result(idx)
        integer(int32), dimension(:), intent(in) :: bit_mask
            !! chunked mask to mark active paralogs
        integer(int32) :: idx
            !! index of last active paralog

        integer(int32) :: i_mask_chunk

        idx = size(bit_mask) * 32
        do i_mask_chunk = size(bit_mask), 1, -1
            idx = idx - leadz(bit_mask(i_mask_chunk))
            if (mod(idx, 32) /= 0) exit
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

subroutine angle_between_c(v1, v2, n_dims, angle) bind(C, name="angle_between_c")
    use tox_paralog_analysis, only: angle_between
    use, intrinsic :: iso_c_binding, only: c_double, c_int, c_ptr, c_f_pointer
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none

    integer(c_int), intent(in), value :: n_dims
        !! second vector for angle calculation
    type(c_ptr), intent(in) :: v1
        !! number of elements in `v1` and `v2`
    type(c_ptr), intent(in) :: v2
        !! first vector for angle calculation
    real(c_double), intent(out) :: angle
        !! will hold calculated angle

    real(real64), pointer :: v1_f(:), v2_f(:)

    call c_f_pointer(v1, v1_f, [n_dims])
    call c_f_pointer(v2, v2_f, [n_dims])

    call angle_between(v1_f, v2_f, n_dims, angle)
end subroutine angle_between_c

subroutine angle_between_r(v1, v2, n_dims, angle)
    use tox_paralog_analysis, only: angle_between
    use, intrinsic :: iso_fortran_env, only: real64, int32
    implicit none

    integer(int32), intent(in) :: n_dims
        !! number of elements in `v1` and `v2`
    real(real64), dimension(n_dims), intent(in) :: v1
        !! first vector for angle calculation
    real(real64), dimension(n_dims), intent(in) :: v2
        !! second vector for angle calculation
    real(real64), intent(out) :: angle
        !! will hold calculated angle

    call angle_between(v1, v2, n_dims, angle)
end subroutine angle_between_r

subroutine detect_dosage_effect_c(ancestor, paralogs, n_paralogs, n_dims, filtered_paralogs_mask, n_mask_chunks, &
                                n_results, max_subset_size, work_arr_paralog_subsets, n_paralog_subsets, &
                                active_mask, temp_paralog_vector, ierr, max_angle, gain_gamma) bind(C, name="detect_dosage_effect_c")

    use tox_paralog_analysis, only: detect_dosage_effect
    use, intrinsic :: iso_c_binding, only: c_double, c_int, c_ptr, c_f_pointer
    use, intrinsic :: iso_fortran_env, only: real64, int32
    implicit none

    ! Inputs from C
    integer(c_int), intent(in), value :: n_dims
        !! size of `ancestor` vector and vectors in `paralogs`
    integer(c_int), intent(in), value :: n_paralogs
        !! number of vectors in `paralogs`
    integer(c_int), intent(in), value :: n_mask_chunks
        !! number of 32 bit chunks a mask needs to encode `n_paralogs` paralogs. Use subroutine `mask_chunk_count` for calculation
    type(c_ptr), intent(in) :: ancestor
        !! expression vector of ancestral ortholog
    type(c_ptr), intent(in) :: paralogs
        !! expression vectors of paralogs
    integer(c_int), intent(in), value :: max_subset_size
        !! maximum subset size of checked paralog subsets. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
    integer(c_int), intent(in), value :: n_paralog_subsets
        !! number of paralog subsets that can be stored in `work_arr_paralog_subsets`. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
    type(c_ptr), intent(in) :: filtered_paralogs_mask
        !! bit mask with paralogs' indices kept by pattern set to 1, else 0. Use `filter_paralogs_by_pattern` for its calculation
    real(c_double), intent(in), value :: gain_gamma
        !! in dosage mode required true positive magnitude gain for dosage, default 0.1
    real(c_double), intent(in), value :: max_angle
        !! in dosage mode maximum angle in radians that a subset candidate must not exceed, otherwise pruned, default is Pi

    ! Outputs to C
    integer(c_int), intent(out) :: n_results
        !! number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`
    type(c_ptr) :: work_arr_paralog_subsets
        !! working array to hold bitmask encoded subsets for detection.
    type(c_ptr) :: active_mask
        !! working array to hold the extended subsets
    type(c_ptr) :: temp_paralog_vector
        !! vector used for pruning subsets
    integer(c_int), intent(out) :: ierr
        !! error code

    real(real64), dimension(:), pointer :: ancestor_f, temp_paralog_vector_f
    real(real64), dimension(:,:), pointer :: paralogs_f
    integer(int32), dimension(:), pointer :: filtered_paralogs_mask_f, active_mask_f
    integer(int32), dimension(:,:), pointer :: work_arr_paralog_subsets_f

    ! Allocate shaped arrays using converted dimensions
    call c_f_pointer(ancestor, ancestor_f, [n_dims])
    call c_f_pointer(temp_paralog_vector, temp_paralog_vector_f, [n_dims])
    call c_f_pointer(paralogs, paralogs_f, [n_dims, n_paralogs])
    call c_f_pointer(filtered_paralogs_mask, filtered_paralogs_mask_f, [n_mask_chunks])
    call c_f_pointer(active_mask, active_mask_f, [n_mask_chunks])
    call c_f_pointer(work_arr_paralog_subsets, work_arr_paralog_subsets_f, [n_mask_chunks, n_paralog_subsets])

    ! Call original routine
    call detect_dosage_effect(ancestor_f, paralogs_f, n_paralogs, n_dims, filtered_paralogs_mask_f, n_mask_chunks, &
                              n_results, max_subset_size, work_arr_paralog_subsets_f, n_paralog_subsets, &
                              active_mask_f, temp_paralog_vector_f, ierr, max_angle=max_angle, gain_gamma=gain_gamma)
end subroutine detect_dosage_effect_c

subroutine detect_dosage_effect_r(ancestor, paralogs, n_paralogs, n_dims, filtered_paralogs_mask, n_mask_chunks, &
                                n_results, max_subset_size, work_arr_paralog_subsets, n_paralog_subsets, &
                                active_mask, temp_paralog_vector, ierr, max_angle, gain_gamma)
    use tox_paralog_analysis, only: detect_dosage_effect
    use, intrinsic :: iso_fortran_env, only: real64, int32
    implicit none

    ! Inputs from R
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
    integer(int32), intent(in) :: max_subset_size
        !! maximum subset size of checked paralog subsets. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
    integer(int32), intent(in) :: n_paralog_subsets
        !! number of paralog subsets that can be stored in `work_arr_paralog_subsets`. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
    integer(int32), dimension(n_mask_chunks), intent(in) :: filtered_paralogs_mask
        !! bit mask with paralogs' indices kept by pattern set to 1, else 0. Use `filter_paralogs_by_pattern` for its calculation
    real(real64), intent(in) :: gain_gamma
        !! in dosage mode required true positive magnitude gain for dosage, default 0.1
    real(real64), intent(in) :: max_angle
        !! in dosage mode maximum angle in radians that a subset candidate must not exceed, otherwise pruned, default is Pi

    ! Outputs to R
    integer(int32), intent(out) :: n_results
        !! number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`
    integer(int32), dimension(n_mask_chunks, n_paralog_subsets), intent(out) :: work_arr_paralog_subsets
        !! working array to hold bitmask encoded subsets for detection.
    integer(int32), dimension(n_mask_chunks), intent(out) :: active_mask
        !! working array to hold the extended subsets
    real(real64), dimension(n_dims), intent(out) :: temp_paralog_vector
        !! vector used for pruning subsets
    integer(int32), intent(out) :: ierr
        !! error code

    call detect_dosage_effect(ancestor, paralogs, n_paralogs, n_dims, filtered_paralogs_mask, n_mask_chunks, &
                              n_results, max_subset_size, work_arr_paralog_subsets, n_paralog_subsets, &
                              active_mask, temp_paralog_vector, ierr, max_angle=max_angle, gain_gamma=gain_gamma)
end subroutine detect_dosage_effect_r

subroutine detect_subfunctionalization_c(ancestor, paralogs, n_paralogs, n_dims, rdi_threshold, filtered_paralogs_mask, &
                                       n_mask_chunks, n_results, max_subset_size, work_arr_paralog_subsets, n_paralog_subsets, &
                                       active_mask, temp_paralog_vector, paralog_norms, sorted_paralog_norms_perm, temp_work_array, ierr) &
                                       bind(C, name="detect_subfunctionalization_c")
    use tox_paralog_analysis, only: detect_subfunctionalization
    use, intrinsic :: iso_c_binding, only: c_double, c_int, c_ptr, c_f_pointer
    use, intrinsic :: iso_fortran_env, only: real64, int32
    implicit none

    ! Inputs from C
    integer(c_int), intent(in), value :: n_paralogs
        !! number of vectors in `paralogs`
    integer(c_int), intent(in), value :: n_dims
        !! size of `ancestor` vector and vectors in `paralogs`
    type(c_ptr), intent(in) :: ancestor
        !! expression vector of ancestral ortholog
    type(c_ptr), intent(in) :: paralogs
        !! expression vectors of paralogs, column-major layout
    real(c_double), intent(in), value :: rdi_threshold
        !! max allowed residual distance from `ancestor`
    integer(c_int), intent(in), value :: n_mask_chunks
        !! number of 32 bit chunks a mask needs to encode `n_paralogs` paralogs
    type(c_ptr), intent(in) :: filtered_paralogs_mask
        !! bit mask with paralogs' indices kept by pattern set to 1, else 0
    integer(c_int), intent(in), value :: max_subset_size
        !! maximum subset size of checked paralog subsets
    integer(c_int), intent(in), value :: n_paralog_subsets
        !! number of paralog subsets that can be stored in `work_arr_paralog_subsets`
    type(c_ptr), intent(in) :: paralog_norms
        !! euclidean norms of paralogs
    type(c_ptr), intent(in) :: sorted_paralog_norms_perm
        !! permutation of paralogs sorted by norm

    ! Outputs to C
    integer(c_int), intent(out) :: n_results
        !! number of resulting subsets
    type(c_ptr) :: work_arr_paralog_subsets
        !! working array to hold bitmask encoded subsets
    type(c_ptr) :: active_mask
        !! working array to hold the extended subsets
    type(c_ptr) :: temp_paralog_vector
        !! vector used for pruning subsets
    type(c_ptr) :: temp_work_array
        !! used for efficient check of minimum value after a certain index
    integer(c_int), intent(out) :: ierr
        !! error code

    real(real64), pointer :: ancestor_f(:), temp_paralog_vector_f(:)
    real(real64), pointer :: paralogs_f(:,:)
    integer(int32), pointer :: filtered_paralogs_mask_f(:), active_mask_f(:)
    integer(int32), pointer :: work_arr_paralog_subsets_f(:,:)
    real(real64), pointer :: paralog_norms_f(:), temp_work_array_f(:)
    integer(int32), pointer :: sorted_paralog_norms_perm_f(:)

    ! Allocate shaped arrays using converted dimensions
    call c_f_pointer(ancestor, ancestor_f, [n_dims])
    call c_f_pointer(temp_paralog_vector, temp_paralog_vector_f, [n_dims])
    call c_f_pointer(paralogs, paralogs_f, [n_dims, n_paralogs])
    call c_f_pointer(filtered_paralogs_mask, filtered_paralogs_mask_f, [n_mask_chunks])
    call c_f_pointer(active_mask, active_mask_f, [n_mask_chunks])
    call c_f_pointer(work_arr_paralog_subsets, work_arr_paralog_subsets_f, [n_mask_chunks, n_paralog_subsets])
    call c_f_pointer(paralog_norms, paralog_norms_f, [n_paralogs])
    call c_f_pointer(temp_work_array, temp_work_array_f, [n_paralogs])
    call c_f_pointer(sorted_paralog_norms_perm, sorted_paralog_norms_perm_f, [n_paralogs])

    ! Call original routine
    call detect_subfunctionalization(ancestor_f, paralogs_f, n_paralogs, n_dims, rdi_threshold, filtered_paralogs_mask_f, &
                                     n_mask_chunks, n_results, max_subset_size, work_arr_paralog_subsets_f, n_paralog_subsets, &
                                     active_mask_f, temp_paralog_vector_f, paralog_norms_f, sorted_paralog_norms_perm_f, temp_work_array_f, ierr)
end subroutine detect_subfunctionalization_c

subroutine detect_subfunctionalization_r(ancestor, paralogs, n_paralogs, n_dims, rdi_threshold, filtered_paralogs_mask, &
                                         n_mask_chunks, n_results, max_subset_size, work_arr_paralog_subsets, n_paralog_subsets, &
                                         active_mask, temp_paralog_vector, paralog_norms, sorted_paralog_norms_perm, temp_work_array, ierr)

    use tox_paralog_analysis, only: detect_subfunctionalization
    use, intrinsic :: iso_fortran_env, only: real64, int32
    implicit none

    ! Inputs from R
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
    integer(int32), dimension(n_mask_chunks), intent(in) :: filtered_paralogs_mask
      !! bit mask with paralogs' indices kept by pattern set to 1, else 0. Use `filter_paralogs_by_pattern` for its calculation
    integer(int32), intent(in) :: max_subset_size
      !! maximum subset size of checked paralog subsets. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
    integer(int32), intent(in) :: n_paralog_subsets
      !! number of paralog subsets that can be stored in `work_arr_paralog_subsets`. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
    real(real64), dimension(n_paralogs), intent(in) :: paralog_norms
      !! in subfunctionalization mode needed for subset pruning, holds the euclidean norms of paralogs (you can use the `norm` function from `f42_utils` function for this)
    integer(int32), dimension(n_paralogs), intent(in) :: sorted_paralog_norms_perm
      !! in subfunctionalization mode needed for subset pruning, as the minimum norm of the paralogs that could extend a subset should not be lower than the subset angle to the ancestor

    ! Outputs to R
    integer(int32), intent(out) :: n_results
      !! number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`
    integer(int32), dimension(n_mask_chunks, n_paralog_subsets), intent(out) :: work_arr_paralog_subsets
      !! working array to hold bitmask encoded subsets for detection.
    integer(int32), dimension(n_mask_chunks), intent(out) :: active_mask
      !! working array to hold the extended subsets
    real(real64), dimension(n_dims), intent(out) :: temp_paralog_vector
      !! vector used for pruning subsets
    real(real64), dimension(n_paralogs), intent(out) :: temp_work_array
      !! in subfunctionalization mode needed for efficient check of minimum value after a certain index
    integer(int32), intent(out) :: ierr
      !! error code

    call detect_subfunctionalization(ancestor, paralogs, n_paralogs, n_dims, rdi_threshold, filtered_paralogs_mask, &
                                   n_mask_chunks, n_results, max_subset_size, work_arr_paralog_subsets, n_paralog_subsets, &
                                   active_mask, temp_paralog_vector, paralog_norms, sorted_paralog_norms_perm, temp_work_array, ierr)
end subroutine detect_subfunctionalization_r

subroutine filter_paralogs_by_pattern_subfunctionalization_c(paralog_angles, threshold, n_paralogs, mask, n_mask_chunks, ierr) bind(C, name="filter_paralogs_by_pattern_subfunctionalization_c")
    use tox_paralog_analysis, only: filter_paralogs_by_pattern_subfunctionalization
    use, intrinsic :: iso_c_binding, only: c_double, c_int, c_ptr, c_f_pointer
    use, intrinsic :: iso_fortran_env, only: real64, int32
    implicit none

    integer(c_int), intent(in), value :: n_paralogs
        !! number of vectors in `paralogs`
    integer(c_int), intent(in), value :: n_mask_chunks
        !! number of 32 bit chunks a mask needs to encode `n_paralogs` paralogs
    type(c_ptr), intent(in) :: paralog_angles
        !! vector, holding the ascending sorted angles between ancestor and paralogs. Needed for filtering paralogs
    real(c_double), intent(in), value :: threshold
        !! filter threshold
    type(c_ptr) :: mask
        !! bit mask that will have indices of paralogs kept by pattern set to 1, else 0
    integer(c_int), intent(out) :: ierr
        !! error code

    real(real64), dimension(:), pointer :: paralog_angles_f
    integer(int32), dimension(:), pointer :: mask_f

    call c_f_pointer(paralog_angles, paralog_angles_f, [n_paralogs])
    call c_f_pointer(mask, mask_f, [n_mask_chunks])

    call filter_paralogs_by_pattern_subfunctionalization(paralog_angles_f, threshold, n_paralogs, mask_f, n_mask_chunks, ierr)
end subroutine filter_paralogs_by_pattern_subfunctionalization_c

subroutine filter_paralogs_by_pattern_subfunctionalization_r(paralog_angles, threshold, n_paralogs, mask, n_mask_chunks, ierr)
    use tox_paralog_analysis, only: filter_paralogs_by_pattern_subfunctionalization
    use, intrinsic :: iso_fortran_env, only: real64, int32
    implicit none

    integer(int32), intent(in) :: n_paralogs
        !! number of vectors in `paralogs`
    integer(int32), intent(in) :: n_mask_chunks
        !! number of 32 bit chunks a mask needs to encode `n_paralogs` paralogs
    real(real64), dimension(n_paralogs), intent(in) :: paralog_angles
        !! vector, holding the ascending sorted angles between ancestor and paralogs. Needed for filtering paralogs
    real(real64), intent(in) :: threshold
        !! filter threshold
    integer(int32), dimension(n_mask_chunks), intent(out) :: mask
        !! bit mask that will have indices of paralogs kept by pattern set to 1, else 0
    integer(int32), intent(out) :: ierr
        !! error code

    call filter_paralogs_by_pattern_subfunctionalization(paralog_angles, threshold, n_paralogs, mask, n_mask_chunks, ierr)
end subroutine filter_paralogs_by_pattern_subfunctionalization_r

subroutine filter_paralogs_by_pattern_dosage_effect_c(paralog_angles, threshold, n_paralogs, mask, n_mask_chunks, ierr) bind(C, name="filter_paralogs_by_pattern_dosage_effect")
    use tox_paralog_analysis, only: filter_paralogs_by_pattern_dosage_effect
    use, intrinsic :: iso_c_binding, only: c_double, c_int, c_ptr, c_f_pointer
    use, intrinsic :: iso_fortran_env, only: real64, int32
    implicit none

    integer(c_int), intent(in), value :: n_paralogs
        !! number of vectors in `paralogs`
    integer(c_int), intent(in), value :: n_mask_chunks
        !! number of 32 bit chunks a mask needs to encode `n_paralogs` paralogs
    type(c_ptr), intent(in) :: paralog_angles
        !! vector, holding the ascending sorted angles between ancestor and paralogs. Needed for filtering paralogs
    real(c_double), intent(in), value :: threshold
        !! filter threshold
    type(c_ptr) :: mask
        !! bit mask that will have indices of paralogs kept by pattern set to 1, else 0
    integer(c_int), intent(out) :: ierr
        !! error code

    real(real64), dimension(:), pointer :: paralog_angles_f
    integer(int32), dimension(:), pointer :: mask_f

    call c_f_pointer(paralog_angles, paralog_angles_f, [n_paralogs])
    call c_f_pointer(mask, mask_f, [n_mask_chunks])

    call filter_paralogs_by_pattern_dosage_effect(paralog_angles_f, threshold, n_paralogs, mask_f, n_mask_chunks, ierr)
end subroutine filter_paralogs_by_pattern_dosage_effect_c

subroutine filter_paralogs_by_pattern_dosage_effect_r(paralog_angles, threshold, n_paralogs, mask, n_mask_chunks, ierr)
    use tox_paralog_analysis, only: filter_paralogs_by_pattern_dosage_effect
    use, intrinsic :: iso_fortran_env, only: real64, int32
    implicit none

    integer(int32), intent(in) :: n_paralogs
        !! number of vectors in `paralogs`
    integer(int32), intent(in) :: n_mask_chunks
        !! number of 32 bit chunks a mask needs to encode `n_paralogs` paralogs
    real(real64), dimension(n_paralogs), intent(in) :: paralog_angles
        !! vector, holding the ascending sorted angles between ancestor and paralogs. Needed for filtering paralogs
    real(real64), intent(in) :: threshold
        !! filter threshold
    integer(int32), dimension(n_mask_chunks), intent(out) :: mask
        !! bit mask that will have indices of paralogs kept by pattern set to 1, else 0
    integer(int32), intent(out) :: ierr
        !! error code

    call filter_paralogs_by_pattern_dosage_effect(paralog_angles, threshold, n_paralogs, mask, n_mask_chunks, ierr)
end subroutine filter_paralogs_by_pattern_dosage_effect_r

subroutine calc_work_arr_paralog_subsets_size_c(max_subset_size, n_paralogs, work_array_size, filtered_paralogs_mask, n_mask_chunks, ierr) bind(C, name="calc_work_arr_paralog_subsets_size")
    use tox_paralog_analysis, only: calc_work_arr_paralog_subsets_size
    use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_f_pointer
    use, intrinsic :: iso_fortran_env, only: int32
    implicit none

    integer(c_int), intent(in), value :: n_paralogs
    integer(c_int), intent(in), value :: n_mask_chunks
        !! number of 32 bit chunks a mask needs to encode `n_paralogs` paralogs
    integer(c_int), intent(inout) :: max_subset_size
    integer(c_int), intent(out) :: work_array_size
    type(c_ptr), intent(in) :: filtered_paralogs_mask
    integer(c_int), intent(out) :: ierr

    integer(int32), dimension(:), pointer :: filtered_paralogs_mask_f
    call c_f_pointer(filtered_paralogs_mask, filtered_paralogs_mask_f, [n_mask_chunks])
    call calc_work_arr_paralog_subsets_size(max_subset_size, n_paralogs, work_array_size, filtered_paralogs_mask_f, n_mask_chunks, ierr)
end subroutine calc_work_arr_paralog_subsets_size_c

subroutine calc_work_arr_paralog_subsets_size_r(max_subset_size, n_paralogs, work_array_size, filtered_paralogs_mask, n_mask_chunks, ierr)
    use tox_paralog_analysis, only: calc_work_arr_paralog_subsets_size
    use, intrinsic :: iso_fortran_env, only: real64, int32
    implicit none

    integer(int32), intent(in) :: n_paralogs
    integer(int32), intent(in) :: n_mask_chunks
        !! number of 32 bit chunks a mask needs to encode `n_paralogs` paralogs
    integer(int32), intent(inout) :: max_subset_size
    integer(int32), intent(out) :: work_array_size
    integer(int32), dimension(n_mask_chunks), intent(in) :: filtered_paralogs_mask
    integer(int32), intent(out) :: ierr

    call calc_work_arr_paralog_subsets_size(max_subset_size, n_paralogs, work_array_size, filtered_paralogs_mask, n_mask_chunks, ierr)
end subroutine calc_work_arr_paralog_subsets_size_r

subroutine mask_chunk_count_c(n_paralogs, count) bind(C, name="mask_chunk_count_c")
    use tox_paralog_analysis, only: mask_chunk_count
    use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_f_pointer
    use, intrinsic :: iso_fortran_env, only: int32
    implicit none

    integer(c_int), intent(in), value :: n_paralogs
        !! number of vectors in `paralogs`
    integer(c_int), intent(out) :: count
        !! number of 32 bit chunks a mask needs to encode `n_paralogs` paralogs


    call mask_chunk_count(n_paralogs, count)
end subroutine mask_chunk_count_c

subroutine mask_chunk_count_r(n_paralogs, count)
    use tox_paralog_analysis, only: mask_chunk_count
    use, intrinsic :: iso_fortran_env, only: int32
    implicit none

    integer(int32), intent(in) :: n_paralogs
        !! number of vectors in `paralogs`
    integer(int32), intent(out) :: count
        !! number of 32 bit chunks a mask needs to encode `n_paralogs` paralogs

    call mask_chunk_count(n_paralogs, count)
end subroutine mask_chunk_count_r

subroutine mask_check_state_c(bit_mask, n_mask_chunks, i_paralog, state, ierr) bind(C, name="mask_check_state_c")
    use tox_paralog_analysis, only: mask_check_state
    use tox_conversions, only: logical_as_c_int
    use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_f_pointer
    use, intrinsic :: iso_fortran_env, only: int32
    implicit none

    integer(c_int), intent(in), value :: n_mask_chunks
        !! number of 32 bit chunks a mask needs to encode `n_paralogs` paralogs
    type(c_ptr), intent(in) :: bit_mask
        !! chunked mask to mark active paralogs
    integer(c_int), intent(in), value :: i_paralog
        !! index of paralog top be marked active
    integer(c_int), intent(out) :: state
        !! check result
    integer(c_int), intent(out) :: ierr
        !! error code

    integer(int32), dimension(:), pointer :: bit_mask_f

    call c_f_pointer(bit_mask, bit_mask_f, [n_mask_chunks])
    call logical_as_c_int(mask_check_state(bit_mask_f, i_paralog), state)
end subroutine mask_check_state_c

subroutine mask_check_state_r(bit_mask, n_mask_chunks, i_paralog, state)
    use tox_paralog_analysis, only: mask_check_state
    use, intrinsic :: iso_fortran_env, only: int32
    implicit none

    integer(int32), intent(in) :: n_mask_chunks
        !! number of 32 bit chunks a mask needs to encode `n_paralogs` paralogs
    integer(int32), dimension(n_mask_chunks), intent(in) :: bit_mask
        !! chunked mask to mark active paralogs
    integer(int32), intent(in) :: i_paralog
        !! index of paralog top be marked active
    logical, intent(out) :: state
        !! check result

    state = mask_check_state(bit_mask, i_paralog)
end subroutine mask_check_state_r