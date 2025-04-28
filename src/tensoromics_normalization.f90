module tensoromics_normalization
    use iso_c_binding
    implicit none
contains

    !-----------------------------------------------------------------------
    ! Subroutine to normalize each gene's expression values across tissues
    ! by the square root of the mean of squares (sqrt(mean(x^2))).
    !
    ! Arguments:
    ! n_genes        : Number of genes (rows)
    ! n_tissues      : Number of tissues (columns)
    ! input_matrix   : Flattened input matrix (column-major)
    ! output_matrix  : Flattened output matrix (normalized values)
    !
    ! This normalization method differs slightly from classical std dev:
    ! it normalizes based on sqrt(mean(x^2)) instead of sqrt(mean((x-mean)^2)).
    !-----------------------------------------------------------------------
    pure subroutine normalize_by_std_dev(n_genes, n_tissues, input_matrix, output_matrix) bind(c, name="normalize_by_std_dev")
        use iso_c_binding
        implicit none
        
        ! === Variable declarations ===
        integer(c_int), intent(in) :: n_genes, n_tissues
        real(c_double), intent(in) :: input_matrix(*)
        real(c_double), intent(out) :: output_matrix(*)
        
        integer :: i, j
        real(c_double) :: mean_value, std_dev, temp_sum

        ! === Loop over each gene ===
        do i = 1, n_genes

            ! -- Calculate sqrt(mean(x^2)) for the gene across tissues --
            temp_sum = 0.0d0
            do j = 1, n_tissues
                ! Sum of squares of the expression values
                temp_sum = temp_sum + input_matrix((j-1)*n_genes + i)**2
            end do

            ! Final calculation of the standard deviation-like value
            std_dev = sqrt(temp_sum / dble(n_tissues))

            ! -- Avoid division by zero --
            if (std_dev == 0.0d0) std_dev = 1.0d0

            ! -- Normalize the values for this gene across all tissues --
            do j = 1, n_tissues
                output_matrix((j-1)*n_genes + i) = input_matrix((j-1)*n_genes + i) / std_dev
            end do

        end do

    end subroutine normalize_by_std_dev


    !-----------------------------------------------------------------------
    ! Perform quantile normalization on the input matrix.
    ! Each column (tissue) is sorted independently, and the average 
    ! expression value across tissues for each rank is computed.
    ! The values are then replaced based on rank averages.
    !-----------------------------------------------------------------------
    subroutine quantile_normalization(n_genes, n_tissues, input_matrix, output_matrix) bind(c, name="quantile_normalization")
        integer(c_int), intent(in) :: n_genes, n_tissues
        real(c_double), intent(in) :: input_matrix(*)
        real(c_double), intent(out) :: output_matrix(*)

        integer :: i, j
        integer, dimension(:), allocatable :: idx
        real(c_double), dimension(:), allocatable :: temp_col, rank_means

        ! Allocate temporary arrays
        allocate(temp_col(n_genes))
        allocate(rank_means(n_genes))
        allocate(idx(n_genes))

        ! Initialize rank means to zero
        rank_means = 0.0d0

        ! === First pass: accumulate values by rank across tissues ===
        do j = 1, n_tissues
            ! Copy the column into temp_col and initialize indices
            do i = 1, n_genes
                temp_col(i) = input_matrix((j-1)*n_genes + i)
                idx(i) = i
            end do

            ! Sort temp_col and update idx
            call quicksort_iterative(temp_col, idx, n_genes)

            ! Accumulate values for each rank
            do i = 1, n_genes
                rank_means(i) = rank_means(i) + temp_col(i)
            end do
        end do

        ! Compute average rank values across tissues
        do i = 1, n_genes
            rank_means(i) = rank_means(i) / dble(n_tissues)
        end do

        ! === Second pass: replace original values by rank averages ===
        do j = 1, n_tissues
            ! Copy the column into temp_col and initialize indices again
            do i = 1, n_genes
                temp_col(i) = input_matrix((j-1)*n_genes + i)
                idx(i) = i
            end do

            ! Sort temp_col and update idx
            call quicksort_iterative(temp_col, idx, n_genes)

            ! Replace values in output matrix according to the sorted order
            do i = 1, n_genes
                output_matrix((j-1)*n_genes + idx(i)) = rank_means(i)
            end do
        end do

        ! Deallocate temporary arrays
        deallocate(temp_col)
        deallocate(rank_means)
        deallocate(idx)
    end subroutine quantile_normalization

    !-----------------------------------------------------------------------
    ! Sort an array of real numbers and track the original indices using
    ! a non-recursive, iterative quicksort algorithm (stack-based).
    !-----------------------------------------------------------------------
    subroutine quicksort_iterative(arr, idx, n)
        real(c_double), intent(inout) :: arr(*)
        integer, intent(inout) :: idx(*)
        integer, intent(in) :: n

        integer, parameter :: max_stack = 10000
        integer :: left, right, i, j, top
        integer, dimension(max_stack) :: lstack, rstack
        real(c_double) :: pivot

        ! Initialize the stack
        top = 1
        lstack(top) = 1
        rstack(top) = n

        ! Main loop
        do while (top > 0)
            left = lstack(top)
            right = rstack(top)
            top = top - 1

            if (left >= right) cycle

            pivot = arr((left + right) / 2)
            i = left
            j = right

            ! Partition the array
            do
                do while (arr(i) < pivot)
                    i = i + 1
                end do
                do while (pivot < arr(j))
                    j = j - 1
                end do
                if (i <= j) then
                    call swap_real(arr(i), arr(j))  ! Swap array values
                    call swap_int(idx(i), idx(j))  ! Swap index values
                    i = i + 1
                    j = j - 1
                end if
                if (i > j) exit
            end do

            ! Push left and right partitions onto the stack
            if (left < j) then
                top = top + 1
                lstack(top) = left
                rstack(top) = j
            end if

            if (i < right) then
                top = top + 1
                lstack(top) = i
                rstack(top) = right
            end if
        end do
    end subroutine quicksort_iterative

    !-----------------------------------------------------------------------
    ! Swap two real (double precision) numbers
    !-----------------------------------------------------------------------
    subroutine swap_real(a, b)
        real(c_double), intent(inout) :: a, b
        real(c_double) :: temp

        temp = a
        a = b
        b = temp
    end subroutine swap_real

    !-----------------------------------------------------------------------
    ! Swap two integers
    !-----------------------------------------------------------------------
    subroutine swap_int(a, b)
        integer, intent(inout) :: a, b
        integer :: temp

        temp = a
        a = b
        b = temp
    end subroutine swap_int


    !-----------------------------------------------------------------------
    ! Apply log2(x + 1) transformation to each element of the input matrix.
    !
    ! Arguments:
    ! n_genes        : Number of genes (rows)
    ! n_tissues      : Number of tissues (columns)
    ! input_matrix   : Flattened input matrix (column-major)
    ! output_matrix  : Flattened output matrix (after log2 transformation)
    !
    ! Each value is transformed by:
    ! output = log(input + 1) / log(2)
    ! (Equivalent to log base 2, but avoiding intrinsic log2 function
    ! for better WASM compatibility.)
    !-----------------------------------------------------------------------
    subroutine log2_transformation(n_genes, n_tissues, input_matrix, output_matrix) bind(c, name="log2_transformation")
        integer(c_int), intent(in) :: n_genes, n_tissues
        real(c_double), intent(in) :: input_matrix(*)
        real(c_double), intent(out) :: output_matrix(*)
        integer :: i

        ! Loop through all elements in the flattened input matrix
        do i = 1, n_genes * n_tissues
            ! Apply the log2(x + 1) transformation
            output_matrix(i) = log(input_matrix(i) + 1.0d0) / log(2.0d0)
        end do
    end subroutine log2_transformation



    !-----------------------------------------------------------------------
    ! Calculate tissue averages by averaging replicates within each group.
    !
    ! Arguments:
    ! n_gene        : Number of genes (rows)
    ! n_col         : Number of columns (tissues)
    ! n_grps        : Number of groups (e.g., tissue groups after averaging replicates)
    ! group_s       : Array containing the starting column index for each group
    ! group_c       : Array containing the number of columns (replicates) for each group
    ! input_matrix  : Flattened input matrix (column-major)
    ! output_matrix : Flattened output matrix (after averaging)
    !
    ! For each group:
    !   - Sum the replicates for each gene
    !   - Divide by the number of replicates to compute the average
    !-----------------------------------------------------------------------
    subroutine calc_tiss_avg(n_gene, n_col, n_grps, group_s, group_c, input_matrix, output_matrix) bind(C, name="calc_tiss_avg")

        ! === Variable declarations ===
        integer(c_int), intent(in) :: n_gene, n_col, n_grps
        integer(c_int), intent(in) :: group_s(*)
        integer(c_int), intent(in) :: group_c(*)
        real(c_double), intent(in) :: input_matrix(*)
        real(c_double), intent(out) :: output_matrix(*)

        integer :: i, j, g, col
        real(c_double) :: sum_val
        integer :: start_idx, count_cols

        ! === Main loop over each group ===
        do g = 1, n_grps
            start_idx = group_s(g)   ! Start column index for this group
            count_cols = group_c(g)  ! Number of replicates in this group

            ! --- Loop over each gene ---
            do i = 1, n_gene
                sum_val = 0.0d0      ! Initialize the sum for this gene

                ! --- Sum the replicates for this gene across tissues ---
                do j = 0, count_cols - 1
                    col = start_idx + j
                    sum_val = sum_val + input_matrix((col-1)*n_gene + i)
                end do

                ! --- Store the average value in the output matrix ---
                output_matrix((g-1)*n_gene + i) = sum_val / dble(count_cols)
            end do
        end do

    end subroutine calc_tiss_avg


    !-----------------------------------------------------------------------
    ! Calculate the log2 fold change between condition and control columns.
    !
    ! Arguments:
    ! n_genes        : Number of genes (rows)
    ! n_pairs        : Number of condition-control pairs
    ! control_cols   : Indices of control columns
    ! cond_cols      : Indices of condition columns
    ! i_matrix       : Flattened input matrix (column-major)
    ! o_matrix       : Flattened output matrix (one column per pair)
    !
    ! For each pair:
    !   - Compute (condition value - control value) for every gene
    !-----------------------------------------------------------------------
    subroutine calc_fchange(n_genes, n_pairs, control_cols, cond_cols, i_matrix, o_matrix) bind(C, name="calc_fchange")

        ! === Variable declarations ===
        integer(c_int), intent(in) :: n_genes, n_pairs
        integer(c_int), intent(in) :: control_cols(*)
        integer(c_int), intent(in) :: cond_cols(*)
        real(c_double), intent(in) :: i_matrix(*)
        real(c_double), intent(out) :: o_matrix(*)

        integer :: i, p
        integer :: control_col, cond_col

        ! === Loop over all condition-control pairs ===
        do p = 1, n_pairs
            control_col = control_cols(p)  ! Get control column index
            cond_col = cond_cols(p)        ! Get condition column index

            ! --- Loop over each gene ---
            do i = 1, n_genes
                ! Compute fold change: condition - control
                o_matrix((p-1)*n_genes + i) = i_matrix((cond_col-1)*n_genes + i) - i_matrix((control_col-1)*n_genes + i)
            end do
        end do

    end subroutine calc_fchange


end module tensoromics_normalization
