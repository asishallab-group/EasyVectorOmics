! ToDos: 
! Update error handling
! implement tests

module tox_data_tools
    use iso_fortran_env, only: real64, int32
    use tox_errors
    use array_utils, only : get_array_metadata
    use serialize_int
    use int_deserialize_mod
    use serialize_real
    use real_deserialize_mod
    use serialize_char
    use char_deserialize_mod
    implicit none
    private

    public :: read_gene_ids_from_file
    public :: read_expression_vectors
    public :: read_family_file
    public :: get_gene_index
    public :: get_family_index
    public :: get_family_for_gene_index
    public :: get_expression_vector
    public :: get_family_centroid
    public :: get_shift_components
    PUBLIC :: split_string
    public :: filter_unassigned_genes

    character(len=*), parameter :: DELIMS = ', ' // char(9)

contains

! Read tabular files (CSV/TSV)
subroutine read_expression_vectors(file_list, gene_ids, expression_vectors, &
                             n_header_rows, gene_col, value_cols, start_row, ierr, delimiter)
    character(len=*), intent(in) :: file_list(:)
    character(len=*), intent(inout) :: gene_ids(:)
    real(real64), intent(inout) :: expression_vectors(:,:)
    integer(int32), intent(in) :: n_header_rows, gene_col
    integer(int32), intent(in) :: value_cols(:) ! Array of column indices
    integer(int32), intent(in) :: start_row ! New parameter to specify the start row
    integer(int32), intent(out) :: ierr
    character(len=1), intent(in), optional :: delimiter ! Optional delimiter parameter

    integer :: i, j, k, unit, ios, idx, row_count, n_genes, expected_idx, n_value_cols
    integer :: current_sample, n_columns_in_file, n_valid_cols
    character(len=2048) :: line
    character(len=:), allocatable :: fields(:), test_fields(:)
    integer(int32), allocatable :: valid_cols(:)
    real(real64) :: value
    character(len=len(gene_ids)) :: gene
    logical :: order_mismatch
    character(len=1) :: actual_delimiter

    ierr = 0
    n_genes = size(gene_ids)
    n_value_cols = size(value_cols)
    
    ! Set the delimiter (use tab as default if not provided)
    if (present(delimiter)) then
        actual_delimiter = delimiter
    else
        actual_delimiter = char(9) ! Default to tab
    end if
    
    ! Allocate temporary array for valid columns
    allocate(valid_cols(n_value_cols))
    
    current_sample = start_row - 1
    
    do i = 1, size(file_list)
        open(newunit=unit, file=trim(file_list(i)), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            ierr = 1
            deallocate(valid_cols)
            return
        end if

        ! Skip header rows
        do j = 1, n_header_rows
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
        end do

        ! Read first data line to determine number of columns in this file
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) then
            ierr = 3
            close(unit)
            deallocate(valid_cols)
            return
        end if
        
        call split_string(line, test_fields, actual_delimiter)
        n_columns_in_file = size(test_fields)
        
        ! Rewind to beginning of data section
        rewind(unit)
        do j = 1, n_header_rows
            read(unit, '(A)', iostat=ios) line
        end do

        ! Determine valid columns for this file
        n_valid_cols = 0
        do k = 1, n_value_cols
            if (value_cols(k) <= n_columns_in_file) then
                n_valid_cols = n_valid_cols + 1
                valid_cols(n_valid_cols) = value_cols(k)
            end if
        end do

        ! Debug output
        ! write(*,*) 'Processing file: ', trim(file_list(i))
        ! write(*,*) 'Columns in file: ', n_columns_in_file
        ! write(*,*) 'Valid value columns: ', valid_cols(1:n_valid_cols)
        ! write(*,*) 'Current sample start: ', current_sample + 1

        row_count = 0
        order_mismatch = .false.
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            row_count = row_count + 1
            
            if (row_count > n_genes) then
                ierr = 2 ! More genes than allocated space 
                exit
            end if

            call split_string(line, fields, actual_delimiter)
            if (size(fields) < max(gene_col, maxval(valid_cols(1:n_valid_cols)))) cycle

            gene = trim(adjustl(fields(gene_col)))

            ! Fast path: check if gene matches expected position 
            expected_idx = row_count
            if (expected_idx <= n_genes .and. trim(adjustl(gene_ids(expected_idx))) == trim(adjustl(gene))) then
                idx = expected_idx
            else
                ! Slow path: search through all genes
                idx = get_gene_index(gene_ids, gene)
                if (idx == 0) then
                    write(*,*) 'Warning: Gene ', trim(gene), ' not found in master gene list'
                    order_mismatch = .true.
                    cycle
                end if
            end if
            
            ! Read all valid value columns for this gene 
            do k = 1, n_valid_cols
                read(fields(valid_cols(k)), *, iostat=ios) value
                if (ios == 0) then
                    expression_vectors(current_sample + k, idx) = value
                else
                    ! Set to zero if read fails
                    expression_vectors(current_sample + k, idx) = 0.0_real64
                end if
            end do
        end do
        close(unit)
        
        ! Update sample offset for next file 
        current_sample = current_sample + n_valid_cols
        
        ! Report if order mismatch was detected
        if (order_mismatch) then
            write(*,*) 'Note: Gene order mismatch detected in file: ', trim(file_list(i))
        end if
    end do
    
    deallocate(valid_cols)
end subroutine read_expression_vectors

subroutine read_gene_ids_from_file(filename, gene_ids, n_header_rows, gene_col, ierr)
    character(len=*), intent(in) :: filename
    character(len=*), intent(out) :: gene_ids(:)
    integer(int32), intent(in) :: n_header_rows, gene_col
    integer(int32), intent(out) :: ierr

    integer :: unit, ios, j, row_count
    character(len=2048) :: line
    character(len=:), allocatable :: fields(:)

    ierr = 0
    row_count = 0

    open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        ierr = 1
        return
    end if

    ! Skip header rows
    do j = 1, n_header_rows
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) exit
    end do

    ! Read data rows
    do
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        
        row_count = row_count + 1
        if (row_count > size(gene_ids)) then
            ierr = 2  ! More genes than allocated space
            exit
        end if

        call split_string(line, fields, char(9))
        if (size(fields) < gene_col) cycle

        gene_ids(row_count) = trim(adjustl(fields(gene_col)))
    end do

    close(unit)
end subroutine read_gene_ids_from_file

subroutine read_family_file(filename, gene_ids, family_ids, gene_to_fam, ierr)
    use hashmap_module
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: gene_ids(:)
    character(len=*), intent(out) :: family_ids(:)
    integer(int32), intent(out) :: gene_to_fam(:)
    integer(int32), intent(out) :: ierr

    integer(int32) :: unit, ios, i, j, fam_idx, gene_idx, n_families, n_genes
    integer(int32) :: pos, start_pos, end_pos 
    character(len=4096) :: line
    character(len=:), allocatable :: fields(:), genes(:)
    character(len=len(family_ids)) :: current_family
    type(hashmap_type) :: gene_map
    logical :: fill_family_ids
    integer(int32) :: hashmap_size

    ierr = 0
    gene_to_fam = 0
    n_families = size(family_ids)
    n_genes = size(gene_ids)
    
    fill_family_ids = all(family_ids == "")

    ! Hashmap für schnelle Gen-Suche erstellen
    hashmap_size = int(1.3 * n_genes)
    call hashmap_create(gene_map, hashmap_size)
    
    ! Hashmap mit allen Gen-IDs füllen
    !write(*,*) 'Building gene hashmap with ', n_genes, ' genes...'
    do i = 1, n_genes
        call hashmap_put(gene_map, gene_ids(i), i)
    end do
    !write(*,*) 'Hashmap complete.'

    open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
        ierr = 1
        call hashmap_destroy(gene_map)
        return
    end if

    ! Header-Zeile überspringen
    read(unit, '(A)', iostat=ios) line
    if (ios /= 0) then
        ierr = 4
        close(unit)
        call hashmap_destroy(gene_map)
        return
    end if

    !write(*,*) 'Reading family file with hashmap...'
    
    fam_idx = 0
    do
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) exit

        fam_idx = fam_idx + 1
        
        ! if (mod(fam_idx, 1000) == 0) then
        !     write(*,*) 'Processed ', fam_idx, ' families...'
        ! end if

        call split_string(line, fields, char(9))
        if (size(fields) < 2) cycle

        if (fill_family_ids .and. fam_idx > n_families) then
            ierr = 2
            exit
        end if

        current_family = trim(adjustl(fields(1)))
        if (fill_family_ids) then
            family_ids(fam_idx) = current_family
        end if

        ! Process all gene columns mit Hashmap-Lookup
        do i = 2, size(fields)
            call split_string(fields(i), genes, ',')
            do j = 1, size(genes)
                if (len_trim(genes(j)) == 0) cycle
                
                ! Hashmap-Lookup
                gene_idx = hashmap_get(gene_map, trim(adjustl(genes(j))))
                if (gene_idx > 0) then
                    gene_to_fam(gene_idx) = fam_idx
                else
                    write(*,*) 'Gene not found in hashmap: ', trim(adjustl(genes(j)))
                end if
            end do
        end do
    end do
    close(unit)
    
    call hashmap_destroy(gene_map)
    ! write(*,*) 'Finished processing ', fam_idx, ' families'
end subroutine read_family_file

!> Filter out genes that are not assigned to any family.
subroutine filter_unassigned_genes(gene_ids, expression_vectors, gene_to_fam, n_genes_kept, ierr)
    character(len=*), allocatable, intent(inout) :: gene_ids(:)
    real(real64), allocatable, intent(inout) :: expression_vectors(:,:)
    integer(int32), allocatable, intent(inout) :: gene_to_fam(:)
    integer(int32), intent(out) :: n_genes_kept
    integer(int32), intent(out) :: ierr
    
    integer(int32) :: i, j, n_genes_total, n_samples
    integer(int32), allocatable :: valid_indices(:)
    character(len=len(gene_ids)), allocatable :: temp_gene_ids(:)
    real(real64), allocatable :: temp_expression_vectors(:,:)
    integer(int32), allocatable :: temp_gene_to_fam(:)
    
    ierr = ERR_OK
    n_genes_total = size(gene_ids)
    n_samples = size(expression_vectors, 1)
    
    ! Count genes with valid family assignments
    n_genes_kept = count(gene_to_fam >= 1)
    
    if (n_genes_kept == 0) then
        ierr = ERR_EMPTY_INPUT
        return
    end if
    
    ! Allocate temporary arrays for valid genes
    allocate(valid_indices(n_genes_kept))
    allocate(temp_gene_ids(n_genes_kept))
    allocate(temp_expression_vectors(n_samples, n_genes_kept))
    allocate(temp_gene_to_fam(n_genes_kept))
    
    ! Collect indices of genes with valid family assignments
    j = 1
    do i = 1, n_genes_total
        if (gene_to_fam(i) >= 1) then
            valid_indices(j) = i
            j = j + 1
        end if
    end do
    
    ! Copy valid genes to temporary arrays
    do i = 1, n_genes_kept
        temp_gene_ids(i) = gene_ids(valid_indices(i))
        temp_expression_vectors(:, i) = expression_vectors(:, valid_indices(i))
        temp_gene_to_fam(i) = gene_to_fam(valid_indices(i))
    end do
    
    ! Deallocate original arrays and replace with filtered ones
    deallocate(gene_ids, expression_vectors, gene_to_fam)
    
    call move_alloc(temp_gene_ids, gene_ids)
    call move_alloc(temp_expression_vectors, expression_vectors)
    call move_alloc(temp_gene_to_fam, gene_to_fam)
    
    write(*,*) 'Filtered out ', n_genes_total - n_genes_kept, ' unassigned genes'
    write(*,*) 'Kept ', n_genes_kept, ' genes with valid family assignments'
end subroutine filter_unassigned_genes

subroutine count_family_file_stats(filename, n_families, n_genes, ierr)
    character(len=*), intent(in) :: filename
    integer(int32), intent(out) :: n_families, n_genes, ierr
    
    integer(int32) :: unit, ios, i
    character(len=1024) :: line
    character(len=:), allocatable :: fields(:), genes(:)
    
    ierr = 0
    n_families = 0
    n_genes = 0
    
    open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
        ierr = 1
        return
    end if

    ! Header-Zeile überspringen
    read(unit, '(A)', iostat=ios) line
    if (ios /= 0) then
        ierr = 4
        close(unit)
        return
    end if

    do
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) exit

        n_families = n_families + 1
        call split_string(line, fields, char(9))
        
        if (size(fields) >= 2) then
            do i = 2, size(fields)
                call split_string(fields(i), genes, ',')
                n_genes = n_genes + size(genes)
            end do
        end if
    end do
    close(unit)
end subroutine count_family_file_stats

! Helper function to find gene index with early exit optimization
integer function get_gene_index(gene_ids, gene) result(idx)
    character(len=*), intent(in) :: gene_ids(:)
    character(len=*), intent(in) :: gene
    integer :: i
    
    idx = 0
    ! Quick check: if gene is empty, return immediately
    if (len_trim(gene) == 0) return
    
    ! Search through gene list
    do i = 1, size(gene_ids)
        if (trim(adjustl(gene_ids(i))) == trim(adjustl(gene))) then
            idx = i
            return
        end if
    end do
end function get_gene_index

! Helper function to find family index
integer function get_family_index(family_ids, family) result(idx)
    character(len=*), intent(in) :: family_ids(:)
    character(len=*), intent(in) :: family
    integer :: i
    
    idx = 0
    do i = 1, size(family_ids)
        if (trim(adjustl(family_ids(i))) == trim(adjustl(family))) then
            idx = i
            return
        end if
    end do
end function get_family_index

! Helper subroutine to split strings
subroutine split_string(input, output, delimiter)
    character(len=*), intent(in) :: input
    character(len=:), allocatable, intent(out) :: output(:)
    character(len=1), optional, intent(in) :: delimiter
    
    character(len=1) :: delim
    integer :: n, i, start_pos, end_pos, field_count
    integer, allocatable :: field_starts(:), field_ends(:)
    
    if (present(delimiter)) then
        delim = delimiter
    else
        delim = ' '
    end if

    ! First pass: find field boundaries
    allocate(field_starts(len_trim(input)))
    allocate(field_ends(len_trim(input)))
    
    n = 0
    start_pos = 1
    do i = 1, len_trim(input)
        if (input(i:i) == delim) then
            if (start_pos <= i-1) then
                n = n + 1
                field_starts(n) = start_pos
                field_ends(n) = i-1
            end if
            start_pos = i + 1
        end if
    end do
    
    ! Add last field if exists
    if (start_pos <= len_trim(input)) then
        n = n + 1
        field_starts(n) = start_pos
        field_ends(n) = len_trim(input)
    end if

    ! Allocate output array
    allocate(character(len=maxval(field_ends(1:n) - field_starts(1:n) + 1)) :: output(n))
    
    ! Extract fields
    do i = 1, n
        output(i) = trim(adjustl(input(field_starts(i):field_ends(i))))
    end do
    
    deallocate(field_starts, field_ends)
end subroutine split_string

subroutine get_family_for_gene_index(gene_idx, gene_to_fam, out_family_idx, ierr)
    integer(int32), intent(in) :: gene_idx
    integer(int32), intent(in) :: gene_to_fam(:)
    integer(int32), intent(out) :: out_family_idx
    integer(int32), intent(out) :: ierr
    
    ierr = ERR_OK
    
    if (gene_idx < 1 .or. gene_idx > size(gene_to_fam)) then
        ierr = ERR_INVALID_INPUT
        out_family_idx = 0
        return
    end if
    
    out_family_idx = gene_to_fam(gene_idx)
    
    if (out_family_idx < 1) then
        ierr = ERR_INVALID_INPUT
    end if
end subroutine get_family_for_gene_index

subroutine get_expression_vector(gene_idx, expression_vectors, out_vector, ierr)
    integer(int32), intent(in) :: gene_idx
    real(real64), intent(in) :: expression_vectors(:,:)
    real(real64), intent(out) :: out_vector(:)
    integer(int32), intent(out) :: ierr
    
    ierr = ERR_OK
    
    if (gene_idx < 1 .or. gene_idx > size(expression_vectors, 2)) then
        ierr = ERR_INVALID_INPUT
        return
    end if
    
    if (size(out_vector) /= size(expression_vectors, 1)) then
        ierr = ERR_SIZE_MISMATCH
        return
    end if
    
    out_vector = expression_vectors(:, gene_idx)
end subroutine get_expression_vector

subroutine get_family_centroid(family_idx, family_centroids, out_centroid, ierr)
    integer(int32), intent(in) :: family_idx
    real(real64), intent(in) :: family_centroids(:,:)
    real(real64), intent(out) :: out_centroid(:)
    integer(int32), intent(out) :: ierr
    
    ierr = ERR_OK
    
    if (family_idx < 1 .or. family_idx > size(family_centroids, 2)) then
        ierr = ERR_INVALID_INPUT
        return
    end if
    
    if (size(out_centroid) /= size(family_centroids, 1)) then
        ierr = ERR_SIZE_MISMATCH
        return
    end if
    
    out_centroid = family_centroids(:, family_idx)
end subroutine get_family_centroid

subroutine get_shift_components(gene_idx, shift_vectors, d, out_start, out_shift, ierr)
    integer(int32), intent(in) :: gene_idx
    real(real64), intent(in) :: shift_vectors(:,:)
    integer(int32), intent(in) :: d
    real(real64), intent(out) :: out_start(:), out_shift(:)
    integer(int32), intent(out) :: ierr
    
    ierr = ERR_OK
    
    if (gene_idx < 1 .or. gene_idx > size(shift_vectors, 2)) then
        ierr = ERR_INVALID_INPUT
        return
    end if
    
    if (size(shift_vectors, 1) < 2*d) then
        ierr = ERR_SIZE_MISMATCH
        return
    end if
    
    if (size(out_start) /= d .or. size(out_shift) /= d) then
        ierr = ERR_SIZE_MISMATCH
        return
    end if
    
    out_start = shift_vectors(1:d, gene_idx)
    out_shift = shift_vectors(d+1:2*d, gene_idx)
end subroutine get_shift_components

end module tox_data_tools