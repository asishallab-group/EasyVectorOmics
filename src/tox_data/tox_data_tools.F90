module tox_data_tools
    use iso_fortran_env, only: real64, int32
    use tox_errors, only: set_ok, set_err_once, is_err, set_err, check_io_stat
    use tox_errors, only: ERR_INVALID_INPUT, ERR_FILE_OPEN, ERR_READ_DATA
    use array_utils, only :ascii_to_string, string_to_ascii, check_okay_ioerror, ERR_SIZE_MISMATCH
    use config, only: DEBUG
    implicit none
    private

    public :: read_gene_ids_from_file
    public :: read_expression_vectors
    public :: read_family_file
    public :: split_string
    public :: filter_unassigned_genes, get_unassigned_mask

    character(len=*), parameter :: DELIMS = ', ' // char(9)

contains

!> Read tabular expression files
subroutine read_expression_vectors(file_list, gene_ids, expression_vectors, &
                             n_header_rows, gene_col, value_cols, start_row, ierr, delimiter)
    use xxh3_hashmap_module
    character(len=*), intent(in) :: file_list(:)
    !! List of files to read from
    character(len=*), intent(in) :: gene_ids(:)
    !! Array of gene IDS
    real(real64), intent(inout) :: expression_vectors(:,:)
    !! Array of expression vectors
    integer(int32), intent(in) :: n_header_rows
    !! Number of header rows to skip
    integer(int32), intent(in) :: gene_col
    !! Index of column with gene_ids
    integer(int32), intent(in) :: value_cols(:)
    !! Indicies of columns containing values
    integer(int32), intent(in) :: start_row
    !! Row in the expression vectors to start in
    integer(int32), intent(out) :: ierr
    !! Error code
    character(len=1), intent(in), optional :: delimiter
    !! optional delimiter, default is tab

    integer(int32) :: i, j, k, unit, ios, idx, n_genes, expected_idx, n_value_cols
    integer(int32) :: current_sample, n_columns_in_file, n_valid_cols
    character(len=2048) :: line
    character(len=:), allocatable :: fields(:), test_fields(:)
    integer(int32), allocatable :: valid_cols(:)
    real(real64) :: value
    character(len=len(gene_ids)) :: gene
    character(len=1) :: actual_delimiter
    integer(int32) :: current_row

    ! Hashmap for gene lookup
    type(hashmap_type) :: gene_map

    call set_ok(ierr)
    call set_ok(ios)
    n_genes = size(gene_ids)
    n_value_cols = size(value_cols)
    
    ! Set the delimiter (use tab as default if not provided)
    if (present(delimiter)) then
        actual_delimiter = delimiter
    else
        actual_delimiter = char(9)
    end if
    
    if (n_value_cols <= 0) then
        call set_err_once(ierr, ERR_INVALID_INPUT)
        if(DEBUG) write(*,*) 'Error: No value columns specified.'
        return
    end if

    if (size(gene_ids) < n_genes) then
        call set_err_once(ierr, ERR_INVALID_INPUT)
        if(DEBUG) write(*,*) 'Error: gene_ids array is too small.'
        return
    end if

    ! Allocate temporary array for valid columns
    allocate(valid_cols(n_value_cols), stat = ios)
    call check_io_stat(ios, ierr)
    if(is_err(ierr)) return

    ! Build gene hashmap for fast lookup
    call hashmap_create(gene_map, initial_size=n_genes)
    do i = 1, n_genes
        call hashmap_put(gene_map, trim(gene_ids(i)), i)
    end do
    
    current_sample = start_row - 1
    
    do i = 1, size(file_list)
        if(DEBUG) write(*,*) 'Reading file: ', trim(file_list(i))
        open(newunit=unit, file=trim(file_list(i)), status='old', action='read', iostat=ios)
        call check_okay_ioerror(ios, ierr, ERR_FILE_OPEN, unit)
        if(is_err(ierr)) then
            call hashmap_destroy(gene_map)
            return
        end if

        ! Skip header rows
        do j = 1, n_header_rows
            read(unit, '(A)', iostat=ios) line
            call check_okay_ioerror(ios, ierr, ERR_READ_DATA, unit)
            if(is_err(ierr)) then
                call hashmap_destroy(gene_map)
                RETURN
            end if
        end do

        ! Read first data line to determine number of columns in this file
        read(unit, '(A)', iostat=ios) line
        call check_okay_ioerror(ios, ierr, ERR_READ_DATA, unit)
        if (is_err(ierr)) then
            call hashmap_destroy(gene_map)
            return
        end if
        
        call split_string(line, test_fields, ierr, actual_delimiter)
        n_columns_in_file = size(test_fields)
        
        ! Rewind to beginning of data section
        rewind(unit)
        do j = 1, n_header_rows
            read(unit, '(A)', iostat=ios) line
            call check_okay_ioerror(ios, ierr, ERR_READ_DATA, unit)
            if(is_err(ierr)) then
                call hashmap_destroy(gene_map)
                RETURN
            end if
        end do

        ! Determine valid columns for this file
        n_valid_cols = 0
        do k = 1, n_value_cols
            if (value_cols(k) <= n_columns_in_file) then
                n_valid_cols = n_valid_cols + 1
                valid_cols(n_valid_cols) = value_cols(k)
            else 
                call set_err(ierr, ERR_INVALID_INPUT)
                return 
            end if
        end do
        if (n_valid_cols == 0) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            if(DEBUG) write(*,*) 'Error: No valid value columns found in file: ', trim(file_list(i))
            close(unit)
            cycle
        end if

        current_row = 0
        do
            if(current_row > size(gene_ids)) then
                call set_err_once(ierr, ERR_INVALID_INPUT)
                if(DEBUG) write(*,*) 'Provided file contains more lines then expected'
                close(unit)
                return
            end if
            current_row = current_row + 1
            read(unit, '(A)', iostat=ios) line
            if(ios < 0) exit !End of file
            call check_okay_ioerror(ios, ierr, ERR_READ_DATA, unit)
            if(is_err(ierr)) then
                call hashmap_destroy(gene_map)
                RETURN
            end if

            call split_string(line, fields, ierr, actual_delimiter)
            if (size(fields) < max(gene_col, maxval(valid_cols(1:n_valid_cols)))) cycle

            gene = trim(adjustl(fields(gene_col)))

            ! Use hashmap for gene lookup
            idx = hashmap_get(gene_map, gene)
            if (idx == -1) then
                if(DEBUG) write(*,*) 'Warning: Gene ', trim(gene), ' not found in master gene list'
                cycle
            end if
            
            ! Read all valid value columns for this gene 
            do k = 1, n_valid_cols
                read(fields(valid_cols(k)), *, iostat=ios) value
                if (ios == 0) then
                    expression_vectors(current_sample + k, idx) = value
                else
                    if(DEBUG) write(*,*) 'Missing values in row: ', k
                    expression_vectors(current_sample + k, idx) = 0.0_real64
                end if
            end do
        end do
        close(unit)
        
        current_sample = current_sample + n_valid_cols

    end do

    call hashmap_destroy(gene_map)
    deallocate(valid_cols)
end subroutine read_expression_vectors

!> Only read the gene ids from a file
subroutine read_gene_ids_from_file(filename, gene_ids, n_header_rows, gene_col, ierr)
    character(len=*), intent(in) :: filename
        !! Name of the file
    character(len=*), intent(out) :: gene_ids(:)
        !! gene ids array
    integer(int32), intent(in) :: n_header_rows
        !! number of headers to skip
    integer(int32), intent(in) :: gene_col
        !! Index of the column containing gene ids
    integer(int32), intent(out) :: ierr
        !! Error code

    integer(int32) :: unit, ios, j, row_count
    character(len=2048) :: line
    character(len=:), allocatable :: fields(:)

    call set_ok(ierr)
    call set_ok(ios)
    row_count = 0

    open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
    call check_okay_ioerror(ios, ierr, ERR_FILE_OPEN, unit)
    if(is_err(ierr)) return

    ! Skip header rows
    do j = 1, n_header_rows
        read(unit, '(A)', iostat=ios) line
        call check_okay_ioerror(ios, ierr, ERR_READ_DATA, unit)
        if(is_err(ierr)) return
    end do

    ! Read data rows
    do
        read(unit, '(A)', iostat=ios) line
        if (ios < 0) exit
        call check_okay_ioerror(ios, ierr, ERR_READ_DATA, unit)
        if(is_err(ierr)) return
        
        row_count = row_count + 1
        if (row_count > size(gene_ids)) then
            call set_err_once(ierr, ERR_SIZE_MISMATCH)  ! More genes than allocated space
            close(unit)
            return
        end if

        call split_string(line, fields, ierr, char(9))
        if (size(fields) < gene_col) cycle

        gene_ids(row_count) = trim(adjustl(fields(gene_col)))
    end do

    close(unit)
end subroutine read_gene_ids_from_file

!> Read a family file (Orthofinder)
subroutine read_family_file(filename, gene_ids, family_ids, gene_to_fam, ierr)
    use xxh3_hashmap_module
    character(len=*), intent(in) :: filename
        !! Name of the file
    character(len=*), intent(in) :: gene_ids(:)
        !! gene ids array
    character(len=*), intent(out) :: family_ids(:)
        !! family ids array
    integer(int32), intent(out) :: gene_to_fam(:)
        !! gene to family mapping
    integer(int32), intent(out) :: ierr
        !! Error code

    integer(int32) :: unit, ios, i, j, fam_idx, gene_idx, n_families, n_genes
    integer(int32) :: pos, start_pos, end_pos 
    character(len=4096) :: line
    character(len=:), allocatable :: fields(:), genes(:)
    character(len=len(family_ids)) :: current_family
    type(hashmap_type) :: gene_map
    integer(int32) :: hashmap_size

    call set_ok(ierr)
    gene_to_fam = 0
    n_families = size(family_ids)
    n_genes = size(gene_ids)
    
    ! Hashmap for efficient lookups
    call hashmap_create(gene_map, initial_size = n_genes)
    
    ! Fill hashmap
    do i = 1, n_genes
        call hashmap_put(gene_map, gene_ids(i), i)
    end do

    open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
    if (is_err(ios)) then
        if(DEBUG) write(*,*) 'Error opening file: ', trim(filename)
        call set_err_once(ierr, ERR_FILE_OPEN)
        call hashmap_destroy(gene_map)
        return
    end if

    ! skip header
    read(unit, '(A)', iostat=ios) line
    call check_okay_ioerror(ios, ierr, ERR_READ_DATA, unit)
    if(is_err(ierr)) then
        call hashmap_destroy(gene_map)
        RETURN
    end if
    
    fam_idx = 0
    do
        read(unit, '(A)', iostat=ios) line
        if (ios < 0) exit
        call check_okay_ioerror(ios, ierr, ERR_READ_DATA, unit)
        if(is_err(ierr)) then
            call hashmap_destroy(gene_map)
            RETURN
        end if

        fam_idx = fam_idx + 1

        call split_string(line, fields, ierr, char(9))


        if (size(fields) < 2) then
            if(DEBUG) write(*,*) 'Warning: Skipping invalid line (less than 2 fields).'
            cycle
        end if

        current_family = trim(adjustl(fields(1)))

        family_ids(fam_idx) = current_family

        ! Process all gene columns with Hashmap-Lookup
        do i = 2, size(fields)
            call split_string(fields(i), genes, ierr, ',')
            if(is_err(ierr)) return
            do j = 1, size(genes)
                if (len_trim(genes(j)) == 0) cycle
                
                ! Hashmap-Lookup
                gene_idx = hashmap_get(gene_map, trim(adjustl(genes(j))))
                if (gene_idx > 0) then
                    gene_to_fam(gene_idx) = fam_idx
                else
                    if(DEBUG) write(*,*) 'Gene not found in hashmap: ', trim(adjustl(genes(j)))
                end if
            end do
        end do
    end do
    close(unit)
    
    call hashmap_destroy(gene_map)
end subroutine read_family_file

!> Filter out genes that are not assigned to any family.
subroutine filter_unassigned_genes(gene_ids, expression_vectors, gene_to_fam, n_genes_kept, ierr)
    character(len=*), allocatable, intent(inout) :: gene_ids(:)
        !! gene ids array
    real(real64), allocatable, intent(inout) :: expression_vectors(:,:)
        !! expression vectors array
    integer(int32), allocatable, intent(inout) :: gene_to_fam(:)
        !! gene to family mapping
    integer(int32), intent(out) :: n_genes_kept
        !! number of genes kept after filtering
    integer(int32), intent(out) :: ierr
        !! Error code
    
    integer(int32) :: i, j, n_genes_total, n_samples, ios
    integer(int32), allocatable :: valid_indices(:)
    character(len=len(gene_ids)), allocatable :: temp_gene_ids(:)
    real(real64), allocatable :: temp_expression_vectors(:,:)
    integer(int32), allocatable :: temp_gene_to_fam(:)
    
    call set_ok(ierr)
    call set_ok(ios)
    n_genes_total = size(gene_ids)
    n_samples = size(expression_vectors, 1)
    
    ! Count genes with valid family assignments
    n_genes_kept = count(gene_to_fam >= 1)
    
    if (n_genes_kept == 0) then
        call set_err_once(ierr, ERR_INVALID_INPUT)
        return
    end if
    
    ! Allocate temporary arrays for valid genes
    allocate(valid_indices(n_genes_kept), stat = ios)
    call check_io_stat(ios, ierr)
    allocate(temp_gene_ids(n_genes_kept), stat = ios)
    call check_io_stat(ios, ierr)
    allocate(temp_expression_vectors(n_samples, n_genes_kept), stat = ios)
    call check_io_stat(ios, ierr)
    allocate(temp_gene_to_fam(n_genes_kept), stat = ios)
    call check_io_stat(ios, ierr)
    if(is_err(ierr)) return
    
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
    
    if(DEBUG) write(*,*) 'Filtered out ', n_genes_total - n_genes_kept, ' unassigned genes'
    if(DEBUG) write(*,*) 'Kept ', n_genes_kept, ' genes with valid family assignments'
end subroutine filter_unassigned_genes

!> Helper subroutine to split strings
subroutine split_string(input, output, ierr, delimiter)
    character(len=*), intent(in) :: input
        !! Input string to split
    character(len=:), allocatable, intent(out) :: output(:)
        !! Output of parts
    character(len=1), optional, intent(in) :: delimiter
        !! delimiter, default is space
    integer(int32), intent(out) :: ierr
        !! Error code
    integer(int32) :: ios
    
    character(len=1) :: delim
    integer(int32) :: n, i, start_pos, end_pos, field_count
    integer(int32), allocatable :: field_starts(:), field_ends(:)
    
    call set_ok(ierr)
    call set_ok(ios)

    if (present(delimiter)) then
        delim = delimiter
    else
        delim = ' '
    end if

    ! First pass: find field boundaries
    allocate(field_starts(len_trim(input)), stat=ios)
    call check_io_stat(ios, ierr)
    allocate(field_ends(len_trim(input)), stat=ios)
    call check_io_stat(ios, ierr)
    if(is_err(ierr)) return
    
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
    allocate(character(len=maxval(field_ends(1:n) - field_starts(1:n) + 1)) :: output(n), stat=ios)
    call check_io_stat(ios, ierr)
    if(is_err(ierr)) return

    ! Extract fields
    do i = 1, n
        output(i) = trim(adjustl(input(field_starts(i):field_ends(i))))
    end do
    
end subroutine split_string

!> Helper to create a mask of genes that are unassigned
subroutine get_unassigned_mask(gene_to_fam, mask, n_genes_kept)
    use iso_fortran_env, only: int32
    implicit none
    
    integer(int32), intent(in) :: gene_to_fam(:)
        !! gene to family mapping
    logical, intent(out) :: mask(size(gene_to_fam))
        !! mask for mapping
    integer(int32), intent(out) :: n_genes_kept
        !! number of genes kept
    
    integer(int32) :: i
    
    n_genes_kept = 0
    do i = 1, size(gene_to_fam)
        mask(i) = (gene_to_fam(i) >= 1)
        if (mask(i)) n_genes_kept = n_genes_kept + 1
    end do
end subroutine get_unassigned_mask

end module tox_data_tools

!> R binding to read gene IDs from a file.
subroutine read_gene_ids_from_file_R(filename_ascii, fn_len, gene_ids_ascii, gene_ids_len, n_genes, &
                                 n_header_rows, gene_col, ierr)
    use iso_fortran_env, only: int32
    use tox_errors, only: set_ok, is_err
    use tox_data_tools, only: read_gene_ids_from_file
    use array_utils, only: ascii_to_string_padded, string_to_ascii
    implicit none
    integer(int32), intent(in) :: fn_len
        !! Length of the filename
    integer(int32), intent(in) :: filename_ascii(fn_len)
        !! Filename
    integer(int32), intent(in) :: gene_ids_len
        !! Length of the gene ids
    integer(int32), intent(in) :: n_genes
        !! Number of genes
    integer(int32), intent(inout) :: gene_ids_ascii(gene_ids_len, n_genes)
        !! Gene ids array
    integer(int32), intent(in) :: n_header_rows
        !! number of headers to skip
    integer(int32), intent(in) :: gene_col
        !! Column index that contains the gene ids
    integer(int32), intent(out) :: ierr
        !! Error code

    character(len=:), allocatable :: filename
    character(len=gene_ids_len) :: gene_ids(n_genes)
    integer(int32) :: i

    call set_ok(ierr)

    call ascii_to_string_padded(filename_ascii, fn_len, filename)
    call read_gene_ids_from_file(filename, gene_ids, n_header_rows, gene_col, ierr)

    if(is_err(ierr)) return

    do i = 1, n_genes
        call string_to_ascii(gene_ids(i), gene_ids_ascii(:, i))
    end do
end subroutine read_gene_ids_from_file_R

!> R binding to read expression vectors from files.
subroutine read_expression_vectors_R(file_list_ascii, file_list_len, n_files, &
                                 gene_ids_ascii, gene_ids_len, n_genes, &
                                 expression_vectors_flat, n_samples, &
                                 n_header_rows, gene_col, value_cols, &
                                 n_value_cols, ierr, delimiter_ascii, dlen)
    use iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, is_err, check_io_stat 
    use array_utils, only: ascii_to_string_padded, string_to_ascii
    use tox_data_tools, only: read_expression_vectors
    implicit none
    integer(int32), intent(in) :: file_list_len
        !! Length of the filenames
    integer(int32), intent(in) :: n_files
        !! Number of files
    integer(int32), intent(in) :: file_list_ascii(file_list_len, n_files)
        !! File list
    integer(int32), intent(in) :: gene_ids_len
        !! Length of the gene ids
    integer(int32), intent(in) :: n_genes
        !! Number of genes
    integer(int32), intent(inout) :: gene_ids_ascii(gene_ids_len, n_genes)
        !! Gene ids array
    integer(int32), intent(in) :: n_samples
        !! Number of samples    
    real(real64), intent(inout) :: expression_vectors_flat(n_samples * n_genes)
        !! Expression vectors array
    integer(int32), intent(in) :: n_header_rows
        !! Number of header rows
    integer(int32), intent(in) :: gene_col
        !! Index of the gene column
    integer(int32), intent(in) :: n_value_cols
        !! Number of value columns    
    integer(int32), intent(in) :: value_cols(n_value_cols)
        !! Indicies of columns with values
    integer(int32), intent(out) :: ierr
        !! Error code
    integer(int32), intent(in) :: dlen  
        !! Length of the delimiter
    integer(int32), intent(in) :: delimiter_ascii(dlen)
        !! delimiter

    character(len=file_list_len), allocatable :: file_list(:)
    character(len=gene_ids_len), allocatable :: gene_ids(:)
    character(len=:), allocatable :: delimiter
    character(len=:), allocatable :: tmp_str
    real(real64), allocatable :: expression_vectors(:,:)
    integer(int32) :: i, j, ios, start_row

    call set_ok(ierr)
    call set_ok(ios)

    start_row = 1
    allocate(file_list(n_files), stat = ios)
    call check_io_stat(ios, ierr)
    allocate(gene_ids(n_genes), stat = ios)
    call check_io_stat(ios, ierr)
    allocate(expression_vectors(n_samples, n_genes), stat = ios)
    call check_io_stat(ios, ierr)
    if(is_err(ierr)) return

    ! Convert 2D ASCII arrays to string arrays
    do i = 1, n_files
      call ascii_to_string_padded(file_list_ascii(:, i), file_list_len, tmp_str)
      file_list(i) = trim(tmp_str)
    end do
    
    do i = 1, n_genes
      call ascii_to_string_padded(gene_ids_ascii(:, i), gene_ids_len, tmp_str)
      gene_ids(i) = trim(tmp_str)
    end do

    if (dlen > 0) then
      call ascii_to_string_padded(delimiter_ascii, dlen, delimiter)
    else
      delimiter = char(9)
    end if

    call read_expression_vectors(file_list, gene_ids, expression_vectors, &
                                n_header_rows, gene_col, value_cols, &
                                start_row, ierr, delimiter)

    if(is_err(ierr)) return

    ! Flatten the 2D array back to 1D for R
    do j = 1, n_genes
      do i = 1, n_samples
        expression_vectors_flat((j-1)*n_samples + i) = expression_vectors(i, j)
      end do
    end do

    do i = 1, n_genes
    call string_to_ascii(gene_ids(i), gene_ids_ascii(:, i))
    end do
end subroutine read_expression_vectors_R

!> R Binding to read a family file
subroutine read_family_file_R(filename_ascii, fn_len, gene_ids_ascii, gene_ids_len, n_genes, &
                             family_ids_ascii, family_ids_len, n_families, gene_to_fam, ierr)
    use iso_fortran_env, only: int32
    use tox_errors, only: set_ok, is_err
    use array_utils, only: ascii_to_string, string_to_ascii
    use tox_data_tools, only: read_family_file
    implicit none
    
    integer(int32), intent(in) :: fn_len
        !! Length of the filename
    integer(int32), intent(in) :: filename_ascii(fn_len)
        !! Filename
    integer(int32), intent(in) :: gene_ids_len
        !! Length of the gene ids
    integer(int32), intent(in) :: n_genes
        !! Number of genes
    integer(int32), intent(in) :: gene_ids_ascii(gene_ids_len, n_genes)
        !! Gene ids array
    integer(int32), intent(in) :: family_ids_len
        !! Length of the family ids
    integer(int32), intent(in) :: n_families
        !! Number of families
    integer(int32), intent(out) :: family_ids_ascii(family_ids_len, n_families)
        !! Family ids
    integer(int32), intent(out) :: gene_to_fam(n_genes)
        !! Gene to family mapping
    integer(int32), intent(out) :: ierr
        !! Error code
    
    character(len=fn_len) :: filename
    character(len=gene_ids_len) :: gene_ids(n_genes)
    character(len=family_ids_len) :: family_ids(n_families)
    character(len=:), allocatable :: temp_str
    integer(int32) :: i

    call set_ok(ierr)
    
    ! Initialize family_ids with spaces
    family_ids = ' '
    
    ! Convert filename from ASCII
    call ascii_to_string(filename_ascii, fn_len, temp_str)
    filename = trim(adjustl(temp_str))
    
    ! Convert gene IDs from ASCII to strings
    do i = 1, n_genes
        call ascii_to_string(gene_ids_ascii(:, i), gene_ids_len, temp_str)
        gene_ids(i) = trim(adjustl(temp_str))
    end do
    
    call read_family_file(filename, gene_ids, family_ids, gene_to_fam, ierr)
    if(is_err(ierr)) return
    
    ! Convert family IDs to ASCII
    do i = 1, n_families
        call string_to_ascii(trim(adjustl(family_ids(i))), family_ids_ascii(:, i))
    end do
end subroutine read_family_file_R

!> R binding to filter unassigned genes
subroutine filter_unassigned_genes_R(gene_ids_ascii, gene_ids_len, n_genes, &
                                    expression_vectors_flat, n_samples, &
                                    gene_to_fam, mask, n_genes_kept, ierr)
    use iso_fortran_env, only: real64, int32
    use tox_errors, only: set_ok, set_err_once, ERR_INVALID_INPUT
    use array_utils, only: ascii_to_string, string_to_ascii
    use tox_data_tools, only: get_unassigned_mask
    implicit none
    
    integer(int32), intent(in) :: gene_ids_len
        !! Length of the gene ids
    integer(int32), intent(in) :: n_genes
        !! Number of genes
    integer(int32), intent(in) :: gene_ids_ascii(gene_ids_len, n_genes)
        !! Gene ids array
    integer(int32), intent(in) :: n_samples
        !! Number of samples
    real(real64), intent(in) :: expression_vectors_flat(n_samples * n_genes)
        !! Expression vectors
    integer(int32), intent(in) :: gene_to_fam(n_genes)
        !! gene to family mapping
    logical, intent(out) :: mask(n_genes)
        !! mask for unassigned genes
    integer(int32), intent(out) :: n_genes_kept
        !! Number of genes that are kept
    integer(int32), intent(out) :: ierr
        !! Error code
    
    call set_ok(ierr)

    call get_unassigned_mask(gene_to_fam, mask, n_genes_kept)
    if (n_genes_kept == 0) then
        call set_err_once(ierr, ERR_INVALID_INPUT)
        return
    end if
end subroutine filter_unassigned_genes_R

!> C binding for reading gene IDs from a file
subroutine read_gene_ids_from_file_C(filename_ascii, fn_len, gene_ids_ascii, gene_ids_len, n_genes, &
                                 n_header_rows, gene_col, ierr) bind(C, name="read_gene_ids_from_file_C")
    use iso_c_binding, only: c_int
    use tox_errors, only: set_ok, is_err
    use array_utils, only: ascii_to_string_padded, string_to_ascii
    use tox_data_tools, only: read_gene_ids_from_file
    implicit none

    integer(c_int), intent(in), value :: fn_len       
        !! Length of filename
    integer(c_int), intent(in) :: filename_ascii(fn_len)
        !! Pointer to filename array
    integer(c_int), intent(in), value :: gene_ids_len 
        !! Length of each gene ID string
    integer(c_int), intent(in), value :: n_genes      
        !! Number of genes    
    integer(c_int), intent(out) :: gene_ids_ascii(gene_ids_len, n_genes)
        !! Pointer to gene_ids array
    integer(c_int), intent(in), value :: n_header_rows
        !! Number of header rows to skip
    integer(c_int), intent(in), value :: gene_col
        !! Index of the gene column
    integer(c_int), intent(out) :: ierr
        !! Error code

    character(len=:), allocatable :: filename
    character(len=gene_ids_len) :: gene_ids(n_genes)
    integer(c_int) :: i

    call set_ok(ierr)

    call ascii_to_string_padded(filename_ascii, fn_len, filename)

    call read_gene_ids_from_file(filename, gene_ids, n_header_rows, gene_col, ierr)
    if(is_err(ierr)) return

    do i = 1, n_genes
        call string_to_ascii(gene_ids(i), gene_ids_ascii(:, i))
    end do
end subroutine read_gene_ids_from_file_C

!> C binding for reading expression vectors from files
subroutine read_expression_vectors_C(file_list_ascii, file_list_len, n_files, &
                                 gene_ids_ascii, gene_ids_len, n_genes, &
                                 expression_vectors, n_samples, &
                                 n_header_rows, gene_col, value_cols, &
                                 n_value_cols, ierr, delimiter_ascii, dlen) bind(C, name="read_expression_vectors_C")
    use iso_c_binding, only: c_int, c_double
    use tox_data_tools, only: read_expression_vectors
    use tox_errors, only: set_ok, is_err, check_io_stat
    use array_utils, only: ascii_to_string_padded, string_to_ascii
    implicit none

    integer(c_int), intent(in), value :: file_list_len     
        !! Length of each filename
    integer(c_int), intent(in), value :: n_files           
        !! Number of files
    integer(c_int), intent(in) :: file_list_ascii(file_list_len, n_files)    
        !! Pointer to file_list array
    integer(c_int), intent(in), value :: gene_ids_len      
        !! Length of each gene ID
    integer(c_int), intent(in), value :: n_genes           
        !! Number of genes
    integer(c_int), intent(in) :: gene_ids_ascii(gene_ids_len, n_genes)     
        !! Pointer to gene_ids array
    integer(c_int), intent(in), value :: n_samples
        !! Number of samples    
    real(c_double), intent(out) :: expression_vectors(n_samples, n_genes)
        !! Pointer to expression vectors (flat array)
    integer(c_int), intent(in), value :: n_header_rows
        !! Number of header rows
    integer(c_int), intent(in), value :: gene_col
        !! Index of column containing gene ids
    integer(c_int), intent(in), value :: n_value_cols
        !! Number of cols containing values    
    integer(c_int), intent(in) :: value_cols (n_value_cols)          
        !! Pointer to value_cols array
    integer(c_int), intent(out) :: ierr
        !! Error code
    integer(c_int), intent(in), value :: dlen   
        !! Length of the delimiter
    integer(c_int), intent(in) :: delimiter_ascii(dlen)
        !! Delimiter
    integer(c_int) :: start_row
    
    character(len=file_list_len), allocatable :: file_list(:)
    character(len=gene_ids_len), allocatable :: gene_ids(:)
    character(len=:), allocatable :: delimiter
    character(len=:), allocatable :: tmp_str
    integer(c_int) :: i, j, ios

    start_row = 1
    call set_ok(ierr)

    allocate(file_list(n_files), stat=ios)
    call check_io_stat(ios, ierr)
    allocate(gene_ids(n_genes), stat=ios)
    call check_io_stat(ios, ierr)
    if(is_err(ierr)) return
   
    ! Convert 2D ASCII arrays to string arrays
    do i = 1, n_files
      call ascii_to_string_padded(file_list_ascii(:, i), file_list_len, tmp_str)
      file_list(i) = trim(tmp_str)
    end do
    
    do i = 1, n_genes
      call ascii_to_string_padded(gene_ids_ascii(:, i), gene_ids_len, tmp_str)
      gene_ids(i) = trim(tmp_str)
    end do

    if (dlen > 0) then
      call ascii_to_string_padded(delimiter_ascii, dlen, delimiter)
    else
      delimiter = char(9)
    end if

    call read_expression_vectors(file_list, gene_ids, expression_vectors, &
                                n_header_rows, gene_col, value_cols, &
                                start_row, ierr, delimiter)

    if(is_err(ierr)) return

end subroutine read_expression_vectors_C

!> C binding for reading family file
subroutine read_family_file_C(filename_ascii, fn_len, gene_ids_ascii, gene_ids_len, n_genes, &
                             family_ids_ascii, family_ids_len, n_families, gene_to_fam, ierr) bind(C, name="read_family_file_C")
    use iso_c_binding, only: c_int
    use tox_errors, only: set_ok, is_err
    use tox_data_tools, only: read_family_file
    use array_utils, only: ascii_to_string_padded, string_to_ascii
    implicit none
    integer(c_int), intent(in), value :: fn_len            
        !! Length of filename
    integer(c_int), intent(in) :: filename_ascii(fn_len)       
        !! Pointer to filename array
    integer(c_int), intent(in), value :: gene_ids_len      
        !! Length of each gene ID
    integer(c_int), intent(in), value :: n_genes           
        !! Number of genes    
    integer(c_int), intent(in) :: gene_ids_ascii(gene_ids_len, n_genes)       
        !! Pointer to gene_ids array
    integer(c_int), intent(in), value :: family_ids_len    
        !! Length of each family ID
    integer(c_int), intent(in), value :: n_families        
        !! Number of families
    integer(c_int), intent(out) :: family_ids_ascii(family_ids_len, n_families)     
        !! Pointer to family_ids array
    integer(c_int), intent(out) :: gene_to_fam(n_genes)          
        !! Pointer to gene_to_fam array
    integer(c_int), intent(out) :: ierr
        !! error code
    
    character(len=fn_len) :: filename
    character(len=gene_ids_len) :: gene_ids(n_genes)
    character(len=family_ids_len) :: family_ids(n_families)
    character(len=:), allocatable :: temp_str
    integer(c_int) :: i

    call set_ok(ierr)

    ! Initialize family_ids with spaces
    family_ids = ' '
    
    ! Convert filename from ASCII
    call ascii_to_string_padded(filename_ascii, fn_len, temp_str)
    filename = trim(adjustl(temp_str))
    
    ! Convert gene IDs from ASCII to strings
    do i = 1, n_genes
        call ascii_to_string_padded(gene_ids_ascii(:, i), gene_ids_len, temp_str)
        gene_ids(i) = trim(adjustl(temp_str))
    end do
    
    call read_family_file(filename, gene_ids, family_ids, gene_to_fam, ierr)
    if(is_err(ierr)) return
    
    ! Convert family IDs to ASCII
    do i = 1, n_families
        call string_to_ascii(trim(adjustl(family_ids(i))), family_ids_ascii(:, i))
    end do
end subroutine read_family_file_C

!> C binding for filtering unassigned genes
subroutine filter_unassigned_genes_C(gene_ids_ascii, gene_ids_len, n_genes, &
                                    gene_to_fam, mask, n_genes_kept, ierr) bind(C, name="filter_unassigned_genes_C")
    use iso_c_binding, only: c_int, c_double
    use iso_fortran_env, only: int32
    use tox_errors, only: set_ok, set_err_once, ERR_INVALID_INPUT, is_err, check_io_stat
    use tox_conversions, only: c_int_as_logical, logical_as_c_int
    use tox_data_tools, only: get_unassigned_mask
    implicit none

    integer(c_int), value :: gene_ids_len      
        !! Length of each gene ID
    integer(c_int), value :: n_genes           
        !! Number of genes
    integer(c_int), intent(in) :: gene_ids_ascii(gene_ids_len, n_genes)       
        !! Pointer to gene_ids array
    integer(c_int), intent(in):: gene_to_fam(n_genes)          
        !! Pointer to gene_to_fam array
    integer(c_int), intent(out) :: mask(n_genes)                 
        !! Pointer to mask array
    integer(c_int), intent(out) :: n_genes_kept
        !! number of genes kept
    integer(c_int), intent(out) :: ierr
        !! Error code
    
    integer(int32) :: i, ios
    logical, allocatable :: mask_logical(:)

    call set_ok(ierr)
    call set_ok(ios)

    allocate(mask_logical(n_genes), stat = ios)
    call check_io_stat(ios, ierr)
    if(is_err(ierr)) return

    call get_unassigned_mask(gene_to_fam, mask_logical, n_genes_kept)

    if(n_genes_kept == 0) then
        call set_err_once(ierr, ERR_INVALID_INPUT)
        return
    end if

    do i = 1, n_genes
        call logical_as_c_int(mask_logical(i), mask(i))
    end do
end subroutine filter_unassigned_genes_C