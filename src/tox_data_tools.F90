! ToDos: 
! Update error handling
! implement tests

module tox_data_tools
    use iso_fortran_env, only: real64, int32
    use tox_errors
    implicit none
    private

    public :: read_gene_ids_from_file
    public :: read_tabular_files
    public :: read_family_file
    public :: get_gene_index
    public :: get_family_index
    public :: get_family_for_gene_index
    public :: get_expression_vector
    public :: get_family_centroid
    public :: get_shift_components
    PUBLIC :: split_string

    character(len=*), parameter :: DELIMS = ', ' // char(9)

contains

! Read tabular files (CSV/TSV)
! Read tabular files (CSV/TSV)
subroutine read_tabular_files(file_list, gene_ids, expression_vectors, &
                             n_header_rows, gene_col, value_cols, start_row, ierr)

    character(len=*), intent(in) :: file_list(:)
    character(len=*), intent(inout) :: gene_ids(:)
    real(real64), intent(inout) :: expression_vectors(:,:)
    integer(int32), intent(in) :: n_header_rows, gene_col
    integer(int32), intent(in) :: value_cols(:) ! Array of column indices [cite: 19, 20]
    integer(int32), intent(in) :: start_row ! New parameter to specify the start row
    integer(int32), intent(out) :: ierr

    integer :: i, j, k, unit, ios, idx, row_count, n_genes, expected_idx, n_value_cols
    integer :: current_sample
    character(len=2048) :: line
    character(len=:), allocatable :: fields(:)
    real(real64) :: value
    character(len=len(gene_ids)) :: gene
    logical :: order_mismatch

    ierr = 0
    n_genes = size(gene_ids)
    n_value_cols = size(value_cols)
    
    current_sample = start_row - 1
    
    do i = 1, size(file_list)
        open(newunit=unit, file=trim(file_list(i)), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            ierr = 1
            return
        end if

        ! Skip header rows [cite: 23, 42]
        do j = 1, n_header_rows
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
        end do

        row_count = 0
        order_mismatch = .false.
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            row_count = row_count + 1
            
            if (row_count > n_genes) then
                ierr = 2 ! More genes than allocated space [cite: 24, 25]
                exit
            end if

            call split_string(line, fields, char(9))
            if (size(fields) < max(gene_col, maxval(value_cols))) cycle

            gene = trim(adjustl(fields(gene_col)))

            ! Fast path: check if gene matches expected position [cite: 30, 31]
            expected_idx = row_count
            if (expected_idx <= n_genes .and. trim(adjustl(gene_ids(expected_idx))) == trim(adjustl(gene))) then
                idx = expected_idx
            else
                ! Slow path: search through all genes [cite: 32]
                idx = get_gene_index(gene_ids, gene)
                if (idx == 0) then
                    write(*,*) 'Warning: Gene ', trim(gene), ' not found in master gene list'
                    order_mismatch = .true.
                    cycle
                end if
            end if
            
            ! Read all value columns for this gene [cite: 35]
            do k = 1, n_value_cols
                read(fields(value_cols(k)), *, iostat=ios) value
                if (ios == 0) then
                    expression_vectors(current_sample + k, idx) = value
                end if
            end do
        end do
        close(unit)
        
        ! Update sample offset for next file [cite: 37]
        current_sample = current_sample + n_value_cols
        
        ! Report if order mismatch was detected [cite: 38]
        if (order_mismatch) then
            write(*,*) 'Note: Gene order mismatch detected in file: ', trim(file_list(i))
        end if
    end do
end subroutine read_tabular_files

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
    character(len=*), intent(in) :: gene_ids(:)  ! Bereits vorhandene Gene-IDs
    character(len=*), intent(out) :: family_ids(:)
    integer(int32), intent(out) :: gene_to_fam(:)
    integer(int32), intent(out) :: ierr

    integer(int32) :: unit, ios, i, j, fam_idx, gene_idx, n_families, n_genes
    character(len=1024) :: line
    character(len=:), allocatable :: fields(:), genes(:)
    character(len=len(family_ids)) :: current_family
    type(hashmap_type) :: gene_map
    logical :: fill_family_ids

    ierr = 0
    gene_to_fam = 0
    n_families = size(family_ids)
    n_genes = size(gene_ids)
    
    fill_family_ids = all(family_ids == "")

    ! Hashmap für schnelle Gen-Suche erstellen
    call hashmap_create(gene_map, n_genes)
    
    ! Hashmap mit allen Gen-IDs füllen
    write(*,*) 'Building gene hashmap with ', n_genes, ' genes...'
    do i = 1, n_genes
        call hashmap_put(gene_map, gene_ids(i), i)
        if (mod(i, 10000) == 0) then
            write(*,*) '  Processed ', i, ' genes...'
        end if
    end do
    write(*,*) 'Hashmap complete.'

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

    write(*,*) 'Reading family file with hashmap...'
    
    fam_idx = 0
    do
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) exit

        fam_idx = fam_idx + 1
        if (mod(fam_idx, 1000) == 0) then
            write(*,*) 'Processed ', fam_idx, ' families...'
        end if

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
    write(*,*) 'Finished processing ', fam_idx, ' families'
end subroutine read_family_file

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
    integer :: n, i, start_pos, field_count
    logical :: in_field
    
    if (present(delimiter)) then
        delim = delimiter
    else
        delim = ' '
    end if

    ! Count fields
    n = 0
    in_field = .false.
    do i = 1, len_trim(input)
        if (input(i:i) == delim) then
            in_field = .false.
        else if (.not. in_field) then
            n = n + 1
            in_field = .true.
        end if
    end do

    allocate(character(len=len(input)) :: output(n))
    
    field_count = 0
    start_pos = 1
    do i = 1, len_trim(input)
        if (input(i:i) == delim) then
            if (start_pos <= i-1) then
                field_count = field_count + 1
                output(field_count) = trim(adjustl(input(start_pos:i-1)))
            end if
            start_pos = i + 1
        end if
    end do
    
    ! Add last field
    if (start_pos <= len_trim(input)) then
        field_count = field_count + 1
        output(field_count) = trim(adjustl(input(start_pos:)))
    end if
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



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! READ/WRITE UTILS -> ignore if not needed
! subroutine read_int_1d(filename, array, ierr)
!     character(len=*), intent(in)  :: filename
!     integer(int32), intent(out) :: array(:)
!     integer(int32), intent(out)   :: ierr

!     character(len=:), pointer :: tmp(:)

!     ! Call deserialize to get the pointer view
!     call deserialize_int_1d(tmp, filename, ierr)
!     if (.not. is_ok(ierr)) return

!     ! Check that the size matches
!     if (size(array) /= size(tmp)) then
!         ierr = ERR_DIM_MISMATCH
!         if (associated(tmp)) nullify(tmp)
!         return
!     end if

!     ! Copy data into user’s array
!     array = tmp

!     ! Clean up pointer
!     if (associated(tmp)) nullify(tmp)
! end subroutine

! subroutine read_char_1d(filename, array, ierr)
!     character(len=*), intent(in)  :: filename
!     character(len=*), intent(out) :: array(:)
!     integer(int32), intent(out)   :: ierr

!     character(len=:), pointer :: tmp(:)

!     ! Call deserialize to get the pointer view
!     call deserialize_char_1d(tmp, filename, ierr)
!     if (.not. is_ok(ierr)) return

!     ! Check that the size matches
!     if (size(array) /= size(tmp)) then
!         ierr = ERR_DIM_MISMATCH
!         if (associated(tmp)) nullify(tmp)
!         return
!     end if

!     ! Copy data into user’s array
!     array = tmp

!     ! Clean up pointer
!     if (associated(tmp)) nullify(tmp)

! end subroutine

! subroutine read_real_2d(filename, array, ierr)
!     character(len=*), intent(in) :: filename
!     real(real64), intent(out)    :: array(:,:)
!     integer(int32), intent(out)  :: ierr

!     real(real64), pointer :: tmp(:,:)

!     ! First, deserialize into a pointer
!     call deserialize_real_2d(tmp, filename, ierr)
!     if (.not. is_ok(ierr)) return

!     ! Check dimensions match user array
!     if (any(shape(array) /= shape(tmp))) then
!         call set_err_once(ierr, ERR_DIM_MISMATCH)
!         return
!     end if

!     ! Copy into user’s array
!     array = tmp

!     ! Clean up
!     if (associated(tmp)) deallocate(tmp)
! end subroutine

! subroutine save_gene_ids(filename, gene_ids, ierr)
!     CHARACTER(len=*), INTENT(IN) :: gene_ids(:)
!     CHARACTER(len=*), INTENT(IN) :: filename
!     integer, intent(out) :: ierr

!     call serialize_char_1D(filename, gene_ids, ierr)
! end subroutine

! subroutine save_expression_vectors(filename, expression_vectors, ierr)
!     CHARACTER(len=*), INTENT(IN) :: filename
!     real(real64), INTENT(IN) :: expression_vectors(:,:)
!     integer, intent(out) :: ierr

!     call serialize_real_2D(filename, expression_vectors, ierr)
! end subroutine

! subroutine save_gene_to_family(filename, gene_to_fam, ierr)
!     CHARACTER(len=*), INTENT(IN) :: filename
!     integer, intent(in) :: gene_to_fam(:)
!     integer, intent(out) :: ierr

!     call serialize_int_1D(filename, gene_to_fam, ierr)
! end subroutine

! subroutine save_gene_family_ids(filename, gene_family_ids, ierr)
!     CHARACTER(len=*), INTENT(IN) :: filename
!     CHARACTER(len=*), INTENT(IN) :: gene_family_ids(:)
!     integer(int32), INTENT(OUT) :: ierr

!     call serialize_char_1D(filename, gene_family_ids, ierr)
! end subroutine

! subroutine save_family_centroids(filename, family_centroids, ierr)
!     CHARACTER(len=*), INTENT(IN) :: filename
!     real(real64), INTENT(IN) :: family_centroids(:,:)
!     integer(int32), INTENT(OUT) :: ierr

!     call serialize_real_2D(filename, family_centroids, ierr)
! end subroutine

! subroutine save_shift_vectors(filename, shift_vectors, ierr)
!     CHARACTER(len=*), INTENT(IN) :: filename
!     real(real64), INTENT(IN) :: shift_vectors(:,:)
!     integer(int32), INTENT(OUT) :: ierr

!     call serialize_real_2D(filename, shift_vectors, ierr)
! end subroutine

! subroutine load_gene_ids(filename, gene_ids, ierr)
!     character(len=*), intent(in)  :: filename
!     character(len=*), intent(out) :: gene_ids(:)
!     integer(int32), intent(out)   :: ierr

!     call read_char_1d(filename, gene_ids, ierr)
! end subroutine load_gene_ids

! subroutine load_expression_vectors(filename, expression_vectors, ierr)
!     character(len=*), intent(in) :: filename
!     real(real64), intent(out)    :: expression_vectors(:,:)
!     integer(int32), intent(out)  :: ierr

!     call read_real_2d(filename, expression_vectors, ierr)
! end subroutine load_expression_vectors

! subroutine load_gene_to_family(filename, gene_to_fam, ierr)
!     character(len=*), INTENT(IN) :: filename
!     integer(int32), INTENT(OUT) :: gene_to_fam
!     integer(int32), INTENT(OUT) :: ierr

!     call read_int_1d(filename, gene_to_fam, ierr)
! end subroutine

! subroutine load_family_ids(filename, family_ids, ierr)
!     character(len=*), intent(in) :: filename
!     CHARACTER(len=*), intent(out) :: family_ids
!     integer(int32), intent(out) :: ierr

!     call read_char_1d(filename, family_ids, ierr)
! end subroutine

! subroutine load_family_centroids(filename, family_centroids, ierr)
!     character(len=*), intent(in) :: filename
!     real(real64), intent(out)    :: family_centroids(:,:)
!     integer(int32), intent(out)  :: ierr

!     call read_real_2d(filename, family_centroids, ierr)
! end subroutine

! subroutine load_shift_vectors(filename, shift_vectors, ierr)
!     character(len=*), intent(in) :: filename
!     real(real64), intent(out)    :: shift_vectors(:,:)
!     integer(int32), intent(out)  :: ierr

!     call read_real_2d(filename, shift_vectors, ierr)
! end subroutine

end module tox_data_tools