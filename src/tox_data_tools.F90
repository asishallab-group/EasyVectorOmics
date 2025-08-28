! ToDos: 
! Update error handling
! implement tests

module tox_data_tools
    use iso_fortran_env, only: real64, int32
    use tox_errors
    implicit none
    private

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
subroutine read_tabular_files(file_list, gene_ids, expression_vectors, &
                             n_header_rows, gene_col, value_cols, ierr)
    character(len=*), intent(in) :: file_list(:)
    character(len=*), intent(inout) :: gene_ids(:)
    real(real64), intent(out) :: expression_vectors(:,:)
    integer(int32), intent(in) :: n_header_rows, gene_col
    integer(int32), intent(in) :: value_cols(:)  ! Array of column indices
    integer(int32), intent(out) :: ierr

    integer :: i, j, k, unit, ios, idx, row_count, n_genes, expected_idx, n_value_cols
    integer :: sample_offset, current_sample
    character(len=1024) :: line
    character(len=:), allocatable :: fields(:)
    real(real64) :: value
    character(len=len(gene_ids)) :: gene
    logical :: fill_gene_ids, order_mismatch

    ierr = 0
    expression_vectors = 0.0_real64
    n_genes = size(gene_ids)
    n_value_cols = size(value_cols)
    
    ! Check if we need to fill gene_ids (only for first file)
    fill_gene_ids = all(gene_ids == "")
    
    ! Initialize sample offset
    sample_offset = 0

    do i = 1, size(file_list)
        open(newunit=unit, file=trim(file_list(i)), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            ierr = 1
            return
        end if

        ! Skip header rows
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
                ierr = 2  ! More genes than allocated space
                exit
            end if

            call split_string(line, fields)
            if (size(fields) < max(gene_col, maxval(value_cols))) cycle

            gene = trim(adjustl(fields(gene_col)))

            if (fill_gene_ids .and. i == 1) then
                ! Fill gene_ids from first file
                gene_ids(row_count) = gene
                
                ! Read all value columns for this gene
                do k = 1, n_value_cols
                    read(fields(value_cols(k)), *, iostat=ios) value
                    if (ios == 0) then
                        expression_vectors(sample_offset + k, row_count) = value
                    end if
                end do
            else
                ! Hybrid approach: first check if gene matches expected position
                expected_idx = row_count
                
                if (expected_idx <= n_genes .and. trim(adjustl(gene_ids(expected_idx))) == trim(adjustl(gene))) then
                    ! Fast path: gene is in expected position
                    idx = expected_idx
                else
                    ! Slow path: search through all genes
                    idx = get_gene_index(gene_ids, gene)
                    if (idx == 0) then
                        ! Gene not found - this shouldn't happen if files have same genes
                        write(*,*) 'Warning: Gene ', trim(gene), ' not found in master gene list'
                        order_mismatch = .true.
                        cycle
                    end if
                end if
                
                ! Read all value columns for this gene
                do k = 1, n_value_cols
                    read(fields(value_cols(k)), *, iostat=ios) value
                    if (ios == 0) then
                        expression_vectors(sample_offset + k, idx) = value
                    end if
                end do
            end if
        end do
        close(unit)
        
        ! Update sample offset for next file
        sample_offset = sample_offset + n_value_cols
        
        ! Report if order mismatch was detected
        if (order_mismatch) then
            write(*,*) 'Note: Gene order mismatch detected in file: ', trim(file_list(i))
        end if
        
        ! After first file, we no longer need to fill gene_ids
        if (i == 1) fill_gene_ids = .false.
    end do
end subroutine read_tabular_files

! Read family file (OrthoFinder format)
subroutine read_family_file(filename, gene_ids, family_ids, gene_to_fam, ierr)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: gene_ids(:)
    character(len=*), intent(inout) :: family_ids(:)
    integer(int32), intent(out) :: gene_to_fam(:)
    integer(int32), intent(out) :: ierr

    integer :: unit, ios, i, j, fam_idx, gene_idx, n_families
    character(len=1024) :: line
    character(len=:), allocatable :: fields(:), genes(:)
    character(len=len(family_ids)) :: current_family
    logical :: fill_family_ids

    ierr = 0
    gene_to_fam = 0
    n_families = size(family_ids)
    
    ! Check if we need to fill family_ids
    fill_family_ids = all(family_ids == "")

    open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
        ierr = 1
        return
    end if

    fam_idx = 0
    do
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) exit

        call split_string(line, fields, char(9))
        if (size(fields) < 2) cycle

        fam_idx = fam_idx + 1
        if (fill_family_ids .and. fam_idx > n_families) then
            ierr = 2  ! More families than allocated space
            exit
        end if

        current_family = trim(adjustl(fields(1)))
        if (fill_family_ids) then
            family_ids(fam_idx) = current_family
        else
            ! Check if family exists in pre-filled family_ids
            if (get_family_index(family_ids, current_family) == 0) then
                ierr = 3  ! Family not found in pre-filled list
                exit
            end if
        end if

        ! Process all gene columns
        do i = 2, size(fields)
            call split_string(fields(i), genes, ',')
            do j = 1, size(genes)
                if (trim(genes(j)) == "") cycle
                
                ! Find gene index in the gene_ids array
                gene_idx = get_gene_index(gene_ids, trim(adjustl(genes(j))))
                if (gene_idx > 0) then
                    if (fill_family_ids) then
                        gene_to_fam(gene_idx) = fam_idx
                    else
                        gene_to_fam(gene_idx) = get_family_index(family_ids, current_family)
                    end if
                end if
            end do
        end do
    end do
    close(unit)
end subroutine read_family_file

! Helper function to find gene index
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