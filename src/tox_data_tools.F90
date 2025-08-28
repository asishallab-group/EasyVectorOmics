! ToDos: 
! Update error handling
! implement tests

module tox_data_tools
    use iso_fortran_env, only: real64, int32
    use tox_errors, only: ERR_OK, ERR_INVALID_INPUT, ERR_SIZE_MISMATCH, ERR_FILE_OPEN, ERR_DIM_MISMATCH, set_ok, set_err_once
    implicit none
    private

    ! Public procedures
    public :: read_kallisto_files
    public :: read_bowtie_files
    public :: get_gene_index
    public :: get_family_index
    public :: get_family_for_gene_index
    public :: get_expression_vector
    public :: get_family_centroid
    public :: get_shift_components

contains

subroutine read_tsv_files(file_list, n_files, gene_ids, n_genes, &
                          expression_vectors, n_header_rows, &
                          gene_col, value_col, ierr)
    use iso_fortran_env, only: int32, real64
    implicit none
    character(len=*), intent(in) :: file_list(:)
    integer(int32), intent(in) :: n_files, n_genes
    character(len=*), intent(inout) :: gene_ids(n_genes)   ! may be empty on entry
    real(real64), intent(out) :: expression_vectors(n_files, n_genes)
    integer(int32), intent(in) :: n_header_rows, gene_col, value_col
    integer(int32), intent(out) :: ierr

    integer :: i, j, ios, unit, idx, row_number, lenline
    character(len=1024) :: line
    character(len=256) :: fields(200)
    integer :: nfields
    character(len=256) :: gene
    real(real64) :: value
    logical :: gene_ids_empty

    ierr = 0
    expression_vectors = 0.0d0

    ! Check if gene_ids are empty (= not initialized by caller)
    gene_ids_empty = all([ (trim(gene_ids(j)) == '', j=1,n_genes) ])

    do i = 1, n_files
        open(newunit=unit, file=file_list(i), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,*) 'Error: cannot open file ', trim(file_list(i))
            ierr = 1
            return
        end if

        ! Skip header rows
        do j = 1, n_header_rows
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
        end do

        row_number = 1
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit

            call split_line_tsv(line, fields, nfields)
            if (nfields < max(gene_col, value_col)) cycle

            gene  = trim(fields(gene_col))
            read(fields(value_col), *, iostat=ios) value
            if (ios /= 0) value = 0.0d0

            if (i == 1 .and. gene_ids_empty) then
                ! First file, first fill gene_ids
                if (row_number <= n_genes) gene_ids(row_number) = gene
                idx = row_number
            else
                ! Otherwise look up gene in existing IDs
                idx = 0
                do j = 1, n_genes
                    if (trim(gene_ids(j)) == gene) then
                        idx = j
                        exit
                    end if
                end do
                if (idx == 0) then
                    write(*,*) 'Warning: gene ', trim(gene), ' not found in provided gene_ids'
                end if
            end if

            if (idx > 0 .and. idx <= n_genes) then
                expression_vectors(i, idx) = value
            end if

            row_number = row_number + 1
        end do

        close(unit)
    end do
end subroutine read_tsv_files

subroutine parse_tsv_line(line, gene_col, tpm_col, gene, tpm)
    use iso_fortran_env, only: real64
    implicit none
    character(len=*), intent(in) :: line
    integer, intent(in)          :: gene_col, tpm_col
    character(len=*), intent(out):: gene
    real(real64), intent(out)    :: tpm

    character(len=256) :: col(100)
    integer :: ncol, i, start, lenline

    lenline = len_trim(line)
    ncol = 0
    start = 1

    ! Split into tab-delimited fields
    do i = 1, lenline
        if (line(i:i) == char(9)) then
            ncol = ncol + 1
            if (ncol <= size(col)) col(ncol) = adjustl(line(start:i-1))
            start = i + 1
        end if
    end do
    ! Add last field
    if (start <= lenline) then
        ncol = ncol + 1
        if (ncol <= size(col)) col(ncol) = adjustl(line(start:lenline))
    end if

    if (ncol < max(gene_col, tpm_col)) then
        gene = ''
        tpm  = 0.0d0
        return
    end if

    gene = trim(col(gene_col))
    read(col(tpm_col), *, err=100) tpm
    return

100 continue
    write(*,*) 'Warning: failed to parse TPM from line: ', trim(line)
    tpm = 0.0d0
end subroutine parse_tsv_line

subroutine read_gene_ids(filename, gene_ids, n_genes, ierr)
    use iso_fortran_env, only: int32
    implicit none
    character(len=*), intent(in) :: filename
    character(len=*), intent(out) :: gene_ids(:)
    integer(int32), intent(out) :: n_genes, ierr

    integer :: unit, ios, count
    character(len=1024) :: line

    ierr = 0
    n_genes = 0

    open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
        ierr = 1
        return
    end if

    count = 0
    do
        read(unit,'(A)',iostat=ios) line
        if (ios /= 0) exit
        count = count + 1
        if (count <= size(gene_ids)) gene_ids(count) = trim(line)
    end do
    n_genes = count

    close(unit)
end subroutine read_gene_ids

subroutine read_family_file(filename, gene_ids, n_genes, &
                            gene_family_ids, n_families, gene_to_fam, ierr)
    implicit none
    character(len=*), intent(in) :: filename
    character(len=:), allocatable, intent(inout) :: gene_ids(:)
    integer(int32), intent(out) :: n_genes
    character(len=:), allocatable, intent(out) :: gene_family_ids(:)
    integer(int32), intent(out) :: n_families
    integer(int32), allocatable, intent(out) :: gene_to_fam(:)
    integer(int32), intent(out) :: ierr

    integer :: ios, unit, row_number, ncols, i, j
    character(len=1024) :: line
    character(len=256), allocatable :: cols(:)
    integer :: idx_gene, idx_fam, n_genes_current, n_fams_current

    ierr = 0
    n_genes = 0
    n_families = 0
    gene_to_fam = [ ]  ! clear array

    ! open file
    open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: cannot open file ', trim(filename)
        ierr = 1
        return
    end if

    ! initialize family array
    allocate(gene_family_ids(0))

    do
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (trim(line) == '') cycle

        ! split line by tabs
        call split_line_tab(line, cols, ncols)
        if (ncols < 2) cycle

        ! process family name
        idx_fam = 0
        do i = 1, n_families
            if (trim(gene_family_ids(i)) == trim(cols(1))) then
                idx_fam = i
                exit
            end if
        end do
        if (idx_fam == 0) then
            ! new family
            n_families = n_families + 1
            call array_append(gene_family_ids, trim(cols(1)))
            idx_fam = n_families
        end if

        ! process genes
        do i = 2, ncols
            ! check if gene_ids already populated
            if (size(gene_ids) > 0) then
                ! lookup gene in existing master list
                idx_gene = 0
                do j = 1, n_genes
                    if (trim(gene_ids(j)) == trim(cols(i))) then
                        idx_gene = j
                        exit
                    end if
                end do
                if (idx_gene == 0) then
                    write(*,*) 'Warning: gene ', trim(cols(i)), ' not in master gene_ids'
                    cycle
                end if
            else
                ! master gene_ids empty: append
                idx_gene = size(gene_ids) + 1
                call array_append(gene_ids, trim(cols(i)))
                n_genes = idx_gene
            end if

            ! append mapping
            call array_append(gene_to_fam, idx_fam)
        end do
    end do

    close(unit)

    ! final number of genes
    if (size(gene_ids) > 0) n_genes = size(gene_ids)

end subroutine read_family_file

!---------------------------
subroutine parse_gene_family_line(line, fam_name, gene_list, ierr)
    use, intrinsic :: iso_fortran_env, only: int32
    implicit none
    character(len=*), intent(in) :: line
    character(len=:), allocatable, intent(out) :: fam_name
    character(len=:), allocatable :: gene_list(:)
    integer(int32), intent(out) :: ierr

    character(len=1024) :: family_line
    character(len=256), allocatable :: tokens(:)
    integer :: n_tokens, i
    character(len=:), allocatable :: gene_line
    character(len=256) :: token
    integer :: pos

    call set_ok(ierr)
    family_line = trim(line)

    ! Split first column (family name) from the rest
    pos = index(family_line, char(9))  ! tab separator
    if (pos == 0) then
        call set_err_once(ierr, 1)  ! malformed line
        return
    end if

    fam_name = adjustl(family_line(1:pos-1))
    gene_line = adjustl(family_line(pos+1:))

    ! Split the gene_line using regex for multiple separators: tab, comma, space
    ! The regex pattern: '[\t, ]+'
    call regex_split(gene_line, '[\t, ]+', tokens, n_tokens, ierr)
    if (.not. is_ok(ierr)) return

    ! Allocate and copy gene list
    allocate(gene_list(n_tokens))
    do i = 1, n_tokens
        gene_list(i) = adjustl(tokens(i))
    end do

end subroutine parse_gene_family_line


!---------------------------
subroutine array_append(arr, value)
    implicit none
    character(len=:), allocatable, intent(inout) :: arr(:)
    character(len=*), intent(in) :: value
    integer :: n

    n = size(arr)
    if (n == 0) then
        allocate(arr(1))
        arr(1) = value
    else
        allocate(arr(n+1))
        arr(n+1) = value
    end if
end subroutine array_append


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Accessors (i assume we need those)
subroutine get_gene_index(gene_id, gene_ids, out_idx, ierr)
    character(len=*), intent(in) :: gene_id
    character(len=*), intent(in) :: gene_ids(:)
    integer(int32), intent(out) :: out_idx
    integer(int32), intent(out) :: ierr

    integer :: i

    ierr = ERR_OK
    out_idx = 0

    do i = 1, size(gene_ids)
        if (trim(gene_ids(i)) == trim(gene_id)) then
            out_idx = i
            return
        end if
    end do

    ierr = ERR_INVALID_INPUT
end subroutine get_gene_index

subroutine get_family_index(family_id, family_ids, out_idx, ierr)
    character(len=*), intent(in) :: family_id
    character(len=*), intent(in) :: family_ids(:)
    integer(int32), intent(out) :: out_idx
    integer(int32), intent(out) :: ierr

    integer :: i

    ierr = ERR_OK
    out_idx = 0

    do i = 1, size(family_ids)
        if (trim(family_ids(i)) == trim(family_id)) then
            out_idx = i
            return
        end if
    end do

    ierr = ERR_INVALID_INPUT
end subroutine get_family_index

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
subroutine read_int_1d(filename, array, ierr)
    character(len=*), intent(in)  :: filename
    integer(int32), intent(out) :: array(:)
    integer(int32), intent(out)   :: ierr

    character(len=:), pointer :: tmp(:)

    ! Call deserialize to get the pointer view
    call deserialize_int_1d(tmp, filename, ierr)
    if (.not. is_ok(ierr)) return

    ! Check that the size matches
    if (size(array) /= size(tmp)) then
        ierr = ERR_DIM_MISMATCH
        if (associated(tmp)) nullify(tmp)
        return
    end if

    ! Copy data into user’s array
    array = tmp

    ! Clean up pointer
    if (associated(tmp)) nullify(tmp)
end subroutine

subroutine read_char_1d(filename, array, ierr)
    character(len=*), intent(in)  :: filename
    character(len=*), intent(out) :: array(:)
    integer(int32), intent(out)   :: ierr

    character(len=:), pointer :: tmp(:)

    ! Call deserialize to get the pointer view
    call deserialize_char_1d(tmp, filename, ierr)
    if (.not. is_ok(ierr)) return

    ! Check that the size matches
    if (size(array) /= size(tmp)) then
        ierr = ERR_DIM_MISMATCH
        if (associated(tmp)) nullify(tmp)
        return
    end if

    ! Copy data into user’s array
    array = tmp

    ! Clean up pointer
    if (associated(tmp)) nullify(tmp)

end subroutine

subroutine read_real_2d(filename, array, ierr)
    character(len=*), intent(in) :: filename
    real(real64), intent(out)    :: array(:,:)
    integer(int32), intent(out)  :: ierr

    real(real64), pointer :: tmp(:,:)

    ! First, deserialize into a pointer
    call deserialize_real_2d(tmp, filename, ierr)
    if (.not. is_ok(ierr)) return

    ! Check dimensions match user array
    if (any(shape(array) /= shape(tmp))) then
        call set_err_once(ierr, ERR_DIMS_MISMATCH)
        return
    end if

    ! Copy into user’s array
    array = tmp

    ! Clean up
    if (associated(tmp)) deallocate(tmp)
end subroutine

subroutine save_gene_ids(filename, gene_ids, ierr)
    CHARACTER(len=*), INTENT(IN) :: gene_ids(:)
    CHARACTER(len=*), INTENT(IN) :: filename
    integer, intent(out) :: ierr

    call serialize_char_1D(filename, gene_ids, ierr)
end subroutine

subroutine save_expression_vectors(filename, expression_vectors, ierr)
    CHARACTER(len=*), INTENT(IN) :: filename
    real(real64), INTENT(IN) :: expression_vectors(:,:)
    integer, intent(out) :: ierr

    call serialize_real_2D(filename, expression_vectors, ierr)
end subroutine

subroutine save_gene_to_family(filename, gene_to_fam, ierr)
    CHARACTER(len=*), INTENT(IN) :: filename
    integer, intent(in) :: gene_to_fam(:)
    integer, intent(out) :: ierr

    call serialize_int_1D(filename, gene_to_fam, ierr)
end subroutine

subroutine save_gene_family_ids(filename, gene_family_ids, ierr)
    CHARACTER(len=*), INTENT(IN) :: filename
    CHARACTER(len=*), INTENT(IN) :: gene_family_ids(:)
    integer(int32), INTENT(OUT) :: ierr

    call serialize_char_1D(filename, gene_family_ids, ierr)
end subroutine

subroutine save_family_centroids(filename, family_centroids, ierr)
    CHARACTER(len=*), INTENT(IN) :: filename
    real(real64), INTENT(IN) :: family_centroids(:,:)
    integer(int32), INTENT(OUT) :: ierr

    call serialize_real_2D(filename, family_centroids, ierr)
end subroutine

subroutine save_shift_vectors(filename, shift_vectors, ierr)
    CHARACTER(len=*), INTENT(IN) :: filename
    real(real64), INTENT(IN) :: shift_vectors(:,:)
    integer(int32), INTENT(OUT) :: ierr

    call serialize_real_2D(filename, shift_vectors, ierr)
end subroutine

subroutine load_gene_ids(filename, gene_ids, ierr)
    character(len=*), intent(in)  :: filename
    character(len=*), intent(out) :: gene_ids(:)
    integer(int32), intent(out)   :: ierr

    call read_char_1d(filename, gene_ids, ierr)
end subroutine read_gene_ids

subroutine load_expression_vectors(filename, expression_vectors, ierr)
    character(len=*), intent(in) :: filename
    real(real64), intent(out)    :: expression_vectors(:,:)
    integer(int32), intent(out)  :: ierr

    call read_real_2d(filename, expression_vectors, ierr)
end subroutine read_expression_vectors

subroutine load_gene_to_family(filename, gene_to_fam, ierr)
    character(len=*), INTENT(IN) :: filename
    integer(int32), INTENT(OUT) :: gene_to_fam
    integer(int32), INTENT(OUT) :: ierr

    call read_int_1d(filename, gene_to_fam, ierr)
end subroutine

subroutine load_family_ids(filename, family_ids, ierr)
    character(len=*), intent(in) :: filename
    CHARACTER(len=*), intent(out) :: family_ids
    integer(int32), intent(out) :: ierr

    call read_char_1d(filename, family_ids, ierr)
end subroutine

subroutine load_family_centroids(filename, family_centroids, ierr)
    character(len=*), intent(in) :: filename
    real(real64), intent(out)    :: family_centroids(:,:)
    integer(int32), intent(out)  :: ierr

    call read_real_2d(filename, family_centroids, ierr)
end subroutine

subroutine load_shift_vectors(filename, shift_vectors, ierr)
    character(len=*), intent(in) :: filename
    real(real64), intent(out)    :: shift_vectors(:,:)
    integer(int32), intent(out)  :: ierr

    call read_real_2d(filename, shift_vectors, ierr)
end subroutine

end module tox_data_tools