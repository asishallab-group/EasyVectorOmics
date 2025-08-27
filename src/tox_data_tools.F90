module tox_data_tools
use iso_fortran_env, only: real64, int32
use tox_errors, only: ERR_OK, ERR_INVALID_INPUT, ERR_SIZE_MISMATCH, ERR_FILE_OPEN, set_ok
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

subroutine read_kallisto_files(file_list, n_files, gene_ids, n_genes, &
                                                          expression_vectors, n_header_rows, ierr)
    implicit none
    character(len=*), intent(in) :: file_list(:)
    integer(int32), intent(in) :: n_files
    integer(int32), intent(in) :: n_genes
    character(len=*), intent(inout) :: gene_ids(n_genes)    ! empty for first file   
    real(real64), intent(out) :: expression_vectors(n_files, n_genes)
    integer(int32), intent(in) :: n_header_rows
    integer(int32), intent(out) :: ierr

    integer :: i, ios, unit, idx, row_number
    character(len=1024) :: line
    character(len=256) :: gene
    real(real64) :: tpm
    integer :: j

    ierr = 0
    expression_vectors = 0.0d0

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

            call parse_kallisto_line(line, gene, tpm)

            if (i == 1) then
                ! first file, fill gene_ids array
                gene_ids(row_number) = gene
                idx = row_number
            else
                ! subsequent files, find correct column
                idx = 0
                do idx = 1, n_genes
                    if (trim(gene_ids(idx)) == gene) exit
                end do
                if (idx > n_genes) then
                    idx = 0
                    write(*,*) 'Warning: gene ', trim(gene), ' not found in master gene_ids'
                end if
            end if

            if (idx > 0) then
                expression_vectors(i, idx) = tpm
            end if

            row_number = row_number + 1
        end do

        close(unit)
    end do

end subroutine read_kallisto_files

subroutine parse_kallisto_line(line, gene, tpm)
    implicit none
    character(len=*), intent(in) :: line
    character(len=*), intent(out) :: gene
    real(real64), intent(out) :: tpm

    integer :: pos1, pos2, i

    ! Find first tab -> gene
    pos1 = index(line, char(9))
    if (pos1 == 0) then
        gene = adjustl(trim(line))
        tpm = 0.0d0
        return
    end if
    gene = adjustl(trim(line(1:pos1-1)))

    ! Find last tab -> tpm is after last tab
    pos2 = 0
    do i = len(line), 1, -1
        if (line(i:i) == char(9)) then
            pos2 = i
            exit
        end if
    end do
    if (pos2 == 0) then
        tpm = 0.0d0
        return
    end if

    read(line(pos2+1:), *) tpm
end subroutine parse_kallisto_line



subroutine read_bowtie_files(file_list, n_files, gene_ids, n_genes, expression_vectors, ierr)
    character(len=*), intent(in) :: file_list(:)
    integer(int32), intent(in) :: n_files
    integer(int32), intent(in) :: n_genes
    character(len=*), intent(inout) :: gene_ids(n_genes)  ! empty for first file
    real(real64), intent(out) :: expression_vectors(n_files, n_genes)
    integer(int32), intent(out) :: ierr

    integer :: i, ios, unit, idx, row_number
    character(len=2048) :: line
    character(len=256) :: gene
    real(real64) :: tpm

    ierr = ERR_OK
    expression_vectors = 0.0d0

    do i = 1, n_files
        open(newunit=unit, file=file_list(i), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,*) 'Error: cannot open file ', trim(file_list(i))
            ierr = ERR_FILE_OPEN
            return
        end if

        row_number = 1
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit

            read(line, *) gene, tpm
            gene = adjustl(trim(gene))

            if (i == 1) then
                gene_ids(row_number) = gene
                idx = row_number
            else
                idx = 0
                do idx = 1, n_genes
                    if (trim(gene_ids(idx)) == gene) exit
                end do
                if (idx > n_genes) then
                    idx = 0
                    write(*,*) 'Warning: gene ', trim(gene), ' not found in master gene_ids'
                end if
            end if

            if (idx > 0) then
                expression_vectors(i, idx) = tpm
            end if

            row_number = row_number + 1
        end do

        close(unit)
    end do
end subroutine read_bowtie_files


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

end module tox_data_tools