program test_expression_readers
    use iso_fortran_env, only: real64, int32
    use tox_data_tools
    implicit none

    integer, parameter :: n_genes = 11599   ! adjust to your number of genes
    integer, parameter :: n_header_rows = 1
    ! integer, parameter :: n_kallisto_files = 3
    integer, parameter :: n_bowtie_files = 3
    ! character(len=256) :: kallisto_files(n_kallisto_files)
    character(len=256) :: bowtie_files(n_bowtie_files)
    character(len=256) :: gene_ids(n_genes)
    ! real(real64) :: kallisto_expr(n_kallisto_files, n_genes)
    real(real64) :: bowtie_expr(n_bowtie_files, n_genes)
    integer(int32) :: ierr
    integer :: i

    integer :: target_column
    real(real64) :: target_expression(n_bowtie_files)

    ! Example file list (replace with your real file paths)
    ! kallisto_files = [ 'material/control_rep1.tsv', &
    !                    'material/control_rep2.tsv', &
    !                    'material/control_rep3.tsv' ]

    bowtie_files = [ 'material/gene_abundance_63.tsv', &
                     'material/gene_abundance_64.tsv', &
                     'material/gene_abundance_65.tsv' ]

    ! ! Read Kallisto files
    ! call read_kallisto_files(kallisto_files, n_kallisto_files, gene_ids, n_genes, kallisto_expr, n_header_rows, ierr)
    ! if (ierr /= 0) stop 'Error reading Kallisto files'

    ! write(*,*) 'Kallisto expression matrix (first 10 genes only):'
    ! do i = 1, n_kallisto_files
    !     write(*, '(10(F8.2,1X))') kallisto_expr(i, 1:min(10, n_genes))
    ! end do


    ! Read Bowtie files
    call read_bowtie_files(bowtie_files, n_bowtie_files, gene_ids, n_genes, bowtie_expr, n_header_rows, ierr)
    if (ierr /= 0) stop 'Error reading Bowtie files'

    write(*,*) 'Bowtie expression matrix:'
    do i = 1, n_bowtie_files
        write(*, '(100(F8.2,1X))') bowtie_expr(i,1:min(10, n_genes))
    end do

    write(*,*) 'Index for YAL068C_mRNA:'
    call get_gene_index('YAL068C_mRNA', gene_ids, target_column, ierr)
    print *, target_column

    call get_expression_vector(target_column, bowtie_expr, target_expression, ierr)

    write(*,*) 'Expression vector for YAL068C_mRNA:'
    do i = 1, n_bowtie_files
        write(*, '(100(F8.2,1X))') target_expression(i)
    end do

end program test_expression_readers
