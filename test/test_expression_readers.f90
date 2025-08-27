program test_expression_readers
    use iso_fortran_env, only: real64, int32
    use tox_data_tools
    implicit none

    integer, parameter :: n_genes = 11599   ! adjust to your number of genes
    integer, parameter :: n_kallisto_files = 3
    integer, parameter :: n_bowtie_files = 3
    character(len=256) :: kallisto_files(n_kallisto_files)
    character(len=256) :: gene_ids(n_genes)
    real(real64) :: kallisto_expr(n_kallisto_files, n_genes)
    integer(int32) :: ierr
    integer :: i
    integer, parameter :: n_header_rows = 1

    ! Example file list (replace with your real file paths)
    kallisto_files = [ 'material/control_rep1.tsv', &
                       'material/control_rep2.tsv', &
                       'material/control_rep3.tsv' ]

    ! Read Kallisto files
    call read_kallisto_files(kallisto_files, n_kallisto_files, gene_ids, n_genes, kallisto_expr, n_header_rows, ierr)
    if (ierr /= 0) stop 'Error reading Kallisto files'

    write(*,*) 'Kallisto expression matrix (first 10 genes only):'
    do i = 1, n_kallisto_files
        write(*, '(10(F8.2,1X))') kallisto_expr(i, 1:min(10, n_genes))
    end do

end program test_expression_readers
