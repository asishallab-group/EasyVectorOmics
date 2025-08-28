program test_expression_readers
    use iso_fortran_env, only: real64, int32
    use tox_data_tools
    implicit none

    integer(int32), allocatable :: gene_to_fam(:)
    character(len=100), allocatable :: gene_family_ids(:), gene_ids(:)
    real(real64), allocatable :: kallisto_expr(:,:), temp_expr(:,:)
    integer(int32) :: n_genes, n_families, ierr, total_samples, sample_offset
    integer, parameter :: n_header_rows = 1

    ! File lists
    character(len=256), allocatable :: files_6_replicates(:)
    character(len=256) :: file_8_replicates
    integer :: i

    ierr = 0
    
    !-------------------------------
    ! Initialize file lists
    allocate(files_6_replicates(12))
    files_6_replicates = [ &
        'material/kallisto_sex_data_Adipose.tsv  ', &
        'material/kallisto_sex_data_Adrenal.tsv  ', &
        'material/kallisto_sex_data_Brain.tsv    ', &
        'material/kallisto_sex_data_Colon.tsv    ', &
        'material/kallisto_sex_data_Heart.tsv    ', &
        'material/kallisto_sex_data_Liver.tsv    ', &
        'material/kallisto_sex_data_Lung.tsv     ', &
        'material/kallisto_sex_data_Muscle.tsv   ', &
        'material/kallisto_sex_data_Pituitary.tsv', &
        'material/kallisto_sex_data_Skin.tsv     ', &
        'material/kallisto_sex_data_Spleen.tsv   ', &
        'material/kallisto_sex_data_Thyroid.tsv  '  &
    ]
    
    file_8_replicates = 'material/kallisto_sex_data_Testis.tsv'

    !-------------------------------
    ! Get array sizes from user knowledge
    n_genes = 88327      ! From your output
    n_families = 15513   ! From your output

    write(*,*) 'Number of genes: ', n_genes
    write(*,*) 'Number of families: ', n_families

    ! Allocate arrays with the correct sizes
    allocate(gene_ids(n_genes))
    allocate(gene_family_ids(n_families))
    allocate(gene_to_fam(n_genes))
    
    ! Initialize arrays
    gene_ids = ""
    gene_family_ids = ""
    gene_to_fam = 0

    !-------------------------------
    ! Step 1: Read gene IDs from first file (6 replicates)
    ! call read_gene_ids_from_file(files_6_replicates(1), gene_ids, n_header_rows, 1, ierr)
    ! if (ierr /= 0) then
    !     write(*,*) 'Error reading gene IDs from first file'
    !     stop
    ! end if

    !-------------------------------
    ! Step 2: Process files with 6 replicates (columns 2-7)
    write(*,*) 'Processing 12 files with 6 replicates each...'
    
    ! Allocate temporary storage for 12 files × 6 replicates = 72 samples
    allocate(temp_expr(72, n_genes))
    temp_expr = 0.0_real64
    
    call read_tabular_files(files_6_replicates, gene_ids, temp_expr, &
                           n_header_rows, 1, [2, 3, 4, 5, 6, 7], ierr)
    if (ierr /= 0) then
        write(*,*) 'Error reading files with 6 replicates'
        stop
    end if

    ! Allocate main expression matrix and copy first 72 samples
    total_samples = 72
    allocate(kallisto_expr(total_samples, n_genes))
    kallisto_expr = temp_expr
    deallocate(temp_expr)
    
    write(*,*) 'Processed ', total_samples, ' samples from 6-replicate files'

    !-------------------------------
    ! Step 3: Process file with 8 replicates (columns 2-9)
    write(*,*) 'Processing 1 file with 8 replicates...'
    
    ! Allocate temporary storage for 1 file × 8 replicates = 8 samples
    allocate(temp_expr(8, n_genes))
    temp_expr = 0.0_real64
    
    call read_tabular_files([file_8_replicates], gene_ids, temp_expr, &
                           n_header_rows, 1, [2, 3, 4, 5, 6, 7, 8, 9], ierr)
    if (ierr /= 0) then
        write(*,*) 'Error reading file with 8 replicates'
        stop
    end if

    ! Resize main expression matrix to accommodate additional 8 samples
    sample_offset = total_samples
    total_samples = total_samples + 8
    
    ! Resize kallisto_expr (this creates a temporary array)
    kallisto_expr = reshape(kallisto_expr, [total_samples, n_genes])
    
    ! Copy the new samples
    kallisto_expr(sample_offset+1:sample_offset+8, :) = temp_expr
    deallocate(temp_expr)
    
    write(*,*) 'Processed additional 8 samples from 8-replicate file'

    !-------------------------------
    ! Step 4: Read orthogroups (family mapping)
    write(*,*) 'Reading orthogroups file...'
    call read_family_file('material/Orthogroups.tsv', gene_ids, gene_family_ids, gene_to_fam, ierr)
    if (ierr /= 0) then
        write(*,*) 'Error reading orthofinder family file'
        stop
    end if

    !-------------------------------
    ! Final summary
    write(*,*) '========================================'
    write(*,*) 'PROCESSING COMPLETE'
    write(*,*) '========================================'
    write(*,*) 'Total genes: ', n_genes
    write(*,*) 'Total families: ', n_families
    write(*,*) 'Total samples: ', total_samples
    write(*,*) 'Expression matrix shape: ', shape(kallisto_expr)
    write(*,*) 'Genes with family assignments: ', count(gene_to_fam > 0)
    write(*,*) 'Genes without family assignments: ', count(gene_to_fam == 0)
    write(*,*) '========================================'

    ! Optional: Show a small preview
    write(*,*) 'First 5 genes, first 5 samples:'
    do i = 1, min(5, total_samples)
        write(*, '(5(F8.2,1X))') kallisto_expr(i, 1:min(5,n_genes))
    end do

end program test_expression_readers