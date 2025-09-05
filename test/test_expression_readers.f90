program test_expression_readers
    use iso_fortran_env, only: real64, int32
    use tox_data_tools
    use tox_data_validation
    use tox_data_read_write
    use array_utils, only: get_array_metadata
    use tox_gene_centroids
    implicit none

    integer(int32), allocatable :: gene_to_fam(:)
    character(len=256), allocatable :: gene_family_ids(:), gene_ids(:)
    real(real64), allocatable :: kallisto_expr(:,:), kallisto_expr_verify(:,:)
    integer(int32) :: n_genes, n_families, ierr, total_samples
    integer, parameter :: n_header_rows = 1
    integer(int32) :: ndims, dims(2)
    real(real64), ALLOCATABLE :: family_centroids(:,:)
    logical, ALLOCATABLE :: ortholog_mask(:)
    integer(int32), ALLOCATABLE :: selected_indices(:)
    integer(int32) :: no_family_entries
    integer(int32), allocatable :: indice_set(:)
    integer(int32) :: current_idx 
    real(real64) :: sample_sum
    real(real64), allocatable :: test_centroid(:)
    integer(int32) :: n_family_genes

    ! File lists
    character(len=256), allocatable :: files_6_replicates(:)
    character(len=256), allocatable :: files_7_replicates(:)
    character(len=256), allocatable :: files_5_replicates(:)
    character(len=256), allocatable :: files_4_replicates(:)

    integer :: i, j

    ierr = 0
    
    !-------------------------------
    ! Initialize file lists
    allocate(files_6_replicates(12))
    allocate(files_4_replicates(1))
    allocate(files_5_replicates(1))
    ALLOCATE(files_7_replicates(1))
    files_6_replicates = [ &
        'material/kallisto_sex_data_Adipose.tsv  ', &
        'material/kallisto_sex_data_Adrenal.tsv  ', &
        'material/kallisto_sex_data_Colon.tsv    ', &
        'material/kallisto_sex_data_Heart.tsv    ', &
        'material/kallisto_sex_data_Liver.tsv    ', &
        'material/kallisto_sex_data_Lung.tsv     ', &
        'material/kallisto_sex_data_Muscle.tsv   ', &   
        'material/kallisto_sex_data_Skin.tsv     ', &
        'material/kallisto_sex_data_Spleen.tsv   ', &
        'material/kallisto_sex_data_Thyroid.tsv  '  &
    ]

    files_5_replicates = [ &
        'material/kallisto_sex_data_Brain.tsv']

    files_4_replicates = [ &
        'material/kallisto_sex_data_Pituitary.tsv']
    
    files_7_replicates = [ &
        'material/kallisto_sex_data_Testis.tsv']

    !-------------------------------
    ! Calculate total samples in advance
    total_samples = 10 * 6 + 7 + 5 + 4 

    !-------------------------------
    ! Step 1: Read gene IDs from first Kallisto file
    n_genes = 88327  ! Known size from Kallisto files
    allocate(gene_ids(n_genes))
    gene_ids = ""
    
    call read_gene_ids_from_file(files_6_replicates(1), gene_ids, n_header_rows, 1, ierr)
    if (ierr /= 0) then
        write(*,*) 'Error reading gene IDs from first file'
        stop
    end if

    !-------------------------------
    ! Step 2: Read expression data from all Kallisto files
    allocate(kallisto_expr(total_samples, n_genes))
    kallisto_expr = 0.0_real64

    ! Process files with 6 replicates (columns 2-7)
    write(*,*) 'Processing 12 files with 6 replicates each...'
    
    call read_tabular_files(files_6_replicates, gene_ids, kallisto_expr, &
                           n_header_rows, 1, [2, 3, 4, 5, 6, 7], 1, ierr)
    if (ierr /= 0) then
        write(*,*) 'Error reading files with 6 replicates'
        stop
    end if

    write(*,*) 'Processed 60 samples from 6-replicate files'

    !-------------------------------
    ! Step 3: Process file with 7 replicates (columns 2-8)
    write(*,*) 'Processing 1 file with 7 replicates...'
    
    call read_tabular_files([files_7_replicates], gene_ids, kallisto_expr, &
                           n_header_rows, 1, [2, 3, 4, 5, 6, 7, 8], 61, ierr)
    if (ierr /= 0) then
        write(*,*) 'Error reading file with 7 replicates'
        stop
    end if

    write(*,*) 'Processed additional 5 samples from 5-replicate file'

    call read_tabular_files([files_5_replicates], gene_ids, kallisto_expr, &
                           n_header_rows, 1, [2, 3, 4, 5, 6], 69, ierr)
    if (ierr /= 0) then
        write(*,*) 'Error reading file with 5 replicates'
        stop
    end if

    write(*,*) 'Processed additional 4 samples from 4-replicate file'

    call read_tabular_files([files_4_replicates], gene_ids, kallisto_expr, &
                           n_header_rows, 1, [2, 3, 4, 5], 73, ierr)
    if (ierr /= 0) then
        write(*,*) 'Error reading file with 4 replicates'
        stop
    end if

    write(*,*) 'Processed additional 4 samples from 4-replicate file'

    !-------------------------------
        ! Step 4: Read orthogroups (family mapping)
    n_families = 15512  ! From your known output
    allocate(gene_family_ids(n_families))
    allocate(gene_to_fam(n_genes))
    gene_family_ids = ""
    gene_to_fam = 0

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
    write(*,*) 'First 5 genes, first 40 samples:'
    do i = 1, min(40, total_samples)
        write(*, '(10(F8.2,1X))') kallisto_expr(i, 1:min(5,n_genes))
    end do

    write(*,*) 'First 5 gene IDs:'
    do i = 1, min(5, n_genes)
        write(*, *) trim(gene_ids(i))
    end do

    write(*,*) 'Verifying data...'
    call validate_expression_data(kallisto_expr, .true., ierr)
    if (ierr /= 0) write(*,*) 'Expression data could not be verified: ', ierr
    
    call validate_gene_to_family_mapping(gene_to_fam, n_families, ierr)
    if (ierr/=0) write(*,*) 'Gene to family mapping could not be verified: ', ierr

    call validate_gene_ids_uniqueness(gene_ids, ierr)
    if (ierr/=0) write(*,*) 'Gene IDs contain duplicates: ', ierr

    allocate(family_centroids(total_samples, n_families))
    ALLOCATE(ortholog_mask(n_genes))
    allocate(selected_indices(n_genes))

    no_family_entries = 0

    do i = 1, n_genes
        if (gene_to_fam(i) < 1 .or. gene_to_fam(i) > n_families) then
          no_family_entries = no_family_entries + 1
        end if
    end do

    write(*,*) 'Genes without families: ', no_family_entries

    do i = 1, n_genes
        ortholog_mask(i) = .true.
    end do

    call group_centroid(kallisto_expr, total_samples, n_genes, gene_to_fam, n_families, family_centroids, &
     .true., ortholog_mask, selected_indices, ierr)
    
    if(ierr/=0) write (*,*) 'Failed to compute centroids: ', ierr

    write(*,*) 'First 5 family centroids:'
    do i = 1, 20
        write(*, '(10(F8.2,1X))') family_centroids(i, 1:min(5,n_families))
    end do

    ! First, count genes in family 2
    n_family_genes = 0
    do i = 1, n_genes
        if (gene_to_fam(i) == 2) then
            n_family_genes = n_family_genes + 1
        end if
    end do

    ! Allocate arrays
    allocate(indice_set(n_family_genes))
    allocate(test_centroid(total_samples))

    ! Fill indice_set with indices of genes in family 2
    current_idx = 1
    do i = 1, n_genes
        if (gene_to_fam(i) == 2) then
            indice_set(current_idx) = i
            current_idx = current_idx + 1
        end if
    end do

    ! Calculate centroid (mean across genes for each sample)
    do j = 1, total_samples
        sample_sum = 0.0_real64
        do i = 1, n_family_genes
            sample_sum = sample_sum + kallisto_expr(j, indice_set(i))  ! Note: j is sample, indice_set(i) is gene
        end do
        test_centroid(j) = sample_sum / n_family_genes
    end do

    ! Output the centroid
    write(*,*) 'Test centroid for family 2:'
    do j = 1, total_samples
        write(*, '(F8.2,1X)', advance='no') test_centroid(j)
        if (mod(j, 10) == 0) write(*,*)  ! New line every 10 values
    end do
    write(*,*) 
    
    call save_expression_vectors(kallisto_expr, 'kallisto_sex_data.bin', ierr)
    if (ierr/=0) then 
        write(*,*) 'Expression data could not be saved: ', ierr
    else 
        write(*,*) 'Expression data saved'
    end if

    call get_array_metadata('kallisto_sex_data.bin', dims, 2, ndims, ierr)
    if(ierr/=0) then 
        write(*,*) 'Metadata could not be read'
    else 
        write(*,*) 'Metadata read'
    end if

    allocate(kallisto_expr_verify(dims(1), dims(2)))

    call load_expression_vectors(kallisto_expr_verify, 'kallisto_sex_data.bin', ierr)
    if(ierr/=0) then
        write(*,*) 'Expression vectors could not be loaded', ierr
    else 
        write(*,*) 'Expression vectors read'
    end if 

    if (any(kallisto_expr/=kallisto_expr_verify)) write(*,*) 'Data not equal'

end program test_expression_readers