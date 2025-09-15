!> Unit test suite for expression readers and data processing routines.
module mod_test_tox_data
  use asserts
  use iso_fortran_env, only: real64, int32
  use tox_data_tools
  use tox_data_validation
  use tox_data_accessors
  use tox_data_read_write
  use array_utils, only: get_array_metadata
  use tox_gene_centroids
  use tox_shift_vectors
  use tox_errors
  implicit none
  public

  abstract interface
    subroutine test_interface()
    end subroutine test_interface
  end interface

  type :: test_case
    character(len=64) :: name
    procedure(test_interface), pointer, nopass :: test_proc => null()
  end type test_case

  ! Global test data
  integer(int32), allocatable :: gene_to_fam(:)
  character(len=256), allocatable :: gene_family_ids(:), gene_ids(:)
  real(real64), allocatable :: kallisto_expr(:,:), kallisto_expr_verify(:,:)
  integer(int32) :: n_genes, n_families, total_samples
  real(real64), allocatable :: family_centroids(:,:)
  real(real64), allocatable :: shift_vectors(:,:)  ! Added for shift vectors

contains

  !> Get array of all available tests.
  function get_all_tests() result(all_tests)
    type(test_case) :: all_tests(13)  ! Increased from 6 to 8
    all_tests(1) = test_case("test_read_gene_ids", test_read_gene_ids)
    all_tests(2) = test_case("test_read_expression_data", test_read_expression_data)
    all_tests(3) = test_case("test_read_family_mapping", test_read_family_mapping)
    all_tests(4) = test_case("test_validate_data", test_validate_data)
    all_tests(5) = test_case("test_compute_centroids", test_compute_centroids)
    all_tests(6) = test_case("test_write_read_expression_data", test_write_read_expression_data)
    all_tests(7) = test_case("test_compute_shift_vectors", test_compute_shift_vectors)
    all_tests(8) = test_case("test_write_read_shift_vectors", test_write_read_shift_vectors)
    all_tests(9) = test_case("test_read_write_gene_ids", test_read_write_gene_ids)
    all_tests(10) = test_case("test_read_write_gene_to_fam", test_read_write_gene_to_fam)
    all_tests(11) = test_case("test_read_write_family_ids", test_read_write_family_ids)
    all_tests(12) = test_case("test_read_write_family_centroids", test_read_write_centroids)
    all_tests(13) = test_case("test_data_accessors", test_data_accessors)
  end function get_all_tests

  !> Setup global test data
  subroutine setup_global_data()
    character(len=256), allocatable :: files_6_replicates(:)
    character(len=256), allocatable :: files_7_replicates(:)
    character(len=256), allocatable :: files_5_replicates(:)
    character(len=256), allocatable :: files_4_replicates(:)
    integer(int32) :: ierr, i, n_genes_kept
    logical, allocatable :: ortholog_mask(:)
    integer(int32), allocatable :: selected_indices(:)

    ! Initialize file lists
    allocate(files_6_replicates(10))
    allocate(files_4_replicates(1))
    allocate(files_5_replicates(1))
    allocate(files_7_replicates(1))
    
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

    ! Calculate total samples
    total_samples = 10 * 6 + 7 + 5 + 4 
    n_genes = 88327  ! Original number of genes
    n_families = 15512

    ! Allocate arrays with original size
    allocate(gene_ids(n_genes))
    allocate(kallisto_expr(total_samples, n_genes))
    allocate(gene_family_ids(n_families))
    allocate(gene_to_fam(n_genes))

    ! Read gene IDs
    call read_gene_ids_from_file(files_6_replicates(1), gene_ids, 1, 1, ierr)
    call assert_equal_int(ierr, 0, "Reading gene IDs should succeed")

    ! Read expression data
    kallisto_expr = 0.0_real64
    call read_expression_vectors(files_6_replicates, gene_ids, kallisto_expr, &
                          1, 1, [2, 3, 4, 5, 6, 7], 1, ierr, char(9))
    
    call assert_equal_int(ierr, 0, "Reading 6-replicate files should succeed")

    call read_expression_vectors([files_7_replicates], gene_ids, kallisto_expr, &
                          1, 1, [2, 3, 4, 5, 6, 7, 8], 61, ierr)
    call assert_equal_int(ierr, 0, "Reading 7-replicate file should succeed")

    call read_expression_vectors([files_5_replicates], gene_ids, kallisto_expr, &
                          1, 1, [2, 3, 4, 5, 6], 69, ierr)
    call assert_equal_int(ierr, 0, "Reading 5-replicate file should succeed")

    call read_expression_vectors([files_4_replicates], gene_ids, kallisto_expr, &
                          1, 1, [2, 3, 4, 5], 73, ierr)
    call assert_equal_int(ierr, 0, "Reading 4-replicate file should succeed")

    ! Read family mapping
    call read_family_file('material/Orthogroups.tsv', gene_ids, gene_family_ids, gene_to_fam, ierr)
    call assert_equal_int(ierr, 0, "Reading family file should succeed")

    ! Filter out genes without family assignments
    call filter_unassigned_genes(gene_ids, kallisto_expr, gene_to_fam, n_genes_kept, ierr)
    call assert_equal_int(ierr, 0, "Filtering unassigned genes should succeed")
    n_genes = n_genes_kept  ! Update n_genes to reflect the filtered count

    ! Print first 10 family IDs
    write(*,*) 'First 10 family IDs:'
    do i = 1, min(10, n_families)
      write(*,*) trim(gene_family_ids(i))
    end do

    ! Compute centroids
    allocate(family_centroids(total_samples, n_families))
    
    ! Create explicit arrays instead of array constructors
    allocate(ortholog_mask(n_genes))
    allocate(selected_indices(n_genes))
    
    do i = 1, n_genes
      ortholog_mask(i) = .true.
      selected_indices(i) = i
    end do

    write(*,*) 'Size of family_centroids: ', size(family_centroids, 1), size(family_centroids, 2)
    
    call group_centroid(kallisto_expr, total_samples, n_genes, gene_to_fam, &
                      n_families, family_centroids, 'all', selected_indices, ierr)
    call assert_equal_int(ierr, 0, "Computing centroids should succeed")
    
    ! Compute shift vectors
    allocate(shift_vectors(2*total_samples, n_genes))
    call compute_shift_vector_field(total_samples, n_genes, n_families, kallisto_expr, &
                                  family_centroids, gene_to_fam, shift_vectors, ierr)
    call assert_equal_int(ierr, 0, "Computing shift vectors should succeed")
    
    ! Clean up temporary arrays
    deallocate(ortholog_mask, selected_indices)
  end subroutine setup_global_data

  !> Run all expression reader tests.
  subroutine run_all_tests_tox_data()
    type(test_case) :: all_tests(13)  ! Updated
    integer(int32) :: i
    
    ! Setup global data first
    call setup_global_data()
    
    all_tests = get_all_tests()
    do i = 1, size(all_tests)
      call all_tests(i)%test_proc()
      print *, trim(all_tests(i)%name), " passed."
    end do
    print *, "All tox data tests passed successfully."
  end subroutine run_all_tests_tox_data

  !> Run specific expression reader tests by name.
  subroutine run_named_tests_tox_data(test_names)
    character(len=*), intent(in) :: test_names(:)
    type(test_case) :: all_tests(13)  ! Updated
    integer(int32) :: i, j
    logical :: found
    
    ! Setup global data first
    call setup_global_data()
    
    all_tests = get_all_tests()
    do i = 1, size(test_names)
      found = .false.
      do j = 1, size(all_tests)
        if (trim(test_names(i)) == trim(all_tests(j)%name)) then
          call all_tests(j)%test_proc()
          print *, trim(test_names(i)), " passed."
          found = .true.
          exit
        end if
      end do
      if (.not. found) then
        print *, "Unknown test: ", trim(test_names(i))
      end if
    end do
  end subroutine run_named_tests_tox_data

  !> Test reading gene IDs
  subroutine test_read_gene_ids()
    call assert_true(allocated(gene_ids), "Gene IDs should be allocated")
    call assert_equal_int(size(gene_ids), n_genes, "Number of gene IDs should match")
    call assert_true(len_trim(gene_ids(1)) > 0, "First gene ID should not be empty")
  end subroutine test_read_gene_ids

  !> Test reading expression data
  subroutine test_read_expression_data()
    call assert_true(allocated(kallisto_expr), "Expression data should be allocated")
    call assert_equal_int(size(kallisto_expr, 1), total_samples, "Number of samples should match")
    call assert_equal_int(size(kallisto_expr, 2), n_genes, "Number of genes should match")
    call assert_true(all(kallisto_expr >= 0.0_real64), "All expression values should be non-negative")
  end subroutine test_read_expression_data

  !> Test reading family mapping
  subroutine test_read_family_mapping()
    call assert_true(allocated(gene_to_fam), "Gene to family mapping should be allocated")
    call assert_true(allocated(gene_family_ids), "Family IDs should be allocated")
    call assert_equal_int(size(gene_to_fam), n_genes, "Gene to family mapping size should match")
    call assert_equal_int(size(gene_family_ids), n_families, "Number of family IDs should match")
  end subroutine test_read_family_mapping

  !> Test data validation
  subroutine test_validate_data()
    integer(int32) :: ierr
    write(*,*) 'Testing data validation... This might take a while'
    
    call validate_all_data(n_genes, n_families, total_samples, total_samples, gene_ids, gene_family_ids, &
                           gene_to_fam, kallisto_expr, family_centroids, shift_vectors, ierr, .true., .true.)
  end subroutine test_validate_data

  !> Test centroid computation
  subroutine test_compute_centroids()
    integer(int32) :: i, j, n_family_genes, current_idx
    integer(int32), allocatable :: indice_set(:)
    real(real64), allocatable :: test_centroid(:)
    real(real64) :: sample_sum
    
    call assert_true(allocated(family_centroids), "Family centroids should be allocated")
    call assert_equal_int(size(family_centroids, 1), total_samples, "Centroid samples should match")
    call assert_equal_int(size(family_centroids, 2), n_families, "Number of centroid families should match")
    
    ! Test specific family (family 2)
    n_family_genes = count(gene_to_fam == 2)
    allocate(indice_set(n_family_genes))
    allocate(test_centroid(total_samples))
    
    current_idx = 1
    do i = 1, n_genes
      if (gene_to_fam(i) == 2) then
        indice_set(current_idx) = i
        current_idx = current_idx + 1
      end if
    end do

    do j = 1, total_samples
      sample_sum = 0.0_real64
      do i = 1, n_family_genes
        sample_sum = sample_sum + kallisto_expr(j, indice_set(i))
      end do
      test_centroid(j) = sample_sum / n_family_genes
    end do

    ! Fixed call to match the assert_equal_array_real interface
    call assert_equal_array_real(test_centroid, family_centroids(:, 2), &
                                total_samples, 1e-12_real64, &
                                "Family 2 centroid should match manual calculation")
  end subroutine test_compute_centroids

  !> Test binary write/read operations
  subroutine test_write_read_expression_data()
    integer(int32) :: ierr, ndims, dims(2)
    integer(int32) :: total_elements
    
    call save_expression_vectors(kallisto_expr, 'test_kallisto_data.bin', ierr)
    call assert_equal_int(ierr, 0, "Saving expression data should succeed")
    
    call get_array_metadata('test_kallisto_data.bin', dims, 2, ndims, ierr)
    call assert_equal_int(ierr, 0, "Reading metadata should succeed")
    call assert_equal_int(ndims, 2, "Array should have 2 dimensions")
    call assert_equal_int(dims(1), total_samples, "First dimension should match sample count")
    call assert_equal_int(dims(2), n_genes, "Second dimension should match gene count")
    
    allocate(kallisto_expr_verify(dims(1), dims(2)))
    call load_expression_vectors(kallisto_expr_verify, 'test_kallisto_data.bin', ierr)
    call assert_equal_int(ierr, 0, "Loading expression data should succeed")
    
    ! Fixed call to match the assert_equal_array_real interface
    total_elements = size(kallisto_expr)
    call assert_equal_array_real(reshape(kallisto_expr, [total_elements]), &
                                reshape(kallisto_expr_verify, [total_elements]), &
                                total_elements, 1e-12_real64, &
                                "Original and loaded data should be identical")
    
    ! Clean up
    deallocate(kallisto_expr_verify)
  end subroutine test_write_read_expression_data

  !> Test computation of shift vectors
  subroutine test_compute_shift_vectors()
    integer(int32) :: ierr, i, j
    real(real64), allocatable :: test_shift_vectors(:,:)
    
    call assert_true(allocated(shift_vectors), "Shift vectors should be allocated")
    call assert_equal_int(size(shift_vectors, 1), 2*total_samples, "Shift vectors should have 2*d rows")
    call assert_equal_int(size(shift_vectors, 2), n_genes, "Shift vectors should have n_genes columns")
    
    ! Verify centroids are correctly placed in first d rows
    do i = 1, n_genes
      if (gene_to_fam(i) >= 1 .and. gene_to_fam(i) <= n_families) then
        call assert_equal_array_real(shift_vectors(1:total_samples, i), &
                                    family_centroids(:, gene_to_fam(i)), &
                                    total_samples, 1e-12_real64, &
                                    "Centroid should match family centroid")
      end if
    end do
    
    ! Verify shift vectors are correctly computed
    do i = 1, n_genes
      if (gene_to_fam(i) >= 1 .and. gene_to_fam(i) <= n_families) then
        call assert_equal_array_real(shift_vectors(total_samples+1:2*total_samples, i), &
                                    kallisto_expr(:, i) - family_centroids(:, gene_to_fam(i)), &
                                    total_samples, 1e-12_real64, &
                                    "Shift vector should be expression minus centroid")
      end if
    end do
  end subroutine test_compute_shift_vectors

  !> Test writing and reading shift vectors
  subroutine test_write_read_shift_vectors()
    integer(int32) :: ierr, ndims, dims(2)
    integer(int32) :: total_elements
    real(real64), allocatable :: shift_vectors_loaded(:,:)
    
    call save_expression_vectors(shift_vectors, 'test_shift_vectors.bin', ierr)
    call assert_equal_int(ierr, 0, "Saving shift vectors should succeed")
    
    call get_array_metadata('test_shift_vectors.bin', dims, 2, ndims, ierr)
    call assert_equal_int(ierr, 0, "Reading metadata should succeed")
    call assert_equal_int(ndims, 2, "Array should have 2 dimensions")
    call assert_equal_int(dims(1), 2*total_samples, "First dimension should match 2*d")
    call assert_equal_int(dims(2), n_genes, "Second dimension should match gene count")
    
    allocate(shift_vectors_loaded(dims(1), dims(2)))
    call load_expression_vectors(shift_vectors_loaded, 'test_shift_vectors.bin', ierr)
    call assert_equal_int(ierr, 0, "Loading shift vectors should succeed")
    
    ! Verify loaded data matches original
    total_elements = size(shift_vectors)
    call assert_equal_array_real(reshape(shift_vectors, [total_elements]), &
                                reshape(shift_vectors_loaded, [total_elements]), &
                                total_elements, 1e-12_real64, &
                                "Original and loaded shift vectors should be identical")
    
    ! Clean up
    deallocate(shift_vectors_loaded)
  end subroutine test_write_read_shift_vectors

  subroutine test_read_write_gene_ids()
    integer(int32) :: ierr, n_loaded_genes, ndims, dims(1)
    character(len=256), allocatable :: loaded_gene_ids(:)
    call save_gene_ids(gene_ids, 'test_gene_ids.bin', ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'Failed to save gene IDs: ', ierr
      error stop
    end if
    call get_array_metadata('test_gene_ids.bin', dims, 1, ndims, ierr)
    if(.not. is_ok(ierr)) then
      write(*,*) 'Failed to get metadata: ', ierr
      error stop
    end if
    call assert_equal_int(ndims, 1, "Gene IDs should be 1D array")
    allocate(loaded_gene_ids(dims(1)))
    call load_gene_ids(loaded_gene_ids, 'test_gene_ids.bin', ierr)
    if(.not. is_ok(ierr)) then
      write(*,*) 'Failed to load gene IDs: ', ierr
      error stop
    end if
    call assert_equal_array_char(loaded_gene_ids, gene_ids, dims(1), &
                                "Loaded gene IDs should match original")
  end subroutine test_read_write_gene_ids

  subroutine test_read_write_gene_to_fam()
    integer(int32) :: ierr, n_loaded_genes, ndims, dims(1)
    integer(int32), allocatable :: loaded_gene_to_fam(:)
    call save_gene_to_family(gene_to_fam, 'test_gene_to_fam.bin', ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'Failed to save gene to family mapping: ', ierr
      error stop
    end if
    call get_array_metadata('test_gene_to_fam.bin', dims, 1, ndims, ierr)
    if(.not. is_ok(ierr)) then
      write(*,*) 'Failed to get metadata: ', ierr
      error stop
    end if
    call assert_equal_int(ndims, 1, "Gene to family mapping should be 1D array")
    allocate(loaded_gene_to_fam(dims(1)))
    call load_gene_to_family(loaded_gene_to_fam, 'test_gene_to_fam.bin', ierr)
    if(.not. is_ok(ierr)) then
      write(*,*) 'Failed to load gene to family mapping: ', ierr
      error stop
    end if
    call assert_equal_array_int(loaded_gene_to_fam, gene_to_fam, dims(1), &
                               "Loaded gene to family mapping should match original")
  end subroutine test_read_write_gene_to_fam

  subroutine test_read_write_family_ids()
    integer(int32) :: ierr, n_loaded_families, ndims, dims(1)
    character(len=256), allocatable :: loaded_family_ids(:)
    call save_family_ids(gene_family_ids, 'test_family_ids.bin', ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'Failed to save family IDs: ', ierr
      error stop
    end if
    call get_array_metadata('test_family_ids.bin', dims, 1, ndims, ierr)
    if(.not. is_ok(ierr)) then
      write(*,*) 'Failed to get metadata: ', ierr
      error stop
    end if
    call assert_equal_int(ndims, 1, "Family IDs should be 1D array")
    allocate(loaded_family_ids(dims(1)))
    call load_family_ids(loaded_family_ids, 'test_family_ids.bin', ierr)
    if(.not. is_ok(ierr)) then
      write(*,*) 'Failed to load family IDs: ', ierr
      error stop
    end if
    call assert_equal_array_char(loaded_family_ids, gene_family_ids, dims(1), &
                                "Loaded family IDs should match original")
  end subroutine test_read_write_family_ids

  subroutine test_read_write_centroids()
    integer(int32) :: ierr, n_loaded_families, ndims, dims(2)
    real(real64), allocatable :: loaded_centroids(:,:)
    call save_family_centroids(family_centroids, 'test_family_centroids.bin', ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'Failed to save family centroids: ', ierr
      error stop
    end if
    call get_array_metadata('test_family_centroids.bin', dims, 2, ndims, ierr)
    if(.not. is_ok(ierr)) then
      write(*,*) 'Failed to get metadata: ', ierr
      error stop
    end if
    call assert_equal_int(ndims, 2, "Family centroids should be 2D array")
    allocate(loaded_centroids(dims(1), dims(2)))
    call load_family_centroids(loaded_centroids, 'test_family_centroids.bin', ierr)
    if(.not. is_ok(ierr)) then
      write(*,*) 'Failed to load family centroids: ', ierr
      error stop
    end if
    call assert_equal_array_real(reshape(loaded_centroids, [size(loaded_centroids)]), &
                                reshape(family_centroids, [size(family_centroids)]), &
                                size(family_centroids), 1e-12_real64, &
                                "Loaded family centroids should match original")
  end subroutine test_read_write_centroids

  subroutine test_data_accessors()
    real(real64), allocatable :: returned_centroid(:)
    integer(int32) :: family_idx, ierr, gene_idx

    allocate(returned_centroid(total_samples))

    gene_idx = get_gene_index(gene_ids, 'NP_001000001.1')
    call assert_equal_int(gene_idx, 2, "Gene index for NP_001000001.1 should be 2")
    call get_family_for_gene_index(gene_idx, gene_to_fam, family_idx, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'Failed to get family for gene index: ', ierr
      error stop
    end if
    call assert_equal_int(family_idx, 48, "Family index for gene index 2")
    call assert_string_equal(gene_family_ids(family_idx), 'OG0000047', &
                        "Family ID for family index 48 should be OG0000047")

    family_idx = get_family_index(gene_family_ids, 'OG0000001')
    call assert_equal_int(family_idx, 2, "Family index for OG0000001 should be 2")
    call get_family_centroid(family_idx, family_centroids, returned_centroid, ierr)

    if (.not. is_ok(ierr)) then
      write(*,*) 'Failed to get family centroid: ', ierr
      error stop
    end if
    call assert_equal_array_real(returned_centroid, family_centroids(:, family_idx), &
                               size(family_centroids, 1), 1e-12_real64, &
                               "Retrieved centroid should match original")
  end subroutine test_data_accessors
end module mod_test_tox_data