!> Unit test suite for expression readers and data processing routines.
module mod_test_tox_data
  use asserts
  use iso_fortran_env, only: real64, int32
  use tox_data_tools
  use tox_data_validation
  use tox_data_read_write
  use array_utils, only: get_array_metadata
  use tox_gene_centroids
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

contains

  !> Get array of all available tests.
  function get_all_tests() result(all_tests)
    type(test_case) :: all_tests(6)
    all_tests(1) = test_case("test_read_gene_ids", test_read_gene_ids)
    all_tests(2) = test_case("test_read_expression_data", test_read_expression_data)
    all_tests(3) = test_case("test_read_family_mapping", test_read_family_mapping)
    all_tests(4) = test_case("test_validate_data", test_validate_data)
    all_tests(5) = test_case("test_compute_centroids", test_compute_centroids)
    all_tests(6) = test_case("test_write_read_binary", test_write_read_binary)
  end function get_all_tests

  !> Setup global test data
  subroutine setup_global_data()
    character(len=256), allocatable :: files_6_replicates(:)
    character(len=256), allocatable :: files_7_replicates(:)
    character(len=256), allocatable :: files_5_replicates(:)
    character(len=256), allocatable :: files_4_replicates(:)
    integer(int32) :: ierr, i
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
    n_genes = 88327
    n_families = 15512

    ! Allocate arrays
    allocate(gene_ids(n_genes))
    allocate(kallisto_expr(total_samples, n_genes))
    allocate(gene_family_ids(n_families))
    allocate(gene_to_fam(n_genes))

    ! Read gene IDs
    call read_gene_ids_from_file(files_6_replicates(1), gene_ids, 1, 1, ierr)
    call assert_equal_int(ierr, 0, "Reading gene IDs should succeed")

    ! Read expression data
    kallisto_expr = 0.0_real64
    call read_tabular_files(files_6_replicates, gene_ids, kallisto_expr, &
                           1, 1, [2, 3, 4, 5, 6, 7], 1, ierr)
    call assert_equal_int(ierr, 0, "Reading 6-replicate files should succeed")

    call read_tabular_files([files_7_replicates], gene_ids, kallisto_expr, &
                           1, 1, [2, 3, 4, 5, 6, 7, 8], 61, ierr)
    call assert_equal_int(ierr, 0, "Reading 7-replicate file should succeed")

    call read_tabular_files([files_5_replicates], gene_ids, kallisto_expr, &
                           1, 1, [2, 3, 4, 5, 6], 69, ierr)
    call assert_equal_int(ierr, 0, "Reading 5-replicate file should succeed")

    call read_tabular_files([files_4_replicates], gene_ids, kallisto_expr, &
                           1, 1, [2, 3, 4, 5], 73, ierr)
    call assert_equal_int(ierr, 0, "Reading 4-replicate file should succeed")

    ! Read family mapping
    call read_family_file('material/Orthogroups.tsv', gene_ids, gene_family_ids, gene_to_fam, ierr)
    call assert_equal_int(ierr, 0, "Reading family file should succeed")

    ! Compute centroids
    allocate(family_centroids(total_samples, n_families))
    
    ! Create explicit arrays instead of array constructors
    allocate(ortholog_mask(n_genes))
    allocate(selected_indices(n_genes))
    
    do i = 1, n_genes
      ortholog_mask(i) = .true.
      selected_indices(i) = i
    end do
    
    call group_centroid(kallisto_expr, total_samples, n_genes, gene_to_fam, &
                       n_families, family_centroids, .true., &
                       ortholog_mask, selected_indices, ierr)
    call assert_equal_int(ierr, 0, "Computing centroids should succeed")
    
    ! Clean up temporary arrays
    deallocate(ortholog_mask, selected_indices)
  end subroutine setup_global_data

  !> Run all expression reader tests.
  subroutine run_all_tests_tox_data()
    type(test_case) :: all_tests(6)
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
    type(test_case) :: all_tests(6)
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
    
    call validate_expression_data(kallisto_expr, .true., ierr)
    call assert_equal_int(ierr, 0, "Expression data validation should pass")
    
    call validate_gene_to_family_mapping(gene_to_fam, n_families, ierr)
    call assert_equal_int(ierr, 0, "Gene to family mapping validation should pass")
    
    call validate_gene_ids_uniqueness(gene_ids, ierr)
    call assert_equal_int(ierr, 0, "Gene IDs should be unique")
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
  subroutine test_write_read_binary()
    integer(int32) :: ierr, ndims, dims(2)
    integer(int32) :: total_elements
    
    call save_expression_vectors(kallisto_expr, ' test_kallisto_data.bin', ierr)
    call assert_equal_int(ierr, 0, "Saving expression data should succeed")
    
    call get_array_metadata(' test_kallisto_data.bin', dims, 2, ndims, ierr)
    call assert_equal_int(ierr, 0, "Reading metadata should succeed")
    call assert_equal_int(ndims, 2, "Array should have 2 dimensions")
    call assert_equal_int(dims(1), total_samples, "First dimension should match sample count")
    call assert_equal_int(dims(2), n_genes, "Second dimension should match gene count")
    
    allocate(kallisto_expr_verify(dims(1), dims(2)))
    call load_expression_vectors(kallisto_expr_verify, ' test_kallisto_data.bin', ierr)
    call assert_equal_int(ierr, 0, "Loading expression data should succeed")
    
    ! Fixed call to match the assert_equal_array_real interface
    total_elements = size(kallisto_expr)
    call assert_equal_array_real(reshape(kallisto_expr, [total_elements]), &
                                reshape(kallisto_expr_verify, [total_elements]), &
                                total_elements, 1e-12_real64, &
                                "Original and loaded data should be identical")
    
    ! Clean up
    deallocate(kallisto_expr_verify)
  end subroutine test_write_read_binary

end module mod_test_tox_data