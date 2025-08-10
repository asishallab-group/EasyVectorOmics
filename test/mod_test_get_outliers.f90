!> Unit test suite for tox_get_outliers routines.
module mod_test_get_outliers
  use asserts
  use tox_get_outliers
  use, intrinsic :: iso_fortran_env, only: real64, int32
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

contains

  !> Get array of all available tests (as subroutine, not function).
  subroutine get_all_tests(all_tests)
    type(test_case), intent(out) :: all_tests(31)
    ! Original tests
    all_tests(1) = test_case("test_scaling_basic", test_scaling_basic)
    all_tests(2) = test_case("test_loess_fallback", test_loess_fallback)
    all_tests(3) = test_case("test_rdi_basic", test_rdi_basic)
    all_tests(4) = test_case("test_identify_outliers_basic", test_identify_outliers_basic)
    all_tests(5) = test_case("test_detect_outliers_basic", test_detect_outliers_basic)
    all_tests(6) = test_case("test_detect_outliers_invalid_indices", test_detect_outliers_invalid_indices)
    all_tests(7) = test_case("test_invalid_indices", test_invalid_indices)
    all_tests(8) = test_case("test_all_outliers", test_all_outliers)
    all_tests(9) = test_case("test_no_outliers", test_no_outliers)
    all_tests(10) = test_case("test_single_gene_family", test_single_gene_family)
    all_tests(11) = test_case("test_all_zero_distances", test_all_zero_distances)
    all_tests(12) = test_case("test_all_genes_one_family", test_all_genes_one_family)
    all_tests(13) = test_case("test_negative_nan_distances", test_negative_nan_distances)
    all_tests(14) = test_case("test_identical_rdi_at_threshold", test_identical_rdi_at_threshold)
    all_tests(15) = test_case("test_single_gene_family_scaling", test_single_gene_family_scaling)
    all_tests(16) = test_case("test_all_negative_rdi", test_all_negative_rdi)
    ! Extended comprehensive tests
    all_tests(17) = test_case("test_median_calculation_odd", test_median_calculation_odd)
    all_tests(18) = test_case("test_median_calculation_even", test_median_calculation_even)
    all_tests(19) = test_case("test_median_two_elements", test_median_two_elements)
    all_tests(20) = test_case("test_scaling_with_extreme_outliers", test_scaling_with_extreme_outliers)
    all_tests(21) = test_case("test_mixed_family_sizes", test_mixed_family_sizes)
    all_tests(22) = test_case("test_loess_with_identical_medians", test_loess_with_identical_medians)
    all_tests(23) = test_case("test_rdi_calculation_precision", test_rdi_calculation_precision)
    all_tests(24) = test_case("test_large_dataset_performance", test_large_dataset_performance)
    all_tests(25) = test_case("test_zero_scaling_factors", test_zero_scaling_factors)
    all_tests(26) = test_case("test_identical_distances_all_genes", test_identical_distances_all_genes)
    all_tests(27) = test_case("test_high_variance_families", test_high_variance_families)
    all_tests(28) = test_case("test_loess_extrapolation", test_loess_extrapolation)
    all_tests(29) = test_case("test_rdi_with_tiny_scaling", test_rdi_with_tiny_scaling)
    all_tests(30) = test_case("test_outlier_detection_stability", test_outlier_detection_stability)
    all_tests(31) = test_case("test_detect_outliers_default_percentile", test_detect_outliers_default_percentile)
  end subroutine get_all_tests

  !> Run specific get_outliers tests by name.
  subroutine run_named_tests_get_outliers(test_names)
    character(len=*), intent(in) :: test_names(:)
    type(test_case) :: all_tests(31)
    integer(int32) :: i, j
    logical :: found

    call get_all_tests(all_tests)

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
  end subroutine run_named_tests_get_outliers

  !> Run all tox_get_outliers tests.
  subroutine run_all_tests_get_outliers()
    type(test_case) :: all_tests(31)
    integer(int32) :: i
    call get_all_tests(all_tests)
    do i = 1, size(all_tests)
      call all_tests(i)%test_proc()
      print *, trim(all_tests(i)%name), " passed."
    end do
    print *, "All tox_get_outliers tests passed successfully."
  end subroutine run_all_tests_get_outliers

  !> Test: Basic scaling with orthologs and non-orthologs.
  !> This test checks that the scaling for each family is set to the corresponding value using loess.
  subroutine test_scaling_basic()
    integer(int32), parameter :: n_genes=4, n_families=2
    real(real64) :: distances(n_genes) = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1,1,2,2]
    real(real64) :: dscale(n_families)
    integer(int32) :: error_code
    real(real64) :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: indices_used(n_families)
    integer(int32) :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    real(real64) :: family_distances(n_genes)
    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, dscale, &
      loess_x, loess_y, indices_used, perm, stack_left, stack_right, family_distances, error_code)
    call assert_equal_int(error_code, 0, 'Error code 0 (test_scaling_basic)')
    call assert_true(all(dscale > 0.0_real64), 'LOESS scaling should be positive for nontrivial families')
  end subroutine test_scaling_basic

  !> Test: Basic RDI calculation.
  !> This test checks that the RDI (relative distance index) is computed correctly for a simple case with two families.
  subroutine test_rdi_basic()
    integer(int32), parameter :: n_genes=3, n_families=2
    real(real64) :: distances(n_genes) = [2.0_real64, 4.0_real64, 6.0_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1,2,2]
    real(real64) :: dscale(n_families) = [2.0_real64, 4.0_real64]
    real(real64) :: rdi(n_genes), sorted_rdi(n_genes)
    integer(int32) :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    integer(int32) :: i
    perm = [(i, i=1,n_genes)]
    stack_left = 0
    stack_right = 0
    call compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi, sorted_rdi, perm, stack_left, stack_right)
    call assert_equal_real(rdi(1), 1.0_real64, 1e-12_real64, "RDI gene 1")
    call assert_equal_real(rdi(2), 1.0_real64, 1e-12_real64, "RDI gene 2")
    call assert_equal_real(rdi(3), 1.5_real64, 1e-12_real64, "RDI gene 3")
  end subroutine test_rdi_basic

  !> Test: Outlier detection in a simple RDI vector.
  !> This test checks that the identify_outliers routine correctly identifies the highest and second highest RDI values as outliers at the 75th percentile.
  subroutine test_identify_outliers_basic()
    integer(int32), parameter :: n_genes=4
    real(real64) :: rdi(n_genes), sorted_rdi(n_genes)
    logical :: is_outlier(n_genes)
    real(real64) :: threshold
    integer(int32) :: i
    rdi = [0.1_real64, 0.2_real64, 0.3_real64, 0.4_real64]
    ! rdi = [0.1, 0.2, 0.3, 0.4] already sorted and no negatives
    sorted_rdi = [0.1_real64, 0.2_real64, 0.3_real64, 0.4_real64]
    call identify_outliers(n_genes, rdi, sorted_rdi, is_outlier, threshold, 75.0_real64)
    call assert_true(is_outlier(4), "Highest RDI is outlier")
    call assert_true(is_outlier(3), "Second highest RDI is outlier")
    call assert_false(is_outlier(1), "Lowest RDI is not outlier")
  end subroutine test_identify_outliers_basic

  !> Test: Full outlier detection workflow.
  !> This test checks the complete workflow of outlier detection, including scaling, RDI calculation, and outlier identification, for a small example with two families.
  subroutine test_detect_outliers_basic()
    integer(int32), parameter :: n_genes=5, n_families=2
    real(real64) :: distances(n_genes)
    integer(int32) :: gene_to_fam(n_genes)
    real(real64) :: work_array(n_genes)
    integer(int32) :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    logical :: is_outlier(n_genes)
    integer(int32) :: error_code
    real(real64) :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: loess_n(n_families)
    distances = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64]
    gene_to_fam = [1,1,2,2,2]
    call detect_outliers(n_genes, n_families, distances, gene_to_fam, work_array, perm, stack_left, stack_right, &
      is_outlier, loess_x, loess_y, loess_n, error_code, &
      60.0_real64)
    call assert_true(any(is_outlier), "At least one outlier detected")
  end subroutine test_detect_outliers_basic

  !> Test: Invalid family indices.
  !> This test checks that if gene_to_fam contains invalid indices (e.g., 0 or out of range), the error code is set to -2 and dscale falls back to -1.0 for all families.
  subroutine test_invalid_indices()
    integer(int32), parameter :: n_genes=3, n_families=2
    real(real64) :: distances(n_genes) = [1.0_real64, 2.0_real64, 3.0_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1,3,0] ! 3 and 0 invalid
    real(real64) :: dscale(n_families)
    integer(int32) :: error_code
    real(real64) :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: indices_used(n_families), loess_n(n_families)
    integer(int32) :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    real(real64) :: family_distances(n_genes)
    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, dscale, &
      loess_x, loess_y, indices_used, perm, stack_left, stack_right, family_distances, error_code)
    call assert_equal_int(error_code, -2, 'Error code -2 for invalid family indices')
    call assert_true(all(dscale == -1.0_real64), 'dscale must be -1.0 on error')
  end subroutine test_invalid_indices

  !> Test: All genes are outliers (percentile 0).
  !> This test checks that if the outlier percentile is set to 0, all genes are identified as outliers.
  subroutine test_all_outliers()
    integer(int32), parameter :: n_genes=4
    real(real64) :: rdi(n_genes) = [10.0_real64, 11.0_real64, 12.0_real64, 13.0_real64]
    real(real64) :: sorted_rdi(n_genes)
    logical :: is_outlier(n_genes)
    real(real64) :: threshold
    ! rdi = [10, 11, 12, 13] already sorted and no negatives
    sorted_rdi = [10.0_real64, 11.0_real64, 12.0_real64, 13.0_real64]
    call identify_outliers(n_genes, rdi, sorted_rdi, is_outlier, threshold, 0.0_real64)
    call assert_true(all(is_outlier), "All are outliers at 0 percentile")
  end subroutine test_all_outliers

  !> Test: Only highest gene is outlier (percentile 100).
  !> This test checks that if the outlier percentile is set to 100, only the highest gene is identified as outlier.
  subroutine test_no_outliers()
    integer(int32), parameter :: n_genes=4
    real(real64) :: rdi(n_genes) = [0.1_real64, 0.2_real64, 0.3_real64, 0.4_real64]
    real(real64) :: sorted_rdi(n_genes)
    logical :: is_outlier(n_genes)
    real(real64) :: threshold
    ! rdi = [0.1, 0.2, 0.3, 0.4] already sorted and no negatives
    sorted_rdi = [0.1_real64, 0.2_real64, 0.3_real64, 0.4_real64]
    call identify_outliers(n_genes, rdi, sorted_rdi, is_outlier, threshold, 100.0_real64)
    call assert_true(is_outlier(4), "Highest value is outlier at 100 percentile")
    call assert_true(.not. any(is_outlier(1:3)), "Others are not outliers at 100 percentile")
  end subroutine test_no_outliers

  !> Test: Single gene and family case.
  !> This test checks that for a single gene in a single family, scaling and RDI are both 0.0, and the gene is not an outlier.
  subroutine test_single_gene_family()
    integer(int32), parameter :: n_genes=1, n_families=1
    real(real64) :: distances(1) = [0.0_real64]
    integer(int32) :: gene_to_fam(1) = [1]
    real(real64) :: dscale(1), rdi(1), sorted_rdi(1), work_array(1)
    integer(int32) :: perm(1), stack_left(1), stack_right(1)
    logical :: is_outlier(1)
    integer(int32) :: i, error_code
    real(real64) :: loess_x(1), loess_y(1)
    integer(int32) :: loess_n(1)
    real(real64) :: family_distances(n_genes)
    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, dscale, &
      loess_x, loess_y, loess_n, perm, stack_left, stack_right, family_distances, error_code)
    perm = [(i, i=1,n_genes)]
    stack_left = 0
    stack_right = 0
    call compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi, sorted_rdi, perm, stack_left, stack_right)
    perm = [(i, i=1,n_genes)]
    call detect_outliers(n_genes, n_families, distances, gene_to_fam, work_array, perm, stack_left, stack_right, &
      is_outlier, loess_x, loess_y, loess_n, error_code, &
      95.0_real64)
    call assert_equal_real(dscale(1), 0.0_real64, 1e-12_real64, &
         "Single gene scaling (should be 0.0)")
    call assert_equal_real(rdi(1), 0.0_real64, 1e-12_real64, &
         "Single gene RDI")
    call assert_false(is_outlier(1), "Single gene not outlier")
  end subroutine test_single_gene_family

  !> Test: All distances are zero.
  !> This test checks that if all distances are zero, scaling is 0.0 and all RDI values are 0.0.
  subroutine test_all_zero_distances()
    integer(int32), parameter :: n_genes=3, n_families=1
    real(real64) :: distances(n_genes) = [0.0_real64, 0.0_real64, 0.0_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1,1,1]
    logical :: is_ortholog(n_genes) = [.true., .true., .true.]
    real(real64) :: dscale(1), rdi(n_genes), sorted_rdi(n_genes)
    integer(int32) :: error_code
    real(real64) :: max_distance_bw_orths(1) = [0.0_real64]
    real(real64) :: loess_x(1), loess_y(1), loess_n(1)
    integer(int32) :: indices_used(1)
    integer(int32) :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    real(real64) :: family_distances(n_genes)
    integer(int32) :: i
    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, dscale, &
      loess_x, loess_y, indices_used, perm, stack_left, stack_right, family_distances, error_code)
    perm = [(i, i=1,n_genes)]
    stack_left = 0
    stack_right = 0
    call compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi, sorted_rdi, perm, stack_left, stack_right)
    call assert_equal_real(dscale(1), 0.0_real64, 1e-12_real64, &
         "All zero distances scaling (should be 0.0)")
    call assert_true(all(rdi == 0.0_real64), "All RDI zero for zero distances")
  end subroutine test_all_zero_distances

  !> Test: All genes in a single family.
  !> This test checks that if all genes belong to a single family, scaling is set to the value in max_distance_bw_orths if orthologs are present.
  subroutine test_all_genes_one_family()
    integer(int32), parameter :: n_genes=4, n_families=1
    real(real64) :: distances(n_genes) = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1,1,1,1]
    real(real64) :: dscale(1)
    integer(int32) :: error_code
    real(real64) :: loess_x(1), loess_y(1), loess_n(1)
    integer(int32) :: indices_used(1)
    integer(int32) :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    real(real64) :: family_distances(n_genes)
    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, dscale, &
      loess_x, loess_y, indices_used, perm, stack_left, stack_right, family_distances, error_code)
    call assert_true(dscale(1) > 0.0_real64, "All genes one family scaling (LOESS, should be positive)")
  end subroutine test_all_genes_one_family

  !> Test: Negative distances only.
  !> This test checks that negative distances with zero scaling result in RDI=0.0 and are not outliers. NaN test is skipped for portability.
  subroutine test_negative_nan_distances()
    integer(int32), parameter :: n_genes=3, n_families=1
    real(real64) :: distances(n_genes)
    integer(int32) :: gene_to_fam(n_genes) = [1,1,1]
    real(real64) :: dscale(1), rdi(n_genes), sorted_rdi(n_genes)
    integer(int32) :: error_code
    real(real64) :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: indices_used(n_families), loess_n(n_families)
    integer(int32) :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    real(real64) :: family_distances(n_genes)
    integer(int32) :: i
    distances = [0.0_real64, 0.0_real64, 0.0_real64]
    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, dscale, &
      loess_x, loess_y, indices_used, perm, stack_left, stack_right, family_distances, error_code)
    distances(1) = -1.0_real64
    perm = [(i, i=1,n_genes)]
    stack_left = 0
    stack_right = 0
    call compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi, sorted_rdi, perm, stack_left, stack_right)
    call assert_true(rdi(1) == 0.0_real64, "Negative distance with zero scaling gives RDI=0.0 (not outlier)")
    ! NaN test skipped for portability
  end subroutine test_negative_nan_distances

  !> Test: Multiple genes with identical RDI at threshold.
  !> This test checks that if multiple genes have identical RDI at the threshold, they are all identified as outliers.
  subroutine test_identical_rdi_at_threshold()
    integer(int32), parameter :: n_genes=4
    real(real64) :: rdi(n_genes) = [1.0_real64, 2.0_real64, 2.0_real64, 3.0_real64]
    real(real64) :: sorted_rdi(n_genes)
    logical :: is_outlier(n_genes)
    real(real64) :: threshold
    integer(int32) :: i
    ! rdi = [1, 2, 2, 3] already sorted and no negatives
    sorted_rdi = [1.0_real64, 2.0_real64, 2.0_real64, 3.0_real64]
    call identify_outliers(n_genes, rdi, sorted_rdi, is_outlier, threshold, 50.0_real64)
    call assert_true(is_outlier(3), "Identical RDI at threshold is outlier")
    call assert_true(is_outlier(2), "Identical RDI at threshold is outlier")
  end subroutine test_identical_rdi_at_threshold

  !> Exhaustive: single-gene family (should scale to 0.0).
  !> This test ensures that for families with only one gene, the scaling is set to 0.0 and error code is 0.
  subroutine test_single_gene_family_scaling()
    integer(int32), parameter :: n_genes=2, n_families=2
    real(real64) :: distances(n_genes) = [1.0_real64, 2.0_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1,2]
    real(real64) :: dscale(n_families)
    integer(int32) :: error_code
    real(real64) :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: indices_used(n_families), loess_n(n_families)
    integer(int32) :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    real(real64) :: family_distances(n_genes)
    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, dscale, &
      loess_x, loess_y, indices_used, perm, stack_left, stack_right, family_distances, error_code)
    call assert_equal_int(error_code, 0, 'Error code 0 for single-gene families')
    call assert_equal_real(dscale(1), 0.0_real64, 1e-12_real64, 'Single-gene family scaling 0.0')
    call assert_equal_real(dscale(2), 0.0_real64, 1e-12_real64, 'Single-gene family scaling 0.0')
  end subroutine test_single_gene_family_scaling

  !> Exhaustive: detect_outliers error propagation (invalid indices).
  !> This test checks that detect_outliers propagates error_code -2 if gene_to_fam contains invalid indices.
  subroutine test_detect_outliers_invalid_indices()
    integer(int32), parameter :: n_genes=3, n_families=2
    real(real64) :: distances(n_genes) = [1.0_real64, 2.0_real64, 3.0_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1,3,0]
    real(real64) :: work_array(n_genes)
    integer(int32) :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    logical :: is_outlier(n_genes)
    integer(int32) :: error_code
    logical :: is_ortholog(n_genes) = [.false., .false., .false.]
    real(real64) :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: indices_used(n_families), loess_n(n_families)
    call detect_outliers(n_genes, n_families, distances, gene_to_fam, work_array, perm, stack_left, stack_right, &
      is_outlier, loess_x, loess_y, loess_n, error_code)
    call assert_equal_int(error_code, -2, 'detect_outliers propagates error_code -2 for invalid indices')
  end subroutine test_detect_outliers_invalid_indices

  !> Exhaustive: LOESS fallback for scaling (workspace arrays required).
  subroutine test_loess_fallback()
    integer(int32), parameter :: n_genes=4, n_families=2
    real(real64) :: distances(n_genes) = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1,1,2,2]
    real(real64) :: dscale(n_families)
    integer(int32) :: error_code
    real(real64) :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: indices_used(n_families),loess_n(n_families)
    integer(int32) :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    real(real64) :: family_distances(n_genes)
    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, dscale, &
      loess_x, loess_y, indices_used, perm, stack_left, stack_right, family_distances, error_code)
    call assert_equal_int(error_code, 0, 'Error code 0 for LOESS fallback')
    call assert_true(any(dscale > 0.0_real64), 'Scaling must be positive if LOESS fallback is successful')
  end subroutine test_loess_fallback

  !> Test: identify_outliers uses default percentile (95) if not provided.
  subroutine test_identify_outliers_default_percentile()
    integer(int32), parameter :: n_genes=4
    real(real64) :: rdi(n_genes) = [0.1_real64, 0.2_real64, 0.3_real64, 0.4_real64]
    real(real64) :: sorted_rdi(n_genes)
    logical :: is_outlier(n_genes)
    real(real64) :: threshold
    ! rdi = [0.1, 0.2, 0.3, 0.4] ya está ordenado y no hay negativos
    sorted_rdi = [0.1_real64, 0.2_real64, 0.3_real64, 0.4_real64]
    call identify_outliers(n_genes, rdi, sorted_rdi, is_outlier, threshold)
    ! At 95th percentile, only the highest value should be outlier
    call assert_true(is_outlier(4), "Highest RDI is outlier (default percentile)")
    call assert_false(any(is_outlier(1:3)), "Others are not outliers (default percentile)")
  end subroutine test_identify_outliers_default_percentile

  !> Test: detect_outliers uses default percentile (95) if not provided.
  subroutine test_detect_outliers_default_percentile()
    integer(int32), parameter :: n_genes=4, n_families=1
    real(real64) :: distances(n_genes) = [0.1_real64, 0.2_real64, 0.3_real64, 0.4_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1,1,1,1]
    real(real64) :: work_array(n_genes)
    integer(int32) :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    logical :: is_outlier(n_genes)
    integer(int32) :: error_code
    real(real64) :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: loess_n(n_families)
    call detect_outliers(n_genes, n_families, distances, gene_to_fam, work_array, perm, stack_left, stack_right, &
      is_outlier, loess_x, loess_y, loess_n, error_code)
    ! At 95th percentile, only the highest value should be outlier
    call assert_true(is_outlier(4), "Highest distance is outlier (default percentile)")
    call assert_false(any(is_outlier(1:3)), "Others are not outliers (default percentile)")
    call assert_equal_int(error_code, 0, "No error for default percentile")
  end subroutine test_detect_outliers_default_percentile

  !> Test: All RDI negative (should be ignored for outlier detection).
  !> This test checks that if all RDI values are negative, none are marked as outliers.
  subroutine test_all_negative_rdi()
    integer(int32), parameter :: n_genes=5
    real(real64) :: rdi(n_genes) = [-0.1_real64, -0.2_real64, -0.3_real64, -0.4_real64, -0.5_real64]
    real(real64) :: sorted_rdi(n_genes)
    logical :: is_outlier(n_genes)
    real(real64) :: threshold
    integer(int32) :: i
    ! rdi = [-0.1, -0.2, -0.3, -0.4, -0.5] -> sorted_rdi = [0,0,0,0,0] (all negatives)
    sorted_rdi = [0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64]
    call identify_outliers(n_genes, rdi, sorted_rdi, is_outlier, threshold, 80.0_real64)
    call assert_true(.not. any(is_outlier), "All-negative RDI should not be outliers")
  end subroutine test_all_negative_rdi

  !> Test median calculation with odd number of elements
  subroutine test_median_calculation_odd()
    integer(int32), parameter :: n_genes = 5, n_families = 1
    real(real64) :: distances(n_genes) = [1.0_real64, 3.0_real64, 2.0_real64, 5.0_real64, 4.0_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1, 1, 1, 1, 1]
    real(real64) :: dscale(n_families)
    real(real64) :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: indices_used(n_families)
    integer(int32) :: perm_tmp(n_genes), stack_left_tmp(n_genes), stack_right_tmp(n_genes)
    real(real64) :: family_distances(n_genes)
    integer(int32) :: error_code

    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, dscale, &
                                loess_x, loess_y, indices_used, perm_tmp, stack_left_tmp, stack_right_tmp, family_distances, error_code)

    ! Sorted distances should be [1, 2, 3, 4, 5], median = 3.0
    call assert_equal_real(loess_x(1), 3.0_real64, 1e-10_real64, "Median calculation for odd n")
    call assert_equal_int(error_code, 0, "No error for median odd test")
  end subroutine test_median_calculation_odd

  !> Test median calculation with even number of elements  
  subroutine test_median_calculation_even()
    integer(int32), parameter :: n_genes = 4, n_families = 1
    real(real64) :: distances(n_genes) = [1.0_real64, 4.0_real64, 2.0_real64, 3.0_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1, 1, 1, 1]
    real(real64) :: dscale(n_families)
    real(real64) :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: indices_used(n_families)
    integer(int32) :: perm_tmp(n_genes), stack_left_tmp(n_genes), stack_right_tmp(n_genes)
    real(real64) :: family_distances(n_genes)
    integer(int32) :: error_code

    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, dscale, &
                                loess_x, loess_y, indices_used, perm_tmp, stack_left_tmp, stack_right_tmp, family_distances, error_code)

    ! Sorted distances should be [1, 2, 3, 4], median = (2+3)/2 = 2.5
    call assert_equal_real(loess_x(1), 2.5_real64, 1e-10_real64, "Median calculation for even n")
    call assert_equal_int(error_code, 0, "No error for median even test")
  end subroutine test_median_calculation_even

  !> Test median with exactly two elements
  subroutine test_median_two_elements()
    integer(int32), parameter :: n_genes = 2, n_families = 1
    real(real64) :: distances(n_genes) = [10.0_real64, 2.0_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1, 1]
    real(real64) :: dscale(n_families)
    real(real64) :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: indices_used(n_families)
    integer(int32) :: perm_tmp(n_genes), stack_left_tmp(n_genes), stack_right_tmp(n_genes)
    real(real64) :: family_distances(n_genes)
    integer(int32) :: error_code

    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, dscale, &
                                loess_x, loess_y, indices_used, perm_tmp, stack_left_tmp, stack_right_tmp, family_distances, error_code)

    ! Median = (2.0 + 10.0)/2 = 6.0
    call assert_equal_real(loess_x(1), 6.0_real64, 1e-10_real64, "Median calculation for two elements")
    call assert_equal_int(error_code, 0, "No error for two elements test")
  end subroutine test_median_two_elements

  !> Test scaling with extreme outliers in family
  subroutine test_scaling_with_extreme_outliers()
    integer(int32), parameter :: n_genes = 6, n_families = 1
    ! Family with normal distances and one extreme outlier
    real(real64) :: distances(n_genes) = [1.0_real64, 1.1_real64, 0.9_real64, 1.05_real64, 100.0_real64, 1.2_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1, 1, 1, 1, 1, 1]
    real(real64) :: dscale(n_families)
    real(real64) :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: indices_used(n_families)
    integer(int32) :: perm_tmp(n_genes), stack_left_tmp(n_genes), stack_right_tmp(n_genes)
    real(real64) :: family_distances(n_genes)
    integer(int32) :: error_code

    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, dscale, &
                                loess_x, loess_y, indices_used, perm_tmp, stack_left_tmp, stack_right_tmp, family_distances, error_code)

    ! Should not crash and should handle extreme outlier gracefully
    call assert_equal_int(error_code, 0, "No error with extreme outliers")
    call assert_true(dscale(1) > 0.0_real64, "Positive scaling factor with outliers")
    call assert_true(loess_y(1) > 10.0_real64, "High standard deviation detected")
  end subroutine test_scaling_with_extreme_outliers

  !> Test with mixed family sizes (some large, some small)
  subroutine test_mixed_family_sizes()
    integer(int32), parameter :: n_genes = 10, n_families = 3
    real(real64) :: distances(n_genes) = [1.0_real64, 1.1_real64, 1.2_real64, 1.3_real64, 1.4_real64, &
                                          2.0_real64, 2.1_real64, 3.0_real64, 4.0_real64, 5.0_real64]
    ! Family 1: 5 genes, Family 2: 2 genes, Family 3: 3 genes
    integer(int32) :: gene_to_fam(n_genes) = [1, 1, 1, 1, 1, 2, 2, 3, 3, 3]
    real(real64) :: dscale(n_families)
    real(real64) :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: indices_used(n_families)
    integer(int32) :: perm_tmp(n_genes), stack_left_tmp(n_genes), stack_right_tmp(n_genes)
    real(real64) :: family_distances(n_genes)
    integer(int32) :: error_code

    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, dscale, &
                                loess_x, loess_y, indices_used, perm_tmp, stack_left_tmp, stack_right_tmp, family_distances, error_code)

    call assert_equal_int(error_code, 0, "No error with mixed family sizes")
    ! All families should have positive scaling factors
    call assert_true(all(dscale > 0.0_real64), "All scaling factors positive")
    ! Family 1 median should be 1.2 (middle of 1.0, 1.1, 1.2, 1.3, 1.4)
    call assert_equal_real(loess_x(1), 1.2_real64, 1e-10_real64, "Family 1 median correct")
    ! Family 2 median should be 2.05 (average of 2.0, 2.1)
    call assert_equal_real(loess_x(2), 2.05_real64, 1e-10_real64, "Family 2 median correct")
  end subroutine test_mixed_family_sizes

  !> Test LOESS behavior when all families have identical medians
  subroutine test_loess_with_identical_medians()
    integer(int32), parameter :: n_genes = 6, n_families = 3
    real(real64) :: distances(n_genes) = [5.0_real64, 5.0_real64, 5.0_real64, 5.0_real64, 5.0_real64, 5.0_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1, 1, 2, 2, 3, 3]
    real(real64) :: dscale(n_families)
    real(real64) :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: indices_used(n_families)
    integer(int32) :: perm_tmp(n_genes), stack_left_tmp(n_genes), stack_right_tmp(n_genes)
    real(real64) :: family_distances(n_genes)
    integer(int32) :: error_code

    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, dscale, &
                                loess_x, loess_y, indices_used, perm_tmp, stack_left_tmp, stack_right_tmp, family_distances, error_code)

    call assert_equal_int(error_code, 0, "No error with identical medians")
    ! All medians should be 5.0
    call assert_true(all(abs(loess_x - 5.0_real64) < 1e-10_real64), "All medians equal to 5.0")
    ! All stddevs should be 0.0
    call assert_true(all(abs(loess_y) < 1e-10_real64), "All stddevs equal to 0.0")
  end subroutine test_loess_with_identical_medians

  !> Test RDI calculation precision with known values
  subroutine test_rdi_calculation_precision()
    integer(int32), parameter :: n_genes = 3
    real(real64) :: distances(n_genes) = [2.0_real64, 4.0_real64, 6.0_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1, 1, 1]
    real(real64) :: dscale(1) = [2.0_real64]  ! Known scaling factor
    real(real64) :: rdi(n_genes), sorted_rdi(n_genes)
    integer(int32) :: perm(n_genes) = [1, 2, 3]
    integer(int32) :: stack_left(n_genes), stack_right(n_genes)

    call compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi, sorted_rdi, perm, stack_left, stack_right)

    ! Expected RDI: [2.0/2.0, 4.0/2.0, 6.0/2.0] = [1.0, 2.0, 3.0]
    call assert_equal_real(rdi(1), 1.0_real64, 1e-10_real64, "RDI calculation gene 1")
    call assert_equal_real(rdi(2), 2.0_real64, 1e-10_real64, "RDI calculation gene 2")
    call assert_equal_real(rdi(3), 3.0_real64, 1e-10_real64, "RDI calculation gene 3")
  end subroutine test_rdi_calculation_precision

  !> Test percentile edge cases (0%, 50%, 100%)
  subroutine test_percentile_edge_cases()
    integer(int32), parameter :: n_genes = 5
    real(real64) :: rdi(n_genes) = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64]
    real(real64) :: sorted_rdi(n_genes) = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64]
    logical :: is_outlier(n_genes)
    real(real64) :: threshold

    ! Test 0th percentile (all outliers)
    call identify_outliers(n_genes, rdi, sorted_rdi, is_outlier, threshold, 0.0_real64)
    call assert_true(all(is_outlier), "0th percentile: all outliers")

    ! Test 100th percentile (only highest is outlier)
    call identify_outliers(n_genes, rdi, sorted_rdi, is_outlier, threshold, 100.0_real64)
    call assert_true(is_outlier(5), "100th percentile: highest value is outlier")
    call assert_true(.not. any(is_outlier(1:4)), "100th percentile: others are not outliers")

    ! Test 50th percentile (median split)
    call identify_outliers(n_genes, rdi, sorted_rdi, is_outlier, threshold, 50.0_real64)
    call assert_true(is_outlier(4) .and. is_outlier(5), "50th percentile: top half outliers")
    call assert_true(.not. is_outlier(1) .and. .not. is_outlier(2), "50th percentile: bottom half not outliers")
  end subroutine test_percentile_edge_cases

  !> Test performance with larger dataset
  subroutine test_large_dataset_performance()
    integer(int32), parameter :: n_genes = 1000, n_families = 50
    real(real64) :: distances(n_genes)
    integer(int32) :: gene_to_fam(n_genes)
    real(real64) :: work_array(n_genes)
    integer(int32) :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    logical :: is_outlier(n_genes)
    real(real64) :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: loess_n(n_families)
    integer(int32) :: error_code
    integer(int32) :: i

    ! Generate synthetic data: 20 genes per family, with some noise
    do i = 1, n_genes
      gene_to_fam(i) = (i-1) / 20 + 1  ! Assign to families
      distances(i) = real(gene_to_fam(i), real64) + 0.1_real64 * real(mod(i, 10), real64)
    end do

    call detect_outliers(n_genes, n_families, distances, gene_to_fam, &
                        work_array, perm, stack_left, stack_right, &
                        is_outlier, loess_x, loess_y, loess_n, error_code, 95.0_real64)

    call assert_equal_int(error_code, 0, "Large dataset: no error")
    call assert_true(any(is_outlier), "Large dataset: some outliers detected")
    call assert_true(.not. all(is_outlier), "Large dataset: not all are outliers")
  end subroutine test_large_dataset_performance

  !> Test behavior when scaling factors are zero
  subroutine test_zero_scaling_factors()
    integer(int32), parameter :: n_genes = 3
    real(real64) :: distances(n_genes) = [1.0_real64, 2.0_real64, 3.0_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1, 1, 1]
    real(real64) :: dscale(1) = [0.0_real64]  ! Zero scaling factor
    real(real64) :: rdi(n_genes), sorted_rdi(n_genes)
    integer(int32) :: perm(n_genes) = [1, 2, 3]
    integer(int32) :: stack_left(n_genes), stack_right(n_genes)

    call compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi, sorted_rdi, perm, stack_left, stack_right)

    ! With zero scaling, all RDI should be 0.0 (not outliers)
    call assert_true(all(abs(rdi) < 1e-10_real64), "Zero scaling: all RDI are zero")
  end subroutine test_zero_scaling_factors

  !> Test all genes having identical distances
  subroutine test_identical_distances_all_genes()
    integer(int32), parameter :: n_genes = 6, n_families = 2
    real(real64) :: distances(n_genes) = [5.0_real64, 5.0_real64, 5.0_real64, 5.0_real64, 5.0_real64, 5.0_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1, 1, 1, 2, 2, 2]
    real(real64) :: work_array(n_genes)
    integer(int32) :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    logical :: is_outlier(n_genes)
    real(real64) :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: loess_n(n_families)
    integer(int32) :: error_code

    call detect_outliers(n_genes, n_families, distances, gene_to_fam, &
                        work_array, perm, stack_left, stack_right, &
                        is_outlier, loess_x, loess_y, loess_n, error_code, 95.0_real64)

    call assert_equal_int(error_code, 0, "Identical distances: no error")
    ! With identical distances and zero variance, no outliers should be detected
    call assert_true(.not. any(is_outlier), "Identical distances: no outliers")
  end subroutine test_identical_distances_all_genes

  !> Test families with very high variance
  subroutine test_high_variance_families()
    integer(int32), parameter :: n_genes = 4, n_families = 1
    real(real64) :: distances(n_genes) = [0.1_real64, 1000.0_real64, 0.2_real64, 999.0_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1, 1, 1, 1]
    real(real64) :: dscale(n_families)
    real(real64) :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: indices_used(n_families)
    integer(int32) :: perm_tmp(n_genes), stack_left_tmp(n_genes), stack_right_tmp(n_genes)
    real(real64) :: family_distances(n_genes)
    integer(int32) :: error_code

    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, dscale, &
                                loess_x, loess_y, indices_used, perm_tmp, stack_left_tmp, stack_right_tmp, family_distances, error_code)

    call assert_equal_int(error_code, 0, "High variance: no error")
    call assert_true(loess_y(1) > 100.0_real64, "High variance detected")
    call assert_true(dscale(1) > 0.0_real64, "Positive scaling with high variance")
  end subroutine test_high_variance_families

  !> Test LOESS extrapolation beyond reference range
  subroutine test_loess_extrapolation()
    integer(int32), parameter :: n_genes = 4, n_families = 2
    ! Create two families with very different median distances
    real(real64) :: distances(n_genes) = [1.0_real64, 1.1_real64, 100.0_real64, 101.0_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1, 1, 2, 2]
    real(real64) :: dscale(n_families)
    real(real64) :: loess_x(n_families), loess_y(n_families)
    integer(int32) :: indices_used(n_families)
    integer(int32) :: perm_tmp(n_genes), stack_left_tmp(n_genes), stack_right_tmp(n_genes)
    real(real64) :: family_distances(n_genes)
    integer(int32) :: error_code

    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, dscale, &
                                loess_x, loess_y, indices_used, perm_tmp, stack_left_tmp, stack_right_tmp, family_distances, error_code)

    call assert_equal_int(error_code, 0, "LOESS extrapolation: no error")
    call assert_true(all(dscale > 0.0_real64), "LOESS extrapolation: positive scaling factors")
    ! Both families should get reasonable scaling factors despite wide separation
    call assert_true(abs(loess_x(1) - 1.05_real64) < 0.1_real64, "Family 1 median reasonable")
    call assert_true(abs(loess_x(2) - 100.5_real64) < 0.1_real64, "Family 2 median reasonable")
  end subroutine test_loess_extrapolation

  !> Test RDI with very small scaling factors (near machine epsilon)
  subroutine test_rdi_with_tiny_scaling()
    integer(int32), parameter :: n_genes = 2
    real(real64) :: distances(n_genes) = [1.0_real64, 2.0_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1, 1]
    real(real64) :: dscale(1) = [1e-17_real64]  ! Smaller than machine epsilon
    real(real64) :: rdi(n_genes), sorted_rdi(n_genes)
    integer(int32) :: perm(n_genes) = [1, 2]
    integer(int32) :: stack_left(n_genes), stack_right(n_genes)

    call compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi, sorted_rdi, perm, stack_left, stack_right)

    ! Should handle tiny scaling factors gracefully, setting RDI to 0 
    call assert_true(all(abs(rdi) < 1e-10_real64), "Tiny scaling: RDI set to zero")
  end subroutine test_rdi_with_tiny_scaling

  !> Test stability of outlier detection with repeated calls
  subroutine test_outlier_detection_stability()
    integer(int32), parameter :: n_genes = 10, n_families = 2
    real(real64) :: distances(n_genes) = [1.0_real64, 1.1_real64, 1.2_real64, 1.3_real64, 1.4_real64, &
                                          10.0_real64, 10.1_real64, 10.2_real64, 10.3_real64, 100.0_real64]
    integer(int32) :: gene_to_fam(n_genes) = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2]
    real(real64) :: work_array1(n_genes), work_array2(n_genes)
    integer(int32) :: perm1(n_genes), perm2(n_genes)
    integer(int32) :: stack_left1(n_genes), stack_right1(n_genes)
    integer(int32) :: stack_left2(n_genes), stack_right2(n_genes)
    logical :: is_outlier1(n_genes), is_outlier2(n_genes)
    real(real64) :: loess_x1(n_families), loess_y1(n_families)
    real(real64) :: loess_x2(n_families), loess_y2(n_families)
    integer(int32) :: loess_n1(n_families), loess_n2(n_families)
    integer(int32) :: error_code1, error_code2

    ! Run detection twice with same data
    call detect_outliers(n_genes, n_families, distances, gene_to_fam, &
                        work_array1, perm1, stack_left1, stack_right1, &
                        is_outlier1, loess_x1, loess_y1, loess_n1, error_code1, 90.0_real64)

    call detect_outliers(n_genes, n_families, distances, gene_to_fam, &
                        work_array2, perm2, stack_left2, stack_right2, &
                        is_outlier2, loess_x2, loess_y2, loess_n2, error_code2, 90.0_real64)

    ! Results should be identical
    call assert_equal_int(error_code1, error_code2, "Stability: same error codes")
    call assert_true(all(is_outlier1 .eqv. is_outlier2), "Stability: identical outlier detection")
    call assert_true(all(abs(loess_x1 - loess_x2) < 1e-10_real64), "Stability: identical LOESS x")
    call assert_true(all(abs(loess_y1 - loess_y2) < 1e-10_real64), "Stability: identical LOESS y")
  end subroutine test_outlier_detection_stability


end module mod_test_get_outliers


