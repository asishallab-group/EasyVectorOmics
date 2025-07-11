!> @brief Unit test suite for tox_get_outliers routines.
module mod_test_get_outliers
  use asserts
  use tox_get_outliers
  use, intrinsic :: iso_fortran_env, only: real64, int64
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

  !> @brief Get array of all available tests (as subroutine, not function).
  subroutine get_all_tests(all_tests)
    type(test_case), intent(out) :: all_tests(16)
    all_tests(1) = test_case("test_scaling_basic", test_scaling_basic)
    all_tests(2) = test_case("test_scaling_no_orthologs", test_scaling_no_orthologs)
    all_tests(3) = test_case("test_rdi_basic", test_rdi_basic)
    all_tests(4) = test_case("test_identify_outliers_basic", test_identify_outliers_basic)
    all_tests(5) = test_case("test_detect_outliers_basic", test_detect_outliers_basic)
    all_tests(6) = test_case("test_all_orthologs", test_all_orthologs)
    all_tests(7) = test_case("test_no_orthologs", test_no_orthologs)
    all_tests(8) = test_case("test_invalid_indices", test_invalid_indices)
    all_tests(9) = test_case("test_all_outliers", test_all_outliers)
    all_tests(10) = test_case("test_no_outliers", test_no_outliers)
    all_tests(11) = test_case("test_single_gene_family", test_single_gene_family)
    all_tests(12) = test_case("test_all_zero_distances", test_all_zero_distances)
    all_tests(13) = test_case("test_all_genes_one_family", test_all_genes_one_family)
    all_tests(14) = test_case("test_percentile_boundaries", test_percentile_boundaries)
    all_tests(15) = test_case("test_negative_nan_distances", test_negative_nan_distances)
    all_tests(16) = test_case("test_identical_rdi_at_threshold", test_identical_rdi_at_threshold)
  end subroutine get_all_tests

  !> @brief Run specific get_outliers tests by name.
  subroutine run_named_tests_get_outliers(test_names)
    character(len=*), intent(in) :: test_names(:)
    type(test_case) :: all_tests(16)
    integer :: i, j
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

  !> @brief Run all tox_get_outliers tests.
  subroutine run_all_tests_get_outliers()
    type(test_case) :: all_tests(16)
    integer :: i
    call get_all_tests(all_tests)
    do i = 1, size(all_tests)
      call all_tests(i)%test_proc()
      print *, trim(all_tests(i)%name), " passed."
    end do
    print *, "All tox_get_outliers tests passed successfully."
  end subroutine run_all_tests_get_outliers

  !> Test: Basic scaling with orthologs and non-orthologs.
  subroutine test_scaling_basic()
    integer, parameter :: n_genes=4, n_families=2
    real(real64) :: distances(n_genes) = [1.0, 2.0, 3.0, 4.0]
    integer :: gene_to_fam(n_genes) = [1,1,2,2]
    logical :: is_ortholog(n_genes) = [.true., .false., .true., .false.]
    real(real64) :: dscale(n_families)
    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, is_ortholog, dscale)
    call assert_equal_real(dscale(1), 1.0_real64, 1e-12_real64, "Family 1 scaling")
    call assert_equal_real(dscale(2), 3.0_real64, 1e-12_real64, "Family 2 scaling")
  end subroutine test_scaling_basic

  !> Test: No orthologs, should use fallback scaling 1.0.
  subroutine test_scaling_no_orthologs()
    integer, parameter :: n_genes=3, n_families=2
    real(real64) :: distances(n_genes) = [1.0, 2.0, 3.0]
    integer :: gene_to_fam(n_genes) = [1,2,2]
    logical :: is_ortholog(n_genes) = [.false., .false., .false.]
    real(real64) :: dscale(n_families)
    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, is_ortholog, dscale)
    call assert_equal_real(dscale(1), 1.0_real64, 1e-12_real64, "No orthologs family 1 fallback")
    call assert_equal_real(dscale(2), 1.0_real64, 1e-12_real64, "No orthologs family 2 fallback")
  end subroutine test_scaling_no_orthologs

  !> Test: Basic RDI calculation.
  subroutine test_rdi_basic()
    integer, parameter :: n_genes=3, n_families=2
    real(real64) :: distances(n_genes) = [2.0, 4.0, 6.0]
    integer :: gene_to_fam(n_genes) = [1,2,2]
    real(real64) :: dscale(n_families) = [2.0, 4.0]
    real(real64) :: rdi(n_genes)
    call compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi)
    call assert_equal_real(rdi(1), 1.0_real64, 1e-12_real64, "RDI gene 1")
    call assert_equal_real(rdi(2), 1.0_real64, 1e-12_real64, "RDI gene 2")
    call assert_equal_real(rdi(3), 1.5_real64, 1e-12_real64, "RDI gene 3")
  end subroutine test_rdi_basic

  !> Test: Outlier detection in a simple RDI vector.
  subroutine test_identify_outliers_basic()
    integer, parameter :: n_genes=4
    real(real64) :: rdi(n_genes)
    real(real64) :: work_array(n_genes)
    integer :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    logical :: is_outlier(n_genes)
    real(real64) :: threshold
    integer :: i
    rdi = [0.1_real64, 0.2_real64, 0.3_real64, 0.4_real64]
    perm = [(i, i=1,n_genes)]
    call identify_outliers(n_genes, rdi, 75.0_real64, work_array, perm, stack_left, stack_right, &
         is_outlier, threshold)
    call assert_true(is_outlier(4), "Highest RDI is outlier")
    call assert_true(is_outlier(3), "Second highest RDI is outlier")
    call assert_false(is_outlier(1), "Lowest RDI is not outlier")
  end subroutine test_identify_outliers_basic

  !> Test: Full outlier detection workflow.
  subroutine test_detect_outliers_basic()
    integer, parameter :: n_genes=5, n_families=2
    real(real64) :: distances(n_genes)
    integer :: gene_to_fam(n_genes)
    logical :: is_ortholog(n_genes)
    real(real64) :: work_array(n_genes)
    integer :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    logical :: is_outlier(n_genes)
    real(real64) :: rdi(n_genes), rdi_threshold
    distances = [1.0, 2.0, 3.0, 4.0, 5.0]
    gene_to_fam = [1,1,2,2,2]
    is_ortholog = [.true., .false., .true., .false., .true.]
    call detect_outliers(n_genes, n_families, distances, gene_to_fam, is_ortholog, 60.0_real64, &
         work_array, perm, stack_left, stack_right, is_outlier, rdi, rdi_threshold)
    call assert_true(any(is_outlier), "At least one outlier detected")
    call assert_true(rdi_threshold > 0.0, "Threshold positive")
  end subroutine test_detect_outliers_basic

  !> Test: All genes are orthologs.
  subroutine test_all_orthologs()
    integer, parameter :: n_genes=3, n_families=1
    real(real64) :: distances(n_genes) = [1.0, 2.0, 3.0]
    integer :: gene_to_fam(n_genes) = [1,1,1]
    logical :: is_ortholog(n_genes) = [.true., .true., .true.]
    real(real64) :: dscale(n_families)
    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, is_ortholog, dscale)
    call assert_equal_real(dscale(1), 3.0_real64, 1e-12_real64, "All orthologs scaling")
  end subroutine test_all_orthologs

  !> Test: No genes are orthologs.
  subroutine test_no_orthologs()
    integer, parameter :: n_genes=2, n_families=1
    real(real64) :: distances(n_genes) = [1.0, 2.0]
    integer :: gene_to_fam(n_genes) = [1,1]
    logical :: is_ortholog(n_genes) = [.false., .false.]
    real(real64) :: dscale(n_families)
    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, is_ortholog, dscale)
    call assert_equal_real(dscale(1), 1.0_real64, 1e-12_real64, "No orthologs fallback scaling")
  end subroutine test_no_orthologs

  !> Test: Invalid family indices.
  subroutine test_invalid_indices()
    integer, parameter :: n_genes=3, n_families=2
    real(real64) :: distances(n_genes) = [1.0, 2.0, 3.0]
    integer :: gene_to_fam(n_genes) = [1,3,0] ! 3 and 0 invalid
    logical :: is_ortholog(n_genes) = [.true., .true., .true.]
    real(real64) :: dscale(n_families)
    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, is_ortholog, dscale)
    call assert_equal_real(dscale(1), 1.0_real64, 1e-12_real64, "Valid index scaling")
    call assert_equal_real(dscale(2), 1.0_real64, 1e-12_real64, "Invalid index fallback")
  end subroutine test_invalid_indices

  !> Test: All genes are outliers (percentile 0).
  subroutine test_all_outliers()
    integer, parameter :: n_genes=4
    real(real64) :: rdi(n_genes) = [10.0_real64, 11.0_real64, 12.0_real64, 13.0_real64]
    real(real64) :: work_array(n_genes)
    integer :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    logical :: is_outlier(n_genes)
    real(real64) :: threshold
    integer :: i
    perm = [(i, i=1,n_genes)]
    call identify_outliers(n_genes, rdi, 0.0_real64, work_array, perm, stack_left, stack_right, is_outlier, threshold)
    call assert_true(all(is_outlier), "All are outliers at 0 percentile")
  end subroutine test_all_outliers

  !> Test: No genes are outliers (percentile 100).
  subroutine test_no_outliers()
    integer, parameter :: n_genes=4
    real(real64) :: rdi(n_genes) = [0.1_real64, 0.2_real64, 0.3_real64, 0.4_real64]
    real(real64) :: work_array(n_genes)
    integer :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    logical :: is_outlier(n_genes)
    real(real64) :: threshold
    integer :: i
    perm = [(i, i=1,n_genes)]
    call identify_outliers(n_genes, rdi, 100.0_real64, work_array, perm, stack_left, stack_right, is_outlier, threshold)
    call assert_true(.not. any(is_outlier(1:3)), "No outliers below 100 percentile")
  end subroutine test_no_outliers

  !> Test: Single gene and family case.
  subroutine test_single_gene_family()
    integer, parameter :: n_genes=1, n_families=1
    real(real64) :: distances(1) = [0.0_real64]
    integer :: gene_to_fam(1) = [1]
    logical :: is_ortholog(1) = [.true.]
    real(real64) :: dscale(1), rdi(1), work_array(1)
    integer :: perm(1), stack_left(1), stack_right(1)
    logical :: is_outlier(1)
    real(real64) :: rdi_threshold
    integer :: i
    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, is_ortholog, dscale)
    call compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi)
    perm = [(i, i=1,n_genes)]
    call detect_outliers(n_genes, n_families, distances, gene_to_fam, is_ortholog, 95.0_real64, &
                        work_array, perm, stack_left, stack_right, is_outlier, rdi, rdi_threshold)
    call assert_equal_real(dscale(1), 0.0_real64, 1e-12_real64, "Single gene scaling")
    call assert_equal_real(rdi(1), 0.0_real64, 1e-12_real64, "Single gene RDI")
    call assert_false(is_outlier(1), "Single gene not outlier")
  end subroutine test_single_gene_family

  !> Test: All distances are zero.
  subroutine test_all_zero_distances()
    integer, parameter :: n_genes=3, n_families=1
    real(real64) :: distances(n_genes) = [0.0_real64, 0.0_real64, 0.0_real64]
    integer :: gene_to_fam(n_genes) = [1,1,1]
    logical :: is_ortholog(n_genes) = [.true., .true., .true.]
    real(real64) :: dscale(1), rdi(n_genes)
    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, is_ortholog, dscale)
    call compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi)
    call assert_equal_real(dscale(1), 0.0_real64, 1e-12_real64, "All zero distances scaling")
    call assert_true(all(rdi == 0.0_real64), "All RDI zero for zero distances")
  end subroutine test_all_zero_distances

  !> Test: All genes in a single family.
  subroutine test_all_genes_one_family()
    integer, parameter :: n_genes=4, n_families=1
    real(real64) :: distances(n_genes) = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
    integer :: gene_to_fam(n_genes) = [1,1,1,1]
    logical :: is_ortholog(n_genes) = [.true., .false., .true., .false.]
    real(real64) :: dscale(1)
    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, is_ortholog, dscale)
    call assert_equal_real(dscale(1), 3.0_real64, 1e-12_real64, "All genes one family scaling")
  end subroutine test_all_genes_one_family

  !> Test: Percentile at boundaries (0 and 100).
  subroutine test_percentile_boundaries()
    integer, parameter :: n_genes=4
    real(real64) :: rdi(n_genes) = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
    real(real64) :: work_array(n_genes)
    integer :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    logical :: is_outlier(n_genes)
    real(real64) :: threshold
    integer :: i
    perm = [(i, i=1,n_genes)]
    call identify_outliers(n_genes, rdi, 0.0_real64, work_array, perm, stack_left, stack_right, is_outlier, threshold)
    call assert_true(all(is_outlier), "All outliers at 0 percentile")
    call identify_outliers(n_genes, rdi, 100.0_real64, work_array, perm, stack_left, stack_right, is_outlier, threshold)
    call assert_true(.not. any(is_outlier(1:3)), "No outliers below 100 percentile")
  end subroutine test_percentile_boundaries

  !> Test: Negative distances only.
  subroutine test_negative_nan_distances()
    integer, parameter :: n_genes=3, n_families=1
    real(real64) :: distances(n_genes)
    integer :: gene_to_fam(n_genes) = [1,1,1]
    logical :: is_ortholog(n_genes) = [.true., .true., .true.]
    real(real64) :: dscale(1), rdi(n_genes)
    ! Case: all orthologs with zero distance, scaling should be zero
    distances = [0.0_real64, 0.0_real64, 0.0_real64]
    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, is_ortholog, dscale)
    ! Now force a negative distance for the first gene
    distances(1) = -1.0_real64
    call compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi)
    call assert_true(rdi(1) == 0.0_real64, "Negative distance with zero scaling gives RDI=0.0 (not outlier)")
    ! NaN test skipped for portability
  end subroutine test_negative_nan_distances

  !> Test: Multiple genes with identical RDI at threshold.
  subroutine test_identical_rdi_at_threshold()
    integer, parameter :: n_genes=4
    real(real64) :: rdi(n_genes) = [1.0_real64, 2.0_real64, 2.0_real64, 3.0_real64]
    real(real64) :: work_array(n_genes)
    integer :: perm(n_genes), stack_left(n_genes), stack_right(n_genes)
    logical :: is_outlier(n_genes)
    real(real64) :: threshold
    integer :: i
    perm = [(i, i=1,n_genes)]
    call identify_outliers(n_genes, rdi, 50.0_real64, work_array, perm, stack_left, stack_right, is_outlier, threshold)
    call assert_true(is_outlier(3), "Identical RDI at threshold is outlier")
    call assert_true(is_outlier(2), "Identical RDI at threshold is outlier")
  end subroutine test_identical_rdi_at_threshold

end module mod_test_get_outliers
