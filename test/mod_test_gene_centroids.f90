! Exhaustive unit test suite for gene_centroid_module, compatible with the project's test framework.
MODULE mod_test_gene_centroids
    USE, INTRINSIC :: iso_fortran_env, ONLY: INT32, REAL64
    USE gene_centroid_module, ONLY: group_centroid
    USE asserts, ONLY: assert_allclose_array_real

    IMPLICIT NONE
    PRIVATE
    PUBLIC :: run_all_tests_gene_centroids, run_named_tests_gene_centroids

    ABSTRACT INTERFACE
        SUBROUTINE test_interface()
        END SUBROUTINE test_interface
    END INTERFACE

    TYPE :: test_case
        CHARACTER(LEN=64) :: name
        PROCEDURE(test_interface), POINTER, NOPASS :: test_proc => NULL()
    END TYPE test_case

CONTAINS

    ! --------------------------------------------------------------------------
    ! Test Registration and Runners
    ! --------------------------------------------------------------------------

    ! Get array of all available tests in this suite.
    FUNCTION get_all_tests() RESULT(all_tests)
        TYPE(test_case) :: all_tests(8)
        all_tests(1) = test_case("test_basic_all_mode", test_basic_all_mode)
        all_tests(2) = test_case("test_basic_ortho_mode", test_basic_ortho_mode)
        all_tests(3) = test_case("test_empty_family", test_empty_family)
        all_tests(4) = test_case("test_no_matching_orthologs", test_no_matching_orthologs)
        all_tests(5) = test_case("test_single_gene_family", test_single_gene_family)
        all_tests(6) = test_case("test_extreme_values", test_extreme_values)
        all_tests(7) = test_case("test_higher_dimensions", test_higher_dimensions)
        all_tests(8) = test_case("test_gene_order_invariance", test_gene_order_invariance)
    END FUNCTION get_all_tests

    ! Run all tests in this module.
    SUBROUTINE run_all_tests_gene_centroids()
        TYPE(test_case) :: all_tests(8)
        INTEGER :: i
        all_tests = get_all_tests()
        WRITE(*, '(A)') "--- Running Suite: gene_centroids ---"
        DO i = 1, SIZE(all_tests)
            WRITE(*, '(A, A, A)', ADVANCE='NO') "  Running test: ", TRIM(all_tests(i)%name), "..."
            CALL all_tests(i)%test_proc()
            WRITE(*, '(A)') " PASSED"
        END DO
        WRITE(*, '(A)') "--- Suite PASSED: gene_centroids ---"
        WRITE(*,*)
    END SUBROUTINE run_all_tests_gene_centroids

    ! Run specific tests by name.
    SUBROUTINE run_named_tests_gene_centroids(test_names)
        CHARACTER(LEN=*), INTENT(IN) :: test_names(:)
        TYPE(test_case) :: all_tests(8)
        INTEGER :: i, j
        LOGICAL :: found
        all_tests = get_all_tests()
        DO i = 1, SIZE(test_names)
            found = .FALSE.
            DO j = 1, SIZE(all_tests)
                IF (TRIM(test_names(i)) == TRIM(all_tests(j)%name)) THEN
                    WRITE(*, '(A, A, A)', ADVANCE='NO') "  Running test: ", TRIM(all_tests(j)%name), "..."
                    CALL all_tests(j)%test_proc()
                    WRITE(*, '(A)') " PASSED"
                    found = .TRUE.
                    EXIT
                END IF
            END DO
            IF (.NOT. found) THEN
                WRITE(*,*) "Unknown test in suite 'gene_centroids': ", TRIM(test_names(i))
            END IF
        END DO
    END SUBROUTINE run_named_tests_gene_centroids

    ! --------------------------------------------------------------------------
    ! Individual Test Cases
    ! --------------------------------------------------------------------------

    ! Test case 1: Basic functionality in 'all' mode.
    SUBROUTINE test_basic_all_mode()
        INTEGER, PARAMETER :: d=2, n_genes=5, n_families=2
        REAL(REAL64) :: vectors(d, n_genes), centroids(d, n_families)
        INTEGER(INT32) :: gene_to_family(n_genes)
        LOGICAL :: ortholog_set(n_genes)
        REAL(REAL64) :: expected(d, n_families)

        vectors = reshape([1.0, 1.0, 3.0, 3.0, 10.0, 10.0, 20.0, 20.0, 5.0, 5.0], [d, n_genes])
        gene_to_family = [1, 1, 2, 2, 1]
        ortholog_set = .false.
        expected = reshape([3.0, 3.0, 15.0, 15.0], [d, n_families])
        
        CALL group_centroid(vectors, d, n_genes, gene_to_family, n_families, centroids, .true., ortholog_set)
        CALL assert_allclose_array_real(centroids, expected, d*n_families, 0.0_REAL64, 1e-9_REAL64, "test_basic_all_mode")
    END SUBROUTINE test_basic_all_mode

    ! Test case 2: Basic functionality in 'ortho' mode.
    SUBROUTINE test_basic_ortho_mode()
        INTEGER, PARAMETER :: d=2, n_genes=5, n_families=2
        REAL(REAL64) :: vectors(d, n_genes), centroids(d, n_families)
        INTEGER(INT32) :: gene_to_family(n_genes)
        LOGICAL :: ortholog_set(n_genes)
        REAL(REAL64) :: expected(d, n_families)

        vectors = reshape([1.0, 1.0, 3.0, 3.0, 10.0, 10.0, 20.0, 20.0, 5.0, 5.0], [d, n_genes])
        gene_to_family = [1, 1, 2, 2, 1]
        ortholog_set = [.true., .false., .true., .true., .true.]
        expected = reshape([3.0, 3.0, 15.0, 15.0], [d, n_families])

        CALL group_centroid(vectors, d, n_genes, gene_to_family, n_families, centroids, .false., ortholog_set)
        CALL assert_allclose_array_real(centroids, expected, d*n_families, 0.0_REAL64, 1e-9_REAL64, "test_basic_ortho_mode")
    END SUBROUTINE test_basic_ortho_mode

    ! Test case 3: A family exists but has no genes assigned to it.
    SUBROUTINE test_empty_family()
        INTEGER, PARAMETER :: d=3, n_genes=4, n_families=2
        REAL(REAL64) :: vectors(d, n_genes), centroids(d, n_families)
        INTEGER(INT32) :: gene_to_family(n_genes)
        LOGICAL :: ortholog_set(n_genes)
        REAL(REAL64) :: expected(d, n_families)

        vectors = 1.0
        gene_to_family = [1, 1, 1, 1] ! All genes in family 1, none in family 2
        ortholog_set = .true.
        expected = 0.0
        expected(:, 1) = 1.0 ! Centroid of family 1 is just the vector itself
        
        CALL group_centroid(vectors, d, n_genes, gene_to_family, n_families, centroids, .true., ortholog_set)
        CALL assert_allclose_array_real(centroids, expected, d*n_families, 0.0_REAL64, 1e-9_REAL64, "test_empty_family")
    END SUBROUTINE test_empty_family

    ! Test case 4: 'ortho' mode is selected, but a family has no orthologs.
    SUBROUTINE test_no_matching_orthologs()
        INTEGER, PARAMETER :: d=2, n_genes=3, n_families=1
        REAL(REAL64) :: vectors(d, n_genes), centroids(d, n_families)
        INTEGER(INT32) :: gene_to_family(n_genes)
        LOGICAL :: ortholog_set(n_genes)
        REAL(REAL64) :: expected(d, n_families)

        vectors = reshape([10.0, 10.0, 20.0, 20.0, 30.0, 30.0], [d, n_genes])
        gene_to_family = [1, 1, 1]
        ortholog_set = .false. ! No genes are orthologs
        expected = 0.0 ! Expect a zero vector
        
        CALL group_centroid(vectors, d, n_genes, gene_to_family, n_families, centroids, .false., ortholog_set)
        CALL assert_allclose_array_real(centroids, expected, d*n_families, 0.0_REAL64, 1e-9_REAL64, "test_no_matching_orthologs")
    END SUBROUTINE test_no_matching_orthologs

    ! Test case 5: A family contains only a single gene.
    SUBROUTINE test_single_gene_family()
        INTEGER, PARAMETER :: d=3, n_genes=1, n_families=1
        REAL(REAL64) :: vectors(d, n_genes), centroids(d, n_families)
        INTEGER(INT32) :: gene_to_family(n_genes)
        LOGICAL :: ortholog_set(n_genes)
        
        vectors(:,1) = [12.3, -4.5, 6.7]
        gene_to_family = [1]
        ortholog_set = .true.
        
        CALL group_centroid(vectors, d, n_genes, gene_to_family, n_families, centroids, .true., ortholog_set)
        CALL assert_allclose_array_real(centroids, vectors, d*n_families, 0.0_REAL64, 1e-9_REAL64, "test_single_gene_family")
    END SUBROUTINE test_single_gene_family

    ! Test case 6: Input vectors with extreme values.
    SUBROUTINE test_extreme_values()
        INTEGER, PARAMETER :: d=2, n_genes=4, n_families=1
        REAL(REAL64) :: vectors(d, n_genes), centroids(d, n_families)
        INTEGER(INT32) :: gene_to_family(n_genes)
        LOGICAL :: ortholog_set(n_genes)
        REAL(REAL64) :: expected(d, n_families)

        vectors(:,1) = [1.0e12, -1.0e-12]
        vectors(:,2) = [-1.0e12, 1.0e-12]
        vectors(:,3) = [0.0, 5.0]
        vectors(:,4) = [0.0, -5.0]
        gene_to_family = [1, 1, 1, 1]
        ortholog_set = .true.
        expected(:,1) = [0.0, 0.0] ! (1e12 - 1e12 + 0 + 0)/4 = 0, (-1e-12 + 1e-12 + 5 - 5)/4 = 0

        CALL group_centroid(vectors, d, n_genes, gene_to_family, n_families, centroids, .true., ortholog_set)
        CALL assert_allclose_array_real(centroids, expected, d*n_families, 0.0_REAL64, 1e-9_REAL64, "test_extreme_values")
    END SUBROUTINE test_extreme_values

    ! Test case 7: Higher dimensional data.
    SUBROUTINE test_higher_dimensions()
        INTEGER, PARAMETER :: d=10, n_genes=100, n_families=5
        REAL(REAL64) :: vectors(d, n_genes), centroids(d, n_families)
        INTEGER(INT32) :: gene_to_family(n_genes)
        LOGICAL :: ortholog_set(n_genes)
        INTEGER :: i
        REAL(REAL64) :: expected(d, n_families)

        DO i = 1, n_genes
            vectors(:, i) = real(i, REAL64)
            gene_to_family(i) = mod(i-1, n_families) + 1
        END DO
        ortholog_set = .true.
        
        ! Calculate expected manually for family 1 (genes 1, 6, 11, ...)
        expected = 0.0
        expected(:,1) = sum(vectors(:, 1:n_genes:n_families), dim=2) / real(count(gene_to_family == 1), REAL64)

        CALL group_centroid(vectors, d, n_genes, gene_to_family, n_families, centroids, .true., ortholog_set)
        CALL assert_allclose_array_real(centroids(:,1), expected(:,1), d, 0.0_REAL64, 1e-9_REAL64, "test_higher_dimensions")
    END SUBROUTINE test_higher_dimensions

    ! Test case 8: Ensure the result is invariant to the order of genes.
    SUBROUTINE test_gene_order_invariance()
        INTEGER, PARAMETER :: d=2, n_genes=5, n_families=2
        REAL(REAL64) :: vectors1(d, n_genes), centroids1(d, n_families)
        INTEGER(INT32) :: gene_to_family1(n_genes)
        LOGICAL :: ortholog_set1(n_genes)

        REAL(REAL64) :: vectors2(d, n_genes), centroids2(d, n_families)
        INTEGER(INT32) :: gene_to_family2(n_genes)
        LOGICAL :: ortholog_set2(n_genes)
        
        ! Setup 1: Original order
        vectors1 = reshape([1.0, 1.0, 3.0, 3.0, 10.0, 10.0, 20.0, 20.0, 5.0, 5.0], [d, n_genes])
        gene_to_family1 = [1, 1, 2, 2, 1]
        ortholog_set1 = [.true., .false., .true., .true., .true.]
        
        ! Setup 2: Shuffled order
        vectors2 = reshape([5.0, 5.0, 10.0, 10.0, 1.0, 1.0, 3.0, 3.0, 20.0, 20.0], [d, n_genes])
        gene_to_family2 = [1, 2, 1, 1, 2]
        ortholog_set2 = [.true., .true., .true., .false., .true.]

        CALL group_centroid(vectors1, d, n_genes, gene_to_family1, n_families, centroids1, .false., ortholog_set1)
        CALL group_centroid(vectors2, d, n_genes, gene_to_family2, n_families, centroids2, .false., ortholog_set2)

        CALL assert_allclose_array_real(centroids1, centroids2, d*n_families, 0.0_REAL64, 1e-9_REAL64, "test_gene_order_invariance")
    END SUBROUTINE test_gene_order_invariance

END MODULE mod_test_gene_centroids
