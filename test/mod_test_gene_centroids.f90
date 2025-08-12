!> @brief Unit test suite for gene_centroid_module.
MODULE mod_test_gene_centroids
    USE, INTRINSIC :: iso_fortran_env, ONLY: INT32, REAL64
    USE gene_centroid_module, ONLY: mean_vector, group_centroid
    USE asserts, ONLY: assert_equal_real

    IMPLICIT NONE
    PUBLIC

    ABSTRACT INTERFACE
        SUBROUTINE test_interface()
        END SUBROUTINE test_interface
    END INTERFACE

    TYPE :: test_case
        CHARACTER(LEN=64) :: name
        PROCEDURE(test_interface), POINTER, NOPASS :: test_proc => NULL()
    END TYPE test_case

CONTAINS

    FUNCTION get_all_tests_gene_centroids() RESULT(all_tests)
        TYPE(test_case) :: all_tests(3)
        all_tests(1) = test_case("test_mean_vector_basic", test_mean_vector_basic)
        all_tests(2) = test_case("test_group_centroid_all_mode", test_group_centroid_all_mode)
        all_tests(3) = test_case("test_group_centroid_orthologs_mode", test_group_centroid_orthologs_mode)
    END FUNCTION get_all_tests_gene_centroids

    SUBROUTINE run_all_tests_gene_centroids()
        TYPE(test_case) :: all_tests(3)
        INTEGER :: i
        all_tests = get_all_tests_gene_centroids()
        WRITE(*, '(A)') "--- Running Suite: gene_centroids ---"
        DO i = 1, SIZE(all_tests)
            WRITE(*, '(A, A, A)', ADVANCE='NO') "  Running test: ", TRIM(all_tests(i)%name), "..."
            CALL all_tests(i)%test_proc()
            WRITE(*, '(A)') " PASSED"
        END DO
        WRITE(*, '(A)') "--- Suite PASSED: gene_centroids ---"
        WRITE(*,*)
    END SUBROUTINE run_all_tests_gene_centroids

    SUBROUTINE run_named_tests_gene_centroids(test_names)
        CHARACTER(LEN=*), INTENT(IN) :: test_names(:)
        ! This is a stub for now, can be implemented fully if needed.
    END SUBROUTINE run_named_tests_gene_centroids

    !> @brief Tests the basic functionality of the mean_vector subroutine.
    SUBROUTINE test_mean_vector_basic()
        REAL(REAL64), ALLOCATABLE :: vectors(:,:)
        REAL(REAL64), ALLOCATABLE :: centroid(:)
        INTEGER(INT32) :: gene_indices(2)
        REAL(REAL64) :: expected_centroid(3)
        INTEGER :: d=3, n=4

        ! Arrange: Create a (d x n) = (3 x 4) matrix of vectors
        ALLOCATE(vectors(d, n), centroid(d))
        vectors(:, 1) = [1.0, 2.0, 3.0]
        vectors(:, 2) = [2.0, 4.0, 6.0]
        vectors(:, 3) = [3.0, 6.0, 9.0]
        vectors(:, 4) = [4.0, 8.0, 12.0]
        
        ! We want to average the 2nd and 4th gene vectors.
        gene_indices = [2, 4]
        
        ! Expected: mean([2,4,6] and [4,8,12]) = [3, 6, 9]
        expected_centroid = [3.0, 6.0, 9.0]
        
        ! Act
        CALL mean_vector(vectors, gene_indices, 2, centroid)
        
        ! Assert
        CALL assert_equal_real(centroid(1), expected_centroid(1), 1e-9_REAL64, "mean_vector: dim 1 failed")
        CALL assert_equal_real(centroid(2), expected_centroid(2), 1e-9_REAL64, "mean_vector: dim 2 failed")
        CALL assert_equal_real(centroid(3), expected_centroid(3), 1e-9_REAL64, "mean_vector: dim 3 failed")
        
        DEALLOCATE(vectors, centroid)
    END SUBROUTINE test_mean_vector_basic
    
    !> @brief Tests the group_centroid function in 'all' mode.
    SUBROUTINE test_group_centroid_all_mode()
        REAL(REAL64) :: vectors(2,5), centroids(2,2)
        INTEGER(INT32) :: gene_to_family(5), selected(5)
        LOGICAL :: orthologs(5)
        
        ! Arrange
        vectors(:,1)=[1,1]; vectors(:,2)=[3,3]; vectors(:,3)=[10,10]; vectors(:,4)=[20,20]; vectors(:,5)=[5,5]
        gene_to_family = [1, 1, 2, 2, 1] ! Genes 1,2,5 in Fam 1; Genes 3,4 in Fam 2
        orthologs = .FALSE. ! Not used in "all" mode, can be anything
        
        ! Act
        CALL group_centroid(vectors, 5, gene_to_family, 2, centroids, "all", orthologs, selected)
        
        ! Assert
        ! Fam 1: mean([1,1], [3,3], [5,5]) = [3,3]
        CALL assert_equal_real(centroids(1,1), 3.0_REAL64, 1e-9_REAL64, "all_mode: fam 1, dim 1")
        CALL assert_equal_real(centroids(1,2), 3.0_REAL64, 1e-9_REAL64, "all_mode: fam 1, dim 2")
        ! Fam 2: mean([10,10], [20,20]) = [15,15]
        CALL assert_equal_real(centroids(2,1), 15.0_REAL64, 1e-9_REAL64, "all_mode: fam 2, dim 1")
        CALL assert_equal_real(centroids(2,2), 15.0_REAL64, 1e-9_REAL64, "all_mode: fam 2, dim 2")
    END SUBROUTINE test_group_centroid_all_mode

    !> @brief Tests the group_centroid function in 'orthologs' mode.
    SUBROUTINE test_group_centroid_orthologs_mode()
        REAL(REAL64) :: vectors(2,5), centroids(2,2)
        INTEGER(INT32) :: gene_to_family(5), selected(5)
        LOGICAL :: orthologs(5)
        
        ! Arrange
        vectors(:,1)=[1,1]; vectors(:,2)=[3,3]; vectors(:,3)=[10,10]; vectors(:,4)=[20,20]; vectors(:,5)=[5,5]
        gene_to_family = [1, 1, 2, 2, 1]
        orthologs = [.TRUE., .FALSE., .TRUE., .TRUE., .TRUE.] ! Orthologs are genes 1, 3, 4, 5
        
        ! Act
        CALL group_centroid(vectors, 5, gene_to_family, 2, centroids, "orthologs", orthologs, selected)
        
        ! Assert
        ! Fam 1 orthologs: genes 1, 5 -> mean([1,1], [5,5]) = [3,3]
        CALL assert_equal_real(centroids(1,1), 3.0_REAL64, 1e-9_REAL64, "ortho_mode: fam 1, dim 1")
        CALL assert_equal_real(centroids(1,2), 3.0_REAL64, 1e-9_REAL64, "ortho_mode: fam 1, dim 2")
        ! Fam 2 orthologs: genes 3, 4 -> mean([10,10], [20,20]) = [15,15]
        CALL assert_equal_real(centroids(2,1), 15.0_REAL64, 1e-9_REAL64, "ortho_mode: fam 2, dim 1")
        CALL assert_equal_real(centroids(2,2), 15.0_REAL64, 1e-9_REAL64, "ortho_mode: fam 2, dim 2")
    END SUBROUTINE test_group_centroid_orthologs_mode

END MODULE mod_test_gene_centroids
