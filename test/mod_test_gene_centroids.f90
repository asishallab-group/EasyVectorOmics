! Unit test suite for gene_centroid_module.
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

    ! Tests the basic functionality of the mean_vector subroutine.
    SUBROUTINE test_mean_vector_basic()
        REAL(REAL64), ALLOCATABLE :: v(:,:), c(:)
        LOGICAL :: selection_mask(4)
        REAL(REAL64) :: exp(3)
        
        ALLOCATE(v(3, 4), c(3))
        v(:,1)=[1.0, 2.0, 3.0]
        v(:,2)=[2.0, 4.0, 6.0]
        v(:,3)=[3.0, 6.0, 9.0]
        v(:,4)=[4.0, 8.0, 12.0]
        
        ! Select vectors at index 2 and 4
        selection_mask = [ .false., .true., .false., .true. ]
        exp=[3.0, 6.0, 9.0]
        
        CALL mean_vector(v, selection_mask, c)
        
        CALL assert_equal_real(c(1), exp(1), 1e-9_REAL64, "mean_vector: dim 1")
        CALL assert_equal_real(c(2), exp(2), 1e-9_REAL64, "mean_vector: dim 2")
        CALL assert_equal_real(c(3), exp(3), 1e-9_REAL64, "mean_vector: dim 3")
        
        DEALLOCATE(v,c)
    END SUBROUTINE test_mean_vector_basic
    
    ! Tests the group_centroid function in 'all' mode.
    SUBROUTINE test_group_centroid_all_mode()
        REAL(REAL64) :: vectors(2,5), centroids(2,2)
        INTEGER(INT32) :: gene_to_family(5)
        LOGICAL :: ortholog_set(5)
        
        ! Arrange
        ! Family 1: genes 1, 2, 5 -> vectors [1,1], [3,3], [5,5]. Centroid should be [3,3].
        ! Family 2: genes 3, 4    -> vectors [10,10], [20,20]. Centroid should be [15,15].
        vectors(:,1)=[1,1]; vectors(:,2)=[3,3]; vectors(:,3)=[10,10]; vectors(:,4)=[20,20]; vectors(:,5)=[5,5]
        gene_to_family = [1, 1, 2, 2, 1]
        ortholog_set = .false. ! Not used in "all" mode, pass dummy value.
        
        ! Act
        CALL group_centroid(vectors, 5, gene_to_family, 2, centroids, .true., ortholog_set)
        
        ! Assert (centroids are columns, so centroids(:,1) is for family 1)
        CALL assert_equal_real(centroids(1,1), 3.0_REAL64, 1e-9_REAL64, "all_mode: fam 1, dim 1")
        CALL assert_equal_real(centroids(2,1), 3.0_REAL64, 1e-9_REAL64, "all_mode: fam 1, dim 2")
        CALL assert_equal_real(centroids(1,2), 15.0_REAL64, 1e-9_REAL64, "all_mode: fam 2, dim 1")
        CALL assert_equal_real(centroids(2,2), 15.0_REAL64, 1e-9_REAL64, "all_mode: fam 2, dim 2")
    END SUBROUTINE test_group_centroid_all_mode

    ! Tests the group_centroid function in 'orthologs' mode.
    SUBROUTINE test_group_centroid_orthologs_mode()
        REAL(REAL64) :: vectors(2,5), centroids(2,2)
        INTEGER(INT32) :: gene_to_family(5)
        LOGICAL :: ortholog_set(5)
        
        ! Arrange
        ! Family 1 orthologs: genes 1, 5 -> vectors [1,1], [5,5]. Centroid should be [3,3].
        ! Family 2 orthologs: genes 3, 4 -> vectors [10,10], [20,20]. Centroid should be [15,15].
        vectors(:,1)=[1,1]; vectors(:,2)=[3,3]; vectors(:,3)=[10,10]; vectors(:,4)=[20,20]; vectors(:,5)=[5,5]
        gene_to_family = [1, 1, 2, 2, 1]
        ortholog_set = [.TRUE., .FALSE., .TRUE., .TRUE., .TRUE.]
        
        ! Act
        CALL group_centroid(vectors, 5, gene_to_family, 2, centroids, .false., ortholog_set)
        
        ! Assert (centroids are columns, so centroids(:,1) is for family 1)
        CALL assert_equal_real(centroids(1,1), 3.0_REAL64, 1e-9_REAL64, "ortho_mode: fam 1, dim 1")
        CALL assert_equal_real(centroids(2,1), 3.0_REAL64, 1e-9_REAL64, "ortho_mode: fam 1, dim 2")
        CALL assert_equal_real(centroids(1,2), 15.0_REAL64, 1e-9_REAL64, "ortho_mode: fam 2, dim 1")
        CALL assert_equal_real(centroids(2,2), 15.0_REAL64, 1e-9_REAL64, "ortho_mode: fam 2, dim 2")
    END SUBROUTINE test_group_centroid_orthologs_mode

END MODULE mod_test_gene_centroids
