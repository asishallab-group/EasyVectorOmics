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
        REAL(REAL64), ALLOCATABLE :: v(:,:), c(:); INTEGER(INT32) :: idx(2); REAL(REAL64) :: exp(3)
        ALLOCATE(v(3, 4), c(3)); v(:,1)=[1,2,3]; v(:,2)=[2,4,6]; v(:,3)=[3,6,9]; v(:,4)=[4,8,12]
        idx=[2,4]; exp=[3,6,9]
        CALL mean_vector(v, idx, 2, c)
        CALL assert_equal_real(c(1), exp(1), 1e-9_REAL64, "mean_vector: dim 1"); DEALLOCATE(v,c)
    END SUBROUTINE test_mean_vector_basic
    
    !> @brief Tests the group_centroid function in 'all' mode.
    SUBROUTINE test_group_centroid_all_mode()
        REAL(REAL64) :: vectors(2,5), centroids(2,2)
        INTEGER(INT32) :: gene_to_family(5), selected(5)
        INTEGER(INT32) :: orthologs_int(5)
        CHARACTER(LEN=3) :: mode_str
        INTEGER(INT32) :: mode_ascii(3)
        INTEGER :: i
        
        ! Arrange
        vectors(:,1)=[1,1]; vectors(:,2)=[3,3]; vectors(:,3)=[10,10]; vectors(:,4)=[20,20]; vectors(:,5)=[5,5]
        gene_to_family = [1, 1, 2, 2, 1]
        orthologs_int = 0 ! Not used in "all" mode
        mode_str = "all"
        DO i = 1, LEN(mode_str); mode_ascii(i) = ICHAR(mode_str(i:i)); END DO
        
        ! Act
        CALL group_centroid(vectors, 5, gene_to_family, 2, centroids, mode_ascii, LEN(mode_str), orthologs_int, selected)
        
        ! Assert
        CALL assert_equal_real(centroids(1,1), 3.0_REAL64, 1e-9_REAL64, "all_mode: fam 1, dim 1")
        CALL assert_equal_real(centroids(1,2), 3.0_REAL64, 1e-9_REAL64, "all_mode: fam 1, dim 2")
        CALL assert_equal_real(centroids(2,1), 15.0_REAL64, 1e-9_REAL64, "all_mode: fam 2, dim 1")
    END SUBROUTINE test_group_centroid_all_mode

    !> @brief Tests the group_centroid function in 'orthologs' mode.
    SUBROUTINE test_group_centroid_orthologs_mode()
        REAL(REAL64) :: vectors(2,5), centroids(2,2)
        INTEGER(INT32) :: gene_to_family(5), selected(5)
        LOGICAL :: orthologs_logical(5)
        INTEGER(INT32) :: orthologs_int(5)
        CHARACTER(LEN=9) :: mode_str
        INTEGER(INT32) :: mode_ascii(9)
        INTEGER :: i
        
        ! Arrange
        vectors(:,1)=[1,1]; vectors(:,2)=[3,3]; vectors(:,3)=[10,10]; vectors(:,4)=[20,20]; vectors(:,5)=[5,5]
        gene_to_family = [1, 1, 2, 2, 1]
        orthologs_logical = [.TRUE., .FALSE., .TRUE., .TRUE., .TRUE.]
        DO i = 1, 5
            IF (orthologs_logical(i)) THEN
                orthologs_int(i) = 1
            ELSE
                orthologs_int(i) = 0
            END IF
        END DO
        mode_str = "orthologs"
        DO i = 1, LEN(mode_str); mode_ascii(i) = ICHAR(mode_str(i:i)); END DO
        
        ! Act
        CALL group_centroid(vectors, 5, gene_to_family, 2, centroids, mode_ascii, LEN(mode_str), orthologs_int, selected)
        
        ! Assert
        CALL assert_equal_real(centroids(1,1), 3.0_REAL64, 1e-9_REAL64, "ortho_mode: fam 1, dim 1")
        CALL assert_equal_real(centroids(2,1), 15.0_REAL64, 1e-9_REAL64, "ortho_mode: fam 2, dim 1")
    END SUBROUTINE test_group_centroid_orthologs_mode

END MODULE mod_test_gene_centroids
