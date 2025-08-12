!> @brief A module for computing expression centroids of gene families.
!> @note Adheres to F42 conventions: stateless, allocation-free loops, modular.
MODULE gene_centroid_module
    USE, INTRINSIC :: iso_fortran_env, ONLY: INT32, REAL64
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: mean_vector, group_centroid

CONTAINS

    !> @brief Computes the element-wise mean for a given set of vectors.
    SUBROUTINE mean_vector(vectors, gene_indices, n_genes, centroid)
        REAL(REAL64), INTENT(IN) :: vectors(:,:)
        INTEGER(INT32), INTENT(IN) :: gene_indices(:)
        INTEGER(INT32), INTENT(IN) :: n_genes
        REAL(REAL64), INTENT(OUT) :: centroid(:)

        INTEGER :: d, i, j
        INTEGER(INT32) :: gene_idx
        d = SIZE(vectors, DIM=1)
        centroid = 0.0_REAL64
        IF (n_genes == 0) RETURN
        DO i = 1, n_genes
            gene_idx = gene_indices(i)
            !$OMP SIMD
            DO j = 1, d
                centroid(j) = centroid(j) + vectors(j, gene_idx)
            END DO
            !$OMP END SIMD
        END DO
        centroid = centroid / REAL(n_genes, REAL64)
    END SUBROUTINE mean_vector

    !> @brief Iterates over orthogroups, filters gene indices, and computes centroids.
    SUBROUTINE group_centroid(vectors, num_genes, gene_to_family_map, num_families, &
                              centroid_matrix, mode, ortholog_set, selected_indices)
        REAL(REAL64), INTENT(IN) :: vectors(:,:)
        INTEGER(INT32), INTENT(IN) :: num_genes
        INTEGER(INT32), INTENT(IN) :: gene_to_family_map(num_genes)
        INTEGER(INT32), INTENT(IN) :: num_families
        REAL(REAL64), INTENT(OUT) :: centroid_matrix(:,:)
        CHARACTER(LEN=*), INTENT(IN) :: mode
        LOGICAL, INTENT(IN) :: ortholog_set(num_genes)
        INTEGER(INT32), INTENT(INOUT) :: selected_indices(:) ! Workspace buffer

        INTEGER(INT32) :: i, j, current_count
        
        DO j = 1, num_families
            current_count = 0
            IF (TRIM(mode) == 'all') THEN
                DO i = 1, num_genes
                    IF (gene_to_family_map(i) == j) THEN
                        current_count = current_count + 1
                        selected_indices(current_count) = i
                    END IF
                END DO
            ELSE IF (TRIM(mode) == 'orthologs') THEN
                DO i = 1, num_genes
                    IF (gene_to_family_map(i) == j .AND. ortholog_set(i)) THEN
                        current_count = current_count + 1
                        selected_indices(current_count) = i
                    END IF
                END DO
            END IF
            
            CALL mean_vector(vectors, selected_indices, current_count, centroid_matrix(j, :))
        END DO
    END SUBROUTINE group_centroid

END MODULE gene_centroid_module

! =============================================================================
! C and R Wrapper Subroutines
! =============================================================================

!> @brief C interface for computing group centroids.
SUBROUTINE group_centroid_c(vectors, d, n, gene_to_family_map, num_families, &
                            centroid_matrix, mode_ascii, mode_len, ortholog_set, &
                            selected_indices, selected_indices_len) &
                            BIND(C, NAME='group_centroid_c')
    USE, INTRINSIC :: iso_c_binding
    USE gene_centroid_module
    IMPLICIT NONE

    REAL(C_DOUBLE), INTENT(IN) :: vectors(d, n)
    INTEGER(C_INT), VALUE, INTENT(IN) :: d, n, num_families, mode_len, selected_indices_len
    INTEGER(C_INT), INTENT(IN) :: gene_to_family_map(n)
    REAL(C_DOUBLE), INTENT(OUT) :: centroid_matrix(num_families, d)
    INTEGER(C_INT), INTENT(IN) :: mode_ascii(*)
    LOGICAL(C_BOOL), INTENT(IN) :: ortholog_set(n)
    ! CORRECTED: Use an explicit-shape declaration for the workspace array
    INTEGER(C_INT), INTENT(INOUT) :: selected_indices(selected_indices_len)

    CHARACTER(LEN=:), ALLOCATABLE :: mode_f
    LOGICAL, ALLOCATABLE :: ortholog_set_f(:)
    INTEGER :: i

    ALLOCATE(CHARACTER(LEN=mode_len) :: mode_f)
    DO i = 1, mode_len; mode_f(i:i) = CHAR(mode_ascii(i)); END DO

    ALLOCATE(ortholog_set_f(n)); ortholog_set_f = ortholog_set

    CALL group_centroid(vectors, n, gene_to_family_map, num_families, &
                        centroid_matrix, mode_f, ortholog_set_f, selected_indices)

END SUBROUTINE group_centroid_c

!> @brief R interface for computing group centroids.
SUBROUTINE group_centroid_r(vectors, d, n, gene_to_family_map, num_families, &
                            centroid_matrix, mode_ascii, mode_len, ortholog_set, &
                            selected_indices, selected_indices_len)
    USE, INTRINSIC :: iso_fortran_env, ONLY: INT32, REAL64
    USE gene_centroid_module
    IMPLICIT NONE

    REAL(REAL64), INTENT(IN) :: vectors(d, n)
    INTEGER(INT32), INTENT(IN) :: d, n, num_families, mode_len, selected_indices_len
    INTEGER(INT32), INTENT(IN) :: gene_to_family_map(n)
    REAL(REAL64), INTENT(OUT) :: centroid_matrix(num_families, d)
    INTEGER(INT32), INTENT(IN) :: mode_ascii(mode_len)
    LOGICAL, INTENT(IN) :: ortholog_set(n)
    INTEGER(INT32), INTENT(INOUT) :: selected_indices(selected_indices_len)

    CHARACTER(LEN=:), ALLOCATABLE :: mode_f
    INTEGER :: i

    ALLOCATE(CHARACTER(LEN=mode_len) :: mode_f)
    DO i = 1, mode_len; mode_f(i:i) = CHAR(mode_ascii(i)); END DO

    CALL group_centroid(vectors, n, gene_to_family_map, num_families, &
                        centroid_matrix, mode_f, ortholog_set, selected_indices)

END SUBROUTINE group_centroid_r
