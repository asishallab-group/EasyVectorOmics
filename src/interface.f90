! ==============================================================================
! MODULE: centroid_interfaces
! ==============================================================================
! Description:
!   Provides C and R compatible wrappers for the `gene_centroid_calculator`
!   module.
!
! Revision Notes:
!   - Corrected C wrapper to use a temporary allocatable matrix to satisfy the
!     intent of the core `group_centroid` routine, resolving the compiler error.
!   - Rewrote the R wrapper to only use the public `group_centroid` function.
!   - Added the `TARGET` attribute to C-interfacing arrays to satisfy `C_LOC`.
!   - Corrected all REAL kinds to default REAL and REAL(C_FLOAT) to prevent
!     type mismatches with the core compiled module.
! ==============================================================================
MODULE centroid_interfaces
  USE, INTRINSIC :: iso_c_binding
  USE gene_centroid_calculator, ONLY: orthogroup, group_centroid
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: compute_centroid_for_group_r
  PUBLIC :: compute_centroids_c

CONTAINS

! ==============================================================================
! SUBROUTINE: compute_centroid_for_group_r
! ==============================================================================
! Purpose:
!   R-compatible wrapper to compute the centroid for a SINGLE orthogroup.
! ==============================================================================
  SUBROUTINE compute_centroid_for_group_r(vectors, n_genes, n_tissues, &
                                          gene_indices_for_group, n_genes_in_group, &
                                          centroid)
    ! --- Argument Declaration ---
    INTEGER, INTENT(IN) :: n_genes, n_tissues, n_genes_in_group
    REAL, INTENT(IN) :: vectors(n_genes, n_tissues)
    INTEGER, INTENT(IN) :: gene_indices_for_group(n_genes_in_group)
    REAL, INTENT(OUT) :: centroid(n_tissues)

    ! --- Local Variables for temporary group construction ---
    TYPE(orthogroup), DIMENSION(1) :: temp_group_list
    REAL, DIMENSION(:,:), ALLOCATABLE :: temp_centroid_matrix

    ! --- Implementation ---
    ! 1. Create a temporary, single-element orthogroup list.
    ALLOCATE(temp_group_list(1)%gene_indices(n_genes_in_group))
    temp_group_list(1)%id = 'R_INTERFACE_GROUP'
    temp_group_list(1)%gene_indices = gene_indices_for_group(1:n_genes_in_group)

    ! 2. Call the public group_centroid function.
    CALL group_centroid(vectors, temp_group_list, temp_centroid_matrix, 'all')

    ! 3. Extract the single resulting centroid vector.
    IF (ALLOCATED(temp_centroid_matrix)) THEN
        centroid = temp_centroid_matrix(1, :)
        DEALLOCATE(temp_centroid_matrix)
    ELSE
        centroid = 0.0 ! Should not happen if logic is correct
    END IF

    ! 4. Clean up allocated memory.
    DEALLOCATE(temp_group_list(1)%gene_indices)

  END SUBROUTINE compute_centroid_for_group_r


! ==============================================================================
! SUBROUTINE: compute_centroids_c
! ==============================================================================
! Purpose:
!   C-compatible wrapper for `group_centroid`. Taking only C-compatible arguments
!   (pointers/assumed-size arrays and scalars).
! ==============================================================================
  SUBROUTINE compute_centroids_c(vectors_flat, n_genes, n_tissues, &
                                 group_indices, group_boundaries, n_groups, &
                                 mode, ortholog_set, centroid_matrix_flat) &
                                 BIND(C, name='compute_centroids_c')

    ! --- Argument Declaration (C types) ---
    INTEGER(C_INT), VALUE, INTENT(IN) :: n_genes, n_tissues, n_groups
    CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) :: mode
    REAL(C_FLOAT), DIMENSION(*), INTENT(IN), TARGET :: vectors_flat
    INTEGER(C_INT), DIMENSION(*), INTENT(IN) :: group_indices, group_boundaries
    TYPE(C_PTR), VALUE, INTENT(IN) :: ortholog_set
    REAL(C_FLOAT), DIMENSION(*), INTENT(OUT), TARGET :: centroid_matrix_flat

    ! --- Fortran Pointers for reconstruction ---
    REAL(C_FLOAT), POINTER :: vectors(:,:)
    REAL(C_FLOAT), POINTER :: centroid_matrix_p(:,:)
    LOGICAL, POINTER :: ortholog_set_p(:)

    ! --- Local Variables ---
    CHARACTER(LEN=10) :: f_mode
    TYPE(orthogroup), DIMENSION(n_groups) :: f_group_list
    INTEGER :: i, start_idx
    ! REVISION: Add a temporary allocatable matrix for the core routine call.
    REAL, ALLOCATABLE :: f_centroid_matrix_temp(:,:)

    ! --- Reconstruct Fortran arrays from C pointers (zero-copy) ---
    CALL C_F_POINTER(C_LOC(vectors_flat), vectors, (/n_genes, n_tissues/))
    CALL C_F_POINTER(C_LOC(centroid_matrix_flat), centroid_matrix_p, (/n_groups, n_tissues/))

    ! --- Convert C string to Fortran string ---
    i = 1
    f_mode = ''
    DO WHILE (mode(i) /= C_NULL_CHAR .AND. i <= LEN(f_mode))
      f_mode(i:i) = mode(i)
      i = i + 1
    END DO
    f_mode = TRIM(f_mode)

    ! --- Reconstruct the orthogroup list from index arrays ---
    start_idx = 1
    DO i = 1, n_groups
        f_group_list(i)%id = 'C_INTERFACE_GROUP'
        ALLOCATE(f_group_list(i)%gene_indices(group_boundaries(i) - start_idx + 1))
        f_group_list(i)%gene_indices = group_indices(start_idx : group_boundaries(i))
        start_idx = group_boundaries(i) + 1
    END DO

    ! --- Handle the optional ortholog_set ---
    IF (C_ASSOCIATED(ortholog_set)) THEN
        CALL C_F_POINTER(ortholog_set, ortholog_set_p, (/n_genes/))
        CALL group_centroid(vectors, f_group_list, f_centroid_matrix_temp, f_mode, ortholog_set_p)
    ELSE
        CALL group_centroid(vectors, f_group_list, f_centroid_matrix_temp, f_mode)
    END IF
    
    ! --- Copy results from temp array back to the C buffer and clean up ---
    IF (ALLOCATED(f_centroid_matrix_temp)) THEN
        centroid_matrix_p = f_centroid_matrix_temp
        DEALLOCATE(f_centroid_matrix_temp)
    END IF

    ! --- Clean up allocated group memory ---
    DO i = 1, n_groups
        DEALLOCATE(f_group_list(i)%gene_indices)
    END DO

  END SUBROUTINE compute_centroids_c

END MODULE centroid_interfaces
