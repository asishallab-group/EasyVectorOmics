! ==============================================================================
! MODULE: gene_centroid_calculator
! ==============================================================================
! Description:
!   Provides procedures to compute gene expression centroids for orthogroups.
!   Supports "full" mode (all genes) and "orthologs" mode.
!   Designed following principles of clarity and modularity, avoiding GOTO.
! ==============================================================================
MODULE gene_centroid_calculator

  IMPLICIT NONE

  ! --- Public Interface ---
  PUBLIC :: orthogroup          ! Derived type for orthogroups
  PUBLIC :: mean_vector         ! Function to compute a single centroid
  PUBLIC :: group_centroid      ! Subroutine to compute all centroids

  ! --- Private Data ---
  PRIVATE ! Default to private for encapsulation

  ! --- Derived Type Definition ---
  TYPE :: orthogroup
     CHARACTER(LEN=50) :: id             ! Identifier for the orthogroup
     INTEGER, DIMENSION(:), ALLOCATABLE :: gene_indices ! Indices of member genes
  END TYPE orthogroup

CONTAINS

! ==============================================================================
! FUNCTION: mean_vector
! ==============================================================================
! Purpose:
!   Computes the element-wise mean (centroid) for a set of vectors
!   identified by their indices.
!
! Arguments:
!   vectors      (Input)  : REAL, DIMENSION(:,:), INTENT(IN)
!                           The full expression matrix (genes x tissues).
!   gene_indices (Input)  : INTEGER, DIMENSION(:), INTENT(IN)
!                           An array of row indices specifying which vectors
!                           (genes) to include in the mean calculation.
!
! Returns:
!   centroid     (Output) : REAL, DIMENSION(SIZE(vectors, 2))
!                           A 1D array representing the mean vector.
! ==============================================================================
  FUNCTION mean_vector(vectors, gene_indices) RESULT(centroid)
    ! --- Arguments ---
    REAL, DIMENSION(:, :), INTENT(IN) :: vectors      ! (n_genes, n_tissues)
    INTEGER, DIMENSION(:), INTENT(IN) :: gene_indices ! List of row indices

    ! --- Return Value ---
    REAL, DIMENSION(SIZE(vectors, 2)) :: centroid     ! (n_tissues)

    ! --- Local Variables ---
    INTEGER :: num_genes
    INTEGER :: num_tissues

    ! --- Implementation ---
    num_genes = SIZE(gene_indices)
    num_tissues = SIZE(vectors, 2)

    ! Case where no genes are selected
    IF (num_genes == 0) THEN
        centroid = 0.0  ! Return a zero vector
        RETURN
    END IF

    ! Check if gene_indices has any valid indices (optional, but good practice)
    ! For simplicity, we assume indices are valid.

    ! Compute the sum along dimension 1 (rows) using array slicing
    ! vectors(gene_indices, :) selects the specified rows and all columns.
    centroid = SUM(vectors(gene_indices, :), DIM=1)

    ! Divide by the number of genes to get the mean
    centroid = centroid / REAL(num_genes)

  END FUNCTION mean_vector

! ==============================================================================
! SUBROUTINE: group_centroid
! ==============================================================================
! Purpose:
!   Iterates over all orthogroups, filters genes based on the selected mode,
!   and computes the centroid for each group using 'mean_vector'.
!   The structure allows for potential parallelization over groups, similar
!   to the approach shown in the guidelines (Sources 19-25).
!
! Arguments:
!   vectors          (Input)  : REAL, DIMENSION(:,:), INTENT(IN)
!                               The full expression matrix (genes x tissues).
!   group_list       (Input)  : TYPE(orthogroup), DIMENSION(:), INTENT(IN)
!                               An array mapping orthogroup IDs to gene indices.
!   centroid_matrix  (Output) : REAL, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT)
!                               The resulting matrix (groups x tissues).
!   mode             (Input)  : CHARACTER(LEN=*), INTENT(IN)
!                               Calculation mode: "all" or "orthologs".
!   ortholog_set     (Input)  : LOGICAL, DIMENSION(:), INTENT(IN), OPTIONAL
!                               A logical mask (genes), TRUE if a gene is an
!                               ortholog. Required if mode = "orthologs".
! ==============================================================================
  SUBROUTINE group_centroid(vectors, group_list, centroid_matrix, mode, ortholog_set)
    ! --- Arguments ---
    REAL, DIMENSION(:, :), INTENT(IN) :: vectors          ! (n_genes, n_tissues)
    TYPE(orthogroup), DIMENSION(:), INTENT(IN) :: group_list ! (n_groups)
    REAL, DIMENSION(:, :), ALLOCATABLE, INTENT(OUT) :: centroid_matrix ! (n_groups, n_tissues)
    CHARACTER(LEN=*), INTENT(IN) :: mode
    LOGICAL, DIMENSION(:), INTENT(IN), OPTIONAL :: ortholog_set ! (n_genes)

    ! --- Local Variables ---
    INTEGER :: i_group, i_gene, j_gene
    INTEGER :: n_groups, n_tissues, n_genes_in_group
    INTEGER :: current_gene_index
    INTEGER :: selected_count
    INTEGER, DIMENSION(:), ALLOCATABLE :: selected_indices
    CHARACTER(LEN=10) :: current_mode
    LOGICAL :: use_orthologs

    ! --- Implementation ---

    n_groups = SIZE(group_list)
    n_tissues = SIZE(vectors, 2)

    ! Allocate the output matrix
    ALLOCATE(centroid_matrix(n_groups, n_tissues))
    centroid_matrix = 0.0 ! Initialize to zero

    ! Determine mode and check for 'ortholog_set' if needed
    current_mode = TRIM(ADJUSTL(mode))
    use_orthologs = (current_mode == 'orthologs')

    IF (use_orthologs .AND. .NOT. PRESENT(ortholog_set)) THEN
        WRITE(*,*) 'ERROR: Mode is "orthologs" but ortholog_set was not provided.'
        STOP 'group_centroid - Missing ortholog_set'
    END IF

    ! --- Main Loop: Iterate over each orthogroup ---
    ! This loop is independent for each group and could be parallelized
    ! using OpenMP or Coarrays, as suggested by the guidelines[cite: 19, 21].
    DO i_group = 1, n_groups

        n_genes_in_group = SIZE(group_list(i_group)%gene_indices)

        ! --- Step 1: Count how many genes will be selected ---
        selected_count = 0
        DO i_gene = 1, n_genes_in_group
            current_gene_index = group_list(i_group)%gene_indices(i_gene)

            ! Basic bounds check (optional but recommended)
            IF (current_gene_index < 1 .OR. current_gene_index > SIZE(vectors, 1)) THEN
                WRITE(*,*) 'WARNING: Gene index ', current_gene_index, &
                          ' for group ', TRIM(group_list(i_group)%id), ' is out of bounds.'
                CYCLE
            END IF

            IF (use_orthologs) THEN
                ! Check bounds for ortholog_set as well
                IF (current_gene_index > SIZE(ortholog_set)) THEN
                     WRITE(*,*) 'WARNING: Gene index ', current_gene_index, &
                          ' for group ', TRIM(group_list(i_group)%id), ' is out of bounds for ortholog_set.'
                     CYCLE
                END IF
                ! Check if it's an ortholog
                IF (PRESENT(ortholog_set) .AND. ortholog_set(current_gene_index)) THEN
                    selected_count = selected_count + 1
                END IF
            ELSE
                ! "all" mode: select every gene in the group (that's in bounds)
                selected_count = selected_count + 1
            END IF
        END DO

        ! If no genes were selected for this group, skip to the next
        IF (selected_count == 0) THEN
            CYCLE ! The centroid remains zero as initialized.
        END IF

        ! --- Step 2: Allocate and populate the selected indices array ---
        ALLOCATE(selected_indices(selected_count))
        j_gene = 0 ! Index for selected_indices

        DO i_gene = 1, n_genes_in_group
            current_gene_index = group_list(i_group)%gene_indices(i_gene)

            ! Re-check bounds (or trust Step 1 if bounds checks are robust)
            IF (current_gene_index < 1 .OR. current_gene_index > SIZE(vectors, 1)) CYCLE

            IF (use_orthologs) THEN
                IF (current_gene_index > SIZE(ortholog_set)) CYCLE
                IF (PRESENT(ortholog_set) .AND. ortholog_set(current_gene_index)) THEN
                    j_gene = j_gene + 1
                    selected_indices(j_gene) = current_gene_index
                END IF
            ELSE
                j_gene = j_gene + 1
                selected_indices(j_gene) = current_gene_index
            END IF
        END DO

        ! --- Step 3: Compute the centroid using mean_vector ---
        centroid_matrix(i_group, :) = mean_vector(vectors, selected_indices)

        ! --- Step 4: Deallocate the temporary array ---
        DEALLOCATE(selected_indices)

    END DO ! End of loop over groups

  END SUBROUTINE group_centroid

END MODULE gene_centroid_calculator