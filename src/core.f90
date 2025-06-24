! ==============================================================================
! MODULE: gene_centroid_calculator
! ==============================================================================
! Description:
!   Provides procedures to compute gene expression centroids for orthogroups.
!   Supports "full" mode (all genes) and "orthologs" mode.
!
! Revision Notes:
!   - Corrected variable declaration placement to conform with Fortran standards.
!   - Refactored `mean_vector` to avoid creating temporary array copies from
!     non-contiguous row slicing
!   - Added OpenMP directives for parallelism based on the patterns
!     recommended in the coding guide (top-level `parallel do` and leaf-level `simd`).
! ==============================================================================
MODULE gene_centroid_calculator

  IMPLICIT NONE

  ! --- Public Interface ---
  PUBLIC :: orthogroup      ! Derived type for orthogroups
  PUBLIC :: group_centroid  ! Subroutine to compute all centroids

  ! --- Private Data ---
  PRIVATE ! Default to private for encapsulation

  ! --- Derived Type Definition ---
  TYPE :: orthogroup
     CHARACTER(LEN=50) :: id              ! Identifier for the orthogroup
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
! REVISION: This version was refactored to use an explicit loop. The original
! use of `SUM(vectors(gene_indices, :), DIM=1)` would likely create a
! temporary copy of the data because the rows selected by `gene_indices` are
! not contiguous in memory. This explicit loop avoids hidden copies. 
! ==============================================================================
  FUNCTION mean_vector(vectors, gene_indices) RESULT(centroid)
    ! --- Arguments ---
    REAL, DIMENSION(:, :), INTENT(IN) :: vectors      ! (n_genes, n_tissues)
    INTEGER, DIMENSION(:), INTENT(IN) :: gene_indices ! List of row indices

    ! --- Return Value ---
    REAL, DIMENSION(SIZE(vectors, 2)) :: centroid      ! (n_tissues)

    ! --- Local Variables ---
    INTEGER :: num_genes
    INTEGER :: i, j

    ! --- Implementation ---
    num_genes = SIZE(gene_indices)
    centroid = 0.0 ! Initialize to zero

    ! Case where no genes are selected
    IF (num_genes == 0) THEN
        RETURN
    END IF


    DO i = 1, num_genes
        !$OMP SIMD
        DO j = 1, SIZE(vectors, 2)
            centroid(j) = centroid(j) + vectors(gene_indices(i), j)
        END DO
    END DO

    ! Divide by the number of genes to get the mean
    centroid = centroid / REAL(num_genes)

  END FUNCTION mean_vector

! ==============================================================================
! SUBROUTINE: group_centroid
! ==============================================================================
! Purpose:
!   Iterates over all orthogroups, filters genes based on the selected mode,
!   and computes the centroid for each group using 'mean_vector'.
!   The loop over groups is parallelized with OpenMP.
! ==============================================================================
  SUBROUTINE group_centroid(vectors, group_list, centroid_matrix, mode, ortholog_set)
    ! --- Arguments ---
    REAL, DIMENSION(:, :), INTENT(IN) :: vectors        ! (n_genes, n_tissues)
    TYPE(orthogroup), DIMENSION(:), INTENT(IN) :: group_list ! (n_groups)
    REAL, DIMENSION(:, :), ALLOCATABLE, INTENT(OUT) :: centroid_matrix ! (n_groups, n_tissues)
    CHARACTER(LEN=*), INTENT(IN) :: mode
    LOGICAL, DIMENSION(:), INTENT(IN), OPTIONAL :: ortholog_set ! (n_genes)

    ! --- Local Variables ---
    INTEGER :: i_group
    INTEGER :: n_groups, n_tissues
    CHARACTER(LEN=10) :: current_mode
    LOGICAL :: use_orthologs
    INTEGER, DIMENSION(:), ALLOCATABLE :: selected_indices

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
    ! `selected_indices` must be private to each thread.
    !$OMP PARALLEL DO PRIVATE(i_group, selected_indices) SCHEDULE(STATIC)
    DO i_group = 1, n_groups
        ! --- Step 1: Filter gene indices based on mode ---
        IF (use_orthologs) THEN
            ! "orthologs" mode: Use PACK to select indices where the ortholog mask is true.
            ! The mask expression is passed directly to PACK to avoid a temp variable.
            selected_indices = PACK(group_list(i_group)%gene_indices, ortholog_set(group_list(i_group)%gene_indices))
        ELSE
            ! "all" mode: select every gene in the group by just copying the indices.
            selected_indices = group_list(i_group)%gene_indices
        END IF

        ! --- Step 2: Compute the centroid using mean_vector ---
        ! If selected_indices is empty, mean_vector will correctly return a zero vector.
        centroid_matrix(i_group, :) = mean_vector(vectors, selected_indices)
        
        ! Deallocate the thread-private array. The IF guard is a safe practice.
        IF (ALLOCATED(selected_indices)) THEN
            DEALLOCATE(selected_indices)
        END IF

    END DO !$OMP END PARALLEL DO

  END SUBROUTINE group_centroid

END MODULE gene_centroid_calculator