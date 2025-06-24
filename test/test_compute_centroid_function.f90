! ==============================================================================
! TEST PROGRAM: test_gene_centroid_calculator
! ==============================================================================
! Purpose:
!   Provides a test suite for the `gene_centroid_calculator` module.
!   This program can be compiled and run with the Fortran Package Manager (fpm)
!   by placing it in the `tests/` directory and executing `fpm test`.
!
! Tests Included:
!   - `test_all_mode`: Verifies centroid calculation when all genes in an
!     orthogroup are included.
!   - `test_orthologs_mode`: Verifies centroid calculation when only a subset
!     of genes (marked as orthologs) are included.
! ==============================================================================
PROGRAM test_gene_centroid_calculator
  USE, INTRINSIC :: iso_fortran_env, ONLY: ERROR_UNIT
  IMPLICIT NONE

  WRITE(*, '(A)') 'Running tests for gene_centroid_calculator...'
  
  CALL test_all_mode()
  CALL test_orthologs_mode()
  
  WRITE(*, '(A)') 'All tests passed successfully.'

CONTAINS

! ==============================================================================
! SUBROUTINE: test_all_mode
! ==============================================================================
! Verifies that the `group_centroid` subroutine works correctly in "all" mode.
! It sets up a simple scenario with two orthogroups and calculates the expected
! centroids manually, then compares them to the output of the subroutine.
! ==============================================================================
  SUBROUTINE test_all_mode()
    USE gene_centroid_calculator, ONLY: orthogroup, group_centroid
    IMPLICIT NONE

    ! --- Test Data Setup ---
    REAL, DIMENSION(5, 3) :: expression_vectors
    TYPE(orthogroup), DIMENSION(2) :: groups
    REAL, DIMENSION(:,:), ALLOCATABLE :: calculated_centroids
    REAL, DIMENSION(2, 3) :: expected_centroids
    REAL, PARAMETER :: tolerance = 1E-6

    WRITE(*, '(/,A)') '--- Running Test: "all" mode ---'

    ! Mock expression data (5 genes, 3 tissues)
    ! NOTE: Using legacy (/ ... /) array constructor syntax for broad compiler compatibility.
    expression_vectors = RESHAPE((/ &
      1.0, 2.0, 3.0,  &! Gene 1
      4.0, 5.0, 6.0,  &! Gene 2
      7.0, 8.0, 9.0,  &! Gene 3
      10.0, 11.0, 12.0, &! Gene 4
      13.0, 14.0, 15.0  &! Gene 5
    /), SHAPE(expression_vectors), ORDER=(/2,1/))

    ! Define two orthogroups
    ! Group 1: contains genes 1, 2
    ALLOCATE(groups(1)%gene_indices(2))
    groups(1)%id = 'OG001'
    groups(1)%gene_indices = (/1, 2/)
    
    ! Group 2: contains genes 3, 5
    ALLOCATE(groups(2)%gene_indices(2))
    groups(2)%id = 'OG002'
    groups(2)%gene_indices = (/3, 5/)

    ! --- Call the Subroutine Under Test ---
    CALL group_centroid(expression_vectors, groups, calculated_centroids, mode='all')

    ! --- Verification ---
    ! Expected centroid for Group 1: AVG([1,2,3], [4,5,6]) = [2.5, 3.5, 4.5]
    expected_centroids(1, :) = (/2.5, 3.5, 4.5/)
    
    ! Expected centroid for Group 2: AVG([7,8,9], [13,14,15]) = [10.0, 11.0, 12.0]
    expected_centroids(2, :) = (/10.0, 11.0, 12.0/)

    IF (ANY(ABS(calculated_centroids - expected_centroids) > tolerance)) THEN
      WRITE(ERROR_UNIT, '(A)') 'TEST FAILED: "all" mode produced incorrect centroids.'
      WRITE(ERROR_UNIT, *) 'Expected:', expected_centroids
      WRITE(ERROR_UNIT, *) 'Got:     ', calculated_centroids
      ERROR STOP 1
    END IF

    WRITE(*, '(A)') 'PASSED'
    DEALLOCATE(groups(1)%gene_indices, groups(2)%gene_indices, calculated_centroids)
  END SUBROUTINE test_all_mode

! ==============================================================================
! SUBROUTINE: test_orthologs_mode
! ==============================================================================
! Verifies that the `group_centroid` subroutine works correctly in "orthologs"
! mode. It uses the same base data as the "all" mode test but applies a
! logical mask to filter genes, simulating an ortholog-only analysis.
! ==============================================================================
  SUBROUTINE test_orthologs_mode()
    USE gene_centroid_calculator, ONLY: orthogroup, group_centroid
    IMPLICIT NONE

    ! --- Test Data Setup ---
    REAL, DIMENSION(5, 3) :: expression_vectors
    TYPE(orthogroup), DIMENSION(2) :: groups
    LOGICAL, DIMENSION(5) :: ortholog_mask
    REAL, DIMENSION(:,:), ALLOCATABLE :: calculated_centroids
    REAL, DIMENSION(2, 3) :: expected_centroids
    REAL, PARAMETER :: tolerance = 1E-6
    
    WRITE(*, '(/,A)') '--- Running Test: "orthologs" mode ---'

    ! Mock expression data (same as before)
    expression_vectors = RESHAPE((/ &
      1.0, 2.0, 3.0,  &! Gene 1
      4.0, 5.0, 6.0,  &! Gene 2
      7.0, 8.0, 9.0,  &! Gene 3
      10.0, 11.0, 12.0, &! Gene 4
      13.0, 14.0, 15.0  &! Gene 5
    /), SHAPE(expression_vectors), ORDER=(/2,1/))

    ! Define two orthogroups (same as before)
    ALLOCATE(groups(1)%gene_indices(2))
    groups(1)%id = 'OG001'
    groups(1)%gene_indices = (/1, 2/)
    
    ALLOCATE(groups(2)%gene_indices(2))
    groups(2)%id = 'OG002'
    groups(2)%gene_indices = (/3, 5/)

    ! Define a logical mask for orthologs (genes 1, 3, 5 are orthologs)
    ortholog_mask = (/.TRUE., .FALSE., .TRUE., .FALSE., .TRUE./)

    ! --- Call the Subroutine Under Test ---
    CALL group_centroid(expression_vectors, groups, calculated_centroids, &
                        mode='orthologs', ortholog_set=ortholog_mask)

    ! --- Verification ---
    ! Group 1 genes are [1, 2]. Mask is [T, F]. Only gene 1 is used.
    ! Expected centroid for Group 1: [1.0, 2.0, 3.0]
    expected_centroids(1, :) = (/1.0, 2.0, 3.0/)
    
    ! Group 2 genes are [3, 5]. Mask is [T, T]. Both genes are used.
    ! Expected centroid for Group 2: AVG([7,8,9], [13,14,15]) = [10.0, 11.0, 12.0]
    expected_centroids(2, :) = (/10.0, 11.0, 12.0/)

    IF (ANY(ABS(calculated_centroids - expected_centroids) > tolerance)) THEN
      WRITE(ERROR_UNIT, '(A)') 'TEST FAILED: "orthologs" mode produced incorrect centroids.'
      WRITE(ERROR_UNIT, *) 'Expected:', expected_centroids
      WRITE(ERROR_UNIT, *) 'Got:     ', calculated_centroids
      ERROR STOP 1
    END IF
    
    WRITE(*, '(A)') 'PASSED'
    DEALLOCATE(groups(1)%gene_indices, groups(2)%gene_indices, calculated_centroids)
  END SUBROUTINE test_orthologs_mode

END PROGRAM test_gene_centroid_calculator