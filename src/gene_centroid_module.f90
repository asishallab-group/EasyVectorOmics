! Module for computing expression centroids of gene families.
!
! This module computes the mean expression vector (centroid) for specified
! groups of genes (e.g., families). It can operate in two modes:
!  1. 'all': Computes the centroid using all genes in a family.
!  2. 'ortho': Computes the centroid using only a pre-defined subset of
!             orthologs within each family.
!
! @note Fortran is column-major, so all vector data is stored with genes
!       as columns. `vectors(dimension, num_genes)`.
module gene_centroid_module
    use, intrinsic :: iso_fortran_env, only: int32, real64
    implicit none

    private
    public :: group_centroid, mean_vector

contains

    ! Computes the element-wise mean for a given set of vectors. It sums vectors specified by
    ! a logical mask and writes the resulting mean vector.
    pure subroutine mean_vector(vectors, d, n, selection_mask, centroid_col)
        ! Input matrix of vectors (genes as columns).
        integer(int32), intent(in) :: d, n
        real(real64), intent(in) :: vectors(d, n)
        ! Logical mask indicating which genes to average.
        logical, intent(in) :: selection_mask(n)
        ! Output centroid vector (a single column).
        real(real64), intent(out) :: centroid_col(d)

        integer(int32) :: i, j, num_selected_genes
        real(real64) :: inv_n_genes

        num_selected_genes = count(selection_mask)

        ! Initialize the output centroid column to zero.
        centroid_col = 0.0_real64
        if (num_selected_genes == 0) return

        ! Loop over each dimension of the vectors.
        do j = 1, d
            !$OMP SIMD REDUCTION(+:centroid_col(j))
            do i = 1, n
                if (selection_mask(i)) then
                    centroid_col(j) = centroid_col(j) + vectors(j, i)
                end if
            end do
            !$OMP END SIMD
        end do

        ! Normalize in-place to get the mean, avoiding a temporary array.
        inv_n_genes = 1.0_real64 / real(num_selected_genes, real64)
        do j = 1, d
            centroid_col(j) = centroid_col(j) * inv_n_genes
        end do
    end subroutine mean_vector

    ! Iterates over families, filters gene indices, and computes centroids.
    ! This is the main entry point. It orchestrates the process of
    ! identifying genes for each family and calling the mean_vector kernel.
    ! The loop over families is parallelized using OpenMP.
    subroutine group_centroid(vectors, d, num_genes, gene_to_family_map, num_families, &
                              centroid_matrix, use_all_mode, ortholog_set)
        ! Matrix of all gene vectors (genes as columns).
        integer(int32), intent(in) :: d, num_genes, num_families
        real(real64), intent(in) :: vectors(d, num_genes)
        ! Map from gene index to family index.
        integer(int32), intent(in) :: gene_to_family_map(num_genes)
        ! Output matrix for centroids (centroids as columns).
        real(real64), intent(out) :: centroid_matrix(d, num_families)
        ! Logical flag to select mode: .true. for 'all', .false. for 'ortho'.
        logical, intent(in) :: use_all_mode
        ! Array indicating ortholog set membership.
        logical, intent(in) :: ortholog_set(num_genes)

        integer(int32) :: i, j
        logical :: selection_mask(num_genes)

        ! Parallelize the loop over families. Each family centroid can be
        ! computed independently.
        !$OMP PARALLEL DO PRIVATE(j, i, selection_mask)
        do j = 1, num_families
            ! This loop builds a logical mask for genes belonging to the
            ! current family `j`, based on the selected mode.
            do i = 1, num_genes
                selection_mask(i) = (gene_to_family_map(i) == j .and. (use_all_mode .or. ortholog_set(i)))
            end do

            ! Compute the centroid for the selected genes.
            ! Pass a single column of the output matrix.
            call mean_vector(vectors, d, num_genes, selection_mask, centroid_matrix(:, j))
        end do
        !$OMP END PARALLEL DO
    end subroutine group_centroid

end module gene_centroid_module

! =============================================================================
! C and R Wrapper Subroutines
! =============================================================================

! C interface for computing group centroids.
! This wrapper makes the Fortran routine callable from C-compatible languages.
! It handles the translation between C and Fortran data types and array indexing.
subroutine group_centroid_c(vectors, d, n, gene_to_family_map, num_families, &
                            centroid_matrix, use_all_mode, ortholog_set) &
                            bind(c, name='group_centroid_c')
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_bool
    use gene_centroid_module, only: group_centroid
    implicit none
    ! Arguments passed from C
    integer(c_int), value, intent(in) :: d, n, num_families
    real(c_double), intent(in) :: vectors(d, n)
    integer(c_int), intent(in) :: gene_to_family_map(n)
    logical(c_bool), value, intent(in) :: use_all_mode
    logical(c_bool), intent(in) :: ortholog_set(n)
    ! Output/Workspace arguments
    real(c_double), intent(out) :: centroid_matrix(d, num_families)

    ! Local variables to match Fortran's default logical kind.
    logical :: use_all_mode_fortran
    logical :: ortholog_set_fortran(n)

    use_all_mode_fortran = use_all_mode
    ortholog_set_fortran = ortholog_set

    ! Call the core Fortran implementation with the kind-matched logicals.
    call group_centroid(vectors, d, n, gene_to_family_map, num_families, &
                        centroid_matrix, use_all_mode_fortran, ortholog_set_fortran)
end subroutine group_centroid_c

! R interface for computing group centroids.
subroutine group_centroid_r(vectors, d, n, gene_to_family_map, num_families, &
                            centroid_matrix, use_all_mode, ortholog_set)
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use gene_centroid_module, only: group_centroid
    implicit none
    ! Arguments passed from R
    integer(int32), intent(in) :: d, n, num_families
    real(real64), intent(in) :: vectors(d, n)
    integer(int32), intent(in) :: gene_to_family_map(n)
    logical, intent(in) :: use_all_mode
    logical, intent(in) :: ortholog_set(n)
    ! Output/Workspace arguments
    real(real64), intent(out) :: centroid_matrix(d, num_families)

    ! Call the core Fortran implementation.
    call group_centroid(vectors, d, n, gene_to_family_map, num_families, &
                        centroid_matrix, use_all_mode, ortholog_set)
end subroutine group_centroid_r
