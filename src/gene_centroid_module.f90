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

    ! Computes the element-wise mean for a given set of vectors.
    ! It sums vectors specified by an array of indices.
    pure subroutine mean_vector(vectors, d, gene_indices, num_selected_genes, centroid_col)
        !| Input matrix of vectors (genes as columns).
        integer(int32), intent(in) :: d
        real(real64), intent(in) :: vectors(d, *)
        !| Number of genes to average.
        integer(int32), intent(in) :: num_selected_genes
        !| Indices of the genes to average.
        integer(int32), intent(in) :: gene_indices(num_selected_genes)
        !| Output centroid vector (a single column).
        real(real64), intent(out) :: centroid_col(d)

        integer(int32) :: i, j, gene_idx
        real(real64) :: inv_n_genes

        ! Initialize the output centroid column to zero.
        centroid_col = 0.0_real64
        if (num_selected_genes == 0) return

        ! Loop over each dimension of the vectors.
        do j = 1, d
            !$OMP SIMD REDUCTION(+:centroid_col(j))
            do i = 1, num_selected_genes
                gene_idx = gene_indices(i)
                centroid_col(j) = centroid_col(j) + vectors(j, gene_idx)
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
    subroutine group_centroid(vectors, d, num_genes, gene_to_family_map, num_families, &
                              centroid_matrix, use_all_mode, ortholog_set, selected_indices)
        !| Matrix of all gene vectors (genes as columns).
        real(real64), intent(in) :: vectors(d, num_genes)
        !| The dimension of the expression vectors (number of rows in `vectors`).
        integer(int32), intent(in) :: d
        !| The total number of genes (number of columns in `vectors`).
        integer(int32), intent(in) :: num_genes
        !| A map from each gene index to its corresponding family index.
        integer(int32), intent(in) :: gene_to_family_map(num_genes)
        !| The total number of unique families.
        integer(int32), intent(in) :: num_families
        !| The output matrix where each column will store a computed family centroid.
        real(real64), intent(out) :: centroid_matrix(d, num_families)
        !| A logical flag to select the mode: .true. for 'all' genes, .false. for 'orthologs' only.
        logical, intent(in) :: use_all_mode
        !| A logical array indicating ortholog set membership for each gene.
        logical, intent(in) :: ortholog_set(num_genes)
        !| A workspace buffer used to store the indices of selected genes for each family.
        integer(int32), intent(out) :: selected_indices(num_genes)

        integer(int32) :: i, j, num_selected

        !$OMP PARALLEL DO PRIVATE(j, i, num_selected, selected_indices)
        do j = 1, num_families
            num_selected = 0
            do i = 1, num_genes
                if (gene_to_family_map(i) == j .and. (use_all_mode .or. ortholog_set(i))) then
                    num_selected = num_selected + 1
                    selected_indices(num_selected) = i
                end if
            end do

            call mean_vector(vectors, d, selected_indices, num_selected, centroid_matrix(:, j))
        end do
        !$OMP END PARALLEL DO
    end subroutine group_centroid

end module gene_centroid_module

! =============================================================================
! C and R Wrapper Subroutines
! =============================================================================

subroutine group_centroid_c(vectors, d, n, gene_to_family_map, num_families, &
                            centroid_matrix, use_all_mode, ortholog_set, &
                            selected_indices, selected_indices_len) &
                            bind(c, name='group_centroid_c')
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_bool
    use gene_centroid_module, only: group_centroid
    implicit none
    integer(c_int), value, intent(in) :: d, n, num_families, selected_indices_len
    real(c_double), intent(in) :: vectors(d, n)
    integer(c_int), intent(in) :: gene_to_family_map(n)
    integer(c_int), value, intent(in) :: use_all_mode
    integer(c_int), intent(in) :: ortholog_set(n)
    real(c_double), intent(out) :: centroid_matrix(d, num_families)
    integer(c_int), intent(out) :: selected_indices(selected_indices_len)

    logical :: use_all_mode_fortran
    logical :: ortholog_set_fortran(n)
    integer :: i

    use_all_mode_fortran = (use_all_mode /= 0)
    do i = 1, n
        ortholog_set_fortran(i) = (ortholog_set(i) /= 0)
    end do

    call group_centroid(vectors, d, n, gene_to_family_map, num_families, &
                        centroid_matrix, use_all_mode_fortran, ortholog_set_fortran, selected_indices)
end subroutine group_centroid_c

subroutine group_centroid_r(vectors, d, n, gene_to_family_map, num_families, &
                            centroid_matrix, use_all_mode, ortholog_set, &
                            selected_indices, selected_indices_len)
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use gene_centroid_module, only: group_centroid
    implicit none
    integer(int32), intent(in) :: d, n, num_families, selected_indices_len
    real(real64), intent(in) :: vectors(d, n)
    integer(int32), intent(in) :: gene_to_family_map(n)
    logical, intent(in) :: use_all_mode
    logical, intent(in) :: ortholog_set(n)
    real(real64), intent(out) :: centroid_matrix(d, num_families)
    integer(int32), intent(out) :: selected_indices(selected_indices_len)

    call group_centroid(vectors, d, n, gene_to_family_map, num_families, &
                        centroid_matrix, use_all_mode, ortholog_set, selected_indices)
end subroutine group_centroid_r
