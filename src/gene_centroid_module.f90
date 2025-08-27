! Module for computing expression centroids of gene families.
!
! This module contains the core scientific kernel. The C and R interface
! wrappers are defined outside the module for compatibility.
module gene_centroid_module
    use, intrinsic :: iso_fortran_env, only: int32, real64
    implicit none

    private
    public :: group_centroid, mean_vector

contains

    !> Computes the element-wise mean for a given set of vectors.
    pure subroutine mean_vector(vectors, d, total_num_genes, gene_indices, num_selected_genes, centroid_col)
        implicit none
        !| Dimensionality of the vectors (e.g., number of tissues).
        integer(int32), intent(in) :: d
        !| Total number of genes in the input matrix.
        integer(int32), intent(in) :: total_num_genes
        !| The input matrix of all gene vectors (d x total_num_genes).
        real(real64), intent(in) :: vectors(d, total_num_genes)
        !| The number of genes in the current family to be averaged.
        integer(int32), intent(in) :: num_selected_genes
        !| An array containing the column indices of the selected genes in 'vectors'.
        integer(int32), intent(in) :: gene_indices(num_selected_genes)
        !| The output vector representing the computed centroid.
        real(real64), intent(out) :: centroid_col(d)

        ! Local variables
        integer(int32) :: i, j, gene_idx
        real(real64) :: inv_n_genes
        real(real64) :: sum_val

        centroid_col = 0.0_real64
        if (num_selected_genes == 0) return

        do j = 1, d
            sum_val = 0.0_real64
            do i = 1, num_selected_genes
                gene_idx = gene_indices(i)
                sum_val = sum_val + vectors(j, gene_idx)
            end do
            centroid_col(j) = sum_val
        end do

        inv_n_genes = 1.0_real64 / real(num_selected_genes, real64)
        centroid_col = centroid_col * inv_n_genes
    end subroutine mean_vector

    !> Iterates over families, filters gene indices, and computes centroids.
    pure subroutine group_centroid(vectors, d, num_genes, gene_to_family_map, num_families, &
                              centroid_matrix, use_all_mode, ortholog_set, selected_indices)
        implicit none
        !| Dimensionality of the expression vectors.
        integer(int32), intent(in) :: d
        !| Total number of genes in the 'vectors' matrix.
        integer(int32), intent(in) :: num_genes
        !| Total number of gene families to compute centroids for.
        integer(int32), intent(in) :: num_families
        !| The input matrix of all gene expression vectors (d x num_genes).
        real(real64), intent(in) :: vectors(d, num_genes)
        !| An array mapping each gene (by index) to a family ID.
        integer(int32), intent(in) :: gene_to_family_map(num_genes)
        !| The output matrix (d x num_families) to store the computed centroids.
        real(real64), intent(out) :: centroid_matrix(d, num_families)
        !| A logical flag; if true, all genes in a family are used.
        logical, intent(in) :: use_all_mode
        !| A logical array indicating if a gene is part of a specific subset (e.g., orthologs).
        logical, intent(in) :: ortholog_set(num_genes)
        !| An output array for storing indices.
        integer(int32), intent(out) :: selected_indices(num_genes)

        ! Local variables
        integer(int32) :: i, j, num_selected
        integer(int32) :: local_selected_indices(num_genes)

        selected_indices = 0

        do j = 1, num_families
            num_selected = 0
            do i = 1, num_genes
                if (gene_to_family_map(i) == j .and. (use_all_mode .or. ortholog_set(i))) then
                    num_selected = num_selected + 1
                    local_selected_indices(num_selected) = i
                end if
            end do
            call mean_vector(vectors, d, num_genes, local_selected_indices, num_selected, centroid_matrix(:, j))
        end do
    end subroutine group_centroid

end module gene_centroid_module

! =============================================================================
! C Wrapper Subroutine
! =============================================================================
!> C interface wrapper for group_centroid.
pure subroutine group_centroid_c(vectors, d, n, gene_to_family_map, num_families, &
                                  centroid_matrix, use_all_mode, ortholog_set, &
                                  selected_indices, selected_indices_len) &
                                  bind(c, name='group_centroid_c')
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use gene_centroid_module, only: group_centroid
    implicit none
    !| Dimensionality of the vectors.
    integer(c_int), value, intent(in) :: d
    !| Total number of genes.
    integer(c_int), value, intent(in) :: n
    !| Total number of families.
    integer(c_int), value, intent(in) :: num_families
    !| The allocated length of the 'selected_indices' array.
    integer(c_int), value, intent(in) :: selected_indices_len
    !| Input expression vectors (passed from C).
    real(c_double), intent(in) :: vectors(d, n)
    !| Array mapping gene index to family ID.
    integer(c_int), intent(in) :: gene_to_family_map(n)
    !| Integer flag from C (0=false, non-zero=true) to use all genes.
    integer(c_int), value, intent(in) :: use_all_mode
    !| Integer array from C indicating subset membership.
    integer(c_int), intent(in) :: ortholog_set(n)
    !| Output matrix for centroids.
    real(c_double), intent(out) :: centroid_matrix(d, num_families)
    !| Output array for selected indices.
    integer(c_int), intent(out) :: selected_indices(selected_indices_len)

    ! Local variables
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

! =============================================================================
! R Wrapper Subroutine
! =============================================================================
!> R interface wrapper for group_centroid.
pure subroutine group_centroid_r(vectors, d, n, gene_to_family_map, num_families, &
                              centroid_matrix, use_all_mode, ortholog_set, &
                              selected_indices, selected_indices_len)
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use gene_centroid_module, only: group_centroid
    implicit none
    !| Dimensionality of the expression vectors.
    integer(int32), intent(in) :: d
    !| Total number of genes in the 'vectors' matrix.
    integer(int32), intent(in) :: n
    !| Total number of gene families to compute centroids for.
    integer(int32), intent(in) :: num_families
    !| The allocated length of the 'selected_indices' array.
    integer(int32), intent(in) :: selected_indices_len
    !| The input matrix of all gene expression vectors (d x n).
    real(real64), intent(in) :: vectors(d, n)
    !| An array mapping each gene (by index) to a family ID.
    integer(int32), intent(in) :: gene_to_family_map(n)
    !| A logical flag; if true, all genes in a family are used.
    logical, intent(in) :: use_all_mode
    !| A logical array indicating if a gene is part of a specific subset.
    logical, intent(in) :: ortholog_set(n)
    !| The output matrix (d x num_families) to store the computed centroids.
    real(real64), intent(out) :: centroid_matrix(d, num_families)
    !| An output array for storing selected gene indices.
    integer(int32), intent(out) :: selected_indices(selected_indices_len)

    call group_centroid(vectors, d, n, gene_to_family_map, num_families, &
                        centroid_matrix, use_all_mode, ortholog_set, selected_indices)
end subroutine group_centroid_r