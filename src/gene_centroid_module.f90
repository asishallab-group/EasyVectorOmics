!> Module for computing expression centroids of gene families.
!| Adheres to F42 conventions: stateless, allocation-free loops, modular.
module gene_centroid_module
    use, intrinsic :: iso_fortran_env, only: int32, real64
    implicit none

    private
    public :: mean_vector, group_centroid

contains

    !> Computes the element-wise mean for a given set of vectors.
    pure subroutine mean_vector(vectors, gene_indices, n_genes, centroid)
        !| Matrix of vectors, each column corresponds to a gene
        real(real64), intent(in) :: vectors(:,:)
        !| Indices of the genes to average
        integer(int32), intent(in) :: gene_indices(:)
        !| Number of genes to average
        integer(int32), intent(in) :: n_genes
        !| Output centroid vector
        real(real64), intent(out) :: centroid(:)

        integer :: d, i, j
        integer(int32) :: gene_idx

        d = size(vectors, dim=1)
        centroid = 0.0_real64
        if (n_genes == 0) return

        do i = 1, n_genes
            gene_idx = gene_indices(i)
            !$omp simd
            do j = 1, d
                centroid(j) = centroid(j) + vectors(j, gene_idx)
            end do
            !$omp end simd
        end do
        centroid = centroid / real(n_genes, real64)
    end subroutine mean_vector

    !> Iterates over orthogroups, filters gene indices, and computes centroids.
    !| Allocation-free and processes `mode_ascii` as raw ASCII codes.
    subroutine group_centroid(vectors, num_genes, gene_to_family_map, num_families, &
                              centroid_matrix, mode_ascii, mode_len, ortholog_set_int, selected_indices)
        !| Matrix of vectors, columns correspond to genes
        real(real64), intent(in) :: vectors(:,:)
        !| Total number of genes
        integer(int32), intent(in) :: num_genes
        !| Mapping from gene index to family index
        integer(int32), intent(in) :: gene_to_family_map(num_genes)
        !| Total number of families
        integer(int32), intent(in) :: num_families
        !| Output matrix of centroids, each row is a family centroid
        real(real64), intent(out) :: centroid_matrix(:,:)
        !| Mode string as ASCII integer codes
        integer(int32), intent(in) :: mode_ascii(mode_len)
        !| Length of the mode string
        integer(int32), intent(in) :: mode_len
        !| Array indicating ortholog set membership (1 = in set, 0 = not)
        integer(int32), intent(in) :: ortholog_set_int(num_genes)
        !| Workspace buffer for selected gene indices
        integer(int32), intent(inout) :: selected_indices(:)

        integer(int32) :: i, j, current_count
        logical :: is_all_mode
        
        is_all_mode = .false.
        if (mode_len == 3) then
            if (mode_ascii(1) == ichar('a') .and. mode_ascii(2) == ichar('l') .and. mode_ascii(3) == ichar('l')) then
                is_all_mode = .true.
            end if
        end if
        
        do j = 1, num_families
            current_count = 0
            if (is_all_mode) then
                do i = 1, num_genes
                    if (gene_to_family_map(i) == j) then
                        current_count = current_count + 1
                        selected_indices(current_count) = i
                    end if
                end do
            else
                do i = 1, num_genes
                    if (gene_to_family_map(i) == j .and. ortholog_set_int(i) == 1) then
                        current_count = current_count + 1
                        selected_indices(current_count) = i
                    end if
                end do
            end if
            
            call mean_vector(vectors, selected_indices, current_count, centroid_matrix(j, :))
        end do
    end subroutine group_centroid

end module gene_centroid_module

! =============================================================================
! C and R Wrapper Subroutines
! =============================================================================

!> C interface for computing group centroids.
subroutine group_centroid_c(vectors, d, n, gene_to_family_map, num_families, &
                            centroid_matrix, mode_ascii, mode_len, ortholog_set_int, &
                            selected_indices, selected_indices_len) &
                            bind(c, name='group_centroid_c')
    use, intrinsic :: iso_c_binding
    use gene_centroid_module
    implicit none
    !| Matrix of vectors, columns correspond to genes
    real(c_double), intent(in) :: vectors(d, n)
    !| Vector dimension
    integer(c_int), value, intent(in) :: d
    !| Number of genes
    integer(c_int), value, intent(in) :: n
    !| Number of families
    integer(c_int), value, intent(in) :: num_families
    !| Length of the mode string
    integer(c_int), value, intent(in) :: mode_len
    !| Workspace buffer length
    integer(c_int), value, intent(in) :: selected_indices_len
    !| Mapping from gene index to family index
    integer(c_int), intent(in) :: gene_to_family_map(n)
    !| Output matrix of centroids
    real(c_double), intent(out) :: centroid_matrix(num_families, d)
    !| Mode string as ASCII integer codes
    integer(c_int), intent(in) :: mode_ascii(mode_len)
    !| Ortholog set membership array
    integer(c_int), intent(in) :: ortholog_set_int(n)
    !| Workspace buffer for selected gene indices
    integer(c_int), intent(inout) :: selected_indices(selected_indices_len)

    call group_centroid(vectors, n, gene_to_family_map, num_families, &
                        centroid_matrix, mode_ascii, mode_len, ortholog_set_int, selected_indices)
end subroutine group_centroid_c

!> R interface for computing group centroids.
subroutine group_centroid_r(vectors, d, n, gene_to_family_map, num_families, &
                            centroid_matrix, mode_ascii, mode_len, ortholog_set_int, &
                            selected_indices, selected_indices_len)
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use gene_centroid_module
    implicit none
    !| Matrix of vectors, columns correspond to genes
    real(real64), intent(in) :: vectors(d, n)
    !| Vector dimension
    integer(int32), intent(in) :: d
    !| Number of genes
    integer(int32), intent(in) :: n
    !| Number of families
    integer(int32), intent(in) :: num_families
    !| Length of the mode string
    integer(int32), intent(in) :: mode_len
    !| Workspace buffer length
    integer(int32), intent(in) :: selected_indices_len
    !| Mapping from gene index to family index
    integer(int32), intent(in) :: gene_to_family_map(n)
    !| Output matrix of centroids
    real(real64), intent(out) :: centroid_matrix(num_families, d)
    !| Mode string as ASCII integer codes
    integer(int32), intent(in) :: mode_ascii(mode_len)
    !| Ortholog set membership array
    integer(int32), intent(in) :: ortholog_set_int(n)
    !| Workspace buffer for selected gene indices
    integer(int32), intent(inout) :: selected_indices(selected_indices_len)

    call group_centroid(vectors, n, gene_to_family_map, num_families, &
                        centroid_matrix, mode_ascii, mode_len, ortholog_set_int, selected_indices)
end subroutine group_centroid_r
