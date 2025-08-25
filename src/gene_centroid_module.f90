! Module for computing expression centroids of gene families.
!
! This module contains the core scientific kernel and the C-interface wrapper.
! The R-interface wrapper is defined outside the module for compatibility.
module gene_centroid_module
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    implicit none

    private
    public :: group_centroid, mean_vector, group_centroid_c

contains

    ! Computes the element-wise mean for a given set of vectors.
    ! REMOVED 'PURE' to allow for more complex OpenMP directives.
    subroutine mean_vector(vectors, d, gene_indices, num_selected_genes, centroid_col)
        integer(int32), intent(in) :: d
        real(real64), intent(in) :: vectors(d, *)
        integer(int32), intent(in) :: num_selected_genes
        integer(int32), intent(in) :: gene_indices(num_selected_genes)
        real(real64), intent(out) :: centroid_col(d)

        integer(int32) :: i, j, gene_idx
        real(real64) :: inv_n_genes
        real(real64) :: sum_val ! Temporary scalar for reduction

        centroid_col = 0.0_real64
        if (num_selected_genes == 0) return

        do j = 1, d
            sum_val = 0.0_real64
            ! CORRECTED OMP DIRECTIVE: Reduction is now on a scalar 'sum_val'.
            !$OMP SIMD REDUCTION(+:sum_val)
            do i = 1, num_selected_genes
                gene_idx = gene_indices(i)
                sum_val = sum_val + vectors(j, gene_idx)
            end do
            !$OMP END SIMD
            centroid_col(j) = sum_val
        end do

        inv_n_genes = 1.0_real64 / real(num_selected_genes, real64)
        do j = 1, d
            centroid_col(j) = centroid_col(j) * inv_n_genes
        end do
    end subroutine mean_vector

    ! Iterates over families, filters gene indices, and computes centroids.
    subroutine group_centroid(vectors, d, num_genes, gene_to_family_map, num_families, &
                              centroid_matrix, use_all_mode, ortholog_set, selected_indices)
        real(real64), intent(in) :: vectors(d, num_genes)
        integer(int32), intent(in) :: d, num_genes, num_families
        integer(int32), intent(in) :: gene_to_family_map(num_genes)
        real(real64), intent(out) :: centroid_matrix(d, num_families)
        logical, intent(in) :: use_all_mode
        logical, intent(in) :: ortholog_set(num_genes)
        integer(int32), intent(out) :: selected_indices(num_genes)

        integer(int32) :: i, j, num_selected
        ! NOTE: selected_indices should be private to each thread to avoid race conditions.
        ! A better approach for larger scale would be to use a different parallelization strategy,
        ! but for now, making it private is the most direct fix.
        integer(int32) :: local_selected_indices(num_genes)


        !$OMP PARALLEL DO PRIVATE(j, i, num_selected, local_selected_indices)
        do j = 1, num_families
            num_selected = 0
            do i = 1, num_genes
                if (gene_to_family_map(i) == j .and. (use_all_mode .or. ortholog_set(i))) then
                    num_selected = num_selected + 1
                    local_selected_indices(num_selected) = i
                end if
            end do
            call mean_vector(vectors, d, local_selected_indices, num_selected, centroid_matrix(:, j))
        end do
        !$OMP END PARALLEL DO
    end subroutine group_centroid

    ! C interface for group_centroid
    subroutine group_centroid_c(vectors, d, n, gene_to_family_map, num_families, &
                                centroid_matrix, use_all_mode, ortholog_set, &
                                selected_indices, selected_indices_len) &
                                bind(c, name='group_centroid_c')
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

end module gene_centroid_module

! =============================================================================
! R Wrapper Subroutine (defined outside the module)
! =============================================================================
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
