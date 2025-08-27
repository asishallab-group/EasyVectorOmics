module tox_data_validation
    use iso_fortran_env, only: real64, int32
    use tox_errors, only: ERR_OK, ERR_INVALID_INPUT, ERR_SIZE_MISMATCH
    implicit none
    private

    ! Public procedures
    public :: validate_data_structure
    public :: validate_gene_to_family_mapping
    public :: validate_expression_data
    public :: validate_family_centroids
    public :: validate_shift_vectors
    public :: check_for_nan_inf
contains
    subroutine validate_data_structure(n_genes, n_families, d, gene_ids, gene_family_ids, &
                                 gene_to_fam, expression_vectors, family_centroids, &
                                 shift_vectors, ierr)
        integer(int32), intent(in) :: n_genes, n_families, d
        character(len=*), intent(in) :: gene_ids(n_genes)
        character(len=*), intent(in) :: gene_family_ids(:)
        integer(int32), intent(in) :: gene_to_fam(n_genes)
        real(real64), intent(in) :: expression_vectors(:,:)
        real(real64), intent(in) :: family_centroids(:,:)
        real(real64), intent(in) :: shift_vectors(:,:)
        integer(int32), intent(out) :: ierr
        
        ierr = ERR_OK
        
        ! Check basic dimensions
        if (n_genes < 0 .or. n_families < 0 .or. d < 0) then
            ierr = ERR_INVALID_INPUT
            return
        end if
        
        ! Check gene_ids array
        if (size(gene_ids) /= n_genes) then
            ierr = ERR_SIZE_MISMATCH
            return
        end if
        
        ! Check gene_family_ids array
        if (size(gene_family_ids) /= n_families) then
            ierr = ERR_SIZE_MISMATCH
            return
        end if
        
        ! Check gene_to_fam array
        if (size(gene_to_fam) /= n_genes) then
            ierr = ERR_SIZE_MISMATCH
            return
        end if
        
        ! Check expression_vectors array
        if (size(expression_vectors, 1) /= d .or. size(expression_vectors, 2) /= n_genes) then
            ierr = ERR_SIZE_MISMATCH
            return
        end if
        
        ! Check family_centroids array
        if (n_families > 0) then
            if (size(family_centroids, 1) /= d .or. size(family_centroids, 2) /= n_families) then
                ierr = ERR_SIZE_MISMATCH
                return
            end if
        end if
        
        ! Check shift_vectors array
        if (size(shift_vectors, 1) /= 2*d .or. size(shift_vectors, 2) /= n_genes) then
            ierr = ERR_SIZE_MISMATCH
            return
        end if
        
    end subroutine validate_data_structure

    subroutine validate_gene_to_family_mapping(gene_to_fam, n_families, ierr)
        integer(int32), intent(in) :: gene_to_fam(:)
        integer(int32), intent(in) :: n_families
        integer(int32), intent(out) :: ierr
        
        integer :: i
        
        ierr = ERR_OK
        
        if (n_families == 0) then
            ! No families defined, mapping should be all zeros
            if (any(gene_to_fam /= 0)) then
                ierr = ERR_INVALID_INPUT
                return
            end if
        else
            ! Check that all gene_to_fam values are valid family indices
            do i = 1, size(gene_to_fam)
                if (gene_to_fam(i) < 1 .or. gene_to_fam(i) > n_families) then
                    ierr = ERR_INVALID_INPUT
                    return
                end if
            end do
        end if
        
    end subroutine validate_gene_to_family_mapping

    subroutine validate_expression_data(expression_vectors, check_non_negative, ierr)
        real(real64), intent(in) :: expression_vectors(:,:)
        logical, intent(in) :: check_non_negative
        integer(int32), intent(out) :: ierr
        
        ierr = ERR_OK
        
        ! Check for NaN and Inf values
        call check_for_nan_inf(expression_vectors, "expression_vectors", ierr)
        if (ierr /= ERR_OK) return
        
        ! Check for negative values if requested
        if (check_non_negative) then
            if (any(expression_vectors < 0.0_real64)) then
                ierr = ERR_INVALID_INPUT
                return
            end if
        end if
        
    end subroutine validate_expression_data

    subroutine validate_family_centroids(family_centroids, ierr)
        real(real64), intent(in) :: family_centroids(:,:)
        integer(int32), intent(out) :: ierr
        
        ierr = ERR_OK
        
        ! Check for NaN and Inf values
        call check_for_nan_inf(family_centroids, "family_centroids", ierr)
        if (ierr /= ERR_OK) return
        
    end subroutine validate_family_centroids

    subroutine validate_shift_vectors(shift_vectors, expression_vectors, family_centroids, &
                                    gene_to_fam, d, ierr)
        real(real64), intent(in) :: shift_vectors(:,:)
        real(real64), intent(in) :: expression_vectors(:,:)
        real(real64), intent(in) :: family_centroids(:,:)
        integer(int32), intent(in) :: gene_to_fam(:)
        integer(int32), intent(in) :: d
        integer(int32), intent(out) :: ierr
        
        integer :: i, fam_idx
        real(real64) :: expected_shift
        
        ierr = ERR_OK
        
        ! Check for NaN and Inf values
        call check_for_nan_inf(shift_vectors, "shift_vectors", ierr)
        if (ierr /= ERR_OK) return
        
        ! Verify shift vectors structure: first d rows should be centroids, next d rows should be shifts
        do i = 1, size(shift_vectors, 2)
            fam_idx = gene_to_fam(i)
            
            ! Check that centroid part matches the family centroid
            if (any(abs(shift_vectors(1:d, i) - family_centroids(:, fam_idx)) > 1.0e-10_real64)) then
                ierr = ERR_INVALID_INPUT
                return
            end if
            
            ! Check that shift part matches expression - centroid
            expected_shift = expression_vectors(:, i) - family_centroids(:, fam_idx)
            if (any(abs(shift_vectors(d+1:2*d, i) - expected_shift) > 1.0e-10_real64)) then
                ierr = ERR_INVALID_INPUT
                return
            end if
        end do
        
    end subroutine validate_shift_vectors

    subroutine check_for_nan_inf(array, array_name, ierr)
        real(real64), intent(in) :: array(:,:)
        character(len=*), intent(in) :: array_name
        integer(int32), intent(out) :: ierr
        
        integer :: i, j
        
        ierr = ERR_OK
        
        do j = 1, size(array, 2)
            do i = 1, size(array, 1)
                if (is_nan(array(i, j)) .or. is_inf(array(i, j))) then
                    ierr = ERR_INVALID_INPUT
                    return
                end if
            end do
        end do
        
    contains

        logical function is_nan(x)
            real(real64), intent(in) :: x
            is_nan = (x /= x)  ! NaN is the only value not equal to itself
        end function is_nan
        
        logical function is_inf(x)
            real(real64), intent(in) :: x
            is_inf = (x > huge(x) .or. x < -huge(x))
        end function is_inf
        
    end subroutine check_for_nan_inf

    subroutine validate_gene_ids_uniqueness(gene_ids, ierr)
        character(len=*), intent(in) :: gene_ids(:)
        integer(int32), intent(out) :: ierr
        
        integer :: i, j
        
        ierr = ERR_OK
        
        do i = 1, size(gene_ids) - 1
            do j = i + 1, size(gene_ids)
                if (trim(gene_ids(i)) == trim(gene_ids(j))) then
                    ierr = ERR_INVALID_INPUT
                    return
                end if
            end do
        end do
        
    end subroutine validate_gene_ids_uniqueness

    subroutine validate_family_ids_uniqueness(gene_family_ids, ierr)
        character(len=*), intent(in) :: gene_family_ids(:)
        integer(int32), intent(out) :: ierr
        
        integer :: i, j
        
        ierr = ERR_OK
        
        do i = 1, size(gene_family_ids) - 1
            do j = i + 1, size(gene_family_ids)
                if (trim(gene_family_ids(i)) == trim(gene_family_ids(j))) then
                    ierr = ERR_INVALID_INPUT
                    return
                end if
            end do
        end do
        
    end subroutine validate_family_ids_uniqueness

    subroutine validate_empty_strings(string_array, array_name, ierr)
        character(len=*), intent(in) :: string_array(:)
        character(len=*), intent(in) :: array_name
        integer(int32), intent(out) :: ierr
        
        integer :: i
        
        ierr = ERR_OK
        
        do i = 1, size(string_array)
            if (len_trim(string_array(i)) == 0) then
                ierr = ERR_INVALID_INPUT
                return
            end if
        end do
        
    end subroutine validate_empty_strings
end module tox_data_validation