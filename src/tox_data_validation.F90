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
    public :: validate_gene_ids_uniqueness
    public :: validate_family_ids_uniqueness
    public :: validate_empty_strings
    
    ! Parameters for validation tolerances
    real(real64), parameter :: FLOAT_TOLERANCE = 1.0e-10_real64
    real(real64), parameter :: MISSING_VALUE = -999.0_real64
    
contains

    subroutine validate_data_structure(n_genes, n_families, n_samples, d, gene_ids, gene_family_ids, &
                                     gene_to_fam, expression_vectors, family_centroids, &
                                     shift_vectors, ierr)
        integer(int32), intent(in) :: n_genes, n_families, n_samples, d
        character(len=*), intent(in) :: gene_ids(:)
        character(len=*), intent(in) :: gene_family_ids(:)
        integer(int32), intent(in) :: gene_to_fam(:)
        real(real64), intent(in) :: expression_vectors(:,:)
        real(real64), intent(in) :: family_centroids(:,:)
        real(real64), intent(in) :: shift_vectors(:,:)
        integer(int32), intent(out) :: ierr
        
        ierr = ERR_OK
        
        ! Check basic dimensions
        if (n_genes <= 0 .or. n_families < 0 .or. n_samples <= 0 .or. d <= 0) then
            ierr = ERR_INVALID_INPUT
            write(*,*) 'Error: Invalid dimensions n_genes=', n_genes, ' n_families=', n_families, &
                       ' n_samples=', n_samples, ' d=', d
            return
        end if
        
        ! Check gene_ids array
        if (size(gene_ids) /= n_genes) then
            ierr = ERR_SIZE_MISMATCH
            write(*,*) 'Error: gene_ids size mismatch. Expected:', n_genes, ' Actual:', size(gene_ids)
            return
        end if
        
        ! Check gene_family_ids array
        if (size(gene_family_ids) /= n_families) then
            ierr = ERR_SIZE_MISMATCH
            write(*,*) 'Error: gene_family_ids size mismatch. Expected:', n_families, ' Actual:', size(gene_family_ids)
            return
        end if
        
        ! Check gene_to_fam array
        if (size(gene_to_fam) /= n_genes) then
            ierr = ERR_SIZE_MISMATCH
            write(*,*) 'Error: gene_to_fam size mismatch. Expected:', n_genes, ' Actual:', size(gene_to_fam)
            return
        end if
        
        ! Check expression_vectors array
        if (size(expression_vectors, 1) /= n_samples .or. size(expression_vectors, 2) /= n_genes) then
            ierr = ERR_SIZE_MISMATCH
            write(*,*) 'Error: expression_vectors size mismatch. Expected: (', n_samples, ',', n_genes, &
                       ') Actual: (', size(expression_vectors, 1), ',', size(expression_vectors, 2), ')'
            return
        end if
        
        ! Check family_centroids array if families exist
        if (n_families > 0) then
            if (size(family_centroids, 1) /= d .or. size(family_centroids, 2) /= n_families) then
                ierr = ERR_SIZE_MISMATCH
                write(*,*) 'Error: family_centroids size mismatch. Expected: (', d, ',', n_families, &
                           ') Actual: (', size(family_centroids, 1), ',', size(family_centroids, 2), ')'
                return
            end if
        end if
        
        ! Check shift_vectors array
        if (size(shift_vectors, 1) /= 2*d .or. size(shift_vectors, 2) /= n_genes) then
            ierr = ERR_SIZE_MISMATCH
            write(*,*) 'Error: shift_vectors size mismatch. Expected: (', 2*d, ',', n_genes, &
                       ') Actual: (', size(shift_vectors, 1), ',', size(shift_vectors, 2), ')'
            return
        end if
        
        ! Check for empty strings in gene_ids
        call validate_empty_strings(gene_ids, "gene_ids", ierr)
        if (ierr /= ERR_OK) return
        
        ! Check for empty strings in gene_family_ids if families exist
        if (n_families > 0) then
            call validate_empty_strings(gene_family_ids, "gene_family_ids", ierr)
            if (ierr /= ERR_OK) return
        end if
        
    end subroutine validate_data_structure

    subroutine validate_gene_to_family_mapping(gene_to_fam, n_families, ierr)
        integer(int32), intent(in) :: gene_to_fam(:)
        integer(int32), intent(in) :: n_families
        integer(int32), intent(out) :: ierr
        
        integer :: i, invalid_count
        
        ierr = ERR_OK
        invalid_count = 0
        
        if (n_families == 0) then
            ! No families defined, mapping should be all zeros
            do i = 1, size(gene_to_fam)
                if (gene_to_fam(i) /= 0) then
                    invalid_count = invalid_count + 1
                    if (invalid_count <= 10) then
                        write(*,*) 'Error: gene_to_fam(', i, ') = ', gene_to_fam(i), ' but should be 0 (no families)'
                    end if
                end if
            end do
        else
            ! Check that all gene_to_fam values are valid family indices
            do i = 1, size(gene_to_fam)
                if (gene_to_fam(i) < 0 .or. gene_to_fam(i) > n_families) then
                    invalid_count = invalid_count + 1
                    if (invalid_count <= 10) then
                        write(*,*) 'Error: gene_to_fam(', i, ') = ', gene_to_fam(i), &
                                   ' but valid range is 0 to ', n_families
                    end if
                end if
            end do
        end if
        
        if (invalid_count > 0) then
            ierr = ERR_INVALID_INPUT
            write(*,*) 'Error: Found ', invalid_count, ' invalid gene-to-family mappings'
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
                write(*,*) 'Error: Negative values found in expression data'
                return
            end if
        end if
        
        ! Check for missing values
        if (any(expression_vectors == MISSING_VALUE)) then
            write(*,*) 'Warning: Missing values (', MISSING_VALUE, ') found in expression data'
        end if
        
    end subroutine validate_expression_data

    subroutine validate_family_centroids(family_centroids, ierr)
        real(real64), intent(in) :: family_centroids(:,:)
        integer(int32), intent(out) :: ierr
        
        ierr = ERR_OK
        
        ! Check for NaN and Inf values
        call check_for_nan_inf(family_centroids, "family_centroids", ierr)
        if (ierr /= ERR_OK) return
        
        ! Check for missing values
        if (any(family_centroids == MISSING_VALUE)) then
            write(*,*) 'Warning: Missing values (', MISSING_VALUE, ') found in family centroids'
        end if
        
    end subroutine validate_family_centroids

    subroutine validate_shift_vectors(shift_vectors, expression_vectors, family_centroids, &
                                    gene_to_fam, d, ierr)
        real(real64), intent(in) :: shift_vectors(:,:)
        real(real64), intent(in) :: expression_vectors(:,:)
        real(real64), intent(in) :: family_centroids(:,:)
        integer(int32), intent(in) :: gene_to_fam(:)
        integer(int32), intent(in) :: d
        integer(int32), intent(out) :: ierr
        
        integer :: i, j, k, fam_idx, n_genes, n_samples, error_count
        real(real64) :: expected_shift
        
        ierr = ERR_OK
        n_genes = size(expression_vectors, 2)
        n_samples = size(expression_vectors, 1)
        error_count = 0
        
        ! Check for NaN and Inf values
        call check_for_nan_inf(shift_vectors, "shift_vectors", ierr)
        if (ierr /= ERR_OK) return
        
        ! Check for missing values
        if (any(shift_vectors == MISSING_VALUE)) then
            write(*,*) 'Warning: Missing values (', MISSING_VALUE, ') found in shift vectors'
        end if
        
        ! Verify shift vectors structure: first d rows should be centroids, next d rows should be shifts
        do i = 1, n_genes
            fam_idx = gene_to_fam(i)
            
            if (fam_idx == 0) then
                ! Gene has no family assignment, skip validation
                cycle
            end if
            
            ! Check that centroid part matches the family centroid
            do j = 1, d
                if (abs(shift_vectors(j, i) - family_centroids(j, fam_idx)) > FLOAT_TOLERANCE) then
                    error_count = error_count + 1
                    if (error_count <= 10) then
                        write(*,*) 'Error: Centroid mismatch for gene ', i, &
                                   ' dimension ', j, ' expected ', family_centroids(j, fam_idx), &
                                   ' got ', shift_vectors(j, i)
                    end if
                end if
            end do
            
            ! Check that shift part matches expression - centroid
            do j = 1, d
                expected_shift = 0.0_real64
                ! Calculate average shift across all samples
                do k = 1, n_samples
                    expected_shift = expected_shift + (expression_vectors(k, i) - family_centroids(j, fam_idx))
                end do
                expected_shift = expected_shift / n_samples
                
                if (abs(shift_vectors(d+j, i) - expected_shift) > FLOAT_TOLERANCE) then
                    error_count = error_count + 1
                    if (error_count <= 10) then
                        write(*,*) 'Error: Shift mismatch for gene ', i, &
                                   ' dimension ', j, ' expected ', expected_shift, &
                                   ' got ', shift_vectors(d+j, i)
                    end if
                end if
            end do
        end do
        
        if (error_count > 0) then
            ierr = ERR_INVALID_INPUT
            write(*,*) 'Error: Found ', error_count, ' inconsistencies in shift vectors'
        end if
        
    end subroutine validate_shift_vectors

    subroutine check_for_nan_inf(array, array_name, ierr)
        real(real64), intent(in) :: array(:,:)
        character(len=*), intent(in) :: array_name
        integer(int32), intent(out) :: ierr
        
        integer :: i, j, nan_inf_count
        
        ierr = ERR_OK
        nan_inf_count = 0
        
        do j = 1, size(array, 2)
            do i = 1, size(array, 1)
                if (is_nan(array(i, j)) .or. is_inf(array(i, j))) then
                    nan_inf_count = nan_inf_count + 1
                    if (nan_inf_count <= 10) then
                        write(*,*) 'Error: NaN/Inf found in ', array_name, &
                                   ' at position (', i, ',', j, ')'
                    end if
                end if
            end do
        end do
        
        if (nan_inf_count > 0) then
            ierr = ERR_INVALID_INPUT
            write(*,*) 'Error: Found ', nan_inf_count, ' NaN/Inf values in ', array_name
        end if
        
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
        
        integer :: i, j, duplicate_count
        
        ierr = ERR_OK
        duplicate_count = 0
        
        do i = 1, size(gene_ids) - 1
            do j = i + 1, size(gene_ids)
                if (trim(gene_ids(i)) == trim(gene_ids(j))) then
                    duplicate_count = duplicate_count + 1
                    if (duplicate_count <= 10) then
                        write(*,*) 'Error: Duplicate gene ID found: "', trim(gene_ids(i)), '"'
                    end if
                end if
            end do
        end do
        
        if (duplicate_count > 0) then
            ierr = ERR_INVALID_INPUT
            write(*,*) 'Error: Found ', duplicate_count, ' duplicate gene IDs'
        end if
        
    end subroutine validate_gene_ids_uniqueness

    subroutine validate_family_ids_uniqueness(gene_family_ids, ierr)
        character(len=*), intent(in) :: gene_family_ids(:)
        integer(int32), intent(out) :: ierr
        
        integer :: i, j, duplicate_count
        
        ierr = ERR_OK
        duplicate_count = 0
        
        do i = 1, size(gene_family_ids) - 1
            do j = i + 1, size(gene_family_ids)
                if (trim(gene_family_ids(i)) == trim(gene_family_ids(j))) then
                    duplicate_count = duplicate_count + 1
                    if (duplicate_count <= 10) then
                        write(*,*) 'Error: Duplicate family ID found: "', trim(gene_family_ids(i)), '"'
                    end if
                end if
            end do
        end do
        
        if (duplicate_count > 0) then
            ierr = ERR_INVALID_INPUT
            write(*,*) 'Error: Found ', duplicate_count, ' duplicate family IDs'
        end if
        
    end subroutine validate_family_ids_uniqueness

    subroutine validate_empty_strings(string_array, array_name, ierr)
        character(len=*), intent(in) :: string_array(:)
        character(len=*), intent(in) :: array_name
        integer(int32), intent(out) :: ierr
        
        integer :: i, empty_count
        
        ierr = ERR_OK
        empty_count = 0
        
        do i = 1, size(string_array)
            if (len_trim(string_array(i)) == 0) then
                empty_count = empty_count + 1
                if (empty_count <= 10) then
                    write(*,*) 'Error: Empty string found in ', array_name, ' at index ', i
                end if
            end if
        end do
        
        if (empty_count > 0) then
            ierr = ERR_INVALID_INPUT
            write(*,*) 'Error: Found ', empty_count, ' empty strings in ', array_name
        end if
        
    end subroutine validate_empty_strings

    ! New: Comprehensive validation routine
    subroutine validate_all_data(n_genes, n_families, n_samples, d, gene_ids, gene_family_ids, &
                               gene_to_fam, expression_vectors, family_centroids, &
                               shift_vectors, ierr, check_uniqueness, check_shift_consistency)
        integer(int32), intent(in) :: n_genes, n_families, n_samples, d
        character(len=*), intent(in) :: gene_ids(:)
        character(len=*), intent(in) :: gene_family_ids(:)
        integer(int32), intent(in) :: gene_to_fam(:)
        real(real64), intent(in) :: expression_vectors(:,:)
        real(real64), intent(in) :: family_centroids(:,:)
        real(real64), intent(in) :: shift_vectors(:,:)
        integer(int32), intent(out) :: ierr
        logical, intent(in), optional :: check_uniqueness, check_shift_consistency
        
        logical :: do_check_uniqueness, do_check_shift_consistency
        
        ! Set defaults for optional parameters
        do_check_uniqueness = .true.
        if (present(check_uniqueness)) do_check_uniqueness = check_uniqueness
        
        do_check_shift_consistency = .true.
        if (present(check_shift_consistency)) do_check_shift_consistency = check_shift_consistency
        
        ierr = ERR_OK
        
        write(*,*) 'Starting comprehensive data validation...'
        
        ! 1. Check basic structure
        call validate_data_structure(n_genes, n_families, n_samples, d, gene_ids, gene_family_ids, &
                                   gene_to_fam, expression_vectors, family_centroids, &
                                   shift_vectors, ierr)
        if (ierr /= ERR_OK) return
        
        ! 2. Check gene-to-family mapping
        call validate_gene_to_family_mapping(gene_to_fam, n_families, ierr)
        if (ierr /= ERR_OK) return
        
        ! 3. Check expression data
        call validate_expression_data(expression_vectors, .true., ierr)  ! Check for non-negative
        if (ierr /= ERR_OK) return
        
        ! 4. Check family centroids
        if (n_families > 0) then
            call validate_family_centroids(family_centroids, ierr)
            if (ierr /= ERR_OK) return
        end if
        
        ! 5. Check shift vectors
        call validate_shift_vectors(shift_vectors, expression_vectors, family_centroids, &
                                  gene_to_fam, d, ierr)
        if (ierr /= ERR_OK .and. do_check_shift_consistency) return
        
        ! 6. Check uniqueness (optional, can be slow for large datasets)
        if (do_check_uniqueness) then
            call validate_gene_ids_uniqueness(gene_ids, ierr)
            if (ierr /= ERR_OK) return
            
            if (n_families > 0) then
                call validate_family_ids_uniqueness(gene_family_ids, ierr)
                if (ierr /= ERR_OK) return
            end if
        end if
        
        write(*,*) 'All data validation checks passed!'
        
    end subroutine validate_all_data
    
end module tox_data_validation