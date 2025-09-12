module tox_data_validation
    use iso_fortran_env, only: real64, int32
    use tox_errors, only: set_ok, is_ok, set_err_once, ERR_INVALID_INPUT, ERR_SIZE_MISMATCH
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
    PUBLIC :: validate_all_data
    
    ! Parameters for validation tolerances
    real(real64), parameter :: FLOAT_TOLERANCE = 1.0e-10_real64
    
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
        
        call set_ok(ierr)
        
        ! Check basic dimensions
        if (n_genes <= 0 .or. n_families < 0 .or. n_samples <= 0 .or. d <= 0) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if
        
        ! Check gene_ids array
        if (size(gene_ids) /= n_genes) then
            call set_err_once(ierr, ERR_SIZE_MISMATCH)
            write(*,*) 'Error: gene_ids size mismatch. Expected:', n_genes, ' Actual:', size(gene_ids)
            return
        end if
        
        ! Check gene_family_ids array
        if (size(gene_family_ids) /= n_families) then
            call set_err_once(ierr, ERR_SIZE_MISMATCH)
            write(*,*) 'Error: gene_family_ids size mismatch. Expected:', n_families, ' Actual:', size(gene_family_ids)
            return
        end if
        
        ! Check gene_to_fam array
        if (size(gene_to_fam) /= n_genes) then
            call set_err_once(ierr, ERR_SIZE_MISMATCH)
            write(*,*) 'Error: gene_to_fam size mismatch. Expected:', n_genes, ' Actual:', size(gene_to_fam)
            return
        end if
        
        ! Check expression_vectors array
        if (size(expression_vectors, 1) /= n_samples .or. size(expression_vectors, 2) /= n_genes) then
            call set_err_once(ierr, ERR_SIZE_MISMATCH)
            write(*,*) 'Error: expression_vectors size mismatch. Expected: (', n_samples, ',', n_genes, &
                       ') Actual: (', size(expression_vectors, 1), ',', size(expression_vectors, 2), ')'
            return
        end if
        
        ! Check family_centroids array if families exist
        if (n_families > 0) then
            if (size(family_centroids, 1) /= d .or. size(family_centroids, 2) /= n_families) then
                call set_err_once(ierr, ERR_SIZE_MISMATCH)
                write(*,*) 'Error: family_centroids size mismatch. Expected: (', d, ',', n_families, &
                           ') Actual: (', size(family_centroids, 1), ',', size(family_centroids, 2), ')'
                return
            end if
        else 
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if
        
        ! Check shift_vectors array
        if (size(shift_vectors, 1) /= 2*d .or. size(shift_vectors, 2) /= n_genes) then
            call set_err_once(ierr, ERR_SIZE_MISMATCH)
            write(*,*) 'Error: shift_vectors size mismatch. Expected: (', 2*d, ',', n_genes, &
                       ') Actual: (', size(shift_vectors, 1), ',', size(shift_vectors, 2), ')'
            return
        end if
        
        ! Check for empty strings in gene_ids
        call validate_empty_strings(gene_ids, "gene_ids", ierr)
        if (.not. is_ok(ierr)) return
        
        ! Check for empty strings in gene_family_ids if families exist
        if (n_families > 0) then
            call validate_empty_strings(gene_family_ids, "gene_family_ids", ierr)
            if (.not. is_ok(ierr)) return
        end if
        
    end subroutine validate_data_structure

    subroutine validate_gene_to_family_mapping(gene_to_fam, n_families, ierr)
        integer(int32), intent(in) :: gene_to_fam(:)
        integer(int32), intent(in) :: n_families
        integer(int32), intent(out) :: ierr
        
        integer :: i, invalid_count
        
        call set_ok(ierr)
        invalid_count = 0
        
        if (n_families == 0) then
            ! No families defined, mapping should be all zeros
            do i = 1, size(gene_to_fam)
                if (gene_to_fam(i) /= 0) then
                    invalid_count = invalid_count + 1
                end if
            end do
            if (invalid_count > 0) then
                call set_err_once(ierr, ERR_INVALID_INPUT)
                write(*,*) 'Error: gene_to_fam should be all zeros when no families are defined'
                RETURN
            end if
        else 
            ! Check that all gene_to_fam values are valid family indices
            do i = 1, size(gene_to_fam)
                if (gene_to_fam(i) < 0 .or. gene_to_fam(i) > n_families) then
                    invalid_count = invalid_count + 1
                    if(invalid_count > 0) then
                        write(*,*) 'Error: gene_to_fam(', i, ') = ', gene_to_fam(i), &
                                   ' but valid range is 0 to ', n_families
                        call set_err_once(ierr, ERR_INVALID_INPUT)
                        return
                    end if
                end if
            end do
        end if
        
    end subroutine validate_gene_to_family_mapping

    subroutine validate_expression_data(expression_vectors, check_non_negative, ierr)
        real(real64), intent(in) :: expression_vectors(:,:)
        logical, intent(in) :: check_non_negative
        integer(int32), intent(out) :: ierr
        
        call set_ok(ierr)
        
        ! Check for NaN and Inf values
        call check_for_nan_inf(expression_vectors, "expression_vectors", ierr)
        if (.not. is_ok(ierr)) return
        
        ! Check for negative values if requested
        if (check_non_negative) then
            if (any(expression_vectors < 0.0_real64)) then
                call set_err_once(ierr, ERR_INVALID_INPUT)
                write(*,*) 'Error: Negative values found in expression data'
                return
            end if
        end if
        
    end subroutine validate_expression_data

    subroutine validate_family_centroids(family_centroids, ierr)
        real(real64), intent(in) :: family_centroids(:,:)
        integer(int32), intent(out) :: ierr
        
        call set_ok(ierr)
        
        ! Check for NaN and Inf values
        call check_for_nan_inf(family_centroids, "family_centroids", ierr)
        if (.not. is_ok(ierr)) return

    end subroutine validate_family_centroids

    subroutine validate_shift_vectors(shift_vectors, expression_vectors, family_centroids, &
                                        gene_to_fam, d, ierr)
        use iso_fortran_env, only: real64, int32
        use tox_errors
        implicit none
        
        real(real64), intent(in) :: shift_vectors(:,:)
        real(real64), intent(in) :: expression_vectors(:,:)
        real(real64), intent(in) :: family_centroids(:,:)
        integer(int32), intent(in) :: gene_to_fam(:)
        integer(int32), intent(in) :: d
        integer(int32), intent(out) :: ierr
        real(real64) :: expected_shift

        
        integer(int32) :: i, j, fam_idx, n_genes, error_count
        
        ierr = ERR_OK
        n_genes = size(expression_vectors, 2)
        error_count = 0
        
        ! Check for NaN and Inf values
        call check_for_nan_inf(shift_vectors, "shift_vectors", ierr)
        if (.not. is_ok(ierr)) return
        
        ! Verify shift vectors structure: first d rows should be centroids, next d rows should be shifts
        do i = 1, n_genes
            fam_idx = gene_to_fam(i)
            
            if (fam_idx == 0) then
                call set_err_once(ierr, ERR_INVALID_INPUT)
                RETURN
            end if
            
            ! Check that centroid part (first d rows) matches the family centroid
            do j = 1, d
                if (abs(shift_vectors(j, i) - family_centroids(j, fam_idx)) > FLOAT_TOLERANCE) then
                    error_count = error_count + 1
                    if (error_count <= 10) then
                        call set_err_once(ierr, ERR_INVALID_INPUT)
                        write(*,*) 'Error: Centroid mismatch for gene ', i, &
                                ' dimension ', j, ' expected ', family_centroids(j, fam_idx), &
                                ' got ', shift_vectors(j, i)
                    end if
                end if
            end do
            if(.not. is_ok(ierr)) return
            
            ! Check that shift part (rows d+1 to 2d) matches expression - centroid
            do j = 1, d
                ! The expected shift is simply the difference between the expression vector and the centroid.
                ! This aligns with the compute_shift_vector_field subroutine.
                expected_shift = expression_vectors(j, i) - family_centroids(j, fam_idx)
                
                if (abs(shift_vectors(d+j, i) - expected_shift) > FLOAT_TOLERANCE) then
                    error_count = error_count + 1
                    if (error_count <= 10) then
                        call set_err_once(ierr, ERR_INVALID_INPUT)
                        write(*,*) 'Error: Shift mismatch for gene ', i, &
                                ' dimension ', j, ' expected ', expected_shift, &
                                ' got ', shift_vectors(d+j, i)
                    end if
                end if
            end do
            if(.not. is_ok(ierr)) return
        end do
        
    end subroutine validate_shift_vectors

    subroutine check_for_nan_inf(array, array_name, ierr)
        real(real64), intent(in) :: array(:,:)
        character(len=*), intent(in) :: array_name
        integer(int32), intent(out) :: ierr
        
        integer :: i, j, nan_inf_count
        
        call set_ok(ierr)
        nan_inf_count = 0
        
        do j = 1, size(array, 2)
            do i = 1, size(array, 1)
                if (is_nan(array(i, j)) .or. is_inf(array(i, j))) then
                    nan_inf_count = nan_inf_count + 1
                    if (nan_inf_count <= 10) then
                        call set_err_once(ierr, ERR_INVALID_INPUT)
                        write(*,*) 'Error: NaN/Inf found in ', array_name, &
                                   ' at position (', i, ',', j, ')'
                    end if
                end if
            end do
        end do
        if(.not. is_ok(ierr)) return
        
    contains

        !Possibly replace by existing asserts

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
        
        call set_ok(ierr)
        duplicate_count = 0
        
        do i = 1, size(gene_ids) - 1
            do j = i + 1, size(gene_ids)
                if (trim(gene_ids(i)) == trim(gene_ids(j))) then
                    duplicate_count = duplicate_count + 1
                    if (duplicate_count <= 10) then
                        call set_err_once(ierr, ERR_INVALID_INPUT)
                        write(*,*) 'Error: Duplicate gene ID found: "', trim(gene_ids(i)), '"'
                    end if
                end if
            end do
        end do
        
    end subroutine validate_gene_ids_uniqueness

    subroutine validate_family_ids_uniqueness(gene_family_ids, ierr)
        character(len=*), intent(in) :: gene_family_ids(:)
        integer(int32), intent(out) :: ierr
        
        integer :: i, j, duplicate_count
        
        call set_ok(ierr)
        duplicate_count = 0
        
        do i = 1, size(gene_family_ids) - 1
            do j = i + 1, size(gene_family_ids)
                if (trim(gene_family_ids(i)) == trim(gene_family_ids(j))) then
                    duplicate_count = duplicate_count + 1
                    if (duplicate_count <= 10) then
                        call set_err_once(ierr, ERR_INVALID_INPUT)
                        write(*,*) 'Error: Duplicate family ID found: "', gene_family_ids(i), '"'
                    end if
                end if
            end do
        end do
        
    end subroutine validate_family_ids_uniqueness

    subroutine validate_empty_strings(string_array, array_name, ierr)
        character(len=*), intent(in) :: string_array(:)
        character(len=*), intent(in) :: array_name
        integer(int32), intent(out) :: ierr
        
        integer :: i, empty_count
        
        call set_ok(ierr)
        empty_count = 0
        
        do i = 1, size(string_array)
            if (len_trim(string_array(i)) == 0) then
                empty_count = empty_count + 1
                if (empty_count <= 10) then
                    call set_err_once(ierr, ERR_INVALID_INPUT)
                    write(*,*) 'Error: Empty string found in ', array_name, ' at index ', i
                end if
            end if
        end do
        
    end subroutine validate_empty_strings

    ! Comprehensive validation routine, combining all checks
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
        
        call set_ok(ierr)
        
        write(*,*) 'Starting comprehensive data validation, this might take a while...'
        
        ! 1. Check basic structure
        call validate_data_structure(n_genes, n_families, n_samples, d, gene_ids, gene_family_ids, &
                                   gene_to_fam, expression_vectors, family_centroids, &
                                   shift_vectors, ierr)
        if (.not. is_ok(ierr)) return
        
        ! 2. Check gene-to-family mapping
        call validate_gene_to_family_mapping(gene_to_fam, n_families, ierr)
        if (.not. is_ok(ierr)) return
        
        ! 3. Check expression data
        call validate_expression_data(expression_vectors, .true., ierr)  ! Check for non-negative
        if (.not. is_ok(ierr)) return
        
        ! 4. Check family centroids
        if (n_families > 0) then
            call validate_family_centroids(family_centroids, ierr)
            if (.not. is_ok(ierr)) return
        else 
            write(*,*) 'Warning: No families defined, skipping family centroid checks.'
        end if
        
        ! 5. Check shift vectors
        if (do_check_shift_consistency) then 
            call validate_shift_vectors(shift_vectors, expression_vectors, family_centroids, &
                                  gene_to_fam, d, ierr)
            if (.not. is_ok(ierr)) return
        end if

        ! 6. Check uniqueness (optional, can be slow for large datasets)
        if (do_check_uniqueness) then
            call validate_gene_ids_uniqueness(gene_ids, ierr)
            if (.not. is_ok(ierr)) return
            
            if (n_families > 0) then
                call validate_family_ids_uniqueness(gene_family_ids, ierr)
                if (.not. is_ok(ierr)) return
            end if
        end if
        
        write(*,*) 'All data validation checks passed!'
        
    end subroutine validate_all_data
    
end module tox_data_validation

subroutine validate_gene_to_family_mapping_R(gene_to_fam, n_genes, n_families, ierr)
    use tox_data_validation, only: validate_gene_to_family_mapping
    use tox_errors, only: set_ok
    use iso_fortran_env, only: int32
    integer(int32), intent(in) :: n_genes, n_families
    integer(int32), intent(in) :: gene_to_fam(n_genes)

    integer(int32), intent(out) :: ierr
    call set_ok(ierr)
    call validate_gene_to_family_mapping(gene_to_fam, n_families, ierr)
end subroutine validate_gene_to_family_mapping_R

subroutine validate_expression_data_R(expression_vectors, n_genes, n_samples, check_non_negative, ierr)
    use tox_data_validation, only: validate_expression_data
    use tox_errors, only: set_ok
    use iso_fortran_env, only: int32, real64
    logical, intent(in) :: check_non_negative
    real(real64), intent(in) :: expression_vectors(n_genes, n_samples)
    integer(int32), intent(in) :: n_genes, n_samples
    integer(int32), intent(out) :: ierr
    call set_ok(ierr)
    call validate_expression_data(expression_vectors, check_non_negative, ierr)
end subroutine validate_expression_data_R

subroutine validate_family_centroids_R(family_centroids, n_families, d, ierr)
    use tox_data_validation, only : validate_family_centroids
    use tox_errors, only: set_ok
    use iso_fortran_env, only: int32, real64
    real(real64), intent(in) :: family_centroids(d, n_families)
    integer(int32), intent(in) :: n_families, d
    integer(int32), intent(out) :: ierr
    call set_ok(ierr)
    call validate_family_centroids(family_centroids, ierr)
end subroutine validate_family_centroids_R

subroutine validate_shift_vectors_R(shift_vectors, expression_vectors, family_centroids, gene_to_fam, d, n_genes, n_samples, n_families, ierr)
    use tox_data_validation, only: validate_shift_vectors
    use tox_errors, only: set_ok
    use iso_fortran_env, only: int32, real64
    real(real64), intent(in) :: shift_vectors(2*d, n_genes)
    real(real64), intent(in) :: expression_vectors(n_samples, n_genes)
    real(real64), intent(in) :: family_centroids(d, n_families)
    integer(int32), intent(in) :: gene_to_fam(n_genes)
    integer(int32), intent(in) :: d, n_genes, n_samples, n_families
    integer(int32), intent(out) :: ierr
    call set_ok(ierr)
    
    call validate_shift_vectors(shift_vectors, expression_vectors, family_centroids, gene_to_fam, d, ierr)
end subroutine validate_shift_vectors_R

subroutine validate_gene_ids_uniqueness_R(gene_ids_ascii, gene_ids_len, n_genes, ierr)
    use tox_data_validation, only: validate_gene_ids_uniqueness
    use tox_errors, only: set_ok, is_ok, set_err_once, ERR_ALLOC_FAIL
    use array_utils, only: ascii_to_string_padded
    use iso_fortran_env, only: int32
    integer(int32), intent(in) :: gene_ids_ascii(gene_ids_len, n_genes)
    integer(int32), intent(in) :: gene_ids_len, n_genes
    integer(int32), intent(out) :: ierr
    character(len=:), allocatable :: gene_ids(:)
    character(len=:), allocatable :: temp_str
    integer(int32) :: i, ios

    call set_ok(ierr)
    call set_ok(ios)

    allocate(character(len=gene_ids_len) :: gene_ids(n_genes), stat=ios)
    if(.not. is_ok(ios)) then
        call set_err_once(ierr, ERR_ALLOC_FAIL)
        return
    end if

    do i = 1, n_genes
        call ascii_to_string_padded(gene_ids_ascii(:,i), gene_ids_len, temp_str)
        gene_ids(i) = temp_str
    end do

    call validate_gene_ids_uniqueness(gene_ids, ierr)
end subroutine validate_gene_ids_uniqueness_R

subroutine validate_family_ids_uniqueness_R(family_ids_ascii, fam_len, n_families, ierr)
    use tox_data_validation, only: validate_family_ids_uniqueness
    use tox_errors, only: set_ok, is_ok, set_err_once, ERR_ALLOC_FAIL
    use array_utils, only: ascii_to_string_padded
    use iso_fortran_env, only: int32
    integer(int32), intent(in) :: family_ids_ascii(fam_len, n_families)
    integer(int32), intent(in) :: fam_len, n_families
    integer(int32), intent(out) :: ierr
    character(len=:), allocatable :: family_ids(:)
    character(len=:), ALLOCATABLE :: temp_str
    integer(int32) :: i, ios

    call set_ok(ierr)
    call set_ok(ios)

    allocate(character(len=fam_len) :: family_ids(n_families), stat=ios)
    if(.not. is_ok(ios)) then
        call set_err_once(ierr, ERR_ALLOC_FAIL)
        return
    end if

    do i = 1, n_families
        call ascii_to_string_padded(family_ids_ascii(:,i), fam_len, temp_str)
        family_ids(i) = temp_str
    end do

    call validate_family_ids_uniqueness(family_ids, ierr)
end subroutine validate_family_ids_uniqueness_R

subroutine validate_empty_strings_R(strings_ascii, str_len, n, ierr)
    use tox_data_validation, only: validate_empty_strings
    use tox_errors, only: set_ok
    use array_utils, only: ascii_to_string_padded
    use iso_fortran_env, only: int32
    integer(int32), intent(in) :: strings_ascii(str_len, n)
    integer(int32), intent(in) :: str_len, n
    integer(int32), intent(out) :: ierr
    character(len=:), allocatable :: strings(:)
    character(len=:), allocatable :: temp_str
    integer(int32) :: i

    call set_ok(ierr)

    allocate(character(len=str_len) :: strings(n))
    do i = 1, n
        call ascii_to_string_padded(strings_ascii(:,i), str_len, temp_str)
        strings(i) = temp_str
    end do

    call validate_empty_strings(strings, "strings", ierr)
end subroutine validate_empty_strings_R

subroutine validate_data_structure_R(n_genes, n_families, n_samples, d, &
                                     gene_ids_ascii, gene_ids_len, &
                                     gene_family_ids_ascii, fam_len, &
                                     gene_to_fam, expression_vectors, family_centroids, &
                                     shift_vectors, ierr)
    use tox_data_validation, only: validate_data_structure
    use tox_errors, only: set_ok, is_ok, set_err_once, ERR_ALLOC_FAIL
    use array_utils, only: ascii_to_string_padded
    use iso_fortran_env, only: int32, real64
    integer(int32), intent(in) :: n_genes, n_families, n_samples, d
    integer(int32), intent(in) :: gene_ids_ascii(gene_ids_len, n_genes)
    integer(int32), intent(in) :: gene_ids_len
    integer(int32), intent(in) :: gene_family_ids_ascii(fam_len, n_genes)
    integer(int32), intent(in) :: fam_len
    integer(int32), intent(in) :: gene_to_fam(n_genes)
    real(real64), intent(in) :: expression_vectors(n_genes, n_samples)
    real(real64), intent(in) :: family_centroids(n_families, d)
    real(real64), intent(in) :: shift_vectors(n_families, d)
    integer(int32), intent(out) :: ierr
    character(len=:), allocatable :: temp_str
    character(len=:), allocatable :: gene_ids(:)
    character(len=:), allocatable :: gene_family_ids(:)
    integer(int32) :: i, ios

    call set_ok(ierr)
    call set_ok(ios)

    allocate(character(len=gene_ids_len) :: gene_ids(n_genes), stat=ios)
    allocate(character(len=fam_len) :: gene_family_ids(n_genes), stat=ios)
    if(.not. is_ok(ios)) then
        call set_err_once(ierr, ERR_ALLOC_FAIL)
        return
    end if
    do i = 1, n_genes
        call ascii_to_string_padded(gene_ids_ascii(:,i), gene_ids_len, temp_str)
        gene_ids(i) = temp_str
        call ascii_to_string_padded(gene_family_ids_ascii(:,i), fam_len, temp_str)
        gene_family_ids(i) = temp_str
    end do

    call validate_data_structure(n_genes, n_families, n_samples, d, gene_ids, gene_family_ids, &
                                 gene_to_fam, expression_vectors, family_centroids, shift_vectors, ierr)
end subroutine validate_data_structure_R

subroutine validate_all_data_R(n_genes, n_families, n_samples, d, &
                               gene_ids_ascii, gene_len, &
                               gene_family_ids_ascii, fam_len, &
                               gene_to_fam, expression_vectors, family_centroids, &
                               shift_vectors, ierr)
    use tox_data_validation, only: validate_all_data
    use tox_errors, only: set_ok, is_ok, set_err_once, ERR_ALLOC_FAIL
    use array_utils, only: ascii_to_string_padded
    use iso_fortran_env, only: int32, real64
    integer(int32), intent(in) :: n_genes, n_families, n_samples, d
    integer(int32), intent(in) :: gene_ids_ascii(gene_len, n_genes)
    integer(int32), intent(in) :: gene_len
    integer(int32), intent(in) :: gene_family_ids_ascii(fam_len, n_families)
    integer(int32), intent(in) :: fam_len
    integer(int32), intent(in) :: gene_to_fam(n_genes)
    real(real64), intent(in) :: expression_vectors(n_samples, n_genes)
    real(real64), intent(in) :: family_centroids(d, n_families)
    real(real64), intent(in) :: shift_vectors(2*d, n_genes)
    integer(int32), intent(out) :: ierr
    character(len=:), allocatable :: temp_str
    character(len=:), allocatable :: gene_ids(:)
    character(len=:), allocatable :: gene_family_ids(:)
    integer(int32) :: i, ios

    call set_ok(ierr)
    call set_ok(ios)

    allocate(character(len=gene_len) :: gene_ids(n_genes), stat=ios)
    allocate(character(len=fam_len) :: gene_family_ids(n_families), stat=ios)
    do i = 1, n_genes
        call ascii_to_string_padded(gene_ids_ascii(:,i), gene_len, temp_str)
        gene_ids(i) = temp_str
    end do

    do i=1, n_families
        call ascii_to_string_padded(gene_family_ids_ascii(:,i), fam_len, temp_str)
        gene_family_ids(i) = temp_str
    end do

    call validate_all_data(n_genes, n_families, n_samples, d, gene_ids, gene_family_ids, &
                           gene_to_fam, expression_vectors, family_centroids, shift_vectors, ierr)
end subroutine validate_all_data_R