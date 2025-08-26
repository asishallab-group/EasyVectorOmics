module tox_data_accessors
    use iso_fortran_env, only: real64, int32
    use tox_data_type, only: tox_data_t
    use tox_errors, only: ERR_OK, ERR_INVALID_INPUT, ERR_SIZE_MISMATCH
implicit none

contains

    subroutine get_all_expression_vectors(this, out_vectors, ierr)
        class(tox_data_t), intent(in) :: this
        real(real64), allocatable, intent(out) :: out_vectors(:,:)
        integer(int32), intent(out), optional :: ierr
        
        if (present(ierr)) call set_ok(ierr)
        
        if (.not. allocated(this%expression_vectors)) then
            if (present(ierr)) ierr = ERR_INVALID_INPUT
            allocate(out_vectors(0,0))
            return
        end if
        
        allocate(out_vectors(size(this%expression_vectors, 1), size(this%expression_vectors, 2)))
        out_vectors = this%expression_vectors
    end subroutine get_all_expression_vectors

    subroutine get_single_expression_vector(this, gene_idx, out_vector, ierr)
        class(tox_data_t), intent(in) :: this
        integer(int32), intent(in) :: gene_idx
        real(real64), allocatable, intent(out) :: out_vector(:)
        integer(int32), intent(out), optional :: ierr
        
        if (present(ierr)) call set_ok(ierr)
        
        if (.not. allocated(this%expression_vectors) .or. gene_idx < 1 .or. gene_idx > this%n_genes) then
            if (present(ierr)) ierr = ERR_INVALID_INPUT
            allocate(out_vector(0))
            return
        end if
        
        allocate(out_vector(size(this%expression_vectors, 1)))
        out_vector = this%expression_vectors(:, gene_idx)
    end subroutine get_single_expression_vector

    subroutine get_gene_ids(this, out_ids, ierr)
        class(tox_data_t), intent(in) :: this
        character(len=:), allocatable, intent(out) :: out_ids(:)
        integer(int32), intent(out), optional :: ierr
        integer :: i
        
        if (present(ierr)) call set_ok(ierr)
        
        if (.not. allocated(this%gene_ids)) then
            if (present(ierr)) ierr = ERR_INVALID_INPUT
            allocate(character(len=1) :: out_ids(0))
            return
        end if
        
        allocate(character(len=len_trim(this%gene_ids(1))) :: out_ids(size(this%gene_ids)))
        do i = 1, size(this%gene_ids)
            out_ids(i) = this%gene_ids(i)
        end do
    end subroutine get_gene_ids

    subroutine get_gene_to_families(this, out_mapping, ierr)
        class(tox_data_t), intent(in) :: this
        integer(int32), allocatable, intent(out) :: out_mapping(:)
        integer(int32), intent(out), optional :: ierr
        
        if (present(ierr)) call set_ok(ierr)
        
        if (.not. allocated(this%gene_to_fam)) then
            if (present(ierr)) ierr = ERR_INVALID_INPUT
            allocate(out_mapping(0))
            return
        end if
        
        allocate(out_mapping(size(this%gene_to_fam)))
        out_mapping = this%gene_to_fam
    end subroutine get_gene_to_families

    subroutine get_family_for_gene(this, gene_idx, out_family, ierr)
        class(tox_data_t), intent(in) :: this
        integer(int32), intent(in) :: gene_idx
        character(len=:), allocatable, intent(out) :: out_family
        integer(int32), intent(out), optional :: ierr
        
        if (present(ierr)) call set_ok(ierr)
        
        if (.not. allocated(this%gene_to_fam) .or. .not. allocated(this%gene_family_ids) .or. &
            gene_idx < 1 .or. gene_idx > this%n_genes) then
            if (present(ierr)) ierr = ERR_INVALID_INPUT
            allocate(character(len=1) :: out_family)
            return
        end if
        
        if (this%gene_to_fam(gene_idx) < 1 .or. this%gene_to_fam(gene_idx) > this%n_families) then
            if (present(ierr)) ierr = ERR_INVALID_INPUT
            allocate(character(len=1) :: out_family)
            return
        end if
        
        out_family = this%gene_family_ids(this%gene_to_fam(gene_idx))
    end subroutine get_family_for_gene

    subroutine get_single_gene_id_index(this, gene_id, out_idx, ierr)
        class(tox_data_t), intent(in) :: this
        character(len=*), intent(in) :: gene_id
        integer(int32), intent(out) :: out_idx
        integer(int32), intent(out), optional :: ierr
        integer :: i
        
        if (present(ierr)) call set_ok(ierr)
        
        if (.not. allocated(this%gene_ids)) then
            if (present(ierr)) ierr = ERR_INVALID_INPUT
            out_idx = 0
            return
        end if
        
        out_idx = 0
        do i = 1, size(this%gene_ids)
            if (trim(this%gene_ids(i)) == trim(gene_id)) then
                out_idx = i
                exit
            end if
        end do
        
        if (out_idx == 0 .and. present(ierr)) then
            ierr = ERR_INVALID_INPUT
        end if
    end subroutine get_single_gene_id_index

    subroutine get_family_centroids(this, out_centroids, ierr)
        class(tox_data_t), intent(in) :: this
        real(real64), allocatable, intent(out) :: out_centroids(:,:)
        integer(int32), intent(out), optional :: ierr
        
        if (present(ierr)) call set_ok(ierr)
        
        if (.not. allocated(this%family_centroids)) then
            if (present(ierr)) ierr = ERR_INVALID_INPUT
            allocate(out_centroids(0,0))
            return
        end if
        
        allocate(out_centroids(size(this%family_centroids, 1), size(this%family_centroids, 2)))
        out_centroids = this%family_centroids
    end subroutine get_family_centroids

    subroutine get_single_family_centroid(this, family_idx, out_centroid, ierr)
        class(tox_data_t), intent(in) :: this
        integer(int32), intent(in) :: family_idx
        real(real64), allocatable, intent(out) :: out_centroid(:)
        integer(int32), intent(out), optional :: ierr
        
        if (present(ierr)) call set_ok(ierr)
        
        if (.not. allocated(this%family_centroids) .or. family_idx < 1 .or. family_idx > this%n_families) then
            if (present(ierr)) ierr = ERR_INVALID_INPUT
            allocate(out_centroid(0))
            return
        end if
        
        allocate(out_centroid(size(this%family_centroids, 1)))
        out_centroid = this%family_centroids(:, family_idx)
    end subroutine get_single_family_centroid

    subroutine get_shift_vectors(this, out_vectors, ierr)
        class(tox_data_t), intent(in) :: this
        real(real64), allocatable, intent(out) :: out_vectors(:,:)
        integer(int32), intent(out), optional :: ierr
        
        if (present(ierr)) call set_ok(ierr)
        
        if (.not. allocated(this%shift_vectors)) then
            if (present(ierr)) ierr = ERR_INVALID_INPUT
            allocate(out_vectors(0,0))
            return
        end if
        
        allocate(out_vectors(size(this%shift_vectors, 1), size(this%shift_vectors, 2)))
        out_vectors = this%shift_vectors
    end subroutine get_shift_vectors

    subroutine get_single_shift_vector(this, gene_idx, out_start, out_shift, ierr)
        class(tox_data_t), intent(in) :: this
        integer(int32), intent(in) :: gene_idx
        real(real64), allocatable, intent(out) :: out_start(:), out_shift(:)
        integer(int32), intent(out), optional :: ierr
        
        if (present(ierr)) call set_ok(ierr)
        
        if (.not. allocated(this%shift_vectors) .or. gene_idx < 1 .or. gene_idx > this%n_genes) then
            if (present(ierr)) ierr = ERR_INVALID_INPUT
            allocate(out_start(0), out_shift(0))
            return
        end if
        
        allocate(out_start(this%d), out_shift(this%d))
        out_start = this%shift_vectors(1:this%d, gene_idx)
        out_shift = this%shift_vectors(this%d+1:2*this%d, gene_idx)
    end subroutine get_single_shift_vector
end module 