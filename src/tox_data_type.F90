!> Type definition and basic operations for tox_data_t
module tox_data_type
use iso_fortran_env, only: real64, int32
use tox_errors, only: ERR_OK, ERR_INVALID_INPUT, ERR_ALLOC_FAIL
implicit none
private

public :: tox_data_t
public :: create_tox_data, destroy_tox_data
public :: validate_tox_data

type :: tox_data_t
    integer(int32) :: n_genes = 0
    integer(int32) :: n_families = 0
    integer(int32) :: d = 0

    character(len=:), allocatable :: gene_ids(:)
    character(len=:), allocatable :: gene_family_ids(:)
    integer(int32), allocatable :: gene_to_fam(:)
    real(real64), allocatable :: expression_vectors(:,:)
    real(real64), allocatable :: shift_vectors(:,:)
    real(real64), allocatable :: family_centroids(:,:)

contains
    procedure :: allocate_fields => tox_allocate_fields
    procedure :: free_fields => tox_free_fields
    procedure :: set_gene_ids => tox_set_gene_ids
    procedure :: get_gene_ids => tox_get_gene_ids
    procedure :: set_gene_families => tox_set_gene_families
    procedure :: set_expression_vectors => tox_set_expression_vectors
    procedure :: get_expression_vectors => tox_get_expression_vectors
end type tox_data_t

contains

subroutine create_tox_data(obj, ierr)
    type(tox_data_t), intent(out) :: obj
    integer(int32), intent(out), optional :: ierr
    if (present(ierr)) call set_ok(ierr)
    obj%n_genes = 0
    obj%n_families = 0
    obj%d = 0
end subroutine create_tox_data

subroutine destroy_tox_data(obj, ierr)
    type(tox_data_t), intent(inout) :: obj
    integer(int32), intent(out), optional :: ierr
    if (present(ierr)) call set_ok(ierr)
    call obj%free_fields()
    obj%n_genes = 0
    obj%n_families = 0
    obj%d = 0
end subroutine destroy_tox_data

subroutine tox_allocate_fields(this, n_genes_in, n_families_in, d_in, ierr)
    class(tox_data_t), intent(inout) :: this
    integer(int32), intent(in) :: n_genes_in, n_families_in, d_in
    integer(int32), intent(out), optional :: ierr
    integer :: stat

    if (present(ierr)) call set_ok(ierr)

    if (n_genes_in < 0 .or. n_families_in < 0 .or. d_in < 0) then
    if (present(ierr)) ierr = ERR_INVALID_INPUT
    return
    end if

    call this%free_fields()

    this%n_genes = n_genes_in
    this%n_families = n_families_in
    this%d = d_in

    if (this%n_genes > 0) then
    allocate(character(len=1) :: this%gene_ids(this%n_genes), stat=stat)
    if (stat /= 0) then
        if (present(ierr)) ierr = ERR_ALLOC_FAIL
        return
    end if
    
    allocate(this%gene_to_fam(this%n_genes), stat=stat)
    if (stat /= 0) then
        if (present(ierr)) ierr = ERR_ALLOC_FAIL
        return
    end if
    
    allocate(this%expression_vectors(max(0,this%d), this%n_genes), stat=stat)
    if (stat /= 0) then
        if (present(ierr)) ierr = ERR_ALLOC_FAIL
        return
    end if
    
    allocate(this%shift_vectors(max(0,2*this%d), this%n_genes), stat=stat)
    if (stat /= 0) then
        if (present(ierr)) ierr = ERR_ALLOC_FAIL
        return
    end if
    else
    allocate(character(len=1) :: this%gene_ids(0))
    allocate(this%gene_to_fam(0))
    allocate(this%expression_vectors(0,0))
    allocate(this%shift_vectors(0,0))
    end if

    if (this%n_families > 0) then
    allocate(character(len=1) :: this%gene_family_ids(this%n_families), stat=stat)
    if (stat /= 0) then
        if (present(ierr)) ierr = ERR_ALLOC_FAIL
        return
    end if
    
    allocate(this%family_centroids(max(0,this%d), this%n_families), stat=stat)
    if (stat /= 0) then
        if (present(ierr)) ierr = ERR_ALLOC_FAIL
        return
    end if
    else
    allocate(character(len=1) :: this%gene_family_ids(0))
    allocate(this%family_centroids(0,0))
    end if

end subroutine tox_allocate_fields

subroutine tox_free_fields(this, ierr)
    class(tox_data_t), intent(inout) :: this
    integer(int32), intent(out), optional :: ierr
    if (present(ierr)) call set_ok(ierr)
    if (allocated(this%gene_ids)) deallocate(this%gene_ids)
    if (allocated(this%gene_family_ids)) deallocate(this%gene_family_ids)
    if (allocated(this%gene_to_fam)) deallocate(this%gene_to_fam)
    if (allocated(this%expression_vectors)) deallocate(this%expression_vectors)
    if (allocated(this%shift_vectors)) deallocate(this%shift_vectors)
    if (allocated(this%family_centroids)) deallocate(this%family_centroids)
    this%n_genes = 0
    this%n_families = 0
    this%d = 0
end subroutine tox_free_fields

subroutine tox_set_gene_ids(this, ids, ierr)
    class(tox_data_t), intent(inout) :: this
    character(len=*), dimension(:), intent(in) :: ids
    integer(int32), intent(out), optional :: ierr
    integer :: n, i, stat

    if (present(ierr)) call set_ok(ierr)
    n = size(ids)
    if (n < 0) then
        if (present(ierr)) ierr = ERR_INVALID_INPUT
        return
    end if

    call this%free_fields()
    this%n_genes = n
    this%n_families = 0
    this%d = 0

    allocate(character(len=len_trim(ids(1))) :: this%gene_ids(n), stat=stat)
    if (stat /= 0) then
        if (present(ierr)) ierr = ERR_ALLOC_FAIL
        return
    end if

    do i = 1, n
    this%gene_ids(i) = adjustl(trim(ids(i)))
    end do

    allocate(this%gene_to_fam(n), stat=stat)
    if (stat /= 0) then
        if (present(ierr)) ierr = ERR_ALLOC_FAIL
        return
    end if
    this%gene_to_fam = 0

    allocate(this%expression_vectors(0, n))
    allocate(this%shift_vectors(0, n))

end subroutine tox_set_gene_ids

subroutine tox_get_gene_ids(this, out_ids, ierr)
    class(tox_data_t), intent(in) :: this
    character(len=:), allocatable, intent(out) :: out_ids(:)
    integer(int32), intent(out), optional :: ierr
    integer :: i

    if (present(ierr)) call set_ok(ierr)

    if (.not. allocated(this%gene_ids)) then
        allocate(character(len=1) :: out_ids(0))
        return
    end if

    allocate(character(len=len_trim(this%gene_ids(1))) :: out_ids(size(this%gene_ids)))
    do i = 1, size(this%gene_ids)
    out_ids(i) = this%gene_ids(i)
    end do

end subroutine tox_get_gene_ids

subroutine tox_set_gene_families(this, fam_ids, gene_to_fam_in, ierr)
    class(tox_data_t), intent(inout) :: this
    character(len=*), dimension(:), intent(in) :: fam_ids
    integer(int32), dimension(:), intent(in) :: gene_to_fam_in
    integer(int32), intent(out), optional :: ierr
    integer :: nf, ng, stat

    if (present(ierr)) call set_ok(ierr)

    nf = size(fam_ids)
    ng = size(gene_to_fam_in)
    if (ng /= this%n_genes) then
        if (present(ierr)) ierr = ERR_SIZE_MISMATCH
        return
    end if

    call this%free_fields()
    this%n_genes = ng
    this%n_families = nf
    this%d = 0

    allocate(character(len=len_trim(fam_ids(1))) :: this%gene_family_ids(nf), stat=stat)
    if (stat /= 0) then
        if (present(ierr)) ierr = ERR_ALLOC_FAIL
        return
    end if
    this%gene_family_ids = fam_ids

    allocate(this%gene_to_fam(ng), stat=stat)
    if (stat /= 0) then
        if (present(ierr)) ierr = ERR_ALLOC_FAIL
        return
    end if
    this%gene_to_fam = gene_to_fam_in

    allocate(this%expression_vectors(0, ng))
    allocate(this%shift_vectors(0, ng))
    allocate(this%family_centroids(0, nf))

end subroutine tox_set_gene_families

subroutine tox_set_expression_vectors(this, mat, ierr)
    class(tox_data_t), intent(inout) :: this
    real(real64), dimension(:,:), intent(in) :: mat
    integer(int32), intent(out), optional :: ierr
    integer :: d_in, ng, stat
    text

    if (present(ierr)) call set_ok(ierr)

    d_in = size(mat,1)
    ng = size(mat,2)

    if (this%n_genes /= 0 .and. this%n_genes /= ng) then
        if (present(ierr)) ierr = ERR_SIZE_MISMATCH
        return
    end if

    this%n_genes = ng
    this%d = d_in

    if (allocated(this%expression_vectors)) deallocate(this%expression_vectors)
    allocate(this%expression_vectors(this%d, this%n_genes), stat=stat)
    if (stat /= 0) then
        if (present(ierr)) ierr = ERR_ALLOC_FAIL
        return
    end if
    this%expression_vectors = mat

    if (allocated(this%shift_vectors)) deallocate(this%shift_vectors)
    allocate(this%shift_vectors(2*this%d, this%n_genes), stat=stat)
    if (stat /= 0) then
        if (present(ierr)) ierr = ERR_ALLOC_FAIL
        return
    end if
    this%shift_vectors = 0.0_real64

    if (allocated(this%family_centroids)) then
    if (size(this%family_centroids,1) /= this%d) then
        deallocate(this%family_centroids)
        allocate(this%family_centroids(this%d, this%n_families), stat=stat)
        if (stat /= 0) then
            if (present(ierr)) ierr = ERR_ALLOC_FAIL
            return
        end if
    end if
    else
    allocate(this%family_centroids(this%d, this%n_families), stat=stat)
    if (stat /= 0) then
        if (present(ierr)) ierr = ERR_ALLOC_FAIL
        return
    end if
    end if

end subroutine tox_set_expression_vectors

subroutine tox_get_expression_vectors(this, out_mat, ierr)
    class(tox_data_t), intent(in) :: this
    real(real64), allocatable, intent(out) :: out_mat(:,:)
    integer(int32), intent(out), optional :: ierr
    if (present(ierr)) call set_ok(ierr)

    if (.not. allocated(this%expression_vectors)) then
        allocate(out_mat(0,0))
        return
    end if
    allocate(out_mat(size(this%expression_vectors,1), size(this%expression_vectors,2)))
    out_mat = this%expression_vectors
end subroutine tox_get_expression_vectors

function validate_tox_data(this) result(ok)
    class(tox_data_t), intent(in) :: this
    logical :: ok

    ok = .true.
    if (this%n_genes < 0 .or. this%n_families < 0 .or. this%d < 0) then
    ok = .false.
    return
    end if

    if (allocated(this%gene_ids)) then
    if (size(this%gene_ids) /= this%n_genes) then
        ok = .false.
        return
    end if
    end if

    if (allocated(this%gene_family_ids)) then
    if (size(this%gene_family_ids) /= this%n_families) then
        ok = .false.
        return
    end if
    end if

    if (allocated(this%gene_to_fam)) then
    if (size(this%gene_to_fam) /= this%n_genes) then
        ok = .false.
        return
    end if
    if (any(this%gene_to_fam < 1) .and. this%n_families > 0) then
        ok = .false.
        return
    end if
    end if

    if (allocated(this%expression_vectors)) then
    if (size(this%expression_vectors,1) /= this%d .or. &
        size(this%expression_vectors,2) /= this%n_genes) then
        ok = .false.
        return
    end if
    end if

    if (allocated(this%shift_vectors)) then
    if (size(this%shift_vectors,1) /= 2*this%d .or. &
        size(this%shift_vectors,2) /= this%n_genes) then
        ok = .false.
        return
    end if
    end if

    if (allocated(this%family_centroids)) then
    if (size(this%family_centroids,1) /= this%d .or. &
        size(this%family_centroids,2) /= this%n_families) then
        ok = .false.
        return
    end if
    end if

end function validate_tox_data

end module tox_data_type