module tox_data_accessors
  use iso_fortran_env, only: int32, real64
  use tox_errors
  implicit none

  public :: get_gene_index, get_family_index
  public :: get_family_for_gene_index, get_expression_vector
  public :: get_family_centroid, get_shift_components

contains
subroutine get_family_for_gene_index(gene_idx, gene_to_fam, out_family_idx, ierr)
    integer(int32), intent(in) :: gene_idx
    integer(int32), intent(in) :: gene_to_fam(:)
    integer(int32), intent(out) :: out_family_idx
    integer(int32), intent(out) :: ierr
    
    ierr = ERR_OK
    
    if (gene_idx < 1 .or. gene_idx > size(gene_to_fam)) then
        ierr = ERR_INVALID_INPUT
        out_family_idx = 0
        return
    end if
    
    out_family_idx = gene_to_fam(gene_idx)
    
    if (out_family_idx < 1) then
        ierr = ERR_INVALID_INPUT
    end if
end subroutine get_family_for_gene_index

subroutine get_expression_vector(gene_idx, expression_vectors, out_vector, ierr)
    integer(int32), intent(in) :: gene_idx
    real(real64), intent(in) :: expression_vectors(:,:)
    real(real64), intent(out) :: out_vector(:)
    integer(int32), intent(out) :: ierr
    
    ierr = ERR_OK
    
    if (gene_idx < 1 .or. gene_idx > size(expression_vectors, 2)) then
        ierr = ERR_INVALID_INPUT
        return
    end if
    
    if (size(out_vector) /= size(expression_vectors, 1)) then
        ierr = ERR_SIZE_MISMATCH
        return
    end if
    
    out_vector = expression_vectors(:, gene_idx)
end subroutine get_expression_vector

subroutine get_family_centroid(family_idx, family_centroids, out_centroid, ierr)
    integer(int32), intent(in) :: family_idx
    real(real64), intent(in) :: family_centroids(:,:)
    real(real64), intent(out) :: out_centroid(:)
    integer(int32), intent(out) :: ierr
    
    ierr = ERR_OK
    
    if (family_idx < 1 .or. family_idx > size(family_centroids, 2)) then
        ierr = ERR_INVALID_INPUT
        return
    end if
    
    if (size(out_centroid) /= size(family_centroids, 1)) then
        ierr = ERR_SIZE_MISMATCH
        return
    end if
    
    out_centroid = family_centroids(:, family_idx)
end subroutine get_family_centroid

subroutine get_shift_components(gene_idx, shift_vectors, d, out_start, out_shift, ierr)
    integer(int32), intent(in) :: gene_idx
    real(real64), intent(in) :: shift_vectors(:,:)
    integer(int32), intent(in) :: d
    real(real64), intent(out) :: out_start(:), out_shift(:)
    integer(int32), intent(out) :: ierr
    
    ierr = ERR_OK
    
    if (gene_idx < 1 .or. gene_idx > size(shift_vectors, 2)) then
        ierr = ERR_INVALID_INPUT
        return
    end if
    
    if (size(shift_vectors, 1) < 2*d) then
        ierr = ERR_SIZE_MISMATCH
        return
    end if
    
    if (size(out_start) /= d .or. size(out_shift) /= d) then
        ierr = ERR_SIZE_MISMATCH
        return
    end if
    
    out_start = shift_vectors(1:d, gene_idx)
    out_shift = shift_vectors(d+1:2*d, gene_idx)
end subroutine get_shift_components

! Helper function to find gene index with early exit optimization
integer function get_gene_index(gene_ids, gene) result(idx)
    character(len=*), intent(in) :: gene_ids(:)
    character(len=*), intent(in) :: gene
    integer :: i
    
    idx = 0
    ! Quick check: if gene is empty, return immediately
    if (len_trim(gene) == 0) return
    
    ! Search through gene list
    do i = 1, size(gene_ids)
        if (trim(adjustl(gene_ids(i))) == trim(adjustl(gene))) then
            idx = i
            return
        end if
    end do
end function get_gene_index

! Helper function to find family index
integer function get_family_index(family_ids, family) result(idx)
    character(len=*), intent(in) :: family_ids(:)
    character(len=*), intent(in) :: family
    integer :: i
    
    idx = 0
    do i = 1, size(family_ids)
        if (trim(adjustl(family_ids(i))) == trim(adjustl(family))) then
            idx = i
            return
        end if
    end do
end function get_family_index

end module tox_data_accessors