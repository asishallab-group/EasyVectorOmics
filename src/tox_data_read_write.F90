module tox_data_read_write
    use iso_fortran_env, only: real64, int32
    use tox_errors
    use serialize_int, only: serialize_int_1d
    use int_deserialize_mod, only: deserialize_int_1D
    use serialize_real, only: serialize_real_2D
    use real_deserialize_mod, only: deserialize_real_2D
    use serialize_char, only: serialize_char_1D
    use char_deserialize_mod, only: deserialize_char_1D
    implicit none
    public :: save_expression_vectors, load_expression_vectors

contains
subroutine save_gene_ids(gene_ids, filename, ierr)
    CHARACTER(len=*), INTENT(IN) :: gene_ids(:)
    CHARACTER(len=*), INTENT(IN) :: filename
    integer, intent(out) :: ierr

    call serialize_char_1D(gene_ids, filename, ierr)
end subroutine

subroutine save_expression_vectors(expression_vectors, filename, ierr)
    CHARACTER(len=*), INTENT(IN) :: filename
    real(real64), INTENT(IN) :: expression_vectors(:,:)
    integer, intent(out) :: ierr

    call serialize_real_2D(expression_vectors, filename, ierr)
end subroutine

subroutine save_gene_to_family(gene_to_fam, filename, ierr)
    CHARACTER(len=*), INTENT(IN) :: filename
    integer, intent(in) :: gene_to_fam(:)
    integer, intent(out) :: ierr

    call serialize_int_1D(gene_to_fam, filename, ierr)
end subroutine

subroutine save_gene_family_ids(gene_family_ids, filename, ierr)
    CHARACTER(len=*), INTENT(IN) :: filename
    CHARACTER(len=*), INTENT(IN) :: gene_family_ids(:)
    integer(int32), INTENT(OUT) :: ierr

    call serialize_char_1D(gene_family_ids, filename, ierr)
end subroutine

subroutine save_family_centroids(family_centroids, filename, ierr)
    CHARACTER(len=*), INTENT(IN) :: filename
    real(real64), INTENT(IN) :: family_centroids(:,:)
    integer(int32), INTENT(OUT) :: ierr

    call serialize_real_2D(family_centroids, filename, ierr)
end subroutine

subroutine save_shift_vectors(shift_vectors, filename, ierr)
    CHARACTER(len=*), INTENT(IN) :: filename
    real(real64), INTENT(IN) :: shift_vectors(:,:)
    integer(int32), INTENT(OUT) :: ierr

    call serialize_real_2D(shift_vectors, filename, ierr)
end subroutine

subroutine load_gene_ids(gene_ids, filename, ierr)
    character(len=*), intent(in)  :: filename
    character(len=*), intent(out) :: gene_ids(:)
    integer(int32), intent(out)   :: ierr

    call deserialize_char_1D(gene_ids, filename, ierr)
end subroutine load_gene_ids

subroutine load_expression_vectors(expression_vectors, filename, ierr)
    character(len=*), intent(in) :: filename
    real(real64), intent(out)    :: expression_vectors(:,:)
    integer(int32), intent(out)  :: ierr

    call deserialize_real_2D(expression_vectors, filename, ierr)
end subroutine load_expression_vectors

subroutine load_gene_to_family(gene_to_fam, filename, ierr)
    character(len=*), INTENT(IN) :: filename
    integer(int32), INTENT(OUT) :: gene_to_fam(:)
    integer(int32), INTENT(OUT) :: ierr

    call deserialize_int_1D(gene_to_fam, filename, ierr)
end subroutine

subroutine load_family_ids(family_ids, filename, ierr)
    character(len=*), intent(in) :: filename
    CHARACTER(len=*), intent(out) :: family_ids(:)
    integer(int32), intent(out) :: ierr

    call deserialize_char_1D(family_ids, filename, ierr)
end subroutine

subroutine load_family_centroids(family_centroids, filename, ierr)
    character(len=*), intent(in) :: filename
    real(real64), intent(out)    :: family_centroids(:,:)
    integer(int32), intent(out)  :: ierr

    call deserialize_real_2D(family_centroids, filename, ierr)
end subroutine

subroutine load_shift_vectors(shift_vectors, filename, ierr)
    character(len=*), intent(in) :: filename
    real(real64), intent(out)    :: shift_vectors(:,:)
    integer(int32), intent(out)  :: ierr

    call deserialize_real_2D(shift_vectors, filename, ierr)
end subroutine

end module tox_data_read_write