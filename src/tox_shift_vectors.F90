!> Module for computing the shift vector field for all genes.
module tox_shift_vectors
   use, intrinsic :: iso_fortran_env, only: real64, int32
contains

   !> Compute the shift vector field for all genes.
   !| Computes the shift vectors by substracting the corresponding family centroid from the expression vector.
   pure subroutine compute_shift_vector_field(d, n_genes, n_families, expression_vectors, family_centroids, gene_to_family, family_ids, shift_vectors, ierr)
      implicit none

      !| Expression vector dimension
      integer(int32), intent(in) :: d
      !| Total number of genes
      integer(int32), intent(in) :: n_genes
      !| Total number of families
      integer(int32), intent(in) :: n_families
      !| Gene expression matrix (d × n_genes)
      real(real64), intent(in) :: expression_vectors(d, n_genes)
      !| Family centroid matrix (d × n_families)
      real(real64), intent(in) :: family_centroids(d, n_families)
      !| Mapping from genes to families (family IDs for each gene in expression_vectors)
      integer(int32), intent(in) :: gene_to_family(n_genes)
      !| Mapping from family_centroids to family IDs (family IDs for each column in family_centroids)
      integer(int32), intent(in) :: family_ids(n_families)
      !| Output, real matrix array, size = 2d x n_genes, stores the centroid of the gene's family in rows 1..d and the shift vectors in rows d+1...2d
      real(real64), intent(out) :: shift_vectors(2*d, n_genes)
      !| TODO  Error code: 0 - success, non-zero = error
      integer(int32), intent(out) :: ierr

      !|Local variables
      integer(int32) :: current_gene, current_family_id, current_centroid, i
      real(real64) :: current_family_centroid(d)

      !| For each gene do
      do current_gene = 1, n_genes
         current_family_id = gene_to_family(current_gene)
         current_centroid = -1

         do i = 1, n_families
            if (family_ids(i) == current_family_id) then
               current_centroid = i
               exit
            end if
         end do

         current_family_centroid = family_centroids(:, current_centroid)
         shift_vectors(1:d, current_gene) = current_family_centroid
         shift_vectors(d + 1:2*d, current_gene) = expression_vectors(:, current_gene) - current_family_centroid
      end do

   end subroutine
end module

!quick testing program
program main
   use tox_shift_vectors
   use, intrinsic :: iso_fortran_env, only: real64, int32

   integer(int32) :: d, n_genes, n_families, ierr
   real(real64), allocatable :: expression_vectors(:, :)
   real(real64), allocatable :: family_centroids(:, :)
   integer(int32), allocatable :: gene_to_family(:)
   integer(int32), allocatable :: family_ids(:)
   real(real64), allocatable :: shift_vectors(:, :)

   ! Example sizes
   d = 2
   n_genes = 3
   n_families = 2

   allocate (expression_vectors(d, n_genes))
   allocate (family_centroids(d, n_families))
   allocate (gene_to_family(n_genes))
   allocate (family_ids(n_families))
   allocate (shift_vectors(2*d, n_genes))

   ! Fill with dummy data
   expression_vectors = reshape([1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64, 6.0_real64], [d, n_genes])
   family_centroids = reshape([10.0_real64, 20.0_real64, 30.0_real64, 40.0_real64], [d, n_families])
   gene_to_family = [1_int32, 2_int32, 1_int32]
   family_ids = [1_int32, 2_int32]
   shift_vectors = 0.0_real64
   call compute_shift_vector_field(d, n_genes, n_families, expression_vectors, family_centroids, gene_to_family, family_ids, shift_vectors, ierr)
   print *, shift_vectors

end program main

