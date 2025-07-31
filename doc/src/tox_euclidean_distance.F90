!> Module with Euclidean distance computation routines for tensor omics.
module tox_euclidean_distance
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none

contains

  !> Compute the Euclidean distance between two vectors.
  !> 
  !> Calculates the L2 norm: result = sqrt(sum((vec1_i - vec2_i)^2))
  !>
  !> @param vec1   First expression vector<br>
  !> @param vec2   Second expression vector<br>
  !> @param d      Dimension of both vectors<br>
  !> @param result Output scalar distance<br>
  pure subroutine euclidean_distance(vec1, vec2, d, result)
    implicit none
    
    integer, intent(in) :: d
    real(real64), intent(in) :: vec1(d), vec2(d)
    real(real64), intent(out) :: result
    
    integer :: i
    real(real64) :: sum_squared_diff
    
    sum_squared_diff = 0.0_real64
    
    do i = 1, d
      sum_squared_diff = sum_squared_diff + (vec1(i) - vec2(i))**2
    end do
    
    result = sqrt(sum_squared_diff)
    
  end subroutine euclidean_distance

  !> Compute distance from each gene to its corresponding family centroid.
  !>
  !> For each gene, extracts its expression vector and the centroid of its
  !> assigned family, then computes the Euclidean distance between them.
  !>
  !> @param n_genes     Total number of genes
  !> @param n_families  Total number of gene families  
  !> @param genes       Gene expression matrix (d × n_genes) - column-major optimized
  !> @param centroids   Family centroid matrix (d × n_families) - column-major optimized
  !> @param gene_to_fam Gene-to-family mapping (1-based indexing)
  !> @param distances   Output distances array
  !> @param d           Expression vector dimension
  pure subroutine distance_to_centroid(n_genes, n_families, genes, centroids, &
                                       gene_to_fam, distances, d)
      implicit none

      ! Arguments
      integer, intent(in) :: n_genes, n_families, d
      real(real64), intent(in) :: genes(d, n_genes)
      real(real64), intent(in) :: centroids(d, n_families)
      integer, intent(in) :: gene_to_fam(n_genes)
      real(real64), intent(out) :: distances(n_genes)

      ! Local variables
      integer :: i, family_idx

      ! Iterate over all genes
      do i = 1, n_genes
          ! Get family index for current gene
          family_idx = gene_to_fam(i)

          ! Validate family index (1-based indexing)
          if (family_idx < 1 .or. family_idx > n_families) then
              distances(i) = -1.0_real64  ! Error indicator
              cycle
          end if

          ! Compute Euclidean distance using column-major 
          call euclidean_distance(genes(:, i), centroids(:, family_idx), d, distances(i))
      end do
    
  end subroutine distance_to_centroid




end module tox_euclidean_distance


subroutine euclidean_distance_r(vec1, vec2, d, result)
  use tox_euclidean_distance
    integer, intent(in) :: d
    real(real64), intent(in) :: vec1(d), vec2(d)
    real(real64), intent(out) :: result
  call euclidean_distance(vec1, vec2, d, result)
end subroutine euclidean_distance_r

subroutine euclidean_distance_c(vec1, vec2, d, result) bind(C, name="euclidean_distance_c")
  use iso_c_binding
  use tox_euclidean_distance
  real(c_double), intent(in), target :: vec1(*)
  real(c_double), intent(in), target :: vec2(*)
  integer(c_int), intent(in), value :: d
  real(c_double), intent(out) :: result

  call euclidean_distance(vec1, vec2, d, result)
end subroutine euclidean_distance_c

subroutine distance_to_centroid_r(n_genes, n_families, genes, centroids, &
                                  gene_to_fam, distances, d)
  use tox_euclidean_distance
  integer, intent(in) :: n_genes, n_families, d
  real(real64), intent(in) :: genes(d, n_genes)
  real(real64), intent(in) :: centroids(d, n_families)
  integer, intent(in) :: gene_to_fam(n_genes)
  real(real64), intent(out) :: distances(n_genes)
  call distance_to_centroid(n_genes, n_families, genes, centroids, &
                                  gene_to_fam, distances, d)
end subroutine distance_to_centroid_r

subroutine distance_to_centroid_c(n_genes, n_families, genes, centroids, & 
                                  gene_to_fam, distances, d) bind(C, name="distance_to_centroid_c")
  use iso_c_binding
  use tox_euclidean_distance
  integer(c_int), intent(in), value :: n_genes
  integer(c_int), intent(in), value :: n_families
  real(c_double), intent(in), target :: genes(*)
  real(c_double), intent(in), target :: centroids(*)
  integer(c_int), intent(in), target :: gene_to_fam(*)
  real(c_double), intent(out), target :: distances(*)
  integer(c_int), intent(in), value :: d

  call distance_to_centroid(n_genes, n_families, genes, centroids, &
                            gene_to_fam, distances, d)
end subroutine distance_to_centroid_c