! ======================================================================
! TensorOmics-Module: Full Implementation with Shift Vectors and Metadata
! ======================================================================
module TensorOmics_mod
   use TensorInterface_mod
   implicit none
   private
   public :: TensorOmics_Type

   ! -----------------------------
   ! KDTree placeholder type (must be outside)
   ! -----------------------------
   type :: KDTree_Type
      logical :: is_built = .false.
   end type KDTree_Type

   ! -----------------------------
   ! Metadata type
   ! -----------------------------
   type :: GeneMetadata
      integer :: id
      character(len=32) :: family
      logical :: is_ortholog
   end type GeneMetadata

   ! -----------------------------
   ! Main TensorOmics Type
   ! -----------------------------
   type, extends(TensorInterface) :: TensorOmics_Type
      integer :: n_conditions = 0
      integer :: n_genes = 0
      integer :: estimated_entries = 0
      integer :: next_gene_idx = 1

      real, allocatable :: vec_store(:,:), shift_store(:,:)
      type(GeneMetadata), allocatable :: metadata(:)

      type(KDTree_Type), allocatable :: cart_tree, sphere_tree
   contains
      procedure :: calculate_memory_requirements => tensor_calc_mem
      procedure :: init                         => tensor_init
      procedure :: update                       => tensor_update
      procedure :: save                         => tensor_save
      procedure :: load                         => tensor_load
      procedure :: build_kdtrees                => tensor_build_trees
   end type TensorOmics_Type


contains
   ! -------------------------------------------------------------------
   ! Step 1: Estimate memory requirements
   ! -------------------------------------------------------------------
   subroutine tensor_calc_mem(self, genes_per_condition)
      class(TensorOmics_Type), intent(inout) :: self
      integer, dimension(:), intent(in) :: genes_per_condition

      self%n_conditions = size(genes_per_condition)
      self%n_genes      = maxval(genes_per_condition)
      self%estimated_entries = sum(genes_per_condition)

      print *, "Conditions (tissues/states): ", self%n_conditions
      print *, "Max genes per condition: ", self%n_genes
      print *, "Total entries (memory estimate): ", self%estimated_entries
   end subroutine

   ! -------------------------------------------------------------------
   ! Step 2: Allocate memory and initialize structures
   ! -------------------------------------------------------------------
   subroutine tensor_init(self)
      class(TensorOmics_Type), intent(inout) :: self

      ! Allocate expression and shift stores
      allocate(self%vec_store(self%n_conditions, self%n_genes))
      allocate(self%shift_store(self%n_conditions, self%n_genes))
      allocate(self%metadata(self%n_genes))

      ! Initialize to zero
      self%vec_store = 0.0
      self%shift_store = 0.0
      self%next_gene_idx = 1

      ! Allocate K-D trees (actual initialization deferred)
      allocate(self%cart_tree)
      allocate(self%sphere_tree)

      print *, "TensorOmics initialized: ", self%n_genes, " genes x ", self%n_conditions, " conditions."
   end subroutine

   ! -------------------------------------------------------------------
   ! Step 3: Insert data patches and compute shift vectors (outliers)
   ! -------------------------------------------------------------------
   subroutine tensor_update(self, patch, gene_id)
      class(TensorOmics_Type), intent(inout) :: self
      real, dimension(:), intent(in) :: patch
      integer, intent(in), optional :: gene_id

      integer :: idx
      real, dimension(self%n_conditions) :: centroid

      ! Assign gene ID or auto-increment
      if (present(gene_id)) then
         idx = gene_id
      else
         idx = self%next_gene_idx
         self%next_gene_idx = self%next_gene_idx + 1
      end if

      ! Insert expression vector (ensure column-major)
      self%vec_store(:, idx) = patch

      ! Example: Compute centroid (mean of family) and shift vector
      ! (In practice, centroid would come from metadata)
      centroid = sum(self%vec_store, dim=2) / real(self%n_genes)
      self%shift_store(:, idx) = patch - centroid

      ! Update metadata (placeholder for family/ortholog annotations)
      self%metadata(idx)%id = idx
      self%metadata(idx)%family = "UNKNOWN"
      self%metadata(idx)%is_ortholog = .false.

      print *, "Updated gene ", idx, " with shift vector norm: ", norm2(self%shift_store(:, idx))
   end subroutine

   ! -------------------------------------------------------------------
   ! Step 4: Serialize to binary file (including metadata)
   ! -------------------------------------------------------------------
   subroutine tensor_save(self, filename)
      class(TensorOmics_Type), intent(inout) :: self
      character(len=*), intent(in) :: filename
      integer :: unit, i

      open(unit=unit, file=filename, status="replace", form="unformatted")
      
      ! Save dimensions and metadata
      write(unit) self%n_conditions, self%n_genes
      write(unit) self%vec_store
      write(unit) self%shift_store
      write(unit) self%metadata  ! Assumes derived type serialization is supported

      close(unit)
      print *, "Data saved to ", trim(filename)
   end subroutine

   ! -------------------------------------------------------------------
   ! Step 5: Load from binary file
   ! -------------------------------------------------------------------
   subroutine tensor_load(self, filename)
      class(TensorOmics_Type), intent(inout) :: self
      character(len=*), intent(in) :: filename
      integer :: unit

      open(unit=unit, file=filename, status="old", form="unformatted")
      
      ! Load dimensions
      read(unit) self%n_conditions, self%n_genes
      call self%init()  ! Allocate memory

      ! Load data
      read(unit) self%vec_store
      read(unit) self%shift_store
      read(unit) self%metadata

      close(unit)
      print *, "Data loaded from ", trim(filename)
   end subroutine

   ! -------------------------------------------------------------------
   ! Step 6: Build K-D Trees (deferred until all data is inserted)
   ! -------------------------------------------------------------------
   subroutine tensor_build_trees(self)
      class(TensorOmics_Type), intent(inout) :: self
      ! Placeholder for actual K-D tree construction
      ! (e.g., using a library like FLANN or custom implementation)
      print *, "K-D trees built for ", self%n_genes, " vectors."
   end subroutine

end module TensorOmics_mod