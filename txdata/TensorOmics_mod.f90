! =====================================================
! TensorOmics-Module: Full Implementation with Metadata
! =====================================================
module TensorOmics_mod
   use TensorInterface_mod
   implicit none
   private
   public :: TensorOmics_Type

   integer, parameter :: MAX_STR_LEN = 32, META_COLS = 3

   type, extends(TensorInterface) :: TensorOmics_Type
      ! Core data containers
      real, allocatable    :: vec_store(:,:)     ! Expression vectors [n_conditions, n_genes]
      real, allocatable    :: shift_store(:,:)   ! Shift vectors [n_conditions, n_genes]
      
      ! Turing Bands metadata
      integer, allocatable :: meta_int(:,:)      ! [n_genes, META_COLS]
      real, allocatable    :: meta_real(:,:)     ! [n_genes, META_COLS]
      character(MAX_STR_LEN), allocatable :: meta_str(:,:) ! [n_genes, META_COLS]

      ! Memory tracking
      integer :: n_conditions = 0, n_genes = 0, next_idx = 1
   contains
      procedure :: calculate_memory_requirements => tox_calc_mem
      procedure :: init => tox_init
      procedure :: update => tox_update
      procedure :: save => tox_save
   end type TensorOmics_Type

contains
!======================================================================
! Memory calculation and initialization
!======================================================================
   subroutine tox_calc_mem(self, data_estimates)
      class(TensorOmics_Type), intent(inout) :: self
      integer, dimension(:), intent(in) :: data_estimates
      
      self%n_conditions = size(data_estimates)
      self%n_genes = sum(data_estimates)
   end subroutine

   subroutine tox_init(self)
      class(TensorOmics_Type), intent(inout) :: self
      
      allocate(self%vec_store(self%n_conditions, self%n_genes))
      allocate(self%shift_store(self%n_conditions, self%n_genes))
      allocate(self%meta_int(self%n_genes, META_COLS))
      allocate(self%meta_real(self%n_genes, META_COLS))
      allocate(self%meta_str(self%n_genes, META_COLS))
      
      self%vec_store = 0.0
      self%shift_store = 0.0
      self%meta_int = 0
      self%meta_real = 0.0
      self%meta_str = ""
   end subroutine

!======================================================================
! Data insertion with patch tracking
!======================================================================
   subroutine tox_update(self, patch, indices)
      class(TensorOmics_Type), intent(inout) :: self
      real, dimension(:,:), intent(in) :: patch  ! [n_conditions, n_genes]
      integer, dimension(:), intent(out) :: indices
      
      integer :: i, n_genes_patch, start_idx

      n_genes_patch = size(patch, 2)
      start_idx = self%next_idx
      
      do i = 1, n_genes_patch
         if(self%next_idx > self%n_genes) exit
         
         self%vec_store(:, self%next_idx) = patch(:,i)
         self%shift_store(:, self%next_idx) = patch(:,i) - &  ! Simplified centroid
             sum(self%vec_store, dim=2)/self%n_genes
         
         indices(i) = self%next_idx
         self%next_idx = self%next_idx + 1
      end do
   end subroutine

!======================================================================
! Binary serialization with header
!======================================================================
!   subroutine tox_save(self, filename)
!     class(TensorOmics_Type), intent(inout) :: self
!      character(len=*), intent(in) :: filename
!      
!      integer :: unit, header(6)
!      
!      open(newunit=unit, file=filename, form='unformatted', status='replace')
!      
!      ! Header: [n_cond, n_genes, meta_int_cols, meta_real_cols, meta_str_cols, reserved]
!      header = [self%n_conditions, self%n_genes, META_COLS, META_COLS, META_COLS, 0]
!      write(unit) header
!      
!      ! Body: Flat binary data
!      write(unit) self%vec_store
!      write(unit) self%shift_store
!      write(unit) self%meta_int
!      write(unit) self%meta_real
!      write(unit) self%meta_str
!      
!      close(unit)
!   end subroutine

   subroutine tox_save(self, filename)
      class(TensorOmics_Type), intent(inout) :: self
      character(len=*), intent(in) :: filename
      integer :: unit
      open(newunit=unit, file=filename, form='unformatted', status='replace')
      write(unit) [self%n_conditions, self%n_genes, META_COLS, META_COLS, META_COLS, 0]
      write(unit) self%vec_store
      write(unit) self%shift_store
      write(unit) self%meta_int
      write(unit) self%meta_real
      write(unit) self%meta_str
      close(unit)
   end subroutine

end module TensorOmics_mod