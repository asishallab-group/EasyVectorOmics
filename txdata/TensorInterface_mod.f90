! TensorInterface: abstract baseclass (Interface)
module TensorInterface_mod
   implicit none
   private
   public :: TensorInterface

   type, abstract :: TensorInterface
   contains
      procedure(calculate_memory_sig), deferred :: calculate_memory_requirements
      procedure(init_sig), deferred :: init
      procedure(update_sig), deferred :: update
      procedure(save_sig), deferred :: save
   end type TensorInterface

   abstract interface
      subroutine calculate_memory_sig(self, data_estimates)
         import :: TensorInterface
         class(TensorInterface), intent(inout) :: self
         integer, dimension(:), intent(in) :: data_estimates
      end subroutine

      subroutine init_sig(self)
         import :: TensorInterface
         class(TensorInterface), intent(inout) :: self
      end subroutine

      subroutine update_sig(self, patch, gene_id)
         import :: TensorInterface
         class(TensorInterface), intent(inout) :: self
         real, dimension(:), intent(in) :: patch
         integer, intent(in), optional :: gene_id  ! For metadata tracking
      end subroutine

      subroutine save_sig(self, filename)
         import :: TensorInterface
         class(TensorInterface), intent(inout) :: self
         character(len=*), intent(in) :: filename
      end subroutine
   end interface

end module TensorInterface_mod