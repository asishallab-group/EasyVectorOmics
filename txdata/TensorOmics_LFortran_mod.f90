module TensorOmics_LFortran_mod
  use iso_c_binding
  use TensorOmics_C_Module
  implicit none

  type :: TensorOmics_Type
    integer :: n_conditions
    integer :: n_genes
    integer :: next_idx = 1
    real, allocatable :: vec_store(:,:)
    real, allocatable :: shift_store(:,:)
  end type

contains

  subroutine init_tensoromics_handle(handle, n_cond, n_genes) bind(c, name="init_tensoromics")
    type(c_ptr), intent(out) :: handle
    integer(c_int), intent(in), value :: n_cond, n_genes
    type(TensorOmics_Type), pointer :: tom
    
    allocate(tom)
    tom%n_conditions = n_cond
    tom%n_genes = n_genes
    allocate(tom%vec_store(n_cond, n_genes))
    allocate(tom%shift_store(n_cond, n_genes))
    tom%vec_store = 0.0
    tom%shift_store = 0.0
    handle = c_loc(tom)
  end subroutine

  subroutine calculate_memory(n_cond, n_genes) bind(c, name="calculate_memory_requirements")
    integer(c_int), intent(in), value :: n_cond, n_genes
    print *, "Memory required:", 8*n_cond*n_genes*2, "bytes"
  end subroutine

  subroutine update_tensoromics(handle, patch, n_patch, indices) bind(c, name="update_tensoromics")
    type(c_ptr), intent(in), value :: handle
    real(c_float), intent(in) :: patch(*)
    integer(c_int), intent(in), value :: n_patch
    integer(c_int), intent(out) :: indices(*)
    type(TensorOmics_Type), pointer :: tom
    integer :: i, idx

    call c_f_pointer(handle, tom)
    do i = 1, n_patch
      if (tom%next_idx > tom%n_genes) exit
      tom%vec_store(:, tom%next_idx) = patch((i-1)*tom%n_conditions+1 : i*tom%n_conditions)
      indices(i) = tom%next_idx
      tom%next_idx = tom%next_idx + 1
    end do
  end subroutine

  subroutine transfer_for_save(handle, c_tom) bind(c)
    type(c_ptr), intent(in), value :: handle
    type(TensorOmics_C), intent(out) :: c_tom
    type(TensorOmics_Type), pointer :: tom

    call c_f_pointer(handle, tom)
    c_tom%n_conditions = tom%n_conditions
    c_tom%n_genes = tom%n_genes
    c_tom%vec_ptr = c_loc(tom%vec_store)
    c_tom%shift_ptr = c_loc(tom%shift_store)
  end subroutine

end module TensorOmics_LFortran_mod