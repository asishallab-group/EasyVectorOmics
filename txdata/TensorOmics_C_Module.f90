module TensorOmics_C_Module
  use iso_c_binding
  implicit none

  ! C-compatible data structure
  type, bind(c) :: TensorOmics_C
    integer(c_int) :: n_conditions
    integer(c_int) :: n_genes
    type(c_ptr) :: vec_ptr
    type(c_ptr) :: shift_ptr
  end type TensorOmics_C

end module TensorOmics_C_Module
