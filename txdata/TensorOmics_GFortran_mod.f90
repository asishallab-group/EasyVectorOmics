module TensorOmics_GFortran_mod
  use iso_c_binding
  use TensorOmics_C_Module  ! Use the shared C-compatible module
  implicit none

contains

  subroutine save_tensoromics(c_tom, filename) bind(c)
    type(TensorOmics_C), intent(in) :: c_tom
    character(kind=c_char), intent(in) :: filename(*)
    
    ! Local variables
    real(c_float), pointer :: vec_data(:,:)
    real(c_float), pointer :: shift_data(:,:)
    integer :: unit, i
    character(len=256) :: fname_f

    ! Convert C string to Fortran
    do i = 1, 256
      if (filename(i) == c_null_char) exit
      fname_f(i:i) = filename(i)
    end do

    ! Convert pointers to arrays
    call c_f_pointer(c_tom%vec_ptr, vec_data, [c_tom%n_conditions, c_tom%n_genes])
    call c_f_pointer(c_tom%shift_ptr, shift_data, [c_tom%n_conditions, c_tom%n_genes])

    ! Write data
    open(newunit=unit, file=trim(fname_f), form='unformatted', access='stream', status='replace')
    write(unit) c_tom%n_conditions, c_tom%n_genes
    write(unit) vec_data
    write(unit) shift_data
    close(unit)
  end subroutine

end module TensorOmics_GFortran_mod
