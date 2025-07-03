subroutine serialize_r(arr, filename)
    use iso_c_binding
    use array_ops
    implicit none
    integer (int32), pointer, intent(in) :: arr(:)
    character(len=*), intent(in) :: filename
    call serialize(arr, filename)
end subroutine serialize_r

! Integer Deserialisierung
subroutine deserialize_int_flat_r(flat, dims, filename) bind(C)
  use array_ops
  implicit none
  integer(int32), intent(out) :: flat(*)
  integer(int32), intent(out) :: dims(*)
  character(len=*), intent(in) :: filename
  integer(int32), pointer :: f_flat(:)
  integer(int32), allocatable :: f_dims(:)
  
  call deserialize_int_flat(f_flat, f_dims, filename)
  flat(1:size(f_flat)) = f_flat
  dims(1:size(f_dims)) = f_dims
  deallocate(f_flat, f_dims)
end subroutine

! Real Deserialisierung
subroutine deserialize_real_flat_r(flat, dims, filename) bind(C)
  use array_ops
  implicit none
  real(real64), intent(out) :: flat(*)
  integer(int32), intent(out) :: dims(*)
  character(len=*), intent(in) :: filename
  real(real64), pointer :: f_flat(:)
  integer(int32), allocatable :: f_dims(:)
  
  call deserialize_real_flat(f_flat, f_dims, filename)
  flat(1:size(f_flat)) = f_flat
  dims(1:size(f_dims)) = f_dims
  deallocate(f_flat, f_dims)
end subroutine

! Character Deserialisierung
subroutine deserialize_char_flat_r(flat, dims, clen, filename) bind(C)
  use array_ops
  implicit none
  character(len=*), intent(out) :: flat(*)
  integer(int32), intent(out) :: dims(*)
  integer, intent(out) :: clen
  character(len=*), intent(in) :: filename
  character(len=:), pointer :: f_flat(:)
  integer(int32), allocatable :: f_dims(:)
  integer :: f_clen, i
  
  call deserialize_char_flat(f_flat, f_dims, f_clen, filename)
  clen = f_clen
  do i = 1, size(f_flat)
    flat(i) = f_flat(i)
  end do
  dims(1:size(f_dims)) = f_dims
  deallocate(f_flat, f_dims)
end subroutine

subroutine reshape_real_to_1D_r(flat, arr, dims)
    use iso_c_binding
    use reshape_utils
    implicit none
    real(real64), intent(in), target :: flat(:)
    real(real64), pointer, intent(out) :: arr(:)
    integer(int32), intent(in) :: dims(1)
    real(real64), pointer :: flat_ptr(:)
    call c_f_pointer(c_loc(flat(1)), flat_ptr, [size(flat)])
    call reshape_real_to_1D(flat_ptr, arr, dims)
end subroutine reshape_real_to_1D_r

subroutine reshape_real_to_2D_r(flat, arr, dims)
    use iso_c_binding
    use reshape_utils
    implicit none
    real(real64), intent(in), target :: flat(:)
    real(real64), pointer, intent(out) :: arr(:,:)
    integer(int32), intent(in) :: dims(2)
    real(real64), pointer :: flat_ptr(:)
    call c_f_pointer(c_loc(flat(1)), flat_ptr, [size(flat)])
    call reshape_real_to_2D(flat_ptr, arr, dims)
end subroutine reshape_real_to_2D_r

subroutine reshape_real_to_3D_r(flat, arr, dims)
    use iso_c_binding
    use reshape_utils
    implicit none
    real(real64), intent(in), target :: flat(:)
    real(real64), pointer, intent(out) :: arr(:,:,:)
    integer(int32), intent(in) :: dims(3)
    real(real64), pointer :: flat_ptr(:)
    call c_f_pointer(c_loc(flat(1)), flat_ptr, [size(flat)])
    call reshape_real_to_3D(flat_ptr, arr, dims)
end subroutine reshape_real_to_3D_r

subroutine reshape_real_to_4D_r(flat, arr, dims)
    use iso_c_binding
    use reshape_utils
    implicit none
    real(real64), intent(in), target :: flat(:)
    real(real64), pointer, intent(out) :: arr(:,:,:,:)
    integer(int32), intent(in) :: dims(4)
    real(real64), pointer :: flat_ptr(:)
    call c_f_pointer(c_loc(flat(1)), flat_ptr, [size(flat)])
    call reshape_real_to_4D(flat_ptr, arr, dims)
end subroutine reshape_real_to_4D_r

subroutine reshape_real_to_5D_r(flat, arr, dims)
    use iso_c_binding
    use reshape_utils
    implicit none
    real(real64), intent(in), target :: flat(:)
    real(real64), pointer, intent(out) :: arr(:,:,:,:,:)
    integer(int32), intent(in) :: dims(5)
    real(real64), pointer :: flat_ptr(:)
    call c_f_pointer(c_loc(flat(1)), flat_ptr, [size(flat)])
    call reshape_real_to_5D(flat_ptr, arr, dims)
end subroutine reshape_real_to_5D_r

subroutine reshape_int_to_1D_r(flat, arr, dims)
    use iso_c_binding
    use reshape_utils
    implicit none
    integer(int32), intent(in), target :: flat(:)
    integer(int32), pointer, intent(out) :: arr(:)
    integer(int32), intent(in) :: dims(1)
    integer(int32), pointer :: flat_ptr(:)
    call c_f_pointer(c_loc(flat(1)), flat_ptr, [size(flat)])
    call reshape_int_to_1D(flat_ptr, arr, dims)
end subroutine reshape_int_to_1D_r

subroutine reshape_int_to_2D_r(flat, arr, dims)
    use iso_c_binding
    use reshape_utils
    implicit none
    integer(int32), intent(in), target :: flat(:)
    integer(int32), pointer, intent(out) :: arr(:,:)
    integer(int32), intent(in) :: dims(2)
    integer(int32), pointer :: flat_ptr(:)
    call c_f_pointer(c_loc(flat(1)), flat_ptr, [size(flat)])
    call reshape_int_to_2D(flat_ptr, arr, dims)
end subroutine reshape_int_to_2D_r

subroutine reshape_int_to_3D_r(flat, arr, dims)
    use iso_c_binding
    use reshape_utils
    implicit none
    integer(int32), intent(in), target :: flat(:)
    integer(int32), pointer, intent(out) :: arr(:,:,:)
    integer(int32), intent(in) :: dims(3)
    integer(int32), pointer :: flat_ptr(:)
    call c_f_pointer(c_loc(flat(1)), flat_ptr, [size(flat)])
    call reshape_int_to_3D(flat_ptr, arr, dims)
end subroutine reshape_int_to_3D_r

subroutine reshape_int_to_4D_r(flat, arr, dims)
    use iso_c_binding
    use reshape_utils
    implicit none
    integer(int32), intent(in), target :: flat(:)
    integer(int32), pointer, intent(out) :: arr(:,:,:,:)
    integer(int32), intent(in) :: dims(4)
    integer(int32), pointer :: flat_ptr(:)
    call c_f_pointer(c_loc(flat(1)), flat_ptr, [size(flat)])
    call reshape_int_to_4D(flat_ptr, arr, dims)
end subroutine reshape_int_to_4D_r

subroutine reshape_int_to_5D_r(flat, arr, dims)
    use iso_c_binding
    use reshape_utils
    implicit none
    integer(int32), intent(in), target :: flat(:)
    integer(int32), pointer, intent(out) :: arr(:,:,:,:,:)
    integer(int32), intent(in) :: dims(5)
    integer(int32), pointer :: flat_ptr(:)
    call c_f_pointer(c_loc(flat(1)), flat_ptr, [size(flat)])
    call reshape_int_to_5D(flat_ptr, arr, dims)
end subroutine reshape_int_to_5D_r

! ============================================
! Reshape-Routinen für CHARACTER-Arrays
! ============================================
subroutine reshape_char_to_1D_r(flat, arr, dims, clen)
    use iso_c_binding
    use reshape_utils
    implicit none
    character(len=*), allocatable, intent(in), target :: flat(:)
    character(len=:), pointer, intent(out) :: arr(:)
    integer(int32), intent(in) :: dims(1)
    integer, intent(in) :: clen
    character(len=:), pointer :: flat_ptr(:)
    call c_f_pointer(c_loc(flat(1)), flat_ptr, [size(flat)])
    call reshape_char_to_1D(flat_ptr, arr, dims, clen)
end subroutine reshape_char_to_1D_r

subroutine reshape_char_to_2D_r(flat, arr, dims, clen)
    use iso_c_binding
    use reshape_utils
    implicit none
    character(len=*), allocatable, intent(in), target :: flat(:)
    character(len=:), pointer, intent(out) :: arr(:,:)
    integer(int32), intent(in) :: dims(2)
    integer, intent(in) :: clen
    character(len=:), pointer :: flat_ptr(:)
    call c_f_pointer(c_loc(flat(1)), flat_ptr, [size(flat)])
    call reshape_char_to_2D(flat_ptr, arr, dims, clen)
end subroutine reshape_char_to_2D_r

subroutine reshape_char_to_3D_r(flat, arr, dims, clen)
    use iso_c_binding
    use reshape_utils
    implicit none
    character(len=*), allocatable, intent(in), target :: flat(:)
    character(len=:), pointer, intent(out) :: arr(:,:,:)
    integer(int32), intent(in) :: dims(3)
    integer, intent(in) :: clen
    character(len=:), pointer :: flat_ptr(:)
    call c_f_pointer(c_loc(flat(1)), flat_ptr, [size(flat)])
    call reshape_char_to_3D(flat_ptr, arr, dims, clen)
end subroutine reshape_char_to_3D_r

subroutine reshape_char_to_4D_r(flat, arr, dims, clen)
    use iso_c_binding
    use reshape_utils
    implicit none
    character(len=*), allocatable, intent(in), target :: flat(:)
    character(len=:), pointer, intent(out) :: arr(:,:,:,:)
    integer(int32), intent(in) :: dims(4)
    integer, intent(in) :: clen
    character(len=:), pointer :: flat_ptr(:)
    call c_f_pointer(c_loc(flat(1)), flat_ptr, [size(flat)])
    call reshape_char_to_4D(flat_ptr, arr, dims, clen)
end subroutine reshape_char_to_4D_r

subroutine reshape_char_to_5D_r(flat, arr, dims, clen)
    use iso_c_binding
    use reshape_utils
    implicit none
    character(len=*), allocatable, intent(in), target :: flat(:)
    character(len=:), pointer, intent(out) :: arr(:,:,:,:,:)
    integer(int32), intent(in) :: dims(5)
    integer, intent(in) :: clen
    character(len=:), pointer :: flat_ptr(:)
    call c_f_pointer(c_loc(flat(1)), flat_ptr, [size(flat)])
    call reshape_char_to_5D(flat_ptr, arr, dims, clen)
end subroutine reshape_char_to_5D_r