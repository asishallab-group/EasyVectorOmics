!> Module providing serialization and routines for integer, real, and character arrays
!! of up to 5 dimensions, as well as utility routines for extracting rows, columns, and cells from 2D arrays.
!! Arrays are serialized to a custom binary format with a magic number and type/dimension metadata.
!! @note this is currently not in use
module serialize_mod
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding, only: c_loc
  use serialize_char
  use serialize_int
  use serialize_real
  implicit none

  public::serialize

  interface serialize
    module procedure serialize_int_1d, serialize_int_2d, serialize_int_3d, serialize_int_4d, serialize_int_5d
    module procedure serialize_real_1d, serialize_real_2d, serialize_real_3d, serialize_real_4d, serialize_real_5d
    module procedure serialize_char_1d, serialize_char_2d, serialize_char_3d, serialize_char_4d, serialize_char_5d
  end interface serialize

end module serialize_mod