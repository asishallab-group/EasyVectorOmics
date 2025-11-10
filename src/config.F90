#include "precompiler_constants.F90"

module config
  implicit none
#ifdef DEFAULT_ALIGNMENT
  integer, parameter :: alignment = DEFAULT_ALIGNMENT
#else
  integer, parameter :: alignment = 32  ! fallback
#endif
  logical, parameter :: DEBUG = .true.
  logical, parameter :: debug_hashing = .false.
end module config
