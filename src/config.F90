module config
  implicit none
#ifdef DEFAULT_ALIGNMENT
  integer, parameter :: alignment = DEFAULT_ALIGNMENT
#endif
end module config
