program test_parallel
  implicit none

integer :: k 
integer :: ENE
integer :: start_k , end_k 

ENE = 10

! ----------------------
! START OF PARALLEL LOOP
! ----------------------
#ifdef USE_COARRAY
  start_k  = (this_image() - 1) * ENE / num_images() + 1
  end_k    = min(this_image() * ENE / num_images(), ENE)
  do k  = start_k , end_k 

#elif defined(USE_OPENMP)
  !$omp parallel do private(k )
  do k  = 1, ENE

#else
  do k  = 1, ENE
#endif

! BODY OF THE PARALLEL LOOP
! -----------------------------
   print *, "Running i =", k
! -----------------------------
! END BODY OF THE PARALLEL LOOP

  end do

#ifdef USE_COARRAY
  sync all
#endif

#ifdef USE_OPENMP
  !$omp end parallel do
#endif
! --------------------
! END OF PARALLEL LOOP
! --------------------


end program test_parallel
