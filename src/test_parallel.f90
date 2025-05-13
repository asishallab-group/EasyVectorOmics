program test_parallel
  implicit none

integer :: i
integer :: N
integer :: start_i, end_i

N = 15

! ----------------------
! START OF PARALLEL LOOP
! ----------------------
#ifdef USE_COARRAY
  start_i = (this_image() - 1) * N / num_images() + 1
  end_i   = min(this_image() * N / num_images(), N)
  do i = start_i, end_i

#elif defined(USE_OPENMP)
  !$omp parallel do private(i)
  do i = 1, N

#else
  do i = 1, N
#endif

! BODY OF THE PARALLEL LOOP
! -----------------------------
    print *, "Running i =", i
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
