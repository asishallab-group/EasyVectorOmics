! test_tox.f90
program test_tox
   use TensorOmics_mod
   implicit none
   type(TensorOmics_Type) :: tox
   real :: patch(5, 2)
   integer :: indices(2), header(6), header_read(6)
   integer :: unit
   character(len=*), parameter :: filename = "data.txdata"

   print *, "=== TensorOmics Test Suite ==="

   ! Test 1: Initialization
   print *, "Testing initialization..."
   call tox%calculate_memory_requirements([1, 2, 1])
   call tox%init()
   print *, "Initialized for", tox%n_conditions, "conditions and", tox%n_genes, "genes"

   ! Test 2: Data insertion
   print *, "Testing data insertion..."
   patch = reshape([1.0,2.0,3.0,4.0,5.0, 6.0,7.0,8.0,9.0,10.0], [5,2])
   call tox%update(patch, indices)
   print *, "Inserted at indices:", indices

   ! Test 3: Binary save
   print *, "Testing binary save..."
   call tox%save(filename)
   print *, "Saved to ", trim(filename)

   ! Test 4: Binary read verification
   print *, "Testing binary read..."
   open(newunit=unit, file=filename, form='unformatted', status='old')
   read(unit) header_read
   close(unit)
   
   header = [tox%n_conditions, tox%n_genes, 3, 3, 3, 0]
   print *, "Header read:   ", header_read
   print *, "Header expected:", header
   if(all(header_read == header)) then
      print *, "✔ Header match!"
   else
      print *, "✘ Header mismatch!"
   end if

   print *, "=== Tests completed ==="
end program test_tox