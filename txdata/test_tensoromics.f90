program test_tensoromics
   use TensorOmics_mod
   type(TensorOmics_Type) :: tensor
   real :: patch(3)
   integer :: genes_per_condition(3) = [100, 150, 200]

   ! Initialize
   call tensor%calculate_memory_requirements(genes_per_condition)
   call tensor%init()

   ! Insert data
   patch = [1.0, 2.0, 3.0]
   call tensor%update(patch, gene_id=1)  ! With explicit ID
   call tensor%update([4.0, 5.0, 6.0])  ! Auto-increment ID

   ! Save and build trees
   call tensor%save("test_data.txdata")
   call tensor%build_kdtrees()
end program