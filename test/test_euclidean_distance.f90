program test_euclidean_distance
    use euclidean_distance_module
    implicit none
    
    ! Parameters for the test
    integer, parameter :: n_genes = 6
    integer, parameter :: n_families = 2
    integer, parameter :: d = 3  ! 3D expression vectors
    
    ! Variables
    real(dp) :: genes(n_genes, d)
    real(dp) :: centroids(n_families, d)
    integer :: gene_to_fam(n_genes)
    real(dp) :: distances(n_genes)
    integer :: i, j
    
    write(*,*) '=== Testing Euclidean Distance Module ==='
    write(*,*)
    
    ! Initialize test data
    ! Family 1: genes 1, 2, 3
    ! Family 2: genes 4, 5, 6
    
    ! Gene expression data (example)
    genes(1, :) = [1.0_dp, 2.0_dp, 3.0_dp]    ! Family 1
    genes(2, :) = [1.5_dp, 2.5_dp, 3.5_dp]    ! Family 1
    genes(3, :) = [0.8_dp, 1.8_dp, 2.8_dp]    ! Family 1
    genes(4, :) = [5.0_dp, 6.0_dp, 7.0_dp]    ! Family 2
    genes(5, :) = [5.2_dp, 6.2_dp, 7.2_dp]    ! Family 2
    genes(6, :) = [4.8_dp, 5.8_dp, 6.8_dp]    ! Family 2
    
    ! Family assignments (1-based indexing)
    gene_to_fam = [1, 1, 1, 2, 2, 2]
    
    ! Centroids (pre-computed for this example)
    centroids(1, :) = [1.1_dp, 2.1_dp, 3.1_dp]    ! Centroid of family 1
    centroids(2, :) = [5.0_dp, 6.0_dp, 7.0_dp]    ! Centroid of family 2
    
    ! Print input data
    write(*,*) 'Input Data:'
    write(*,*) '----------'
    write(*,*) 'Number of genes: ', n_genes
    write(*,*) 'Number of families: ', n_families
    write(*,*) 'Expression vector dimension: ', d
    write(*,*)
    
    write(*,*) 'Gene expression vectors:'
    do i = 1, n_genes
        write(*,'(A,I0,A,3F8.2,A,I0)') 'Gene ', i, ': [', &
              (genes(i, j), j=1,d), '] -> Family ', gene_to_fam(i)
    end do
    write(*,*)
    
    write(*,*) 'Family centroids:'
    do i = 1, n_families
        write(*,'(A,I0,A,3F8.2,A)') 'Family ', i, ': [', &
              (centroids(i, j), j=1,d), ']'
    end do
    write(*,*)
    
    ! Compute distances
    call distance_to_centroid(n_genes, n_families, genes, centroids, &
                             gene_to_fam, distances, d)
    
    ! Print results
    write(*,*) 'Results:'
    write(*,*) '--------'
    do i = 1, n_genes
        write(*,'(A,I0,A,F10.6)') 'Distance from gene ', i, &
              ' to its family centroid: ', distances(i)
    end do
    write(*,*)
    
    ! Print statistics
    call print_distance_statistics(distances, n_genes)
    
    ! Test individual euclidean_distance subroutine
    write(*,*)
    write(*,*) '=== Testing individual euclidean_distance subroutine ==='
    block
        real(dp) :: vec1(3), vec2(3), test_distance
        vec1 = [1.0_dp, 2.0_dp, 3.0_dp]
        vec2 = [4.0_dp, 5.0_dp, 6.0_dp]
        
        call euclidean_distance(vec1, vec2, 3, test_distance)
        
        write(*,'(A,3F6.2,A)') 'Vector 1: [', vec1, ']'
        write(*,'(A,3F6.2,A)') 'Vector 2: [', vec2, ']'
        write(*,'(A,F10.6)') 'Euclidean distance: ', test_distance
        write(*,'(A,F10.6)') 'Expected (sqrt(27)): ', sqrt(27.0_dp)
    end block
    
end program test_euclidean_distance
