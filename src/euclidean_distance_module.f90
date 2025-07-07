module euclidean_distance_module
    implicit none
    
    ! Precision for real numbers (double precision)
    integer, parameter :: dp = selected_real_kind(15, 307)
    
contains

    !===========================================================================
    ! Subroutine 1: euclidean_distance
    ! Computes the Euclidean distance between two vectors
    !
    ! Input:
    !   vec1(d)     - First vector
    !   vec2(d)     - Second vector  
    !   d           - Dimension of vectors
    !
    ! Output:
    !   result      - Euclidean distance (scalar)
    !===========================================================================
    subroutine euclidean_distance(vec1, vec2, d, result)
        integer, intent(in) :: d
        real(dp), intent(in) :: vec1(d), vec2(d)
        real(dp), intent(out) :: result
        
        integer :: i
        real(dp) :: sum_squared_diff
        
        ! Initialize sum
        sum_squared_diff = 0.0_dp
        
        ! Calculate sum of squared differences
        do i = 1, d
            sum_squared_diff = sum_squared_diff + (vec1(i) - vec2(i))**2
        end do
        
        ! Take square root to get Euclidean distance
        result = sqrt(sum_squared_diff)
        
    end subroutine euclidean_distance

    !===========================================================================
    ! Subroutine 2: distance_to_centroid
    ! Computes distance from each gene to its corresponding family centroid
    !
    ! Input:
    !   n_genes     - Total number of genes
    !   n_families  - Total number of families
    !   genes       - Matrix of gene expression vectors (n_genes x d)
    !   centroids   - Matrix of centroids (n_families x d)
    !   gene_to_fam - Array mapping each gene to its family (1-based indexing)
    !   d           - Dimension of expression vectors
    !
    ! Output:
    !   distances   - Array of distances from genes to their centroids
    !===========================================================================
    subroutine distance_to_centroid(n_genes, n_families, genes, centroids, &
                                   gene_to_fam, distances, d)
        integer, intent(in) :: n_genes, n_families, d
        real(dp), intent(in) :: genes(n_genes, d)
        real(dp), intent(in) :: centroids(n_families, d)
        integer, intent(in) :: gene_to_fam(n_genes)
        real(dp), intent(out) :: distances(n_genes)
        
        integer :: i, family_idx
        real(dp) :: gene_vector(d), centroid_vector(d)
        
        ! Iterate over all genes
        do i = 1, n_genes
            ! Get the family index for current gene
            family_idx = gene_to_fam(i)
            
            ! Validate family index (1-based indexing)
            if (family_idx < 1 .or. family_idx > n_families) then
                write(*,*) 'Error: Invalid family index for gene ', i, &
                          ': family_idx = ', family_idx
                distances(i) = -1.0_dp  ! Error indicator
                cycle
            end if
            
            ! Extract gene expression vector
            gene_vector(:) = genes(i, :)
            
            ! Extract corresponding centroid vector
            centroid_vector(:) = centroids(family_idx, :)
            
            ! Compute Euclidean distance
            call euclidean_distance(gene_vector, centroid_vector, d, distances(i))
        end do
        
    end subroutine distance_to_centroid

    !===========================================================================
    ! Optional utility subroutine: print_distance_statistics
    ! Prints basic statistics about the computed distances
    !===========================================================================
    subroutine print_distance_statistics(distances, n_genes)
        integer, intent(in) :: n_genes
        real(dp), intent(in) :: distances(n_genes)
        
        real(dp) :: min_dist, max_dist, mean_dist, sum_dist
        integer :: i, valid_count
        
        ! Initialize
        min_dist = huge(1.0_dp)
        max_dist = -huge(1.0_dp)
        sum_dist = 0.0_dp
        valid_count = 0
        
        ! Calculate statistics (excluding error values)
        do i = 1, n_genes
            if (distances(i) >= 0.0_dp) then
                valid_count = valid_count + 1
                sum_dist = sum_dist + distances(i)
                min_dist = min(min_dist, distances(i))
                max_dist = max(max_dist, distances(i))
            end if
        end do
        
        if (valid_count > 0) then
            mean_dist = sum_dist / real(valid_count, dp)
            
            write(*,*) '=== Distance Statistics ==='
            write(*,*) 'Valid distances computed: ', valid_count, ' / ', n_genes
            write(*,*) 'Minimum distance: ', min_dist
            write(*,*) 'Maximum distance: ', max_dist
            write(*,*) 'Mean distance: ', mean_dist
            write(*,*) '=========================='
        else
            write(*,*) 'Error: No valid distances computed!'
        end if
        
    end subroutine print_distance_statistics

end module euclidean_distance_module
