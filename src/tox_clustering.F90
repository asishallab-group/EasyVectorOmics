#include "macros.h"

module tox_clustering
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf, ieee_is_nan
    use tox_errors, only: is_err, set_err, set_ok, ERR_NAN_INF, ERR_INVALID_INPUT, validate_dimension_size, validate_distance_matrix, validate_all_in_range_real, validate_in_range_int
    implicit none

    integer(int32), parameter :: METHOD_UPGMA = 0
    integer(int32), parameter :: METHOD_WPGMA = 1
contains

    !> Performs k-means clustering on factor trajectories, so factor evolution over time
    pure subroutine cluster_factor_trajectories_k_means(n_clusters, trajectories, n_factors, n_samples, n_timepoints, centroids, labels, label_counts, ierr, max_iter)
        integer(int32), intent(in) :: n_clusters
            !! number (`k`) of clusters
        integer(int32), intent(in) :: n_factors
            !! number of factors
        integer(int32), intent(in) :: n_timepoints
            !! number of timepoints
        integer(int32), intent(in) :: n_samples
            !! number of samples
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
            !! matrix with data points to cluster
        real(real64), dimension(n_factors, n_clusters), intent(inout) :: centroids
            !! matrix with initial centroids of the clusters, could be random data or actual points or unassigned garbage.
            !! The centroids should be unique. This is not checked in this routine.
            !!
            !! The final values will be the final centroids of the clusters
        integer(int32), dimension(n_samples * n_timepoints), intent(out) :: labels
            !! array of labels, each index corresponds to the respective point's index, so first label is first point's label.
            !!
            !! each label is the index of its related cluster -> `1<=label<=n_clusters=k`
        integer(int32), dimension(n_clusters), intent(out) :: label_counts
            !! holds the number of points having the respective label assigned
        integer(int32), intent(out) :: ierr
            !! Error code
        integer(int32), intent(in) :: max_iter
            !! number of maximum iterations of the clustering
    
        call k_means_clustering(n_clusters, trajectories, n_samples * n_timepoints, n_factors, centroids, labels, label_counts, ierr, max_iter)
    end subroutine cluster_factor_trajectories_k_means

    !> k-means clustering algorithm:
    !|
    !| 1. Assigns each data point to one of `k` clusters whose centroid is clostest
    !| 2. Recalculates the centroids using the mean of its assigned points
    !| 3. repeat 1-2 until assignment remains unchanged
    pure subroutine k_means_clustering(n_clusters, data_points, n_points, n_dims, centroids, labels, label_counts, ierr, max_iterations)
        integer(int32), intent(in) :: n_clusters
            !! number (`k`) of clusters
        integer(int32), intent(in) :: n_points
            !! number of points to cluster
        integer(int32), intent(in) :: n_dims
            !! number of elements a point has
        real(real64), dimension(n_dims, n_points), intent(in) :: data_points
            !! matrix with data points to cluster
        real(real64), dimension(n_dims, n_clusters), intent(inout) :: centroids
            !! matrix with initial centroids of the clusters, could be random data or actual points or unassigned garbage.
            !! The centroids should be unique. This is not checked in this routine.
            !!
            !! The final values will be the final centroids of the clusters
        integer(int32), dimension(n_points), intent(out) :: labels
            !! array of labels, each index corresponds to the respective point's index, so first label is first point's label.
            !!
            !! each label is the index of its related cluster -> `1<=label<=n_clusters=k`
        integer(int32), dimension(n_clusters), intent(out) :: label_counts
            !! holds the number of points having the respective label assigned
        integer(int32), intent(out) :: ierr
            !! Error code
        integer(int32), intent(in), optional :: max_iterations
            !! number of maximum iterations of the clustering, default 300

        integer(int32) :: label, iteration, i_point, max_iter
        logical :: labels_changed

        call set_ok(ierr)

        M_DEFAULT_VAL(max_iterations, max_iter, 300_int32)

        call validate_dimension_size(n_clusters, ierr)
        call validate_dimension_size(n_points, ierr)
        call validate_dimension_size(n_dims, ierr)
        if (n_clusters > n_points) call set_err(ierr, ERR_INVALID_INPUT)

        call validate_all_in_range_real(data_points, n_points * n_dims, ierr)
        call validate_all_in_range_real(centroids, n_clusters * n_dims, ierr)

        if (is_err(ierr)) return

        iteration = 0_int32
        do while (iteration < max_iter)
            iteration = iteration + 1

            ! assign points to their closest cluster centroid
            labels_changed = .false.
            label_counts = 0
            do i_point = 1, n_points
                call k_means_assign_cluster_helper(i_point, n_clusters, data_points, n_points, n_dims, centroids, label)

                label_counts(label) = label_counts(label) + 1
                if (labels(i_point) /= label) then
                    labels(i_point) = label
                    labels_changed = .true.
                end if
            end do

            ! if the assignments did not change, the clustering is done
            if (is_err(ierr) .or. .not. labels_changed) exit

            call k_means_recompute_cluster_centroids(data_points, n_points, n_dims, centroids, n_clusters, labels, label_counts)
        end do
    end subroutine k_means_clustering

    !> Helper to recompute the centroids in k-means by taking the mean of assigned points
    pure subroutine k_means_recompute_cluster_centroids(data_points, n_points, n_dims, centroids, n_clusters, labels, label_counts)
        integer(int32), intent(in) :: n_clusters
            !! number (`k`) of clusters
        integer(int32), intent(in) :: n_points
            !! number of points to cluster
        integer(int32), intent(in) :: n_dims
            !! number of elements a point has
        real(real64), dimension(n_dims, n_points), intent(in) :: data_points
            !! matrix with data points to cluster
        real(real64), dimension(n_dims, n_clusters), intent(inout) :: centroids
            !! matrix with initial centroids of the clusters, could be random data or actual points or unassigned garbage.
            !! The centroids should be unique. This is not checked in this routine.
            !!
            !! The final values will be the final centroids of the clusters
        integer(int32), dimension(n_points), intent(out) :: labels
            !! array of labels, each index corresponds to the respective point's index, so first label is first point's label.
            !!
            !! each label is the index of its related cluster -> `1<=label<=n_clusters=k`
        integer(int32), dimension(n_clusters), intent(out) :: label_counts
            !! holds the number of points having the respective label assigned
    
        integer(int32) :: i_point, label, i_dim

        centroids = 0.0_real64
        do i_point = 1, n_points
            label = labels(i_point)
            do i_dim = 1, n_dims
                centroids(i_dim, label) = centroids(i_dim, label) + data_points(i_dim, i_point) / real(max(label_counts(label), 1_int32), real64)
            end do
        end do
    end subroutine k_means_recompute_cluster_centroids

    !> Helper routine to assign a cluster's label to a point in k-means
    pure subroutine k_means_assign_cluster_helper(i_point, n_clusters, data_points, n_points, n_dims, centroids, label)
        integer(int32), intent(in) :: i_point
            !! index of the point that should be assigned a cluster
        integer(int32), intent(in) :: n_clusters
            !! number (`k`) of clusters
        integer(int32), intent(in) :: n_points
            !! number of points to cluster
        integer(int32), intent(in) :: n_dims
            !! number of elements a point has
        real(real64), dimension(n_dims, n_points), intent(in) :: data_points
            !! matrix with data points to cluster
        real(real64), dimension(n_dims, n_clusters), intent(in) :: centroids
            !! matrix with initial centroids of the clusters, could be random data or actual points or unassigned garbage.
            !! The centroids should be unique. This is not checked in this routine.
            !!
            !! The final values will be the final centroids of the clusters
        integer(int32), intent(out) :: label
            !! assigned label

        integer(int32) :: i_cluster, i_dim
        real(real64) :: squared_dist_to_centroid, smallest_squared_dist

        label = 0_int32
        smallest_squared_dist = ieee_value(1.0_real64, ieee_positive_inf)

        ! calculate distance to each centroid and assign to closest
        do i_cluster = 1, n_clusters
            squared_dist_to_centroid = 0.0_real64
            do i_dim = 1, n_dims
                squared_dist_to_centroid = squared_dist_to_centroid + (data_points(i_dim, i_point) - centroids(i_dim, i_cluster)) ** 2
            end do
            if (squared_dist_to_centroid <= smallest_squared_dist) then
                smallest_squared_dist = squared_dist_to_centroid
                label = i_cluster
            end if
        end do
    end subroutine k_means_assign_cluster_helper

    pure subroutine upgma(distances, n_points, merge_i, merge_j, heights, cluster_sizes, ierr)
        integer(int32), intent(in) :: n_points
            !! number of points to cluster
        real(real64), dimension(n_points, n_points), intent(inout) :: distances
            !! symmetric distance matrix, holding the positive distances between points. Distance of X->X is always zero.
            !!
            !! @note
            !! This subroutine operates in-place in the bottom triangle of the distance matrix and recovers it using the top triangle once done or on error.
            !! So there is no need to copy an existing distance matrix, just pass the original.
            !! @endnote
        integer(int32), dimension(n_points - 1), intent(out) :: merge_i
            !! holds cluster labels of the merged node pair at iteration k -> positives relate to leafs/data point indices, negatives to inner nodes
        integer(int32), dimension(n_points - 1), intent(out) :: merge_j
            !! holds cluster labels of the merged node pair at iteration k -> positives relate to leafs/data point indices, negatives to inner nodes
        real(real64), dimension(n_points - 1), intent(out) :: heights
            !! height of the shorter branch of the merge, e.g. if (A,B)+(C) merges to ((A,B),C), the branch to (A,B) is shorter
        integer(int32), dimension(n_points - 1), intent(out) :: cluster_sizes
            !! size of cluster at iteration k
        integer(int32), intent(out) :: ierr
            !! Error code
    
        call hierarchical_linkage(distances, n_points, merge_i, merge_j, heights, cluster_sizes, METHOD_UPGMA, ierr)
    end subroutine upgma

    pure subroutine wpgma(distances, n_points, merge_i, merge_j, heights, cluster_sizes, ierr)
        integer(int32), intent(in) :: n_points
            !! number of points to cluster
        real(real64), dimension(n_points, n_points), intent(inout) :: distances
            !! symmetric distance matrix, holding the positive distances between points. Distance of X->X is always zero.
            !!
            !! @note
            !! This subroutine operates in-place in the bottom triangle of the distance matrix and recovers it using the top triangle once done or on error.
            !! So there is no need to copy an existing distance matrix, just pass the original.
            !! @endnote
        integer(int32), dimension(n_points - 1), intent(out) :: merge_i
            !! holds cluster labels of the merged node pair at iteration k -> positives relate to leafs/data point indices, negatives to inner nodes
        integer(int32), dimension(n_points - 1), intent(out) :: merge_j
            !! holds cluster labels of the merged node pair at iteration k -> positives relate to leafs/data point indices, negatives to inner nodes
        real(real64), dimension(n_points - 1), intent(out) :: heights
            !! height of the shorter branch of the merge, e.g. if (A,B)+(C) merges to ((A,B),C), the branch to (A,B) is shorter
        integer(int32), dimension(n_points - 1), intent(out) :: cluster_sizes
            !! size of cluster at iteration k
        integer(int32), intent(out) :: ierr
            !! Error code
    
        call hierarchical_linkage(distances, n_points, merge_i, merge_j, heights, cluster_sizes, METHOD_WPGMA, ierr)
    end subroutine wpgma

    pure subroutine hierarchical_linkage(distances, n_points, merge_i, merge_j, heights, cluster_sizes, method, ierr)
        integer(int32), intent(in) :: n_points
            !! number of points to cluster
        real(real64), dimension(n_points, n_points), intent(inout) :: distances
            !! symmetric distance matrix, holding the positive distances between points. Distance of X->X is always zero.
            !!
            !! @note
            !! This subroutine operates in-place in the bottom triangle of the distance matrix and recovers it using the top triangle once done or on error.
            !! So there is no need to copy an existing distance matrix, just pass the original.
            !! @endnote
        integer(int32), dimension(n_points - 1), intent(out) :: merge_i
            !! holds cluster labels of the merged node pair at iteration k -> positives relate to leafs/data point indices, negatives to inner nodes
        integer(int32), dimension(n_points - 1), intent(out) :: merge_j
            !! holds cluster labels of the merged node pair at iteration k -> positives relate to leafs/data point indices, negatives to inner nodes
        real(real64), dimension(n_points - 1), intent(out) :: heights
            !! height of the shorter branch of the merge, e.g. if (A,B)+(C) merges to ((A,B),C), the branch to (A,B) is shorter
        integer(int32), dimension(n_points - 1), intent(out) :: cluster_sizes
            !! size of cluster at iteration k
        integer(int32), intent(in) :: method
            !! used algorithm
            !!
            !! | Mode  | Value |
            !! |-------|-------|
            !! | UPGMA |   0   |
            !! | WPGMA |   1   |
            !!
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i, row_idx, col_idx, weight_col, weight_row, new_weight, cluster_idx
        real(real64) :: min_dist

        call set_ok(ierr)

        call validate_dimension_size(n_points, ierr)
        call validate_distance_matrix(distances, n_points, ierr)
        call validate_in_range_int(method, ierr, min=0_int32, max=1_int32)
        if (is_err(ierr)) return

        do i = 1, n_points - 1
            call get_min_distance_indices_helper(distances, n_points, row_idx, col_idx, min_dist)

            heights(i) = min_dist

            ! Get Weight and cluster index/label
            call get_cluster_data_hier_linkage_helper(distances, n_points, cluster_sizes, col_idx, weight_col, cluster_idx)
            merge_j(i) = cluster_idx
            call get_cluster_data_hier_linkage_helper(distances, n_points, cluster_sizes, row_idx, weight_row, cluster_idx)
            merge_i(i) = cluster_idx


            new_weight = weight_col + weight_row
            cluster_sizes(i) = new_weight

            select case (method)
                case (METHOD_UPGMA)
                    call merge_distances_hier_linkage_helper(distances, n_points, row_idx, col_idx, real(weight_row, real64), real(weight_col, real64), real(new_weight, real64), i)
                case (METHOD_WPGMA)
                    call merge_distances_hier_linkage_helper(distances, n_points, row_idx, col_idx, 1.0_real64, 1.0_real64, 2.0_real64, i)
            end select
        end do

        call recover_distance_matrix_helper(distances, n_points)
    end subroutine hierarchical_linkage

    pure subroutine recover_distance_matrix_helper(distances, n_points)
        integer(int32), intent(in) :: n_points
            !! number of points to cluster
        real(real64), dimension(n_points, n_points), intent(inout) :: distances
            !! top triangle holds distances, bottom triangle will be recovered

        integer(int32) :: i_row, i_col

        do i_col = 1, n_points
            ! recover self distance
            distances(i_col, i_col) = 0.0_real64
            ! recover triangle
            do i_row = i_col + 1, n_points
                distances(i_row, i_col) = distances(i_col, i_row)
            end do
        end do
    end subroutine recover_distance_matrix_helper

    pure subroutine get_cluster_data_hier_linkage_helper(distances, n_points, cluster_sizes, idx, weight, cluster_idx)
        integer(int32), intent(in) :: n_points
            !! number of points to cluster
        real(real64), dimension(n_points, n_points), intent(in) :: distances
            !! distance matrix, holding the distances between points
        integer(int32), dimension(n_points - 1), intent(in) :: cluster_sizes
            !! size of cluster at iteration k
        integer(int32), intent(in) :: idx
            !! index of cluster in `distances`
        integer(int32), intent(out) :: weight
            !! number of leafs of cluster
        integer(int32), intent(out) :: cluster_idx
            !! index/label of cluster

        integer(int32) :: iteration_k

        iteration_k = int(distances(idx, idx), int32)
        if (iteration_k == 0) then
            weight = 1
            cluster_idx = idx
        else
            weight = cluster_sizes(iteration_k)
            cluster_idx = -iteration_k
        end if   
    end subroutine get_cluster_data_hier_linkage_helper

    pure subroutine merge_distances_hier_linkage_helper(distances, n_points, row_idx, col_idx, weight_row, weight_col, new_weight, iteration_k)
        integer(int32), intent(in) :: n_points
            !! number of points to cluster
        real(real64), dimension(n_points, n_points), intent(inout) :: distances
            !! distance matrix, holding the distances between points
        integer(int32), intent(in) :: row_idx
            !! row index of min distance
        integer(int32), intent(in) :: col_idx
            !! column index of min distance
        real(real64), intent(in) :: weight_row
            !! cluster size of cluster at `row_idx`
        real(real64), intent(in) :: weight_col
            !! cluster size of cluster at `col_idx`
        real(real64), intent(in) :: new_weight
            !! cluster size of new cluster
        integer(int32), intent(in) :: iteration_k
            !! iteration of current merge

        integer(int32) :: i_node
        real(real64) :: new_dist

        ! Merge distances
        ! Update merged nodes with arithmetic mean distances
        ! Use column index for merged distances -> fill row_idx with infinity values
        ! update bottom triangle
        do i_node = 1, n_points
            if (distances(i_node, i_node) >= 0.0_real64) then
                ! i_node < col_idx < row_idx -> fill rows in both
                if (i_node < col_idx) then
                    new_dist = (distances(row_idx, i_node) * weight_row + distances(col_idx, i_node) * weight_col) / new_weight
                    distances(col_idx, i_node) = new_dist
                    distances(row_idx, i_node) = ieee_value(1.0_real64, ieee_positive_inf)
                ! col_idx <= i_node < row_idx -> fill col of col_idx, row of row_idx
                else if (i_node < row_idx) then
                    new_dist = (distances(row_idx, i_node) * weight_row + distances(i_node, col_idx) * weight_col) / new_weight
                    distances(i_node, col_idx) = new_dist
                    distances(row_idx, i_node) = ieee_value(1.0_real64, ieee_positive_inf)
                ! col_idx < row_idx < i_node -> fill columns in both
                else if (i_node > row_idx) then
                    new_dist = (distances(i_node, row_idx) * weight_row + distances(i_node, col_idx) * weight_col) / new_weight
                    distances(i_node, col_idx) = new_dist
                    ! don't fill row_idx with infinity, as whole column will be marked inactive
                end if
            end if
        end do

        ! Update index and cluster size
        distances(col_idx, col_idx) = real(iteration_k, real64)
        distances(row_idx, col_idx) = ieee_value(1.0_real64, ieee_positive_inf) ! self distance of cluster should never be minimum
        distances(row_idx, row_idx) = -1.0_real64 ! mark column inactive
    end subroutine merge_distances_hier_linkage_helper

    !> Helper routine to find the indices of the minimum value a distance matrix.
    !| 
    !| - Expects `n>1`
    !| - Ignores columns with self-distance (i,i)<0
    !| - Searches in bottom triangle (excluding self-distance (i, i))
    !| - bottom triangle -> `row_idx > col_idx`
    pure subroutine get_min_distance_indices_helper(distances, n, row_idx, col_idx, min_dist)
        integer(int32), intent(in) :: n
            !! number of rows and columns of square matrix `distances`
        real(real64), dimension(n, n), intent(in) :: distances
            !! distance matrix, holding the distances between points
        integer(int32), intent(out) :: row_idx
            !! row index of min distance
        integer(int32), intent(out) :: col_idx
            !! column index of min distance
        real(real64), intent(out) :: min_dist
            !! min distance

        integer(int32) :: i_row, i_col

        min_dist = ieee_value(1.0_real64, ieee_positive_inf)

        do i_col = 1, n - 1
            if (distances(i_col, i_col) >= 0.0_real64) then
                do i_row = i_col + 1, n
                    if (distances(i_row, i_col) < min_dist) then
                        min_dist = distances(i_row, i_col)
                        row_idx = i_row
                        col_idx = i_col
                    end if
                end do
            end if
        end do
    end subroutine get_min_distance_indices_helper
end module tox_clustering

!> k-means clustering algorithm:
!|
!| 1. Assigns each data point to one of `k` clusters whose centroid is clostest
!| 2. Recalculates the centroids using the mean of its assigned points
!| 3. repeat 1-2 until assignment remains unchanged
pure subroutine k_means_clustering_c(n_clusters, data_points, n_points, n_dims, centroids, labels, label_counts, ierr, max_iterations) bind(C, name="k_means_clustering_c")
    use tox_clustering, only: k_means_clustering
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    M_USE_NULL_VALIDATION

    integer(c_int), intent(in), target :: n_clusters
        !! number (`k`) of clusters
    integer(c_int), intent(in), target :: n_points
        !! number of points to cluster
    integer(c_int), intent(in), target :: n_dims
        !! number of elements a point has
    real(c_double), dimension(n_dims, n_points), intent(in), target :: data_points
        !! matrix with data points to cluster
    real(c_double), dimension(n_dims, n_clusters), intent(inout), target :: centroids
        !! matrix with initial centroids of the clusters, could be random data or actual points or unassigned garbage.
        !! The centroids should be unique. This is not checked in this routine.
        !!
        !! The final values will be the final centroids of the clusters
    integer(c_int), dimension(n_points), intent(out), target :: labels
        !! array of labels, each index corresponds to the respective point's index, so first label is first point's label.
        !!
        !! each label is the index of its related cluster -> `1<=label<=n_clusters=k`
    integer(c_int), dimension(n_clusters), intent(out), target :: label_counts
        !! holds the number of points having the respective label assigned
    integer(c_int), intent(out), target :: ierr
        !! Error code
    integer(c_int), intent(in), target :: max_iterations
        !! number of maximum iterations of the clustering

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_clusters)
    M_CHECK_NON_NULL(n_points)
    M_CHECK_NON_NULL(n_dims)
    M_CHECK_NON_NULL(data_points)
    M_CHECK_NON_NULL(centroids)
    M_CHECK_NON_NULL(labels)
    M_CHECK_NON_NULL(label_counts)
    M_CHECK_NON_NULL(max_iterations)

    call k_means_clustering(n_clusters, data_points, n_points, n_dims, centroids, labels, label_counts, ierr, max_iterations)
end subroutine k_means_clustering_c
