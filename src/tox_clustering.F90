#include "macros.h"

module tox_clustering
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_positive_inf, ieee_value
    use tox_errors
    implicit none
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

        if (max_iter < 1) call set_err(ierr, ERR_INVALID_INPUT)
        if (is_err(ierr)) return

        iteration = 0_int32
        do while (iteration < max_iter)
            iteration = iteration + 1

            ! assign points to their closest cluster centroid
            labels_changed = .false.
            label_counts = 0
            do i_point = 1, n_points
                
                call k_means_assign_cluster_helper(i_point, n_clusters, data_points, n_points, n_dims, centroids, label)

                if (.not. label > 0) then
                    call set_err(ierr, ERR_NAN_INF)
                    exit
                end if

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

    pure subroutine maximum_linkage_clustering(n_clusters, data_points, n_points, n_dims, distances, sorted_distances_perm, labels, ierr)
        integer(int32), intent(in) :: n_clusters
            !! number of clusters
        integer(int32), intent(in) :: n_points
            !! number of points to cluster
        integer(int32), intent(in) :: n_dims
            !! number of elements a point has
        real(real64), dimension(n_dims, n_points), intent(in) :: data_points
            !! matrix with data points to cluster
        real(real64), dimension(n_points, n_points), intent(in) :: distances
            !! distance matrix, holding the distances between points
        integer(int32), dimension(n_points * n_points), intent(in) :: sorted_distances_perm
            !! permutation vector that sorts distances as flat array
        integer(int32), dimension(n_points), intent(out) :: labels
            !! array of labels, each index corresponds to the respective point's index, so first label is first point's label.
            !!
            !! each label is the index of its related cluster -> `1<=label<=n_clusters=k`
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i_point, i_pair, cluster_count, i_point_a, i_point_b

        call set_ok(ierr)

        call validate_dimension_size(n_clusters, ierr)
        call validate_dimension_size(n_points, ierr)
        call validate_dimension_size(n_dims, ierr)
        if (n_clusters > n_points) call set_err(ierr, ERR_INVALID_INPUT)
        if (is_err(ierr)) return

        labels = 0_int32

        cluster_count = n_points
        i_pair = n_points * n_points
        do while (cluster_count > n_clusters)
            i_point_a = mod(sorted_distances_perm(i_pair) - 1, n_points) + 1
            i_point_b = (sorted_distances_perm(i_pair) - 1) / n_points + 1

            call merge_linkage_clusters_helper(labels, n_points, i_point_a, i_point_b, cluster_count)

            i_pair = i_pair - 1
        end do

        call assign_remaining_linkage_clusters_helper(labels, n_points, n_points - n_clusters)
    end subroutine maximum_linkage_clustering

    !> Helper to merge clusters in linkage hierarchy algorithms
    pure subroutine merge_linkage_clusters_helper(labels, n_points, i_point_a, i_point_b, cluster_count)
        integer(int32), intent(in) :: n_points
            !! number of points to cluster
        integer(int32), dimension(n_points), intent(inout) :: labels
            !! array of labels `0<=label<=n_points`, if `label==0` it is a cluster with only a single point
        integer(int32), intent(in) :: i_point_a
            !! index to label of cluster A in `labels`
        integer(int32), intent(in) :: i_point_b
            !! index to label of cluster B in `labels`
        integer(int32), intent(inout) :: cluster_count
            !! number of existing clusters in labels (zeros are separate clusters)

        integer(int32) :: label_a, label_b, i_point, lower_label, upper_label

        label_a = labels(i_point_a)
        label_b = labels(i_point_b)

        ! if point a is not standalone cluster
        if (label_a > 0) then
            ! b is standalone cluster -> add b to cluster of a
            if (label_b == 0) then
                labels(i_point_b) = label_a
            ! b is not standalone -> merge a and b clusters
            else
                lower_label = min(label_a, label_b)
                upper_label = max(label_a, label_b)
                do i_point = 1, n_points
                    if ((labels(i_point) == label_a) .or. (labels(i_point) == label_b)) then
                        labels(i_point) = lower_label
                    else if (labels(i_point) > upper_label) then
                        labels(i_point) = labels(i_point) - 1
                    end if
                end do
            end if
        ! if point a is a standalone cluster, but b is not, add a to cluster of b
        else if (label_b > 0) then
            labels(i_point_a) = label_b
        ! if both points are standalone clusters, merge them
        else
            labels(i_point_a) = n_points - cluster_count + 1
            labels(i_point_b) = n_points - cluster_count + 1
        end if

        cluster_count = cluster_count - 1
    end subroutine merge_linkage_clusters_helper

    !> Helper to assign labels to remaining points where label is 0 yet (standalone cluster with only one point)
    pure subroutine assign_remaining_linkage_clusters_helper(labels, n_points, max_label)
        integer(int32), intent(in) :: n_points
            !! number of points to cluster
        integer(int32), dimension(n_points), intent(inout) :: labels
            !! array of labels `0<=label<=n_points`, if `label==0` it is a cluster with only a single point
        integer(int32), intent(in) :: max_label
            !! max label value in `labels` -> each 0 gets a value above max

        integer(int32) :: i_point, label

        label = max_label + 1
        do i_point = 1, n_points
            if (labels(i_point) == 0) then
                labels(i_point) = label
                label = label + 1
            end if
        end do
    end subroutine assign_remaining_linkage_clusters_helper

end module tox_clustering