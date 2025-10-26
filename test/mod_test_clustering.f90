!> Unit test suite for tox_clustering routine.
module mod_test_tox_clustering
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use tox_clustering
    use tox_errors
    implicit none

    ! Abstract interface for all test procedures
    abstract interface
        subroutine test_interface()
        end subroutine test_interface
    end interface

    ! Type to hold test name and procedure pointer
    type :: test_case
        character(len=128) :: name
        procedure(test_interface), pointer, nopass :: test_proc => null()
    end type test_case

    real(real64), parameter :: TOL = 1d-12

contains

    !> Get array of all available tests.
    function get_all_tests() result(all_tests)
        type(test_case) :: all_tests(4)

        all_tests(1) = test_case("test_tox_clustering_merge_linkage_clusters_helper", test_merge_linkage_clusters_helper)
        all_tests(2) = test_case("test_tox_clustering_k_means_assign_cluster_helper", test_k_means_assign_cluster_helper)
        all_tests(3) = test_case("test_tox_clustering_k_means_recompute_cluster_centroids", test_k_means_recompute_cluster_centroids)
        all_tests(4) = test_case("test_tox_clustering_k_means_clustering", test_k_means_clustering)
    end function get_all_tests

    subroutine test_k_means_clustering()
        integer(int32), parameter :: n_dims = 2, n_points = 4, n_clusters = 2
        real(real64) :: data_points(n_dims, n_points)
        real(real64) :: centroids(n_dims, n_clusters), expected_centroids(n_dims, n_clusters)
        integer(int32) :: labels(n_points), label_counts(n_clusters), ierr

        ! Define 2D points
        ! Cluster 1: (1,1), (2,2)
        ! Cluster 2: (8,8), (10,10)
        data_points(:,1) = [1.0_real64, 1.0_real64]
        data_points(:,2) = [2.0_real64, 2.0_real64]
        data_points(:,3) = [8.0_real64, 8.0_real64]
        data_points(:,4) = [10.0_real64, 10.0_real64]

        ! Initialize centroids to first and last points
        centroids(:,1) = data_points(:,1)
        centroids(:,2) = data_points(:,2)

        ! Expected final centroids
        expected_centroids(:,1) = [1.5_real64, 1.5_real64]
        expected_centroids(:,2) = [9.0_real64, 9.0_real64]

        ! Call clustering routine
        call k_means_clustering(n_clusters, data_points, n_points, n_dims, centroids, labels, label_counts, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_k_means_clustering: expected OK status")

        ! Validate final centroids
        call assert_equal_array_real(centroids(:,1), expected_centroids(:,1), n_dims, TOL, "test_k_means_clustering: centroid(:,1) mismatch")
        call assert_equal_array_real(centroids(:,2), expected_centroids(:,2), n_dims, TOL, "test_k_means_clustering: centroid(:,2) mismatch")

        ! Validate label assignments
        call assert_equal_int(labels(1), 1_int32, "test_k_means_clustering: labels(1) should be 1")
        call assert_equal_int(labels(2), 1_int32, "test_k_means_clustering: labels(2) should be 1")
        call assert_equal_int(labels(3), 2_int32, "test_k_means_clustering: labels(3) should be 2")
        call assert_equal_int(labels(4), 2_int32, "test_k_means_clustering: labels(4) should be 2")

        ! Validate label counts
        call assert_equal_int(label_counts(1), 2_int32, "test_k_means_clustering: label_counts(1) mismatch")
        call assert_equal_int(label_counts(2), 2_int32, "test_k_means_clustering: label_counts(2) mismatch")

        ! Test max_iterations=1, should be different -> not clusters [1,2] [3,4], instead [1] [2,3,4]
        ! Initialize centroids to first and last points
        centroids(:,1) = data_points(:,1)
        centroids(:,2) = data_points(:,2)

        ! Expected final centroids
        expected_centroids(:,1) = [1.0_real64, 1.0_real64]
        expected_centroids(:,2) = [20.0_real64 / 3, 20.0_real64 / 3]

        ! Call clustering routine
        call k_means_clustering(n_clusters, data_points, n_points, n_dims, centroids, labels, label_counts, ierr, 1_int32)
        call assert_equal_int(ierr, ERR_OK, "test_k_means_clustering: expected OK status")

        ! Validate final centroids
        call assert_equal_array_real(centroids(:,1), expected_centroids(:,1), n_dims, TOL, "test_k_means_clustering: centroid(:,1) mismatch")
        call assert_equal_array_real(centroids(:,2), expected_centroids(:,2), n_dims, TOL, "test_k_means_clustering: centroid(:,2) mismatch")

        ! Validate label assignments
        call assert_equal_int(labels(1), 1, "test_k_means_clustering: labels(1) should be 1")
        call assert_equal_int(labels(2), 2, "test_k_means_clustering: labels(2) should be 1")
        call assert_equal_int(labels(3), 2, "test_k_means_clustering: labels(3) should be 2")
        call assert_equal_int(labels(4), 2, "test_k_means_clustering: labels(4) should be 2")

        ! Validate label counts
        call assert_equal_int(label_counts(1), 1, "test_k_means_clustering: label_counts(1) mismatch")
        call assert_equal_int(label_counts(2), 3, "test_k_means_clustering: label_counts(2) mismatch")
    end subroutine test_k_means_clustering

    subroutine test_k_means_recompute_cluster_centroids()
        integer(int32), parameter :: n_dims = 2, n_points = 4, n_clusters = 3
        real(real64) :: data_points(n_dims, n_points)
        real(real64) :: centroids(n_dims, n_clusters), expected_centroids(n_dims, n_clusters)
        integer(int32) :: labels(n_points), label_counts(n_clusters)

        ! Define 2D points
        ! Cluster 1: (1,1), (2,2)
        ! Cluster 2: (8,8), (10,10)
        ! Cluster 3: no points assigned
        data_points(:,1) = [1.0_real64, 1.0_real64]
        data_points(:,2) = [2.0_real64, 2.0_real64]
        data_points(:,3) = [8.0_real64, 8.0_real64]
        data_points(:,4) = [10.0_real64, 10.0_real64]

        ! Assign labels manually
        labels = [1, 1, 2, 2]
        label_counts = [2, 2, 0]

        ! Expected centroids:
        expected_centroids(:,1) = [1.5_real64, 1.5_real64]  ! mean of (1,1) and (2,2)
        expected_centroids(:,2) = [9.0_real64, 9.0_real64]  ! mean of (8,8) and (10,10)
        expected_centroids(:,3) = [0.0_real64, 0.0_real64]  ! no points assigned

        ! Initialize centroids to sentinel
        centroids = -999.0_real64

        ! Call routine
        call k_means_recompute_cluster_centroids(data_points, n_points, n_dims, centroids, n_clusters, labels, label_counts)

        ! Validate centroids
        call assert_equal_array_real(centroids(:,1), expected_centroids(:,1), n_dims, TOL, "test_k_means_recompute_cluster_centroids: centroid(:,1) mismatch")
        call assert_equal_array_real(centroids(:,2), expected_centroids(:,2), n_dims, TOL, "test_k_means_recompute_cluster_centroids: centroid(:,2) mismatch")
        call assert_equal_array_real(centroids(:,3), expected_centroids(:,3), n_dims, TOL, "test_k_means_recompute_cluster_centroids: centroid(:,3) should be zero")

        ! Validate label counts
        call assert_equal_int(label_counts(1), 2, "test_k_means_recompute_cluster_centroids: label_counts(1) mismatch")
        call assert_equal_int(label_counts(2), 2, "test_k_means_recompute_cluster_centroids: label_counts(2) mismatch")
        call assert_equal_int(label_counts(3), 0, "test_k_means_recompute_cluster_centroids: label_counts(3) should be zero")

    end subroutine test_k_means_recompute_cluster_centroids

    subroutine test_k_means_assign_cluster_helper()
        integer(int32), parameter :: n_dims = 2, n_points = 3, n_clusters = 2
        real(real64) :: data_points(n_dims, n_points)
        real(real64) :: centroids(n_dims, n_clusters)
        integer(int32) :: i_point, label, expected_label

        ! Define 2D points
        ! Point 1: (1.0, 1.0)
        ! Point 2: (5.0, 5.0)
        ! Point 3: (9.0, 9.0)
        data_points(:,1) = [1.0_real64, 1.0_real64]
        data_points(:,2) = [5.0_real64, 5.0_real64]
        data_points(:,3) = [9.0_real64, 9.0_real64]

        ! Define centroids
        ! Cluster 1: (0.0, 0.0)
        ! Cluster 2: (10.0, 10.0)
        centroids(:,1) = [0.0_real64, 0.0_real64]
        centroids(:,2) = [10.0_real64, 10.0_real64]

        ! Point 1 should be closer to Cluster 1
        i_point = 1
        expected_label = 1
        call k_means_assign_cluster_helper(i_point, n_clusters, data_points, n_points, n_dims, centroids, label)
        call assert_equal_int(label, expected_label, "test_k_means_assign_cluster_helper: Point 1 should be assigned to Cluster 1")

        ! Point 2 has same distance to both -> should be assigned to last
        i_point = 2
        expected_label = 2
        call k_means_assign_cluster_helper(i_point, n_clusters, data_points, n_points, n_dims, centroids, label)
        call assert_equal_int(label, expected_label, "test_k_means_assign_cluster_helper: Point 2 should be assigned to Cluster 2")

        ! Point 3 should be closer to Cluster 2
        i_point = 3
        expected_label = 2
        call k_means_assign_cluster_helper(i_point, n_clusters, data_points, n_points, n_dims, centroids, label)
        call assert_equal_int(label, expected_label, "test_k_means_assign_cluster_helper: Point 3 should be assigned to Cluster 2")

    end subroutine test_k_means_assign_cluster_helper

    subroutine test_merge_linkage_clusters_helper()
        integer(int32), parameter :: n_points = 6
        integer(int32) :: labels(n_points), expected(n_points)
        integer(int32) :: point_a, point_b, cluster_count

        ! Merge two standalone points
        labels = 0
        cluster_count = n_points
        point_a = 1
        point_b = 2
        expected = [1, 1, 0, 0, 0, 0]
        call merge_linkage_clusters_helper(labels, n_points, point_a, point_b, cluster_count)
        call assert_equal_array_int(labels, expected, n_points, "test_merge_linkage_clusters_helper: standalone merge failed")
        call assert_equal_int(cluster_count, n_points - 1, "test_merge_linkage_clusters_helper: cluster count incorrect")

        ! Merge standalone into existing cluster
        labels = [1, 1, 0, 0, 0, 0]
        cluster_count = 5
        point_a = 1
        point_b = 3
        expected = [1, 1, 1, 0, 0, 0]
        call merge_linkage_clusters_helper(labels, n_points, point_a, point_b, cluster_count)
        call assert_equal_array_int(labels, expected, n_points, "test_merge_linkage_clusters_helper: standalone into existing merge failed")
        call assert_equal_int(cluster_count, 4_int32, "test_merge_linkage_clusters_helper: cluster count incorrect")

        ! Merge two existing clusters
        labels = [3, 2, 4, 5, 1, 6]
        cluster_count = 6
        point_a = 1
        point_b = 3
        expected = [3, 2, 3, 4, 1, 5]
        call merge_linkage_clusters_helper(labels, n_points, point_a, point_b, cluster_count)
        call assert_equal_array_int(labels, expected, n_points, "test_merge_linkage_clusters_helper: standalone into existing merge failed")
        call assert_equal_int(cluster_count, 5_int32, "test_merge_linkage_clusters_helper: cluster count incorrect")

        ! Merge standalone into merged cluster
        labels = [1, 1, 1, 0, 0, 0]
        cluster_count = 4
        point_a = 1
        point_b = 4
        expected = [1, 1, 1, 1, 0, 0]
        call merge_linkage_clusters_helper(labels, n_points, point_a, point_b, cluster_count)
        call assert_equal_array_int(labels, expected, n_points, "test_merge_linkage_clusters_helper: standalone into existing merge failed")
        call assert_equal_int(cluster_count, 3, "test_merge_linkage_clusters_helper: cluster count incorrect")
    end subroutine test_merge_linkage_clusters_helper


    !> Run all tox_clustering tests.
    subroutine run_all_tests_tox_clustering
        type(test_case), allocatable :: all_tests(:)
        integer(int32) :: i

        all_tests = get_all_tests()

        do i = 1, size(all_tests)
            call all_tests(i)%test_proc()
            print "(' ',A,' passed.')", trim(all_tests(i)%name)
        end do
        print *, "All tox_clustering tests passed successfully."
    end subroutine run_all_tests_tox_clustering

    !> Run specific tox_clustering tests by name.
    subroutine run_named_tests_tox_clustering(test_names)
        character(len=*), intent(in) :: test_names(:)
        type(test_case), allocatable :: all_tests(:)
        integer(int32) :: i, j
        logical :: found

        all_tests = get_all_tests()

        do i = 1, size(test_names)
            found = .false.
            do j = 1, size(all_tests)
                if (trim(test_names(i)) == trim(all_tests(j)%name)) then
                    call all_tests(j)%test_proc()
                    print "(' ',A,' passed.')", trim(test_names(i))
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                print *, "Unknown test: ", trim(test_names(i))
            end if
        end do
    end subroutine run_named_tests_tox_clustering
end module mod_test_tox_clustering
