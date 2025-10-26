"""
Test script for clock hand angle functions
Python equivalent of the R and Fortran clock hand angle tests
"""

import numpy as np
import sys
import os

# Add parent directory to path to import tensoromics_functions
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensoromics_functions import tox_k_means_clustering, tox_k_means_clustering


def test_tox_k_means_clustering():
    """
    Test tox_k_means_clustering using known 2D points and initial centroids.
    """

    # Define 2D points: two clusters
    # Cluster 1: (1,1), (2,2)
    # Cluster 2: (8,8), (10,10)
    data_points = np.array([
        [1.0, 2.0, 8.0, 10.0],  # x
        [1.0, 2.0, 8.0, 10.0]   # y
    ], dtype=np.float64)

    # Initial centroids: first and last points
    centroids_init = np.array([
        [1.0, 10.0],
        [1.0, 10.0]
    ], dtype=np.float64)

    # Expected final centroids
    expected_centroids = np.array([
        [1.5, 9.0],
        [1.5, 9.0]
    ], dtype=np.float64)

    # Expected labels and counts
    expected_labels = np.array([1, 1, 2, 2], dtype=np.int32)
    expected_counts = np.array([2, 2], dtype=np.int32)

    # Run clustering
    centroids, labels, label_counts = tox_k_means_clustering(data_points, centroids_init, max_iter=10).values()

    # Validate centroids
    assert np.allclose(centroids, expected_centroids, atol=1e-12), "tox_k_means_clustering: centroid mismatch"

    # Validate labels
    assert np.array_equal(labels, expected_labels), "tox_k_means_clustering: label assignment mismatch"

    # Validate label counts
    assert np.array_equal(label_counts, expected_counts), "tox_k_means_clustering: label_counts mismatch"

    print("✅ tox_k_means_clustering passed all tests.")


def main():
    print("=================================================")
    print("    CLUSTERING PYTHON INTERFACE TESTS")
    print("=================================================")
    print()

    test_tox_k_means_clustering()


if __name__ == "__main__":
    main()
