"""
Test script for clock hand angle functions
Python equivalent of the R and Fortran clock hand angle tests
"""

import numpy as np
import sys
import os
from math import pi as PI

# Add parent directory to path to import tensoromics_functions
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensoromics_functions import (
    tox_compute_contributions,
    tox_compute_all_contributions
)

# Constants
TOL = 1e-12


def test_compute_contributions():
    # -------------------------------
    # Case 1: RAW baseline
    # -------------------------------
    factor = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64)
    dependent = np.array([2.0, 1.0, 0.0, -1.0], dtype=np.float64)
    local, total = tox_compute_contributions(factor, dependent, mode="raw").values()

    expected_local = factor * dependent
    expected_total = sum(expected_local)
    assert np.allclose(local, expected_local, atol=TOL), "Case 1 local contributions mismatch"
    assert abs(total - expected_total) < TOL, "Case 1 total contribution mismatch"

    # -------------------------------
    # Case 2: MIN baseline
    # -------------------------------
    factor = np.array([3.0, 5.0, 2.0, 4.0], dtype=np.float64)
    dependent = np.array([1.0, 2.0, 0.0, -1.0], dtype=np.float64)
    local, total = tox_compute_contributions(factor, dependent, mode="min").values()

    expected_local = (factor - factor.min()) * (dependent - dependent.min())
    expected_total = sum(expected_local)
    assert np.allclose(local, expected_local, atol=TOL), "Case 2 local contributions mismatch"
    assert abs(total - expected_total) < TOL, "Case 2 total contribution mismatch"

    # -------------------------------
    # Case 3: MEAN baseline
    # -------------------------------
    factor = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64)
    dependent = np.array([4.0, 3.0, 2.0, 1.0], dtype=np.float64)
    local, total = tox_compute_contributions(factor, dependent, mode="mean").values()

    expected_local = (factor - factor.mean()) * (dependent - dependent.mean())
    expected_total = sum(expected_local)
    assert np.allclose(local, expected_local, atol=TOL), "Case 3 local contributions mismatch"
    assert abs(total - expected_total) < TOL, "Case 3 total contribution mismatch"

    print("✅ Compute contributions test passed.")


def test_compute_all_contributions():
    # -------------------------------
    # Case 1: MEAN baseline
    # -------------------------------
    # Factor trajectory: [1,2,3]
    # Dependent trajectory: [4,5,6]
    trajectories = np.empty((2, 1, 3), dtype=np.float64, order="F")
    trajectories[0, 0, :] = [1.0, 2.0, 3.0]
    trajectories[1, 0, :] = [4.0, 5.0, 6.0]

    factor_indices = np.array([1], dtype=np.int32, order="F")
    dependent_indices = np.array([2], dtype=np.int32, order="F")

    result = tox_compute_all_contributions(trajectories, factor_indices, dependent_indices, mode="mean")
    local = result["local_contributions"]
    total = result["total_contributions"]

    # Baselines: mean(factor)=2.0, mean(dependent)=5.0
    expected_local = np.array([1.0, 0.0, 1.0], dtype=np.float64, order="F")
    expected_total = 2.0

    assert np.allclose(local[:, 0, 0, 0], expected_local, atol=TOL), "Case 1 local contributions mismatch"
    assert abs(total[0, 0, 0] - expected_total) < TOL, "Case 1 total contribution mismatch"

    # -------------------------------
    # Case 2: MIN baseline
    # -------------------------------
    # Factor trajectory: [2,4,6]
    # Dependent trajectory: [1,3,5]
    trajectories[0, 0, :] = [2.0, 4.0, 6.0]
    trajectories[1, 0, :] = [1.0, 3.0, 5.0]

    factor_indices = np.array([1], dtype=np.int32, order="F")
    dependent_indices = np.array([2], dtype=np.int32, order="F")

    result = tox_compute_all_contributions(trajectories, factor_indices, dependent_indices, mode="min")
    local = result["local_contributions"]
    total = result["total_contributions"]

    # Baselines: min(factor)=2.0, min(dependent)=1.0
    expected_local = np.array([0.0, 4.0, 16.0], dtype=np.float64, order="F")
    expected_total = 20.0

    assert np.allclose(local[:, 0, 0, 0], expected_local, atol=TOL), "Case 2 local contributions mismatch"
    assert abs(total[0, 0, 0] - expected_total) < TOL, "Case 2 total contribution mismatch"

    print("✅ Compute all contributions test passed.")


def main():
    print("=================================================")
    print("    TRAJECTORY CONTRIBUTION ANALYSIS PYTHON INTERFACE TESTS")
    print("=================================================")
    print()

    test_compute_contributions()
    test_compute_all_contributions()


if __name__ == "__main__":
    main()
