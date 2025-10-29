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
    tox_trajectory_contribution,
    tox_spike_contribution,
    tox_calc_contributions,
    tox_calc_contributions_expert
)

# Constants
TOL = 1e-12
MODE_NORMAL = 1
MODE_RAP = 2


def test_tox_trajectory_contribution():
    """
    Test tox_trajectory_contribution for both modes using known vectors.
    """

    # Aligned vectors
    factor = np.array([1.0, 2.0, 3.0])
    dependent = np.array([2.0, 4.0, 6.0])
    dot_prod = np.dot(factor, dependent)
    magnitude = np.sqrt(np.sum(factor**2)) * np.sqrt(np.sum(dependent**2))

    # MODE_NORMAL
    expected_normal = dot_prod / magnitude
    result_normal = tox_trajectory_contribution(factor, dependent, mode=MODE_NORMAL)
    assert np.isclose(result_normal, expected_normal, atol=1e-12), f"MODE_NORMAL failed: {result_normal} vs {expected_normal}"

    # MODE_RAP
    expected_rap = np.arccos(expected_normal)
    result_rap = tox_trajectory_contribution(factor, dependent, mode=MODE_RAP)
    assert np.isclose(result_rap, expected_rap, atol=1e-12), f"MODE_RAP failed: {result_rap} vs {expected_rap}"

    print("✅ tox_trajectory_contribution passed all tests.")


def test_tox_spike_contribution():
    """
    Test tox_spike_contribution for both modes using known vectors.
    """

    # Aligned vectors
    factor = np.array([1.0, 2.0, 3.0])
    dependent = np.array([2.0, 4.0, 6.0])
    magnitude = np.sqrt(np.sum(factor**2)) * np.sqrt(np.sum(dependent**2))

    # MODE_NORMAL
    expected_normal = (factor * dependent) / magnitude
    result_normal = tox_spike_contribution(factor, dependent, mode=MODE_NORMAL)
    assert np.allclose(result_normal, expected_normal, atol=1e-12), f"MODE_NORMAL failed: {result_normal} vs {expected_normal}"

    # MODE_RAP
    expected_rap = np.arccos(expected_normal)
    result_rap = tox_spike_contribution(factor, dependent, mode=MODE_RAP)
    assert np.allclose(result_rap, expected_rap, atol=1e-12), f"MODE_RAP failed: {result_rap} vs {expected_rap}"

    print("✅ tox_spike_contribution passed all tests.")


def test_tox_calc_contributions():
    """
    Test both tox_calc_contributions and tox_calc_contributions_expert
    using a reproducible trajectory tensor with known alignment.
    """

    # Setup: T[i,j,k] = 100*i + 10*j + k
    n_factors, n_samples, n_timepoints = 2, 3, 4
    T = np.empty((n_factors, n_samples, n_timepoints), dtype=np.float64)
    for i in range(n_factors):
        for j in range(n_samples):
            for k in range(n_timepoints):
                T[i, j, k] = 100.0*(i+1) + 10.0*(j+1) + (k+1)

    # Aligned vectors: factor=1, dependent=0 (0-based)
    i_factor = 1
    dependent_idx = 0
    mode = MODE_NORMAL

    # Reference computation
    expected_spike = np.empty((n_timepoints, n_samples), dtype=np.float64)
    expected_integrated = np.empty(n_samples, dtype=np.float64)
    for j in range(n_samples):
        f = T[i_factor, j, :]
        d = T[dependent_idx, j, :]
        mag = np.sqrt(np.sum(f**2)) * np.sqrt(np.sum(d**2))
        expected_spike[:, j] = (f * d) / mag
        expected_integrated[j] = np.sum(f * d) / mag

    # Test calc_contributions_c
    res = tox_calc_contributions(T, i_factor, dependent_idx, mode)
    assert np.allclose(res["spikes"], expected_spike, atol=1e-12), "tox_calc_contributions: spike mismatch"
    assert np.allclose(res["trajectory"], expected_integrated, atol=1e-12), "tox_calc_contributions: integrated mismatch"

    # Test calc_contributions_expert_c
    temp_factor_vector = np.empty(n_timepoints, dtype=np.float64)
    temp_dependent_vector = np.empty(n_timepoints, dtype=np.float64)
    res_expert = tox_calc_contributions_expert(
        T, i_factor, dependent_idx, mode, temp_factor_vector, temp_dependent_vector
    )

    assert np.allclose(res_expert["spikes"], expected_spike, atol=1e-12), "tox_calc_contributions_expert: spike mismatch"
    assert np.allclose(res_expert["trajectory"], expected_integrated, atol=1e-12), "tox_calc_contributions_expert: integrated mismatch"

    print("✅ Both tox_calc_contributions wrappers passed all tests.")


def main():
    print("=================================================")
    print("    TRAJECTORY CONTRIBUTION ANALYSIS PYTHON INTERFACE TESTS")
    print("=================================================")
    print()

    test_tox_spike_contribution()
    test_tox_trajectory_contribution()
    test_tox_calc_contributions()


if __name__ == "__main__":
    main()
