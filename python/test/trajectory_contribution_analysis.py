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
from tensoromics_functions import tox_trajectory_contribution, tox_spike_contribution

# Constants
TOL = 1e-12


def test_tox_trajectory_contribution():
    """
    Test tox_trajectory_contribution for both modes using known vectors.
    """
    import numpy as np

    # Aligned vectors
    factor = np.array([1.0, 2.0, 3.0])
    dependent = np.array([2.0, 4.0, 6.0])
    dot_prod = np.dot(factor, dependent)
    magnitude = np.sqrt(np.sum(factor**2)) * np.sqrt(np.sum(dependent**2))

    # MODE_NORMAL
    expected_normal = dot_prod / magnitude
    result_normal = tox_trajectory_contribution(factor, dependent, mode=1)
    assert np.isclose(result_normal, expected_normal, atol=1e-12), f"MODE_NORMAL failed: {result_normal} vs {expected_normal}"

    # MODE_RAP
    expected_rap = np.arccos(expected_normal)
    result_rap = tox_trajectory_contribution(factor, dependent, mode=2)
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
    result_normal = tox_spike_contribution(factor, dependent, mode=1)
    assert np.allclose(result_normal, expected_normal, atol=1e-12), f"MODE_NORMAL failed: {result_normal} vs {expected_normal}"

    # MODE_RAP
    expected_rap = np.arccos(expected_normal)
    result_rap = tox_spike_contribution(factor, dependent, mode=2)
    assert np.allclose(result_rap, expected_rap, atol=1e-12), f"MODE_RAP failed: {result_rap} vs {expected_rap}"

    print("✅ tox_spike_contribution passed all tests.")


def main():
    print("=================================================")
    print("    TRAJECTORY CONTRIBUTION ANALYSIS PYTHON INTERFACE TESTS")
    print("=================================================")
    print()

    test_tox_spike_contribution()
    test_tox_trajectory_contribution()


if __name__ == "__main__":
    main()
