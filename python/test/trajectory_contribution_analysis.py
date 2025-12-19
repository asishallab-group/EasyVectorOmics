
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
    tox_compute_all_contributions,
    tox_compute_baselines_factor_dependent,
    tox_compute_velocity_trajectories,
    tox_compute_acceleration_from_velocity,
    tox_compute_velocity_acceleration_contributions,
    tox_compute_velocity_acceleration_contributions_expert,
)


# Constants
MODE_NORMAL = 1
TOL = 1e-12

def _expected_velocity(trajectories: np.ndarray) -> np.ndarray:
    trajectories = np.asarray(trajectories, dtype=np.float64)
    if trajectories.ndim != 3:
        raise ValueError("trajectories must be 3D")

    traj_f = np.transpose(trajectories, (0, 2, 1))
    velocity_f = np.zeros_like(traj_f)
    if traj_f.shape[1] > 1:
        velocity_f[:, 1:, :] = traj_f[:, 1:, :] - traj_f[:, :-1, :]
    return np.transpose(velocity_f, (0, 2, 1))


def _expected_acceleration(velocity: np.ndarray) -> np.ndarray:
    velocity = np.asarray(velocity, dtype=np.float64)
    if velocity.ndim != 3:
        raise ValueError("velocity must be 3D")

    vel_f = np.transpose(velocity, (0, 2, 1))
    acceleration_f = np.zeros_like(vel_f)
    if vel_f.shape[1] > 2:
        acceleration_f[:, 2:, :] = vel_f[:, 2:, :] - vel_f[:, 1:-1, :]
    return np.transpose(acceleration_f, (0, 2, 1))


def test_tox_compute_velocity_trajectories():
    """Test velocity computation wrapper."""
    trajectories = np.array(
        [[[1.0, 2.0, 4.0, 7.0],
          [0.0, -1.0, -1.0,  0.0]]],
        dtype=np.float64,
    )

    result = tox_compute_velocity_trajectories(trajectories)
    velocity = result["velocity"]

    expected = _expected_velocity(trajectories)

    assert velocity.shape == trajectories.shape
    assert np.allclose(velocity, expected, atol=TOL)

    print("✅ tox_compute_velocity_trajectories passed.")


def test_tox_compute_acceleration_from_velocity():
    """Test acceleration computation wrapper."""
    velocity = np.array(
        [[[0.0, 1.0, 2.0, 3.0],
          [0.0, -1.0, 0.0, 1.0]]],
        dtype=np.float64,
    )

    result = tox_compute_acceleration_from_velocity(velocity)
    acceleration = result["acceleration"]

    expected = _expected_acceleration(velocity)

    assert acceleration.shape == velocity.shape
    assert np.allclose(acceleration, expected, atol=TOL)

    print("✅ tox_compute_acceleration_from_velocity passed.")


def test_tox_compute_velocity_acceleration_contributions():
    
    """Test velocity and acceleration contribution computation wrapper."""
    trajectories = np.array(
        [[[1.0, 3.0, 6.0, 10.0],
          [1.0, 2.0, 2.0, 1.0]]],
        dtype=np.float64,
    )
    mode = MODE_NORMAL

    result = tox_compute_velocity_acceleration_contributions(trajectories, mode)

    C_vel = result["C_velocity"]
    series_vel = result["velocity_contribution_series"]
    C_acc = result["C_acceleration"]
    series_acc = result["acceleration_contribution_series"]

    assert C_vel.shape == (1, 2, 2)
    assert series_vel.shape == (1, 2, 2, 4)
    assert C_acc.shape == (1, 2, 2)
    assert series_acc.shape == (1, 2, 2, 4)

    expected_velocity = _expected_velocity(trajectories)
    expected_acceleration = _expected_acceleration(expected_velocity)

    factor_velocity = expected_velocity[0, 0, 1:]
    dependent_velocity = expected_velocity[0, 1, 1:]
    factor_acceleration = expected_acceleration[0, 0, 2:]
    dependent_acceleration = expected_acceleration[0, 1, 2:]

    expected_series_vel = np.zeros(4)
    if factor_velocity.size > 0:
        vel_contribs = factor_velocity * dependent_velocity
        expected_total_vel = vel_contribs.sum()
        expected_series_vel[1:] = vel_contribs
    else:
        expected_total_vel = 0.0

    expected_series_acc = np.zeros(4)
    if factor_acceleration.size > 0:
        acc_contribs = factor_acceleration * dependent_acceleration
        expected_total_acc = acc_contribs.sum()
        expected_series_acc[2:] = acc_contribs
    else:
        expected_total_acc = 0.0

    assert np.isclose(C_vel[0, 0, 1], expected_total_vel, atol=TOL)
    assert np.allclose(series_vel[0, 0, 1, :], expected_series_vel, atol=TOL)
    assert np.isclose(C_acc[0, 0, 1], expected_total_acc, atol=TOL)
    assert np.allclose(series_acc[0, 0, 1, :], expected_series_acc, atol=TOL)

    print("✅ tox_compute_velocity_acceleration_contributions passed.")


def test_tox_compute_velocity_acceleration_contributions_expert():
    trajectories = np.array([
        [[1.0, 2.0, 3.0, 4.0],
         [2.0, 4.0, 6.0, 8.0]]
    ], dtype=np.float64)

    result_base = tox_compute_velocity_acceleration_contributions(trajectories, MODE_NORMAL)
    result_expert = tox_compute_velocity_acceleration_contributions_expert(trajectories, MODE_NORMAL)

    for key in result_base:
        assert key in result_expert, f"Missing key {key} in expert result"
        assert result_base[key].shape == result_expert[key].shape
        assert np.allclose(result_base[key], result_expert[key], atol=TOL)

    print("✅ tox_compute_velocity_acceleration_contributions_expert passed.")

def test_tox_compute_baselines_factor_dependent():
    """Test baseline computation wrapper across all supported modes and error cases."""

    factor = np.array([1.0, 3.0, 2.0, 4.0], dtype=np.float64)
    dependent = np.array([5.0, 7.0, 6.0, 8.0], dtype=np.float64)

    # RAW mode => zero baselines
    res_raw = tox_compute_baselines_factor_dependent(factor, dependent, mode="raw")
    assert np.isclose(res_raw['baseline_factor'], 0.0, atol=TOL)
    assert np.isclose(res_raw['baseline_dependent'], 0.0, atol=TOL)

    # MIN mode => min values
    res_min = tox_compute_baselines_factor_dependent(factor, dependent, mode="min")
    assert np.isclose(res_min['baseline_factor'], np.min(factor), atol=TOL)
    assert np.isclose(res_min['baseline_dependent'], np.min(dependent), atol=TOL)

    # MEAN mode => arithmetic mean
    res_mean = tox_compute_baselines_factor_dependent(factor, dependent, mode="mean")
    assert np.isclose(res_mean['baseline_factor'], np.mean(factor), atol=TOL)
    assert np.isclose(res_mean['baseline_dependent'], np.mean(dependent), atol=TOL)

    # Mismatched lengths should raise ValueError
    try:
        tox_compute_baselines_factor_dependent(factor, dependent[:-1], mode=1)
        raise AssertionError("Expected ValueError for mismatched lengths")
    except ValueError:
        pass

    # Invalid mode should bubble up as RuntimeError from Fortran layer
    try:
        tox_compute_baselines_factor_dependent(factor, dependent, mode="unknown_mode")
        raise AssertionError("Expected RuntimeError for invalid mode")
    except RuntimeError:
        pass

    print("✅ tox_compute_baselines_factor_dependent passed all tests.")


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

    local, total = tox_compute_all_contributions(trajectories, factor_indices, dependent_indices, mode="mean").values()

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
    trajectories[1, 0, :] = [5.0, 3.0, 5.0]

    factor_indices = np.array([1], dtype=np.int32, order="F")
    dependent_indices = np.array([2], dtype=np.int32, order="F")

    local, total = tox_compute_all_contributions(trajectories, factor_indices, dependent_indices, mode="min").values()

    # Baselines: min(factor)=2.0, min(dependent)=1.0
    expected_local = np.array([0.0, 0.0, 8.0], dtype=np.float64, order="F")
    expected_total = 8.0

    assert np.allclose(local[:, 0, 0, 0], expected_local, atol=TOL), "Case 2 local contributions mismatch"
    assert abs(total[0, 0, 0] - expected_total) < TOL, "Case 2 total contribution mismatch"

    print("✅ Compute all contributions test passed.")


def main():
    print("=================================================")
    print("    TRAJECTORY CONTRIBUTION ANALYSIS PYTHON INTERFACE TESTS")
    print("=================================================")
    print()

    test_tox_compute_baselines_factor_dependent()
    test_compute_contributions()
    test_compute_all_contributions()
    test_tox_compute_velocity_trajectories()
    test_tox_compute_acceleration_from_velocity()
    test_tox_compute_velocity_acceleration_contributions()
    test_tox_compute_velocity_acceleration_contributions_expert()

if __name__ == "__main__":
    main()
