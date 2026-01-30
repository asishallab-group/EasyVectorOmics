import numpy as np
import sys
import os

# Add parent directory to path to import tensoromics_functions
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from test_helpers import run_all_tests
from tensoromics_functions import (
    tox_determine_shared_residual_range,
    tox_determine_shared_residual_range_expert,
    tox_build_residual_histograms,
    tox_compute_divergence_per_reference_point,
    tox_compute_weighted_global_divergence
)


TOL = 1e-12


def test_tox_determine_shared_residual_range():

    # ============================================================
    # Test 1 — Basic correctness with simple values
    # ============================================================
    S1 = np.array([
        [ 1,  2,  3],
        [ 4,  5,  6],
        [-7,  8,  9],
        [10, 11, 12],
    ], dtype=np.float64, order="F")

    S2 = np.array([
        [ 2, -4,  6],
        [ 8,  1,  3],
        [ 5,  7,  9],
        [ 0,  1,  2],
    ], dtype=np.float64, order="F")

    R = tox_determine_shared_residual_range(S1, S2, 95.0)
    assert abs(R - 10.85) < TOL, f"Test 1 failed: expected 10.85, got {R}"

    # ============================================================
    # Test 2 — Custom quantile (50%)
    # ============================================================
    R = tox_determine_shared_residual_range(S1, S2, 50.0)
    assert abs(R - 5.0) < TOL, f"Test 2 failed: expected 5.0, got {R}"

    # ============================================================
    # Test 3 — Quantile < 0 → error
    # ============================================================
    try:
        tox_determine_shared_residual_range(S1, S2, -1.0)
        assert False, "Test 3 failed: expected ERR_INVALID_INPUT"
    except Exception as e:
        pass

    # ============================================================
    # Test 4 — Quantile > 100 → error
    # ============================================================
    try:
        tox_determine_shared_residual_range(S1, S2, 150.0)
        assert False, "Test 4 failed: expected ERR_INVALID_INPUT"
    except Exception as e:
        pass

    # ============================================================
    # Test 5 — NaNs must be ignored
    # ============================================================
    S1 = np.array([
        [-1,  2,   3],
        [ 4,  5,   6],
        [-7,  8,  -11],
        [ 9, 10,  12],
    ], dtype=np.float64, order="F")

    S2 = S1.copy(order="F")

    S1[0, 0] = np.nan
    S2[3, 2] = np.nan

    R = tox_determine_shared_residual_range(S1, S2, 95.0)
    assert abs(R - 11.0) < TOL, f"Test 5 failed: expected 11.0, got {R}"

    # ============================================================
    # Test 6 — All zeros
    # ============================================================
    S1 = np.zeros((4, 3), dtype=np.float64, order="F")
    S2 = np.zeros((4, 3), dtype=np.float64, order="F")

    R = tox_determine_shared_residual_range(S1, S2, 95.0)
    assert abs(R - 0.0) < TOL, f"Test 6 failed: expected 0.0, got {R}"

    # ============================================================
    # Test 7 — Single residual (1×1)
    # ============================================================
    S1 = np.array([[3.0]], dtype=np.float64, order="F")
    S2 = np.array([[-4.0]], dtype=np.float64, order="F")

    R = tox_determine_shared_residual_range(S1, S2, 95.0)
    assert abs(R - 3.95) < TOL, f"Test 7 failed: expected 3.95, got {R}"


def test_tox_determine_shared_residual_range_expert():
    # Helper to build pool + permutation
    def make_pool(S1, S2):
        pool = np.abs(np.concatenate([S1.ravel(order="F"), S2.ravel(order="F")]))
        perm = np.argsort(pool) + 1
        return pool, perm

    # ============================================================
    # Test 1 — Basic correctness with simple values
    # ============================================================
    S1 = np.array([
        [ 1,  2,  3],
        [ 4,  5,  6],
        [-7,  8,  9],
        [10, 11, 12],
    ], dtype=np.float64, order="F")

    S2 = np.array([
        [ 2, -4,  6],
        [ 8,  1,  3],
        [ 5,  7,  9],
        [ 0,  1,  2],
    ], dtype=np.float64, order="F")

    pool, perm = make_pool(S1, S2)
    R = tox_determine_shared_residual_range_expert(pool, perm)
    assert abs(R - 10.85) < TOL, f"Test 1 failed: expected 10.85, got {R}"

    # ============================================================
    # Test 2 — Custom quantile (50%)
    # ============================================================
    R = tox_determine_shared_residual_range_expert(pool, perm, 50.0)
    assert abs(R - 5.0) < TOL, f"Test 2 failed: expected 5.0, got {R}"

    # ============================================================
    # Test 3 — Quantile < 0 → error
    # ============================================================
    try:
        tox_determine_shared_residual_range_expert(pool, perm, -1.0)
        assert False, "Test 3 failed: expected ERR_INVALID_INPUT"
    except Exception as e:
        pass

    # ============================================================
    # Test 4 — Quantile > 100 → error
    # ============================================================
    try:
        tox_determine_shared_residual_range_expert(pool, perm, 150.0)
        assert False, "Test 4 failed: expected ERR_INVALID_INPUT"
    except Exception as e:
        pass

    # ============================================================
    # Test 5 — NaNs must be ignored
    # ============================================================
    S1 = np.array([
        [np.nan,  2,   3],
        [ 4,      5,   6],
        [-7,      8, -11],
        [ 9,     10,  12],
    ], dtype=np.float64, order="F")

    S2 = S1.copy(order="F")
    S2[3, 2] = np.nan

    pool, perm = make_pool(S1, S2)
    R = tox_determine_shared_residual_range_expert(pool, perm, 95.0)
    assert abs(R - 11.0) < TOL, f"Test 5 failed: expected 11.0, got {R}"

    # ============================================================
    # Test 6 — All zeros
    # ============================================================
    S1 = np.zeros((4, 3), dtype=np.float64, order="F")
    S2 = np.zeros((4, 3), dtype=np.float64, order="F")

    pool, perm = make_pool(S1, S2)
    R = tox_determine_shared_residual_range_expert(pool, perm, 95.0)
    assert abs(R - 0.0) < TOL, f"Test 6 failed: expected 0.0, got {R}"

    # ============================================================
    # Test 7 — Single residual (1×1)
    # ============================================================
    S1 = np.array([[3.0]], dtype=np.float64, order="F")
    S2 = np.array([[-4.0]], dtype=np.float64, order="F")

    pool, perm = make_pool(S1, S2)
    R = tox_determine_shared_residual_range_expert(pool, perm, 95.0)
    assert abs(R - 3.95) < TOL, f"Test 7 failed: expected 3.95, got {R}"


def test_tox_build_residual_histograms():
    n_residuals = 6
    n_neighbors = 3
    n_bins = 4
    R = 2.0

    # ============================================================
    # Test 1 — Simple symmetric case, no NaNs
    # ============================================================
    E = np.zeros((n_residuals, n_neighbors), dtype=np.float64, order="F")
    E[:, 0] = [-2.0, -0.5, 0.2, 1.7, 0.9, -1.2]
    E[:, 1] = 0.0
    E[:, 2] = [2.5, -3.0, 1.2, 0.4, -0.1, 0.0]

    out = tox_build_residual_histograms(E, R, n_bins)

    counts = out["counts"]
    pmf = out["pmf"]
    included = out["included_n_residuals"]

    expected_counts = np.array([
        [2, 1, 2, 1],
        [0, 0, 6, 0],
        [1, 1, 2, 2],
    ], dtype=np.int32, order="F")

    expected_pmf = np.array([
        [2/6, 1/6, 2/6, 1/6],
        [0/6, 0/6, 6/6, 0/6],
        [1/6, 1/6, 2/6, 2/6],
    ], dtype=np.float64, order="F")

    assert np.array_equal(counts, expected_counts), "Test 1: counts mismatch"
    assert np.allclose(pmf, expected_pmf, atol=TOL), "Test 1: pmf mismatch"
    assert included[0] == 6
    assert included[1] == 6
    assert included[2] == 6

    # ============================================================
    # Test 2 — NaNs must be ignored
    # ============================================================
    E = np.zeros((n_residuals, n_neighbors), dtype=np.float64, order="F")
    E[0, 0] = np.nan
    E[2, 0] = np.nan
    E[1, 1] = np.nan
    E[5, 2] = np.nan

    out = tox_build_residual_histograms(E, R, n_bins)
    counts = out["counts"]
    pmf = out["pmf"]
    included = out["included_n_residuals"]

    expected_counts = np.array([
        [0, 0, 4, 0],
        [0, 0, 5, 0],
        [0, 0, 5, 0],
    ], dtype=np.int32, order="F")

    expected_pmf = np.array([
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
    ], dtype=np.float64, order="F")

    assert np.array_equal(counts, expected_counts), "Test 2: counts mismatch"
    assert np.allclose(pmf, expected_pmf, atol=TOL), "Test 2: pmf mismatch"
    assert included[0] == 4
    assert included[1] == 5
    assert included[2] == 5

    # ============================================================
    # Test 3 — All NaN → pmf = 0, counts = 0, included = 0
    # ============================================================
    E = np.full((n_residuals, n_neighbors), np.nan, dtype=np.float64, order="F")

    out = tox_build_residual_histograms(E, R, n_bins)
    counts = out["counts"]
    pmf = out["pmf"]
    included = out["included_n_residuals"]

    assert np.all(counts == 0), "Test 3: counts mismatch"
    assert np.all(pmf == 0), "Test 3: pmf mismatch"
    assert np.all(included == 0), "Test 3: included mismatch"

    # ============================================================
    # Test 4 — Residuals exactly on boundaries
    # ============================================================
    E = np.array([
        [-2.0, -2.0, -2.0],
        [-1.0, -1.0, -1.0],
        [ 0.0,  0.0,  0.0],
        [ 1.0,  1.0,  1.0],
        [ 2.0,  2.0,  2.0],
        [ 0.0,  0.0,  0.0],
    ], dtype=np.float64, order="F")

    out = tox_build_residual_histograms(E, R, n_bins)
    counts = out["counts"]
    pmf = out["pmf"]
    included = out["included_n_residuals"]

    expected_counts = np.array([
        [1, 1, 2, 2],
        [1, 1, 2, 2],
        [1, 1, 2, 2],
    ], dtype=np.int32, order="F")

    expected_pmf = np.array([
        [1/6, 1/6, 1/3, 1/3],
        [1/6, 1/6, 1/3, 1/3],
        [1/6, 1/6, 1/3, 1/3],
    ], dtype=np.float64, order="F")

    assert np.array_equal(counts, expected_counts), "Test 4: counts mismatch"
    assert np.allclose(pmf, expected_pmf, atol=TOL), "Test 4: pmf mismatch"
    assert np.all(included == 6), "Test 4: included mismatch"


def test_tox_compute_divergence_per_reference_point():
    n_neighbors = 3
    n_bins = 4

    # ============================================================
    # Test 1 — Identical PMFs → JSD = 0
    # ============================================================
    p = np.array([
        [0.1, 0.2, 0.3, 0.4],
        [0.25, 0.25, 0.25, 0.25],
        [1.0, 0.0, 0.0, 0.0],
    ], dtype=np.float64, order="F")

    q = p.copy(order="F")

    jsd = tox_compute_divergence_per_reference_point(p, q)
    expected = np.zeros(n_neighbors, dtype=np.float64)

    assert np.allclose(jsd, expected, atol=TOL), "Test 1 failed: identical PMFs → JSD must be zero"

    # ============================================================
    # Test 2 — Completely disjoint PMFs → JSD = log(2)
    # ============================================================
    p = np.zeros((n_neighbors, n_bins), dtype=np.float64, order="F")
    q = np.zeros((n_neighbors, n_bins), dtype=np.float64, order="F")

    p[0, 0] = 1.0
    q[0, 1] = 1.0

    jsd = tox_compute_divergence_per_reference_point(p, q)

    expected = np.zeros(n_neighbors, dtype=np.float64)
    expected[0] = np.log(2.0)

    assert np.allclose(jsd, expected, atol=TOL), "Test 2 failed: disjoint PMFs → JSD=log(2)"
    assert jsd[1] == 0.0, "Test 2 failed: row 2 must be zero"
    assert jsd[2] == 0.0, "Test 2 failed: row 3 must be zero"

    # ============================================================
    # Test 3 — Partially overlapping PMFs (analytic)
    # ============================================================
    p = np.zeros((n_neighbors, n_bins), dtype=np.float64, order="F")
    q = np.zeros((n_neighbors, n_bins), dtype=np.float64, order="F")

    p[0, :] = [0.5, 0.5, 0.0, 0.0]
    q[0, :] = [0.0, 1.0, 0.0, 0.0]

    jsd = tox_compute_divergence_per_reference_point(p, q)

    expected = np.zeros(n_neighbors, dtype=np.float64)
    expected[0] = 0.5 * (
        0.5 * np.log(2.0) +
        0.5 * np.log(2.0 / 3.0) +
        np.log(1.0 / 0.75)
    )

    assert np.allclose(jsd, expected, atol=TOL), "Test 3 failed: analytic partial-overlap JSD mismatch"

    # ============================================================
    # Test 4 — Zero-probability bins handled correctly
    # ============================================================
    p = np.zeros((n_neighbors, n_bins), dtype=np.float64, order="F")
    q = np.zeros((n_neighbors, n_bins), dtype=np.float64, order="F")

    p[0, 0] = 1.0
    q[0, 2] = 1.0

    jsd = tox_compute_divergence_per_reference_point(p, q)

    expected = np.zeros(n_neighbors, dtype=np.float64)
    expected[0] = np.log(2.0)

    assert np.allclose(jsd, expected, atol=TOL), "Test 4 failed: zero-probability bins not handled correctly"

    # ============================================================
    # Test 5 — Multiple neighbors, mixed patterns
    # ============================================================
    p = np.array([
        [0.2, 0.0, 0.0, 0.25],
        [0.3, 1.0, 0.0, 0.25],
        [0.5, 0.0, 0.25, 0.25],
    ], dtype=np.float64, order="F")

    q = np.array([
        [0.2, 0.0, 0.0, 0.25],
        [0.3, 0.0, 0.0, 0.25],
        [0.5, 1.0, 0.25, 0.25],
    ], dtype=np.float64, order="F")

    jsd = tox_compute_divergence_per_reference_point(p, q)

    expected = np.zeros(n_neighbors, dtype=np.float64)
    expected[1] = 0.5 * (1.0 * np.log(1.0 / 0.5))
    expected[2] = 0.5 * (1.0 * np.log(1.0 / 0.5))

    assert np.allclose(jsd, expected, atol=TOL), "Test 5 failed: mixed patterns JSD mismatch"


def test_tox_compute_weighted_global_divergence():
    # ============================================================
    # Test 1 — Simple case: equal sample counts → uniform weights
    # ============================================================
    jsd = np.array([0.1, 0.2, 0.3, 0.4], dtype=np.float64)
    n1  = np.array([5, 5, 5, 5], dtype=np.int32)
    n2  = np.array([5, 5, 5, 5], dtype=np.int32)

    global_jsd, w = tox_compute_weighted_global_divergence(jsd, n1, n2).values()

    expected_weights = np.array([0.25, 0.25, 0.25, 0.25], dtype=np.float64)

    assert np.allclose(w, expected_weights, atol=TOL), "Test 1 failed: uniform weights mismatch"
    assert abs(global_jsd - 0.25) < TOL, "Test 1 failed: global JSD mismatch"

    # ============================================================
    # Test 2 — Unequal sample counts → weighted average
    # ============================================================
    jsd = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64)
    n1  = np.array([10, 20, 30, 40], dtype=np.int32)
    n2  = np.array([ 0, 10, 10, 10], dtype=np.int32)

    global_jsd, w = tox_compute_weighted_global_divergence(jsd, n1, n2).values()

    expected_weights = np.array([10, 30, 40, 50], dtype=np.float64) / 130.0

    assert np.allclose(w, expected_weights, atol=TOL), "Test 2 failed: weights mismatch"

    expected = (
        1.0 * (10/130) +
        2.0 * (30/130) +
        3.0 * (40/130) +
        4.0 * (50/130)
    )

    assert abs(global_jsd - expected) < TOL, "Test 2 failed: weighted global JSD mismatch"

    # ============================================================
    # Test 3 — Some neighborhoods have zero samples → weight = 0
    # ============================================================
    jsd = np.array([0.5, 1.0, 2.0, 4.0], dtype=np.float64)
    n1  = np.array([0, 10, 0, 5], dtype=np.int32)
    n2  = np.array([0,  0, 0, 5], dtype=np.int32)

    global_jsd, w = tox_compute_weighted_global_divergence(jsd, n1, n2).values()

    expected_weights = np.array([0.0, 0.5, 0.0, 0.5], dtype=np.float64)

    assert np.allclose(w, expected_weights, atol=TOL), "Test 3 failed: zero-sample weights mismatch"
    assert abs(global_jsd - 2.5) < TOL, "Test 3 failed: global JSD mismatch"

    # ============================================================
    # Test 4 — All neighborhoods have zero samples → weights=0, global JSD=0
    # ============================================================
    jsd = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64)
    n1  = np.array([0, 0, 0, 0], dtype=np.int32)
    n2  = np.array([0, 0, 0, 0], dtype=np.int32)

    global_jsd, w = tox_compute_weighted_global_divergence(jsd, n1, n2).values()

    expected_weights = np.zeros(4, dtype=np.float64)

    assert np.allclose(w, expected_weights, atol=TOL), "Test 4 failed: all-zero weights mismatch"
    assert abs(global_jsd - 0.0) < TOL, "Test 4 failed: global JSD mismatch"

    # ============================================================
    # Test 5 — Mixed jsd values, mixed sample counts
    # ============================================================
    jsd = np.array([0.0, 0.5, 1.0, 2.0], dtype=np.float64)
    n1  = np.array([5, 0, 10, 5], dtype=np.int32)
    n2  = np.array([5, 5,  0, 5], dtype=np.int32)

    global_jsd, w = tox_compute_weighted_global_divergence(jsd, n1, n2).values()

    expected_weights = np.array([10, 5, 10, 10], dtype=np.float64) / 35.0

    assert np.allclose(w, expected_weights, atol=TOL), "Test 5 failed: mixed weights mismatch"

    expected = (
        0.0 * (10/35) +
        0.5 * ( 5/35) +
        1.0 * (10/35) +
        2.0 * (10/35)
    )

    assert abs(global_jsd - expected) < TOL, "Test 5 failed: weighted global JSD mismatch"


def main():
    run_all_tests(globals().values(), sys.argv[0])


if __name__ == "__main__":
    main()
