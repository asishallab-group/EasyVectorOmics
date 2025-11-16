import numpy as np
import math
import sys
import os

# Add parent directory to path to import tensoromics_functions
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from tensoromics_functions import (
    tox_mask_check_state,
    tox_filter_paralogs_by_pattern_dosage_effect,
    tox_filter_paralogs_by_pattern_subfunctionalization,
    tox_calc_work_arr_paralog_subsets_size,
    tox_detect_dosage_effect,
    tox_detect_subfunctionalization,
    tox_mask_chunk_count,
    tox_detect_neofunctionalization,
    tox_normalize_unit_length
)


def test_paralog_functions():
    print("=== Testing Mask Logic ===")

    n_paralogs = 5
    i_paralog = 2
    chunk_count = tox_mask_chunk_count(n_paralogs)
    print(f"Chunk count for {n_paralogs} paralogs: {chunk_count}")

    bit_mask = np.zeros(chunk_count, dtype=np.int32)
    bit_mask[i_paralog // 32] = 1 << (i_paralog % 32)
    state = tox_mask_check_state(bit_mask, i_paralog + 1)
    print(f"Paralog {i_paralog + 1} active in mask: {state}")

    print("\n=== Testing Pattern Filtering ===")

    sorted_gene_to_fam_perm = np.arange(1, n_paralogs+1)
    perm_first_paralog_idx = 1
    angles = np.array([0.1, 0.3, 0.5, 0.7, 0.9], dtype=np.float64)
    threshold = 0.6

    dosage_mask = tox_filter_paralogs_by_pattern_dosage_effect(angles, threshold, sorted_gene_to_fam_perm, perm_first_paralog_idx, n_paralogs)
    print("Dosage effect mask:", dosage_mask)

    subfunc_mask = tox_filter_paralogs_by_pattern_subfunctionalization(angles, threshold, sorted_gene_to_fam_perm, perm_first_paralog_idx, n_paralogs)
    print("Subfunctionalization mask:", subfunc_mask)

    print("\n=== Testing Work Array Size Calculation ===")

    max_subset_size = 3
    work_size_info = tox_calc_work_arr_paralog_subsets_size(max_subset_size, n_paralogs, dosage_mask)
    print("Work array size:", work_size_info['work_array_size'])
    print("Adjusted max subset size:", work_size_info['actual_max_subset_size'])

    print("\n=== Testing Dosage Effect Detection ===")

    ancestor = np.array([1.0, 1.0], dtype=np.float64, order="F")
    paralogs = np.array([
        [1.1, 1.2, 1.3, 1.4, 1.5],
        [0.9, 0.8, 0.7, 0.6, 0.5]
    ], dtype=np.float64, order="F")

    dosage_result = tox_detect_dosage_effect(
        ancestor=ancestor,
        genes=paralogs,
        sorted_gene_to_fam_perm=sorted_gene_to_fam_perm,
        perm_first_paralog_idx=perm_first_paralog_idx,
        n_paralogs=n_paralogs,
        filtered_paralogs_mask=dosage_mask,
        max_subset_size=n_paralogs,
        gain_gamma=0.1,
        max_angle=math.pi
    )
    print(f"Dosage effect results: {dosage_result['n_results']} subsets")
    print(dosage_result['results'])

    print("\n=== Testing Subfunctionalization Detection ===")

    norms = np.sqrt(np.sum(paralogs**2, axis=0))
    sorted_perm = np.argsort(norms).astype(np.int32) + 1

    subfunc_result = tox_detect_subfunctionalization(
        ancestor=ancestor,
        genes=paralogs,
        sorted_gene_to_fam_perm=sorted_gene_to_fam_perm,
        perm_first_paralog_idx=perm_first_paralog_idx,
        n_paralogs=n_paralogs,
        rdi_threshold=0.5,
        filtered_paralogs_mask=subfunc_mask,
        max_subset_size=n_paralogs,
        paralog_norms=norms,
        sorted_paralog_norms_perm=sorted_perm
    )
    print(f"Subfunctionalization results: {subfunc_result['n_results']} subsets")
    print(subfunc_result['results'])

    print("\n=== Testing Edge Cases ===")

    try:
        empty_mask = tox_mask_chunk_count(0)
        print("Empty paralog count test: Success")
    except Exception as e:
        print("Empty paralog count test: Error")

    single_mask = tox_mask_chunk_count(1)
    print("Single paralog mask chunk count:", single_mask)

    print("✅ Paralog functions passed.")


def test_detect_neofunctionalization():
    # -------------------------------
    # Case 1: Differences below threshold → all false (all zeros)
    # -------------------------------
    ancestors = np.array([[5, 2],
                          [3, 1]], dtype=np.float64, order="F")

    # Normalize ancestors
    ancestors = np.apply_along_axis(tox_normalize_unit_length, axis=0, arr=ancestors)

    gene_to_fam = np.array([1, 2, 1], dtype=np.int32, order="F")
    thresholds = np.array([0.05, 0.05], dtype=np.float64, order="F")

    # Build genes identical to ancestors for each gene's family
    genes = np.empty((2, 3), dtype=np.float64, order="F")
    for i_gene in range(3):
        genes[:, i_gene] = ancestors[:, gene_to_fam[i_gene] - 1]  # adjust index for Python 0-based

    neofunc = tox_detect_neofunctionalization(ancestors, genes, gene_to_fam, thresholds)
    expected = np.zeros((3, 2), dtype=np.int32, order="F")
    assert np.array_equal(neofunc, expected), "Case 1 output mismatch"

    # -------------------------------
    # Case 2: Differences above threshold → some true (some ones)
    # -------------------------------
    ancestors = np.array([[5, 2],
                          [3, 1]], dtype=np.float64, order="F")

    # Normalize ancestors
    ancestors = np.apply_along_axis(tox_normalize_unit_length, axis=0, arr=ancestors)

    gene_to_fam = np.array([1, 2, 1], dtype=np.int32, order="F")
    thresholds = np.array([0.2, 0.2], dtype=np.float64, order="F")

    # Build genes offset by threshold * family index
    genes = np.empty((2, 3), dtype=np.float64, order="F")
    for i_gene in range(3):
        genes[:, i_gene] = ancestors[:, gene_to_fam[i_gene] - 1] - thresholds * gene_to_fam[i_gene]

    neofunc = tox_detect_neofunctionalization(ancestors, genes, gene_to_fam, thresholds)
    expected = np.array([[False, False],
                         [True, True],
                         [False, False]], dtype=bool, order="F")
    assert np.array_equal(neofunc, expected), "Case 2 output mismatch"

    print("✅ Neofunctionalization passed.")


def main():
    print("=================================================")
    print("    TOX PARALOG ANALYSIS PYTHON INTERFACE TESTS")
    print("=================================================")
    print()

    # Run the tests
    test_paralog_functions()
    test_detect_neofunctionalization()


if __name__ == '__main__':
    main()
