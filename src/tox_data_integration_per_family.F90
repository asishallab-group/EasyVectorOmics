#include "macros.h"

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) JSD Calculation per family
!|
!| This module implements subroutines to obtain the JSD value from neighborhood residuals obtained from [[tox_data_integration_preprocessing(submodule)]]
!| for specific sub-neighborhoods by using the pipeline from [[tox_data_integration_jsd(submodule)]].
module tox_data_integration_per_family
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_data_integration, only: jct_compute_jsd_pipeline_helper
    use tox_errors, only: set_ok, validate_dimension_size, validate_in_range_int, validate_all_in_range_int, is_err, ERR_ALLOC_FAIL, set_err, validate_in_range_real
    implicit none

contains

    !> Computes the family-level compatibility score `global_js_divergence` between two studies for a single gene family (`family_idx`), by reusing the same conditioning-on-mean-expression pipeline as the global gJCT, but restricting residual samples to genes belonging to the specified family
    pure subroutine fjct_compute_jsd_alloc(family_idx, gene_to_family_S1, gene_to_family_S2, n_genes_S1, n_genes_S2, neighborhood_residuals_S1, neighborhood_residuals_S2, &
            neighborhood_genes_S1, neighborhood_genes_S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, resid_S1, resid_S2, n_bins, shared_residual_range, js_divergences, &
            included_n_reps_S1, included_n_reps_S2, global_js_divergence, weights, ierr &
        )
        integer(int32), intent(in) :: n_genes_S1
            !! Number of genes in study 1
        integer(int32), intent(in) :: n_genes_S2
            !! Number of genes in study 2
        integer(int32), intent(in) :: n_reps_S1
            !! Number of replicates in study 1
        integer(int32), intent(in) :: n_reps_S2
            !! Number of replicates in study 2
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors in the studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        integer(int32), intent(in) :: family_idx
            !! Index of the family that should be analyzed
        integer(int32), dimension(n_genes_S1), intent(in) :: gene_to_family_S1
            !! Mapping for study 1: Each index (gene) holds the index of its family
        integer(int32), dimension(n_genes_S2), intent(in) :: gene_to_family_S2
            !! Mapping for study 2: Each index (gene) holds the index of its family
        real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
        real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
        integer(int32), dimension(n_neighbors, n_points), intent(in) :: neighborhood_genes_S1
            !! Indices of selected neighborhood genes, obtained from `neighborhood_indices` of [[tox_data_integration(module):construct_neighborhoods(interface)]]
        integer(int32), dimension(n_neighbors, n_points), intent(in) :: neighborhood_genes_S2
            !! Indices of selected neighborhood genes, obtained from `neighborhood_indices` of [[tox_data_integration(module):construct_neighborhoods(interface)]]
        real(real64), dimension(n_reps_S1, n_genes_S1), intent(in) :: resid_S1
            !! Matrix of signed residuals for study 1, from [[tox_data_integration(module):compute_residuals(interface)]]
        real(real64), dimension(n_reps_S2, n_genes_S2), intent(in) :: resid_S2
            !! Matrix of signed residuals for study 1, from [[tox_data_integration(module):compute_residuals(interface)]]
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins used for the studies in [[tox_data_integration(module):build_residual_histograms(interface)]]
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range for both studies, from [[tox_data_integration(module):determine_shared_residual_range(interface)]]
        real(real64), dimension(n_points), intent(out) :: js_divergences
            !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
        integer(int32), dimension(n_points), intent(out) :: included_n_reps_S1
            !! Count of non-NaN residuals (included ones) in study 1 (obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
        integer(int32), dimension(n_points), intent(out) :: included_n_reps_S2
            !! Count of non-NaN residuals (included ones) in study 2 (obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
        real(real64), intent(out) :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        real(real64), dimension(n_points), intent(out) :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i_point, i_neighbor
        logical, dimension(:, :), allocatable :: neighbor_mask_S1, neighbor_mask_S2
        real(real64), dimension(:, :), allocatable :: pmf_S1, pmf_S2
        integer(int32), dimension(:, :), allocatable :: tmp_counts

        call set_ok(ierr)

        call validate_dimension_size(n_neighbors, ierr)
        call validate_in_range_int(family_idx, ierr, min=1_int32)
        call validate_all_in_range_int(gene_to_family_S1, n_genes_S1, ierr, min=0_int32)
        call validate_all_in_range_int(gene_to_family_S2, n_genes_S2, ierr, min=0_int32)
        call validate_all_in_range_int(neighborhood_genes_S1, size(neighborhood_genes_S1, kind=int32), ierr, min=1_int32, max=n_genes_S1)
        call validate_all_in_range_int(neighborhood_genes_S2, size(neighborhood_genes_S2, kind=int32), ierr, min=1_int32, max=n_genes_S2)

        if (is_err(ierr)) return

        M_ALLOCATE(neighbor_mask_S1(n_neighbors, n_points))
        M_ALLOCATE(neighbor_mask_S2(n_neighbors, n_points))
        M_ALLOCATE(pmf_S1(n_points, n_bins))
        M_ALLOCATE(pmf_S2(n_points, n_bins))
        M_ALLOCATE(tmp_counts(n_points, n_bins))

        ! Set up mask for filtered analysis -> only include neighbors being part of the family
        do concurrent (i_point = 1:n_points)
            do concurrent (i_neighbor = 1:n_neighbors) shared(family_idx, neighbor_mask_S1, neighbor_mask_S2, gene_to_family_S1, neighborhood_genes_S1, gene_to_family_S2, neighborhood_genes_S2)
                neighbor_mask_S1(i_neighbor, i_point) = gene_to_family_S1(neighborhood_genes_S1(i_neighbor, i_point)) == family_idx
                neighbor_mask_S2(i_neighbor, i_point) = gene_to_family_S2(neighborhood_genes_S2(i_neighbor, i_point)) == family_idx
            end do
        end do

        call fjct_compute_jsd(neighborhood_residuals_S1, neighborhood_residuals_S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, neighbor_mask_S1, neighbor_mask_S2, n_bins, shared_residual_range, js_divergences, included_n_reps_S1, included_n_reps_S2, global_js_divergence, weights, pmf_S1, pmf_S2, tmp_counts, ierr)
    end subroutine fjct_compute_jsd_alloc

    !> TODO: #118
    pure subroutine fjct_compute_jsd(neighborhood_residuals_S1, neighborhood_residuals_S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, neighbor_mask_S1, neighbor_mask_S2, n_bins, shared_residual_range, js_divergences, included_n_reps_S1, included_n_reps_S2, global_js_divergence, weights, pmf_S1, pmf_S2, tmp_counts, ierr)
        integer(int32), intent(in) :: n_reps_S1
            !! Number of replicates in study 1
        integer(int32), intent(in) :: n_reps_S2
            !! Number of replicates in study 2
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors in the studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
        real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
        logical, dimension(n_neighbors, n_points), intent(in) :: neighbor_mask_S1
            !! Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
        logical, dimension(n_neighbors, n_points), intent(in) :: neighbor_mask_S2
            !! Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins used for the studies in [[tox_data_integration(module):build_residual_histograms(interface)]]
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range for both studies, from [[tox_data_integration(module):determine_shared_residual_range(interface)]]
        real(real64), dimension(n_points), intent(out) :: js_divergences
            !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
        integer(int32), dimension(n_points), intent(out) :: included_n_reps_S1
            !! Count of non-NaN residuals (included ones) in study 1 (obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
        integer(int32), dimension(n_points), intent(out) :: included_n_reps_S2
            !! Count of non-NaN residuals (included ones) in study 2 (obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
        real(real64), intent(out) :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        real(real64), dimension(n_points), intent(out) :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
        real(real64), dimension(n_points, n_bins), intent(out) :: pmf_S1
            !! Absolute counts of a residual per bin obtained from [[tox_data_integration(module):build_residual_histograms(interface)]]
        real(real64), dimension(n_points, n_bins), intent(out) :: pmf_S2
            !! Absolute counts of a residual per bin obtained from [[tox_data_integration(module):build_residual_histograms(interface)]]
        integer(int32), dimension(n_points, n_bins), intent(out) :: tmp_counts
            !! Working array for [[tox_data_integration(module):build_residual_histograms(interface)]]
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_reps_S1, ierr)
        call validate_dimension_size(n_reps_S2, ierr)
        call validate_dimension_size(n_neighbors, ierr)
        call validate_dimension_size(n_points, ierr)
        call validate_in_range_int(n_bins, ierr, min=1_int32)
        call validate_in_range_real(shared_residual_range, ierr, min=0.0_real64)

        if (is_err(ierr)) return

        call jct_compute_jsd_pipeline_helper(neighborhood_residuals_S1, neighborhood_residuals_S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, n_bins, shared_residual_range, js_divergences, included_n_reps_S1, included_n_reps_S2, global_js_divergence, weights, pmf_S1, pmf_S2, tmp_counts, neighbor_mask_S1, neighbor_mask_S2)
    end subroutine fjct_compute_jsd

end module tox_data_integration_per_family