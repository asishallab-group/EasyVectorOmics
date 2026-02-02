!> In multi-study omics analyses, it is often unclear whether biological replicates originating from different studies can be safely treated as sampling the same biological condition.
!| Even when studies nominally target the same tissue and condition, differences in sample handling, sequencing technologies, preprocessing pipelines, or cohort, 
!| composition can introduce batch effects that are not easily detectable from mean expression levels alone.
!|
!| This ambiguity has direct consequences for downstream analyses in Tensor Omics. Integrating incompatible replicate sets can:
!|  - distort expression spaces,
!|  - affect distance-based analyses,
!|  - bias machine learning models,
!| while unnecessarily separating compatible datasets reduces statistical power.
!|
!| To address this, we introduce a Jensen–Shannon-Divergence based compatibility test (JSD-Comp-Test)
!| that empirically evaluates whether two sets of biological replicates exhibit comparable replicate-level variability.
!| Rather than comparing mean expression values, the method focuses on the distribution of signed residuals (replicate deviations from the gene-wise mean),
!| conditioned on mean expression levels to account for heteroscedasticity, which is a well-known property of omics data.
!|
!| The goal of this issue is to define, implement, and validate this compatibility test as a diagnostic tool that can be applied prior to data integration.
!| The test is intended to support principled decisions
!| on whether replicate sets from different studies should be merged or treated as distinct conditions within Tensor Omics workflows.
module tox_data_integration
    use iso_fortran_env, only: int32, real64
    implicit none

    interface compute_gene_means
        !> Compute per-gene mean expression, ignoring NaN values
        pure module subroutine compute_gene_means(n_genes, n_reps, expr, means, ierr)
            import
            integer(int32), intent(in) :: n_genes
                !! Number of genes in the study
            integer(int32), intent(in) :: n_reps
                !! Number of biological replicates in the study
            real(real64), intent(in) :: expr(n_reps, n_genes)
                !! Expression matrix
            real(real64), intent(out) :: means(n_genes)
                !! Per-gene mean expression values
            integer(int32), intent(out) :: ierr
                !! Error code
        end subroutine compute_gene_means
    end interface compute_gene_means

    interface compute_gene_means_helper
        !> (no input validation) Compute per-gene mean expression, ignoring NaN values
        pure module subroutine compute_gene_means_helper(n_genes, n_reps, expr, means)
            import
            integer(int32), intent(in) :: n_genes
                !! Number of genes in the study
            integer(int32), intent(in) :: n_reps
                !! Number of biological replicates in the study
            real(real64), intent(in) :: expr(n_reps, n_genes)
                !! Expression matrix
            real(real64), intent(out) :: means(n_genes)
                !! Per-gene mean expression values
        end subroutine compute_gene_means_helper
    end interface compute_gene_means_helper

    interface compute_residuals
        !> Compute signed residuals (centering by mean)
        pure module subroutine compute_residuals(n_genes, n_reps, expr, means, resid, ierr)
            import
            integer(int32), intent(in) :: n_genes
                !! Number of genes in the study
            integer(int32), intent(in) :: n_reps
                !! Number of biological replicates in the study
            real(real64), intent(in) :: expr(n_reps, n_genes)
                !! Expression matrix containing
            real(real64), intent(in) :: means(n_genes)
                !! Per-gene mean expression values
            real(real64), intent(out) :: resid(n_reps, n_genes)
                !! Matrix of signed residuals
            integer(int32), intent(out) :: ierr
                !! Error code
        end subroutine compute_residuals
    end interface compute_residuals

    interface compute_residuals_helper
        !> (no input validation) Compute signed residuals (centering by mean)
        pure module subroutine compute_residuals_helper(n_genes, n_reps, expr, means, resid)
            import
            integer(int32), intent(in) :: n_genes
                !! Number of genes in the study
            integer(int32), intent(in) :: n_reps
                !! Number of biological replicates in the study
            real(real64), intent(in) :: expr(n_reps, n_genes)
                !! Expression matrix containing
            real(real64), intent(in) :: means(n_genes)
                !! Per-gene mean expression values
            real(real64), intent(out) :: resid(n_reps, n_genes)
                !! Matrix of signed residuals
        end subroutine compute_residuals_helper
    end interface compute_residuals_helper

    interface pool_means_alloc
        !> Pool per-gene mean expression values across studies
        pure module subroutine pool_means_alloc(n_genes_S1, mean_S1, n_genes_S2, mean_S2, n_points, n_pool, x_star, ierr)
            import
            integer(int32), intent(in) :: n_genes_S1
                !! Number of genes in study S1
            integer(int32), intent(in) :: n_genes_S2
                !! Number of genes in study S2
            integer(int32), intent(in) :: n_points
                !! Number of reference points to define
            real(real64), intent(in) :: mean_S1(n_genes_S1)
                !! Per-gene mean expression values
            real(real64), intent(in) :: mean_S2(n_genes_S2)
                !! Per-gene mean expression values
            integer(int32), intent(out) :: n_pool
                !! Total number of included (non-NaN) pooled mean-expression values
            real(real64), intent(out) :: x_star(n_points)
                !! Mean-expression reference points
            integer(int32), intent(out) :: ierr
                !! Error code
        end subroutine pool_means_alloc
    end interface pool_means_alloc

    interface pool_means
        !> Pool per-gene mean expression values across studies
        pure module subroutine pool_means(pooled_means, pooled_means_perm, pool_size, n_points, n_pool, x_star, ierr)
            import
            integer(int32), intent(in), target :: pool_size
                !! Number of means in the pool, usually `n_genes_S1 + n_genes_S2`
            integer(int32), intent(in) :: n_points
                !! Number of reference points to define
            integer(int32), intent(out) :: n_pool
                !! Total number of included (non-NaN) pooled mean-expression values
            real(real64), intent(in) :: pooled_means(pool_size)
                !! Pooled means
            integer(int32), intent(in) :: pooled_means_perm(pool_size)
                !! Sorting permutation for `pooled_means`
            real(real64), intent(out) :: x_star(n_points)
                !! Mean-expression reference points
            integer(int32), intent(out) :: ierr
                !! Error code
        end subroutine pool_means
    end interface pool_means

    interface pool_means_helper
        !> (no input validation) Pool per-gene mean expression values across studies
        pure module subroutine pool_means_helper(pooled_means, pooled_means_perm, pool_size, n_points, n_pool, x_star)
            import
            integer(int32), intent(in), target :: pool_size
                !! Number of means in the pool, usually `n_genes_S1 + n_genes_S2`
            integer(int32), intent(in) :: n_points
                !! Number of reference points to define
            integer(int32), intent(out) :: n_pool
                !! Total number of included (non-NaN) pooled mean-expression values
            real(real64), intent(in) :: pooled_means(pool_size)
                !! Pooled means
            integer(int32), intent(in) :: pooled_means_perm(pool_size)
                !! Sorting permutation for `pooled_means`
            real(real64), intent(out) :: x_star(n_points)
                !! Mean-expression reference points
        end subroutine pool_means_helper
    end interface pool_means_helper

        !> Calculate the number of neighbors to be used for [[tox_data_integration(module):construct_neighborhoods(interface)]].
        !|
    interface calc_neighborhood_size
        !| The `desired_size` works as upper limit, as the actual neighborhood size might be lower due to few genes with non-NaN mean.
        pure module function calc_neighborhood_size(n_pool, n_points, n_genes_S, mean_S, desired_size) result(n_neighbors)
            import
            integer(int32), intent(in) :: n_pool
                !! Total number of pooled mean-expression values across both studies
            integer(int32), intent(in) :: n_points
                !! Number of reference points
            integer(int32), intent(in) :: n_genes_S
                !! Number of genes in the current study
            real(real64), intent(in) :: mean_S(n_genes_S)
                !! Per-gene mean expression values
            integer(int32), intent(in), optional :: desired_size
                !! Optional desired neighborhood size, default=1000
            integer(int32) :: n_neighbors
                !! Calculated neighborhood size
        end function calc_neighborhood_size
    end interface calc_neighborhood_size

    interface construct_neighborhoods_alloc
        !> Construct neighborhood-based residual sets (kNN)
        pure module subroutine construct_neighborhoods_alloc(n_points, x_star, n_genes_S, mean_S, n_reps_S, resid_S, &
                                                      neighborhood_residuals, neighborhood_indices, n_neighbors, ierr)
            import
            integer(int32), intent(in) :: n_points
                !! Number of reference points
            integer(int32), intent(in) :: n_genes_S
                !! Number of genes in the current study
            integer(int32), intent(in) :: n_reps_S
                !! Number of biological replicates in the study
            real(real64), intent(in) :: x_star(n_points)
                !! Mean-expression reference points
            real(real64), intent(in) :: mean_S(n_genes_S)
                !! Per-gene mean expression values
            real(real64), intent(in) :: resid_S(n_reps_S, n_genes_S)
                !! Matrix of signed residuals
            real(real64), intent(out) :: neighborhood_residuals(n_reps_S, n_neighbors, n_points)
                !! Collection of residual vectors for each neighborhood
            integer(int32), intent(out) :: neighborhood_indices(n_neighbors, n_points)
                !! Indices of selected neighborhood genes
            integer(int32), intent(in) :: n_neighbors
                !! Number of neighbors, **CALCULATE IT WITH [[tox_data_integration(module):calc_neighborhood_size(interface)]]**
            integer(int32), intent(out) :: ierr
                !! Error code
        end subroutine construct_neighborhoods_alloc
    end interface construct_neighborhoods_alloc

    interface construct_neighborhoods
        !> Construct neighborhood-based residual sets (kNN)
        pure module subroutine construct_neighborhoods(n_points, x_star, n_genes_S, mean_S, n_reps_S, resid_S, tmp_distances, tmp_distances_perm, neighborhood_residuals, neighborhood_indices, n_neighbors, ierr)
            import
            integer(int32), intent(in) :: n_points
                !! Number of reference points
            integer(int32), intent(in) :: n_genes_S
                !! Number of genes in the current study
            integer(int32), intent(in) :: n_reps_S
                !! Number of biological replicates in the study
            real(real64), intent(in) :: x_star(n_points)
                !! Mean-expression reference points
            real(real64), intent(in) :: mean_S(n_genes_S)
                !! Per-gene mean expression values
            real(real64), intent(in) :: resid_S(n_reps_S, n_genes_S)
                !! Matrix of signed residuals
            real(real64), intent(out) :: tmp_distances(n_genes_S)
                !! Distances work array
            integer(int32), intent(out) :: tmp_distances_perm(n_genes_S)
                !! Work array for permutation vector to sort `tmp_distances_perm`
            real(real64), intent(out) :: neighborhood_residuals(n_reps_S, n_neighbors, n_points)
                !! Collection of residual vectors for each neighborhood
            integer(int32), intent(out) :: neighborhood_indices(n_neighbors, n_points)
                !! Indices of selected neighborhood genes
            integer(int32), intent(in) :: n_neighbors
                !! Number of neighbors, **CALCULATE IT WITH [[tox_data_integration(module):calc_neighborhood_size(interface)]]**
            integer(int32), intent(out) :: ierr
                !! Error code
        end subroutine construct_neighborhoods
    end interface construct_neighborhoods

    interface construct_neighborhoods_helper
        !> (no input validation) Construct neighborhood-based residual sets (kNN)
        pure module subroutine construct_neighborhoods_helper(n_points, x_star, n_genes_S, mean_S, n_reps_S, resid_S, tmp_distances, tmp_distances_perm, &
                                                       neighborhood_residuals, neighborhood_indices, n_neighbors)
            import
            integer(int32), intent(in) :: n_points
                !! Number of reference points
            integer(int32), intent(in) :: n_genes_S
                !! Number of genes in the current study
            integer(int32), intent(in) :: n_reps_S
                !! Number of biological replicates in the study
            real(real64), intent(in) :: x_star(n_points)
                !! Mean-expression reference points
            real(real64), intent(in) :: mean_S(n_genes_S)
                !! Per-gene mean expression values
            real(real64), intent(in) :: resid_S(n_reps_S, n_genes_S)
                !! Matrix of signed residuals
            real(real64), intent(out) :: tmp_distances(n_genes_S)
                !! Distances work array
            integer(int32), intent(out) :: tmp_distances_perm(n_genes_S)
                !! Work array for permutation vector to sort `tmp_distances_perm`
            real(real64), intent(out) :: neighborhood_residuals(n_reps_S, n_neighbors, n_points)
                !! Collection of residual vectors for each neighborhood
            integer(int32), intent(out) :: neighborhood_indices(n_neighbors, n_points)
                !! Indices of selected neighborhood genes
            integer(int32), intent(in) :: n_neighbors
                !! Number of neighbors, **CALCULATE IT WITH [[tox_data_integration(module):calc_neighborhood_size(interface)]]**
        end subroutine construct_neighborhoods_helper
    end interface construct_neighborhoods_helper

    interface gjct_permutation_test_alloc
        !> Estimates how likely the observed divergence is to occur by chance under the null hypothesis that both studies are exchangeable
        module subroutine gjct_permutation_test_alloc(neighborhood_residuals_S1, neighborhood_residuals_S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, global_jsd_observed, n_bins, shared_residual_range, n_permutations, jsd_null, p_value, ierr, random_seed, neighbor_mask_S1, neighbor_mask_S2)
            import
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
            real(real64), intent(in) :: global_jsd_observed
                !! Observed global JSD value for both studies (from [[tox_data_integration_jsd(module):compute_weighted_global_divergence(subroutine)]])
            integer(int32), intent(in) :: n_bins
                !! Number of equally sized histogram bins used for the studies in [[tox_data_integration_jsd(module):build_residual_histograms(subroutine)]]
            real(real64), intent(in) :: shared_residual_range
                !! Computed residual range for both studies, from [[tox_data_integration_jsd(module):determine_shared_residual_range(subroutine)]]
            integer(int32), intent(in) :: n_permutations
                !! Number of permutations to perform
            real(real64), dimension(n_permutations), intent(out) :: jsd_null
                !! Vector of global divergence values obtained under the null hypothesis
            real(real64), intent(out) :: p_value
                !! Empirical p-value of the permutation test: \( \frac{\text{count}(jsd\_null \ge global\_jsd\_observed) + 1}{n\_permutations} \)
            integer(int32), intent(out) :: ierr
                !! Error code
            integer(int32), intent(in), optional :: random_seed
                !! Seed to use for shuffling
            logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S1
                !! Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
            logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S2
                !! Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)
        end subroutine gjct_permutation_test_alloc
    end interface gjct_permutation_test_alloc

    interface gjct_permutation_test
        !> Estimates how likely the observed divergence is to occur by chance under the null hypothesis that both studies are exchangeable
        module subroutine gjct_permutation_test( &
                neighborhood_residuals_S1_copy, neighborhood_residuals_S2_copy, n_reps_S1, n_reps_S2, n_neighbors, n_points, global_jsd_observed, n_bins, shared_residual_range, n_permutations, jsd_null, p_value, &
                tmp_pool, tmp_pmf_S1, tmp_pmf_S2, tmp_counts, tmp_included_n_reps_S1, tmp_included_n_reps_S2, tmp_js_divergences, tmp_weights, &
                ierr, random_seed, neighbor_mask_S1, neighbor_mask_S2 &
            )
            import
            integer(int32), intent(in) :: n_reps_S1
                !! Number of replicates in study 1
            integer(int32), intent(in) :: n_reps_S2
                !! Number of replicates in study 2
            integer(int32), intent(in) :: n_neighbors
                !! Number of neighbors in the studies
            integer(int32), intent(in) :: n_points
                !! Number of reference points in the studies
            real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(inout) :: neighborhood_residuals_S1_copy
                !! Copy (if wanted) of the computed neighborhood residuals for study 1, will be shuffled in-place
            real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(inout) :: neighborhood_residuals_S2_copy
                !! Copy (if wanted) of the computed neighborhood residuals for study 2, will be shuffled in-place
            real(real64), intent(in) :: global_jsd_observed
                !! Observed global JSD value for both studies (from [[tox_data_integration_jsd(module):compute_weighted_global_divergence(subroutine)]])
            integer(int32), intent(in) :: n_bins
                !! Number of equally sized histogram bins used for the studies in [[tox_data_integration_jsd(module):build_residual_histograms(subroutine)]]
            real(real64), intent(in) :: shared_residual_range
                !! Computed residual range for both studies, from [[tox_data_integration_jsd(module):determine_shared_residual_range(subroutine)]]
            integer(int32), intent(in) :: n_permutations
                !! Number of permutations to perform
            real(real64), dimension(n_permutations), intent(out) :: jsd_null
                !! Vector of global divergence values obtained under the null hypothesis
            real(real64), intent(out) :: p_value
                !! Empirical p-value of the permutation test: \( \frac{\text{count}(jsd\_null \ge global\_jsd\_observed) + 1}{n\_permutations} \)
            real(real64), dimension(n_reps_S1 + n_reps_S2, n_neighbors), intent(out) :: tmp_pool
                !! Working array for shuffling the concatenated residuals from both studies per reference point
            real(real64), dimension(n_points, n_bins), intent(out) :: tmp_pmf_S1
                !! Absolute counts of a residual per bin obtained from [[tox_data_integration_jsd(module):build_residual_histograms(subroutine)]]
            real(real64), dimension(n_points, n_bins), intent(out) :: tmp_pmf_S2
                !! Absolute counts of a residual per bin obtained from [[tox_data_integration_jsd(module):build_residual_histograms(subroutine)]]
            integer(int32), dimension(n_points, n_bins), intent(out) :: tmp_counts
                !! Working array for [[tox_data_integration_jsd(module):build_residual_histograms(subroutine)]]
            integer(int32), dimension(n_points), intent(out) :: tmp_included_n_reps_S1
                !! Working array for [[tox_data_integration_jsd(module):build_residual_histograms(subroutine)]]
            integer(int32), dimension(n_points), intent(out) :: tmp_included_n_reps_S2
                !! Working array for [[tox_data_integration_jsd(module):build_residual_histograms(subroutine)]]
            real(real64), dimension(n_points), intent(out) :: tmp_js_divergences
                !! Working array for [[tox_data_integration_jsd(module):compute_divergence_per_reference_point(subroutine)]]
            real(real64), dimension(n_points), intent(out) :: tmp_weights
                !! Working array for [[tox_data_integration_jsd(module):compute_weighted_global_divergence(subroutine)]]
            integer(int32), intent(out) :: ierr
                !! Error code
            integer(int32), intent(in), optional :: random_seed
                !! Seed to use for shuffling
            logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S1
                !! Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
            logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S2
                !! Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)
        end subroutine gjct_permutation_test
    end interface gjct_permutation_test

    interface gjct_permutation_test_helper
        !> (no input validation) Estimates how likely the observed divergence is to occur by chance under the null hypothesis that both studies are exchangeable
        module subroutine gjct_permutation_test_helper( &
                neighborhood_residuals_S1_copy, neighborhood_residuals_S2_copy, n_reps_S1, n_reps_S2, n_neighbors, n_points, global_jsd_observed, n_bins, shared_residual_range, n_permutations, jsd_null, p_value, &
                tmp_pool, tmp_pmf_S1, tmp_pmf_S2, tmp_counts, tmp_included_n_reps_S1, tmp_included_n_reps_S2, tmp_js_divergences, tmp_weights, &
                random_seed, neighbor_mask_S1, neighbor_mask_S2 &
            )
            import
            integer(int32), intent(in) :: n_reps_S1
                !! Number of replicates in study 1
            integer(int32), intent(in) :: n_reps_S2
                !! Number of replicates in study 2
            integer(int32), intent(in) :: n_neighbors
                !! Number of neighbors in the studies
            integer(int32), intent(in) :: n_points
                !! Number of reference points in the studies
            real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(inout), target :: neighborhood_residuals_S1_copy
                !! Copy (if wanted) of the computed neighborhood residuals for study 1, will be shuffled in-place
            real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(inout), target :: neighborhood_residuals_S2_copy
                !! Copy (if wanted) of the computed neighborhood residuals for study 2, will be shuffled in-place
            real(real64), intent(in) :: global_jsd_observed
                !! Observed global JSD value for both studies (from [[tox_data_integration_jsd(module):compute_weighted_global_divergence(subroutine)]])
            integer(int32), intent(in) :: n_bins
                !! Number of equally sized histogram bins used for the studies in [[tox_data_integration_jsd(module):build_residual_histograms(subroutine)]]
            real(real64), intent(in) :: shared_residual_range
                !! Computed residual range for both studies, from [[tox_data_integration_jsd(module):determine_shared_residual_range(subroutine)]]
            integer(int32), intent(in) :: n_permutations
                !! Number of permutations to perform
            real(real64), dimension(n_permutations), intent(out) :: jsd_null
                !! Vector of global divergence values obtained under the null hypothesis
            real(real64), intent(out) :: p_value
                !! Empirical p-value of the permutation test: \( \frac{\text{count}(jsd\_null \ge global\_jsd\_observed) + 1}{n\_permutations} \)
            real(real64), dimension(n_reps_S1 + n_reps_S2, n_neighbors), intent(out), target :: tmp_pool
                !! Working array for shuffling the concatenated residuals from both studies per reference point
            real(real64), dimension(n_points, n_bins), intent(out) :: tmp_pmf_S1
                !! Absolute counts of a residual per bin obtained from [[tox_data_integration_jsd(module):build_residual_histograms(subroutine)]]
            real(real64), dimension(n_points, n_bins), intent(out) :: tmp_pmf_S2
                !! Absolute counts of a residual per bin obtained from [[tox_data_integration_jsd(module):build_residual_histograms(subroutine)]]
            integer(int32), dimension(n_points, n_bins), intent(out) :: tmp_counts
                !! Working array for [[tox_data_integration_jsd(module):build_residual_histograms(subroutine)]]
            integer(int32), dimension(n_points), intent(out) :: tmp_included_n_reps_S1
                !! Working array for [[tox_data_integration_jsd(module):build_residual_histograms(subroutine)]]
            integer(int32), dimension(n_points), intent(out) :: tmp_included_n_reps_S2
                !! Working array for [[tox_data_integration_jsd(module):build_residual_histograms(subroutine)]]
            real(real64), dimension(n_points), intent(out) :: tmp_js_divergences
                !! Working array for [[tox_data_integration_jsd(module):compute_divergence_per_reference_point(subroutine)]]
            real(real64), dimension(n_points), intent(out) :: tmp_weights
                !! Working array for [[tox_data_integration_jsd(module):compute_weighted_global_divergence(subroutine)]]
            integer(int32), intent(in), optional :: random_seed
                !! Seed to use for shuffling
            logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S1
                !! Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
            logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S2
                !! Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)
        end subroutine gjct_permutation_test_helper
    end interface gjct_permutation_test_helper

    interface shuffle_reference_point_helper
        !> Helper for [[tox_data_integration(module):gjct_permutation_test_helper(interface)]] to shuffle reference points
        module subroutine shuffle_reference_point_helper(reference_point_S1, reference_point_S2, n_reps_S1, n_reps_S2, n_neighbors, pool_flat)
            integer(int32), intent(in) :: n_reps_S1
                !! Number of replicates in study 1
            integer(int32), intent(in) :: n_reps_S2
                !! Number of replicates in study 2
            integer(int32), intent(in) :: n_neighbors
                !! Number of neighbors in the studies
            real(real64), dimension(n_reps_S1 * n_neighbors), intent(inout) :: reference_point_S1
                !! Residuals for one reference point in study 1, will be shuffled in-place
            real(real64), dimension(n_reps_S2 * n_neighbors), intent(inout) :: reference_point_S2
                !! Residuals for one reference point in study 2, will be shuffled in-place
            real(real64), dimension((n_reps_S1 + n_reps_S2) * n_neighbors), intent(out) :: pool_flat
                !! Working array for shuffling the concatenated residuals from both studies per reference point
        end subroutine shuffle_reference_point_helper
    end interface shuffle_reference_point_helper
end module tox_data_integration