!> Module to identify gene outliers based on their distances to family centroids.
module tox_get_outliers
  use, intrinsic :: iso_fortran_env, only: real64
  use tox_sorting, only: sort_array
  use loess_module, only: loess_smooth  ! LOESS smoothing module
  implicit none

contains

  !> Compute family scaling factors (dscale) to normalize distances.
  !> If orthologs are present, uses the maximum distance between orthologs. Otherwise, uses LOESS on the median of intra-family distances or stddev.
  !>
  !> @param n_genes     Total number of genes
  !> @param n_families  Total number of gene families
  !> @param distances   Array of Euclidean distances for each gene
  !> @param gene_to_fam Mapping of each gene to its family (1-based)
  !> @param is_ortholog Boolean array indicating orthologs
  !> @param dscale      Output: array of scaling factors per family
  !> @param loess_x     Work array, dimension n_families (LOESS reference x)
  !> @param loess_y     Work array, dimension n_families (LOESS reference y)
  !> @param indices_used Work array, dimension n_families (indices for LOESS)
  !> @param perm_tmp, stack_left_tmp, stack_right_tmp: work arrays, dimension n_genes
  !> @param workspace_weights, workspace_values: work arrays, dimension n_families (LOESS workspace)
  !> @param max_distance_bw_orths (optional) array of max distances between orthologs per family
  !> @param error_code  Error code: 0=ok, -1=missing max_distance_bw_orths, -2=invalid family indices (required)
  pure subroutine compute_family_scaling_hybrid(n_genes, n_families, distances, gene_to_fam, dscale, &
    loess_x, loess_y, indices_used, perm_tmp, stack_left_tmp, stack_right_tmp, workspace_weights, workspace_values, error_code, &
    is_ortholog, max_distance_bw_orths)
    implicit none
    integer, intent(in) :: n_genes, n_families
    real(real64), intent(in) :: distances(n_genes)
    integer, intent(in) :: gene_to_fam(n_genes)
    real(real64), intent(out) :: dscale(n_families)
    real(real64), intent(inout) :: loess_x(n_families), loess_y(n_families)
    integer, intent(inout) :: indices_used(n_families)
    integer, intent(inout) :: perm_tmp(n_genes), stack_left_tmp(n_genes), stack_right_tmp(n_genes)
    real(real64), intent(inout) :: workspace_weights(n_families)
    real(real64), intent(inout) :: workspace_values(1, n_families)
    integer, intent(out) :: error_code  ! Required
    logical, intent(in), optional :: is_ortholog(:)
    real(real64), intent(in), optional :: max_distance_bw_orths(:)

    integer :: i, family_idx, n_in_family, n_orth_in_fam
    real(real64) :: family_distances(n_genes)
    real(real64) :: median_dist, stddev_dist, mean_dist, sumsq
    real(real64), parameter :: default_sigma = 0.5_real64, default_cutoff = 3.0_real64
    real(real64) :: sigma, cutoff
    integer :: err
    integer :: j, m
    integer :: n_valid
    real(real64) :: loess_pred(1,1)

    dscale = 0.0_real64
    err = 0
    ! Use default values for LOESS
    sigma = default_sigma
    cutoff = default_cutoff

    ! Check for invalid family indices
    do i = 1, n_genes
      if (gene_to_fam(i) < 1 .or. gene_to_fam(i) > n_families) then
        dscale = -1.0_real64  ! Set to -1 to indicate error, do not use if error_code /= 0
        error_code = -2
        return
      end if
    end do

    ! Only require max_distance_bw_orths if is_ortholog is present and at least one ortholog is true
    if (present(is_ortholog)) then
      if (any(is_ortholog)) then
        if (.not. present(max_distance_bw_orths)) then
          dscale = -1.0_real64  ! Set to -1 to indicate error, do not use if error_code /= 0
          error_code = -1
          return
        end if
      end if
    end if

    ! First pass: compute median and stddev for each family, and build indices_used for LOESS
    n_valid = 0
    do family_idx = 1, n_families
      n_in_family = 0
      do i = 1, n_genes
        if (gene_to_fam(i) == family_idx) then
          n_in_family = n_in_family + 1
          family_distances(n_in_family) = abs(distances(i))
        end if
      end do
      if (n_in_family <= 1) then
        loess_x(family_idx) = 0.0_real64
        loess_y(family_idx) = 0.0_real64
        cycle
      end if
      do j = 1, n_in_family
        perm_tmp(j) = j
        stack_left_tmp(j) = 0
        stack_right_tmp(j) = 0
      end do
      call sort_array(family_distances(1:n_in_family), perm_tmp(1:n_in_family), stack_left_tmp(1:n_in_family), &
                      stack_right_tmp(1:n_in_family))
      if (mod(n_in_family,2) == 0) then
        median_dist = 0.5_real64 * (family_distances(n_in_family/2) + family_distances(n_in_family/2+1))
      else
        median_dist = family_distances((n_in_family+1)/2)
      end if
      mean_dist = sum(family_distances(1:n_in_family)) / n_in_family
      sumsq = sum((family_distances(1:n_in_family) - mean_dist)**2)
      stddev_dist = sqrt(sumsq / (n_in_family-1))
      loess_x(family_idx) = median_dist
      loess_y(family_idx) = stddev_dist
      n_valid = n_valid + 1
      indices_used(n_valid) = family_idx
    end do
    ! n_valid is now the number of valid families for LOESS
    ! Only use indices_used(1:n_valid) in LOESS calls

    ! Second pass: assign dscale per family
    do family_idx = 1, n_families
      n_in_family = 0
      n_orth_in_fam = 0
      do i = 1, n_genes
        if (gene_to_fam(i) == family_idx) then
          n_in_family = n_in_family + 1
          family_distances(n_in_family) = abs(distances(i))
          if (present(is_ortholog)) then
            if (is_ortholog(i)) n_orth_in_fam = n_orth_in_fam + 1
          end if
        end if
      end do
      if (n_in_family <= 1) then
        dscale(family_idx) = 0.0_real64
        cycle  ! Skip single-gene or empty families
      end if
      ! If more than one ortholog in family, use max_distance_bw_orths if present
      if (present(is_ortholog) .and. n_orth_in_fam > 1) then
        if (present(max_distance_bw_orths)) then
          dscale(family_idx) = max_distance_bw_orths(family_idx)
        else
          dscale(family_idx) = 0.0_real64
          err = -1
        end if
      else
        ! Otherwise, use LOESS smoothing on median_dist
        if (n_valid > 0) then
          if (mod(n_in_family,2) == 0) then
            median_dist = 0.5_real64 * (family_distances(n_in_family/2) + family_distances(n_in_family/2+1))
          else
            median_dist = family_distances((n_in_family+1)/2)
          end if
          ! Only pass first n_valid elements of LOESS arrays
          call loess_smooth(n_valid, 1, 1, loess_x, reshape(loess_y, [n_valid,1]), indices_used, [median_dist], &
                  sigma, cutoff, loess_pred, workspace_weights, workspace_values)
          dscale(family_idx) = loess_pred(1,1)
        else
          dscale(family_idx) = 0.0_real64
        end if
      end if
    end do
    error_code = err
  end subroutine compute_family_scaling_hybrid

  !> Compute the hybrid RDI for each gene.
  !>
  !> RDI = Euclidean distance / family scaling factor
  !>
  !> @param n_genes     Total number of genes
  !> @param distances   Array of Euclidean distances for each gene to its centroid
  !> @param gene_to_fam Gene-to-family mapping (1-based indexing)
  !> @param dscale      Array of scaling factors for each family
  !> @param rdi         Output array of RDI values for each gene
  pure subroutine compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi)
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    
    integer, intent(in) :: n_genes
    real(real64), intent(in) :: distances(n_genes)
    integer, intent(in) :: gene_to_fam(n_genes)
    real(real64), intent(in) :: dscale(:)
    real(real64), intent(out) :: rdi(n_genes)
    integer :: i, family_idx
    real(real64), parameter :: tol = epsilon(1.0_real64)
    
    ! Calculate RDI for each gene
    do i = 1, n_genes
      family_idx = gene_to_fam(i)
      
      ! Handle invalid family indices
      if (family_idx < 1 .or. family_idx > size(dscale)) then
        rdi(i) = -1.0_real64  ! Error indicator
        cycle
      end if
      
      ! Detect NaN input (portable)
      if (distances(i) /= distances(i)) then
        rdi(i) = distances(i)
      else if (abs(dscale(family_idx)) < tol) then
        rdi(i) = 0.0_real64  ! If scaling is zero, set RDI to zero (not outlier)
      else
        ! Calculate RDI
        rdi(i) = abs(distances(i)) / dscale(family_idx)
      end if
    end do
    
  end subroutine compute_rdi

  !> Identify gene outliers based on the top percentile of RDI values.
  !>
  !> @param n_genes      Total number of genes
  !> @param rdi          Array of RDI values for each gene
  !> @param percentile   (optional) Percentile threshold (default: 95 for top 5%)
  !> @param sorted_rdi   Work array for sorting (dimension n_genes)
  !> @param perm         Permutation array for sorting (dimension n_genes, should be pre-initialized with 1:n_genes)
  !> @param stack_left   Stack array for sorting (dimension n_genes)
  !> @param stack_right  Stack array for sorting (dimension n_genes)
  !> @param is_outlier   Output boolean array indicating outliers
  !> @param threshold    Output threshold value used for detection
  pure subroutine identify_outliers(n_genes, rdi, sorted_rdi, perm, &
                                    stack_left, stack_right, is_outlier, threshold, percentile)
    implicit none
    
    integer, intent(in) :: n_genes
    real(real64), intent(in) :: rdi(n_genes)
    real(real64), intent(inout) :: sorted_rdi(n_genes)
    integer, intent(inout) :: perm(n_genes)
    integer, intent(inout) :: stack_left(n_genes)
    integer, intent(inout) :: stack_right(n_genes)
    logical, intent(out) :: is_outlier(n_genes)
    real(real64), intent(out) :: threshold
    real(real64), intent(in), optional :: percentile
    
    integer :: i, idx
    real(real64) :: perc_pos, percentile_val

    ! Set default percentile if not present
    if (present(percentile)) then
      percentile_val = percentile
    else
      percentile_val = 95.0_real64
    end if

    ! Initialize output
    is_outlier = .false.

    ! Create a copy of RDI for sorting (excluding error values)
    sorted_rdi = rdi

    ! Filter out error values (negative RDIs)
    where (sorted_rdi < 0.0_real64)
      sorted_rdi = 0.0_real64
    end where

    ! Sort RDI values using the tox_sorting module
    call sort_array(sorted_rdi, perm, stack_left, stack_right)

    ! Calculate the position corresponding to the desired percentile
    perc_pos = (n_genes * percentile_val) / 100.0_real64
    idx = ceiling(perc_pos)

    ! Get the threshold value from the sorted array using the permutation vector
    if (idx >= 1 .and. idx <= n_genes) then
      threshold = sorted_rdi(perm(idx))
    else
      threshold = 0.0_real64  ! Default if percentile is invalid
    end if

    ! Mark genes as outliers if their RDI exceeds the threshold
    do i = 1, n_genes
      if (rdi(i) >= threshold .and. rdi(i) > 0.0_real64) then
        is_outlier(i) = .true.
      end if
    end do
  end subroutine identify_outliers

  !> Main hybrid routine to detect outliers using RDI.
  !> If is_ortholog is not present, or if a family has no orthologs, uses LOESS.
  !>
  !> @param n_genes      Total number of genes
  !> @param n_families   Total number of gene families
  !> @param distances    Array of Euclidean distances for each gene to its centroid
  !> @param gene_to_fam  Gene-to-family mapping (1-based indexing)
  !> @param is_ortholog  Boolean array indicating orthologs
  !> @param work_array   Work array for sorting (dimension n_genes)
  !> @param perm         Permutation array for sorting (dimension n_genes)
  !> @param stack_left   Stack array for sorting (dimension n_genes)
  !> @param stack_right  Stack array for sorting (dimension n_genes)
  !> @param is_outlier   Output boolean array indicating outliers
  !> @param percentile   (optional) Percentile threshold for outlier detection (default: 95)
  !> @param max_distance_bw_orths (optional) array of max distances between orthologs per family
  !> @param loess_x, loess_y, loess_n, workspace_weights, workspace_values: work arrays, dimension n_families (for LOESS)
  pure subroutine detect_outliers(n_genes, n_families, distances, gene_to_fam, &
                            work_array, perm, stack_left, stack_right, &
                            is_outlier, loess_x, loess_y, loess_n, workspace_weights, workspace_values, error_code, &
                            is_ortholog, percentile, max_distance_bw_orths)
    implicit none
    integer, intent(in) :: n_genes, n_families
    real(real64), intent(in) :: distances(n_genes)
    integer, intent(in) :: gene_to_fam(n_genes)
    real(real64), intent(inout) :: work_array(n_genes)
    integer, intent(inout) :: perm(n_genes)
    integer, intent(inout) :: stack_left(n_genes)
    integer, intent(inout) :: stack_right(n_genes)
    logical, intent(out) :: is_outlier(n_genes)
    real(real64), intent(inout) :: loess_x(n_families), loess_y(n_families)
    integer, intent(inout) :: loess_n(n_families)
    real(real64), intent(inout) :: workspace_weights(n_families)
    real(real64), intent(inout) :: workspace_values(1, n_families)
    integer, intent(out) :: error_code  ! Required
    logical, intent(in), optional :: is_ortholog(n_genes)
    real(real64), intent(in), optional :: percentile
    real(real64), intent(in), optional :: max_distance_bw_orths(:)
    ! Local variables
    real(real64) :: dscale(n_families)
    real(real64) :: rdi(n_genes)
    real(real64) :: threshold
    integer :: i
    real(real64) :: percentile_val

    ! Set default percentile if not present
    if (present(percentile)) then
      percentile_val = percentile
    else
      percentile_val = 95.0_real64
    end if

    ! Always initialize permutation array
    do i = 1, n_genes
      perm(i) = i
    end do

    if (present(is_ortholog)) then
      if (present(max_distance_bw_orths)) then
        call compute_family_scaling_hybrid(n_genes, n_families, distances, gene_to_fam, dscale, &
                                          loess_x, loess_y, loess_n, perm, stack_left, stack_right, workspace_weights, &
                                          workspace_values, error_code, is_ortholog, max_distance_bw_orths)
      else
        call compute_family_scaling_hybrid(n_genes, n_families, distances, gene_to_fam, dscale, &
                                          loess_x, loess_y, loess_n, perm, stack_left, stack_right, workspace_weights, &
                                          workspace_values, error_code, is_ortholog)
      end if
    else
      if (present(max_distance_bw_orths)) then
        call compute_family_scaling_hybrid(n_genes, n_families, distances, gene_to_fam, dscale, &
                                          loess_x, loess_y, loess_n, perm, stack_left, stack_right, workspace_weights, &
                                          workspace_values, error_code, max_distance_bw_orths=max_distance_bw_orths)
      else
        call compute_family_scaling_hybrid(n_genes, n_families, distances, gene_to_fam, dscale, &
                                          loess_x, loess_y, loess_n, perm, stack_left, stack_right, workspace_weights, &
                                          workspace_values, error_code)
      end if
    end if
    if (error_code /= 0) return
    call compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi)
    call identify_outliers(n_genes, rdi, work_array, perm, stack_left, &
                                stack_right, is_outlier, threshold, percentile_val)
  end subroutine detect_outliers

end module tox_get_outliers

