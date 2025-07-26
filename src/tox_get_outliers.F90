!> Module to identify gene outliers based on their distances to family centroids.
module tox_get_outliers
  use, intrinsic :: iso_fortran_env, only: real64
  use f42_utils, only: loess_smooth_2d,sort_array  
  implicit none

contains

  !> Compute family scaling factors (dscale) to normalize distances.
  !> Now always uses LOESS on the median/stddev of intra-family distances for scaling, regardless of orthologs.
  !>
  !> @param n_genes     Total number of genes
  !> @param n_families  Total number of gene families
  !> @param distances   Array of Euclidean distances for each gene
  !> @param gene_to_fam Mapping of each gene to its family (1-based)
  !> @param dscale      Output: array of scaling factors per family
  !> @param loess_x     Work array, dimension n_families (LOESS reference x)
  !> @param loess_y     Work array, dimension n_families (LOESS reference y)
  !> @param indices_used Work array, dimension n_families (indices for LOESS)
  !> @param perm_tmp, stack_left_tmp, stack_right_tmp: work arrays, dimension n_genes
  !> @param workspace_weights, workspace_values: work arrays, dimension n_families (LOESS workspace)
  !> @param error_code  Error code: 0=ok, -2=invalid family indices (required)
  pure subroutine compute_family_scaling(n_genes, n_families, distances, gene_to_fam, dscale, &
    loess_x, loess_y, indices_used, perm_tmp, stack_left_tmp, stack_right_tmp, workspace_weights, workspace_values, error_code)
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
      do i = 1, n_genes
        if (gene_to_fam(i) == family_idx) then
          n_in_family = n_in_family + 1
          family_distances(n_in_family) = abs(distances(i))
        end if
      end do
      if (n_in_family <= 1) then
        dscale(family_idx) = 0.0_real64
        cycle  ! Skip single-gene or empty families
      end if
      if (n_valid > 0) then
        if (mod(n_in_family,2) == 0) then
          median_dist = 0.5_real64 * (family_distances(n_in_family/2) + family_distances(n_in_family/2+1))
        else
          median_dist = family_distances((n_in_family+1)/2)
        end if
        ! Only pass first n_valid elements of LOESS arrays
        call loess_smooth_2d(n_valid, 1, loess_x, loess_y, indices_used, [median_dist], &
                sigma, cutoff, loess_pred, workspace_weights, workspace_values)
        dscale(family_idx) = loess_pred(1,1)
      else
        dscale(family_idx) = 0.0_real64
      end if
    end do
    error_code = err
  end subroutine compute_family_scaling

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

  !> Main routine to detect outliers using RDI and LOESS-based scaling.
  !>
  !> @param n_genes      Total number of genes
  !> @param n_families   Total number of gene families
  !> @param distances    Array of Euclidean distances for each gene to its centroid
  !> @param gene_to_fam  Gene-to-family mapping (1-based indexing)
  !> @param work_array   Work array for sorting (dimension n_genes)
  !> @param perm         Permutation array for sorting (dimension n_genes)
  !> @param stack_left   Stack array for sorting (dimension n_genes)
  !> @param stack_right  Stack array for sorting (dimension n_genes)
  !> @param is_outlier   Output boolean array indicating outliers
  !> @param percentile   (optional) Percentile threshold for outlier detection (default: 95)
  !> @param loess_x, loess_y, loess_n, workspace_weights, workspace_values: work arrays, dimension n_families (for LOESS)
  pure subroutine detect_outliers(n_genes, n_families, distances, gene_to_fam, &
                            work_array, perm, stack_left, stack_right, &
                            is_outlier, loess_x, loess_y, loess_n, workspace_weights, workspace_values, error_code, &
                            percentile)
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
    real(real64), intent(in), optional :: percentile
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

    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, dscale, &
                                loess_x, loess_y, loess_n, perm, stack_left, stack_right, workspace_weights, &
                                workspace_values, error_code)
    if (error_code /= 0) return
    call compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi)
    call identify_outliers(n_genes, rdi, work_array, perm, stack_left, &
                                stack_right, is_outlier, threshold, percentile_val)
  end subroutine detect_outliers


end module tox_get_outliers


  subroutine compute_family_scaling_r(n_genes, n_families, distances, gene_to_fam, dscale, &
      loess_x, loess_y, indices_used, perm_tmp, stack_left_tmp, stack_right_tmp, workspace_weights, workspace_values, error_code)
  use tox_get_outliers
  use iso_fortran_env, only: real64
  integer, intent(in) :: n_genes, n_families
  real(real64), intent(in) :: distances(n_genes)
  integer, intent(in) :: gene_to_fam(n_genes)
  real(real64), intent(out) :: dscale(n_families)
  real(real64), intent(inout) :: loess_x(n_families), loess_y(n_families)
  integer, intent(inout) :: indices_used(n_families)
  integer, intent(inout) :: perm_tmp(n_genes), stack_left_tmp(n_genes), stack_right_tmp(n_genes)
  real(real64), intent(inout) :: workspace_weights(n_families)
  real(real64), intent(inout) :: workspace_values(1, n_families)
  integer, intent(out) :: error_code
    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, dscale, &
      loess_x, loess_y, indices_used, perm_tmp, stack_left_tmp, stack_right_tmp, workspace_weights, workspace_values, error_code)
end subroutine compute_family_scaling_r

subroutine compute_rdi_r(n_genes, n_families, distances, gene_to_fam, dscale, rdi)
  use tox_get_outliers
  use iso_fortran_env, only: real64
  integer, intent(in) :: n_genes, n_families
  real(real64), intent(in) :: distances(n_genes)
  integer, intent(in) :: gene_to_fam(n_genes)
  real(real64), intent(in) :: dscale(n_families)
  real(real64), intent(out) :: rdi(n_genes)
  integer :: i

  call compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi)

end subroutine compute_rdi_r

subroutine identify_outliers_r(n_genes, rdi, sorted_rdi, perm, stack_left, stack_right, is_outlier, threshold, percentile)
  use tox_get_outliers
  use iso_fortran_env, only: real64
  integer, intent(in) :: n_genes
  real(real64), intent(in) :: rdi(n_genes)
  real(real64), intent(inout) :: sorted_rdi(n_genes)
  integer, intent(inout) :: perm(n_genes)
  integer, intent(inout) :: stack_left(n_genes)
  integer, intent(inout) :: stack_right(n_genes)
  logical, intent(out) :: is_outlier(n_genes)
  real(real64), intent(out) :: threshold
  real(real64), intent(in) :: percentile
  real(real64) :: percentile_val

  call identify_outliers(n_genes, rdi, sorted_rdi, perm, stack_left, stack_right, is_outlier, threshold, percentile)
end subroutine identify_outliers_r

subroutine detect_outliers_r(n_genes, n_families, distances, gene_to_fam, &
                            work_array, perm, stack_left, stack_right, &
                            is_outlier, loess_x, loess_y, loess_n, workspace_weights, workspace_values, error_code, &
                            percentile)
  use tox_get_outliers
  use iso_fortran_env, only: real64
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
  integer, intent(out) :: error_code
  real(real64), intent(in) :: percentile

  call detect_outliers(n_genes, n_families, distances, gene_to_fam, &
                      work_array, perm, stack_left, stack_right, &
                      is_outlier, loess_x, loess_y, loess_n, workspace_weights, workspace_values, error_code, &
                      percentile)
end

! C wrappers for RDI/outlier routines

  subroutine compute_family_scaling_c(n_genes, n_families, distances, gene_to_fam, dscale, &
    loess_x, loess_y, indices_used, perm_tmp, stack_left_tmp, stack_right_tmp, workspace_weights, workspace_values, &
    error_code) bind(C, name="compute_family_scaling_c")
  use iso_c_binding
  use tox_get_outliers
  integer(c_int), intent(in), value :: n_genes, n_families
  real(c_double), intent(in), target :: distances(n_genes)
  integer(c_int), intent(in), target :: gene_to_fam(n_genes)
  real(c_double), intent(out), target :: dscale(n_families)
  real(c_double), intent(inout), target :: loess_x(n_families), loess_y(n_families)
  integer(c_int), intent(inout), target :: indices_used(n_families)
  integer(c_int), intent(inout), target :: perm_tmp(n_genes), stack_left_tmp(n_genes), stack_right_tmp(n_genes)
  real(c_double), intent(inout), target :: workspace_weights(n_families)
  real(c_double), intent(inout), target :: workspace_values(1,n_families )
  integer(c_int), intent(out) :: error_code
    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, dscale, &
      loess_x, loess_y, indices_used, perm_tmp, stack_left_tmp, stack_right_tmp, workspace_weights, workspace_values, error_code)
end subroutine compute_family_scaling_c

subroutine compute_rdi_c(n_genes, n_families, distances, gene_to_fam, dscale, rdi) bind(C, name="compute_rdi_c")
  use iso_c_binding
  use tox_get_outliers
  integer(c_int), intent(in), value :: n_genes, n_families
  real(c_double), intent(in), target :: distances(n_genes)
  integer(c_int), intent(in), target :: gene_to_fam(n_genes)
  real(c_double), intent(in), target :: dscale(n_families)
  real(c_double), intent(out), target :: rdi(n_genes)
  call compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi)
end subroutine compute_rdi_c

subroutine identify_outliers_c(n_genes, rdi, sorted_rdi, perm, stack_left, stack_right, is_outlier_int, threshold, percentile) &
                              bind(C, name="identify_outliers_c")
  use iso_c_binding
  use tox_get_outliers
  integer(c_int), intent(in), value :: n_genes
  real(c_double), intent(in), target :: rdi(n_genes)
  real(c_double), intent(inout), target :: sorted_rdi(n_genes)
  integer(c_int), intent(inout), target :: perm(n_genes)
  integer(c_int), intent(inout), target :: stack_left(n_genes)
  integer(c_int), intent(inout), target :: stack_right(n_genes)
  integer(c_int), intent(out), target :: is_outlier_int(n_genes)
  real(c_double), intent(out) :: threshold
  real(c_double), intent(in), value :: percentile
  logical :: is_outlier(n_genes)
  integer :: i

  ! Convert integer (0/1) to logical (.false./.true.) for is_outlier
  do i = 1, n_genes
    is_outlier(i) = (is_outlier_int(i) /= 0)
  end do

  call identify_outliers(n_genes, rdi, sorted_rdi, perm, stack_left, stack_right, is_outlier, threshold, percentile)
  ! Convert logical (.true./.false.) to integer (1/0)
  do i = 1, n_genes
    if (is_outlier(i)) then
      is_outlier_int(i) = 1
    else
      is_outlier_int(i) = 0
    end if
  end do
end subroutine identify_outliers_c


subroutine detect_outliers_c(n_genes, n_families, distances, gene_to_fam, &
                          work_array, perm, stack_left, stack_right, &
                          is_outlier_int, loess_x, loess_y, loess_n, workspace_weights, workspace_values, error_code, &
                          percentile) bind(C, name="detect_outliers_c")
  use iso_c_binding
  use tox_get_outliers
  integer(c_int), intent(in), value :: n_genes, n_families
  real(c_double), intent(in), target :: distances(n_genes)
  integer(c_int), intent(in), target :: gene_to_fam(n_genes)
  real(c_double), intent(inout), target :: work_array(n_genes)
  integer(c_int), intent(inout), target :: perm(n_genes)
  integer(c_int), intent(inout), target :: stack_left(n_genes)
  integer(c_int), intent(inout), target :: stack_right(n_genes)
  integer(c_int), intent(out), target :: is_outlier_int(n_genes)
  real(c_double), intent(inout), target :: loess_x(n_families), loess_y(n_families)
  integer(c_int), intent(inout), target :: loess_n(n_families)
  real(c_double), intent(inout), target :: workspace_weights(n_families)
  real(c_double), intent(inout), target :: workspace_values(1,n_families)
  integer(c_int), intent(out) :: error_code
  real(c_double), intent(in), value :: percentile
  logical :: is_outlier(n_genes)
  integer :: i
  call detect_outliers(n_genes, n_families, distances, gene_to_fam, &
                    work_array, perm, stack_left, stack_right, &
                    is_outlier, loess_x, loess_y, loess_n, workspace_weights, workspace_values, error_code, &
                    percentile)
  ! Print output for debugging
  write(*, '(A)', advance='no') 'detect_outliers_c: is_outlier = ['
  do i = 1, n_genes
    if (is_outlier(i)) then
      write(*, '(I1)', advance='no') 1
    else
      write(*, '(I1)', advance='no') 0
    end if
    if (i < n_genes) write(*, '(A)', advance='no') ', '
  end do
  write(*, '(A)') ']'
  write(*, '(A, I0)') 'detect_outliers_c: error_code = ', error_code
  ! Convert logical (.true./.false.) to integer (1/0) for is_outlier
  do i = 1, n_genes
    if (is_outlier(i)) then
      is_outlier_int(i) = 1
    else
      is_outlier_int(i) = 0
    end if
  end do
end subroutine detect_outliers_c