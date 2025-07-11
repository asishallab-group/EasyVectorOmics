!> Module to identify gene outliers based on their distances to family centroids.
module tox_get_outliers
  use, intrinsic :: iso_fortran_env, only: real64
  use tox_sorting, only: sort_array
  implicit none

contains

  !> Calculate the family scaling factors (dscale) for normalizing distances.
  !>
  !> For each family, computes the maximum distance between ortholog genes
  !> and the family centroid, which serves as the scaling factor.
  !>
  !> @param n_genes     Total number of genes
  !> @param n_families  Total number of gene families
  !> @param distances   Array of Euclidean distances for each gene to its centroid
  !> @param gene_to_fam Gene-to-family mapping (1-based indexing)
  !> @param is_ortholog Boolean array indicating which genes are orthologous
  !> @param dscale      Output array of scaling factors for each family
  pure subroutine compute_family_scaling(n_genes, n_families, distances, gene_to_fam, is_ortholog, dscale)
    implicit none
    
    integer, intent(in) :: n_genes, n_families
    real(real64), intent(in) :: distances(n_genes)
    integer, intent(in) :: gene_to_fam(n_genes)
    logical, intent(in) :: is_ortholog(n_genes)  ! New parameter
    real(real64), intent(out) :: dscale(n_families)
    
    integer :: i, family_idx
    integer :: ortholog_count(n_families)
    
    ! Initialize scaling factors to 0
    dscale = 0.0_real64
    ortholog_count = 0
    
    ! Find maximum distance for each family, considering only orthologous genes
    do i = 1, n_genes
      ! Skip non-orthologous genes
      if (.not. is_ortholog(i)) cycle
      
      family_idx = gene_to_fam(i)
      
      ! Skip invalid family indices
      if (family_idx < 1 .or. family_idx > n_families) cycle
      
      ! Update maximum distance if current distance is larger
      if (abs(distances(i)) > dscale(family_idx)) then
        dscale(family_idx) = abs(distances(i))
      end if
      ortholog_count(family_idx) = ortholog_count(family_idx) + 1
    end do
    
    ! Handle edge case: if any scaling factor is 0, set it to 1 to avoid division by zero
    do i = 1, n_families
      if (ortholog_count(i) == 0) dscale(i) = 1.0_real64
    end do
    
  end subroutine compute_family_scaling

  !> Calculate the Relative Distance Index (RDI) for each gene.
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
    real(real64), intent(in) :: dscale(:)  ! Dimension should be n_families
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
  !> @param percentile   Percentile threshold (e.g., 95 for top 5%)
  !> @param sorted_rdi   Work array for sorting (dimension n_genes)
  !> @param perm         Permutation array for sorting (dimension n_genes, &
  !>                     should be pre-initialized with 1:n_genes)
  !> @param stack_left   Stack array for sorting (dimension n_genes)
  !> @param stack_right  Stack array for sorting (dimension n_genes)
  !> @param is_outlier   Output boolean array indicating outliers
  !> @param threshold    Output threshold value used for detection
  !> @param initialize_perm Logical flag to control whether to initialize &
  !>                        the permutation array (default: false)
  pure subroutine identify_outliers(n_genes, rdi, percentile, sorted_rdi, perm, &
                                    stack_left, stack_right, is_outlier, threshold, initialize_perm)
    implicit none
    
    integer, intent(in) :: n_genes
    real(real64), intent(in) :: rdi(n_genes)
    real(real64), intent(in) :: percentile
    real(real64), intent(inout) :: sorted_rdi(n_genes)
    integer, intent(inout) :: perm(n_genes)
    integer, intent(inout) :: stack_left(n_genes)
    integer, intent(inout) :: stack_right(n_genes)
    logical, intent(out) :: is_outlier(n_genes)
    real(real64), intent(out) :: threshold
    logical, intent(in), optional :: initialize_perm
    
    integer :: i, idx
    real(real64) :: perc_pos
    logical :: init_perm
    
    ! Determine if we need to initialize the permutation array
    if (present(initialize_perm)) then
      init_perm = initialize_perm
    else
      init_perm = .false.  ! Default is to not re-initialize
    end if
    
    ! Initialize output
    is_outlier = .false.
    
    ! Create a copy of RDI for sorting (excluding error values)
    sorted_rdi = rdi
    
    ! Filter out error values (negative RDIs)
    where (sorted_rdi < 0.0_real64)
      sorted_rdi = 0.0_real64
    end where
    
    ! Initialize permutation array with identity permutation only if requested or first time
    if (init_perm) then
      do i = 1, n_genes
        perm(i) = i
      end do
    end if
    
    ! Sort RDI values using the tox_sorting module
    call sort_array(sorted_rdi, perm, stack_left, stack_right)
    
    ! Calculate the position corresponding to the desired percentile
    perc_pos = (n_genes * percentile) / 100.0_real64
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

  !> Main routine to detect gene outliers using the RDI method.
  !>
  !> This is the primary entry point that combines all steps:
  !> 1. Computes family scaling factors using orthologous genes only
  !> 2. Calculates RDI for each gene
  !> 3. Identifies outliers based on percentile threshold
  !>
  !> @param n_genes      Total number of genes
  !> @param n_families   Total number of gene families
  !> @param distances    Array of Euclidean distances for each gene to its centroid
  !> @param gene_to_fam  Gene-to-family mapping (1-based indexing)
  !> @param is_ortholog  Boolean array indicating which genes are orthologous
  !> @param percentile   Percentile threshold for outlier detection (default: 95)
  !> @param work_array   Work array for sorting (dimension n_genes)
  !> @param perm         Permutation array for sorting (dimension n_genes)
  !> @param stack_left   Stack array for sorting (dimension n_genes)
  !> @param stack_right  Stack array for sorting (dimension n_genes)
  !> @param is_outlier   Output boolean array indicating outliers
  !> @param rdi_values   Output array of RDI values (optional)
  !> @param rdi_threshold Output threshold value used for detection (optional)
  !> @param init_perm    Flag to control permutation array initialization (optional)
  pure subroutine detect_outliers(n_genes, n_families, distances, gene_to_fam, &
                            is_ortholog, percentile, work_array, perm, stack_left, stack_right, &
                            is_outlier, rdi_values, rdi_threshold, init_perm)
    implicit none
    
    integer, intent(in) :: n_genes, n_families
    real(real64), intent(in) :: distances(n_genes)
    integer, intent(in) :: gene_to_fam(n_genes)
    logical, intent(in) :: is_ortholog(n_genes)
    real(real64), intent(in) :: percentile
    real(real64), intent(inout) :: work_array(n_genes)
    integer, intent(inout) :: perm(n_genes)
    integer, intent(inout) :: stack_left(n_genes)
    integer, intent(inout) :: stack_right(n_genes)
    logical, intent(out) :: is_outlier(n_genes)
    real(real64), intent(out), optional :: rdi_values(n_genes)
    real(real64), intent(out), optional :: rdi_threshold
    logical, intent(in), optional :: init_perm
    
    real(real64) :: dscale(n_families)
    real(real64) :: rdi(n_genes)
    real(real64) :: threshold
    logical :: initialize_perm
    integer :: i
    
    ! Determine if we need to initialize the permutation array
    if (present(init_perm)) then
      initialize_perm = init_perm
    else
      ! By default, initialize the permutation array once
      initialize_perm = .true.
    end if
    
    ! Initialize permutation array here if requested, to avoid doing it in every call
    if (initialize_perm) then
      do i = 1, n_genes
        perm(i) = i
      end do
    end if
    
    ! Step 1: Compute family scaling factors (using only orthologous genes)
    call compute_family_scaling(n_genes, n_families, distances, gene_to_fam, is_ortholog, dscale)
    
    ! Step 2: Calculate RDI for each gene
    call compute_rdi(n_genes, distances, gene_to_fam, dscale, rdi)
    
    ! Step 3: Identify outliers based on percentile
    ! Pass false for initialize_perm since we already initialized it above if needed
    call identify_outliers(n_genes, rdi, percentile, work_array, perm, stack_left, &
                                stack_right, is_outlier, threshold, .false.)
    
    ! Return optional outputs if requested
    if (present(rdi_values)) then
      rdi_values = rdi
    end if
    
    if (present(rdi_threshold)) then
      rdi_threshold = threshold
    end if
    
  end subroutine detect_outliers

end module tox_get_outliers

!> Wrapper for R and Fortran usage (R .Fortran)
subroutine detect_outliers_r(n_genes, n_families, distances, gene_to_fam, &
                            is_ortholog, percentile, work_array, perm, stack_left, stack_right, &
                            is_outlier, rdi_values, rdi_threshold, init_perm)
  use tox_get_outliers
  integer, intent(in) :: n_genes, n_families
  real(real64), intent(in) :: distances(n_genes)
  integer, intent(in) :: gene_to_fam(n_genes)
  logical, intent(in) :: is_ortholog(n_genes)
  real(real64), intent(in) :: percentile
  real(real64), intent(inout) :: work_array(n_genes)
  integer, intent(inout) :: perm(n_genes)
  integer, intent(inout) :: stack_left(n_genes)
  integer, intent(inout) :: stack_right(n_genes)
  logical, intent(out) :: is_outlier(n_genes)
  real(real64), intent(out) :: rdi_values(n_genes)
  real(real64), intent(out) :: rdi_threshold
  logical, intent(in), optional :: init_perm
  
  call detect_outliers(n_genes, n_families, distances, gene_to_fam, &
                      is_ortholog, percentile, work_array, perm, stack_left, stack_right, &
                      is_outlier, rdi_values, rdi_threshold, init_perm)
end subroutine detect_outliers_r

!> C/Python interface (bind(C)), that uses integers (0/1) instead of booleans for logical flags. 
!> This simplifies the C wrapper and avoids C-Fortran logical type conversion issues.
!>
!> For each logical array:
!>   0 = false
!>   1 = true (or any non-zero value)
!>
!> @param n_genes      Total number of genes
!> @param n_families   Total number of gene families
!> @param distances    Array of Euclidean distances for each gene to its centroid
!> @param gene_to_fam  Gene-to-family mapping (1-based indexing)
!> @param is_ortholog_int Integer array indicating which genes are orthologous (0=false, 1=true)
!> @param percentile   Percentile threshold for outlier detection
!> @param work_array   Work array for sorting (dimension n_genes)
!> @param perm         Permutation array for sorting (dimension n_genes)
!> @param stack_left   Stack array for sorting (dimension n_genes)
!> @param stack_right  Stack array for sorting (dimension n_genes)
!> @param is_outlier_int Output integer array indicating outliers (0=false, 1=true)
!> @param rdi_values   Output array of RDI values
!> @param rdi_threshold Output threshold value used for detection
!> @param init_perm_int Flag to control permutation array initialization (0=false, 1=true)
!>
subroutine detect_outliers_c_int(n_genes, n_families, distances, gene_to_fam, &
                              is_ortholog_int, percentile, work_array, perm, &
                              stack_left, stack_right, is_outlier_int, &
                              rdi_values, rdi_threshold, init_perm_int) &
                              bind(C, name="detect_outliers_c_int")
  use iso_c_binding
  use tox_get_outliers
  integer(c_int), intent(in), value :: n_genes
  integer(c_int), intent(in), value :: n_families
  real(c_double), intent(in), target :: distances(*)
  integer(c_int), intent(in), target :: gene_to_fam(*)
  integer(c_int), intent(in), target :: is_ortholog_int(*)    ! 0=false, 1=true
  real(c_double), intent(in), value :: percentile
  real(c_double), intent(inout), target :: work_array(*)
  integer(c_int), intent(inout), target :: perm(*)
  integer(c_int), intent(inout), target :: stack_left(*)
  integer(c_int), intent(inout), target :: stack_right(*)
  integer(c_int), intent(out), target :: is_outlier_int(*)    ! 0=false, 1=true
  real(c_double), intent(out), target :: rdi_values(*)
  real(c_double), intent(out) :: rdi_threshold
  integer(c_int), intent(in), value :: init_perm_int          ! 0=false, 1=true
  
  !> Convert integer flags to Fortran logical arrays for internal use.
  logical :: is_ortholog_fortran(n_genes)
  logical :: is_outlier_fortran(n_genes)
  logical :: init_perm_fortran
  integer :: i
  
  ! Convert C integers to Fortran logical
  do i = 1, n_genes
    is_ortholog_fortran(i) = (is_ortholog_int(i) /= 0)
  end do
  
  ! Convert scalar
  init_perm_fortran = (init_perm_int /= 0)
  
  ! Call the main routine with Fortran logical types
  call detect_outliers(n_genes, n_families, distances, gene_to_fam, &
                      is_ortholog_fortran, percentile, work_array, perm, stack_left, stack_right, &
                      is_outlier_fortran, rdi_values, rdi_threshold, init_perm_fortran)
  
  ! Convert Fortran logical back to C integers for output
  do i = 1, n_genes
    is_outlier_int(i) = merge(1, 0, is_outlier_fortran(i))
  end do
  
end subroutine detect_outliers_c_int
