! filepath: test/mod_test_jensen_shannon.f90
!> Unit test suite for Jensen-Shannon divergence and related statistical routines.
module mod_test_jensen_shannon
  use asserts
  use tox_jensen_shannon_test
  use tox_errors
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
  implicit none
  public

  ! Abstract interface for all test procedures
  abstract interface
    subroutine test_interface()
    end subroutine test_interface
  end interface

  ! Type to hold test name and procedure pointer
  type :: test_case
    character(len=64) :: name
    procedure(test_interface), pointer, nopass :: test_proc => null()
  end type test_case

  real(real64), parameter :: TOL = 1d-12

contains

  !> Get array of all available tests.
  function get_all_tests() result(all_tests)
    type(test_case) :: all_tests(14)
    all_tests(1) = test_case("test_compute_gene_means_basic", test_compute_gene_means_basic)
    all_tests(2) = test_case("test_compute_gene_means_with_nan", test_compute_gene_means_with_nan)
    all_tests(3) = test_case("test_compute_gene_means_all_nan", test_compute_gene_means_all_nan)
    all_tests(4) = test_case("test_compute_gene_means_invalid_input", test_compute_gene_means_invalid_input)
    
    all_tests(5) = test_case("test_compute_residuals_basic", test_compute_residuals_basic)
    all_tests(6) = test_case("test_compute_residuals_with_nan", test_compute_residuals_with_nan)
    all_tests(7) = test_case("test_compute_residuals_all_nan", test_compute_residuals_all_nan)
    all_tests(8) = test_case("test_compute_residuals_invalid_input", test_compute_residuals_invalid_input)
    
    all_tests(9) = test_case("test_pool_means_alloc_basic", test_pool_means_alloc_basic)
    all_tests(10) = test_case("test_pool_means_alloc_with_nan", test_pool_means_alloc_with_nan)
    all_tests(11) = test_case("test_pool_means_alloc_single_study", test_pool_means_alloc_single_study)
    all_tests(12) = test_case("test_pool_means_alloc_invalid_input", test_pool_means_alloc_invalid_input)
    
    all_tests(13) = test_case("test_construct_neighborhoods_basic", test_construct_neighborhoods_basic)
    all_tests(14) = test_case("test_construct_neighborhoods_nan_means", test_construct_neighborhoods_nan_means)
  end function get_all_tests

  !> Run all Jensen-Shannon tests.
  subroutine run_all_tests_jensen_shannon()
    type(test_case), allocatable :: all_tests(:)
    integer(int32) :: i
    all_tests = get_all_tests()
    do i = 1, size(all_tests)
      call all_tests(i)%test_proc()
      print *, trim(all_tests(i)%name), " passed."
    end do
    print *, "All Jensen-Shannon tests passed successfully."
  end subroutine run_all_tests_jensen_shannon

  !> Run specific tests by name.
  subroutine run_named_tests_jensen_shannon(test_names)
    character(len=*), intent(in) :: test_names(:)
    type(test_case), allocatable :: all_tests(:)
    integer(int32) :: i, j
    logical :: found
    all_tests = get_all_tests()
    do i = 1, size(test_names)
      found = .false.
      do j = 1, size(all_tests)
        if (trim(test_names(i)) == trim(all_tests(j)%name)) then
          call all_tests(j)%test_proc()
          print *, trim(test_names(i)), " passed."
          found = .true.
          exit
        end if
      end do
      if (.not. found) then
        print *, "Unknown test: ", trim(test_names(i))
      end if
    end do
  end subroutine run_named_tests_jensen_shannon

  ! --------------------------------------------------------------------------
  ! Helper Functions
  ! --------------------------------------------------------------------------

  !> Convert integer to string for error messages.
  function str(i) result(res)
    integer, intent(in) :: i
    character(len=32) :: res
    write(res, *) i
    res = adjustl(res)
  end function str

  ! --------------------------------------------------------------------------
  ! Test Cases for compute_gene_means
  ! --------------------------------------------------------------------------

  ! Test case 1: Basic compute_gene_means functionality.
  subroutine test_compute_gene_means_basic()
    integer, parameter :: n_genes = 4, n_reps = 3
    real(real64) :: expr(n_reps, n_genes), means(n_genes)
    real(real64) :: expected_means(n_genes)
    integer(int32) :: ierr
    
    ! Test data
    expr = reshape([1.0, 2.0, 3.0,    &   ! Gene 1: mean = 2.0
                    4.0, 5.0, 6.0,    &   ! Gene 2: mean = 5.0
                    10.0, 20.0, 30.0, &   ! Gene 3: mean = 20.0
                    0.0, 0.0, 0.0],   &   ! Gene 4: mean = 0.0
                   [n_reps, n_genes])
    
    expected_means = [2.0, 5.0, 20.0, 0.0]
    
    call compute_gene_means(n_genes, n_reps, expr, means, ierr)
    
    call assert_equal_int(ierr, ERR_OK, "test_compute_gene_means_basic: should succeed")
    call assert_allclose_array_real(means, expected_means, n_genes, 0.0_real64, &
                                    TOL, "test_compute_gene_means_basic: means")
  end subroutine test_compute_gene_means_basic

  ! Test case 2: compute_gene_means with NaN values.
  subroutine test_compute_gene_means_with_nan()
    integer, parameter :: n_genes = 3, n_reps = 4
    real(real64) :: expr(n_reps, n_genes), means(n_genes)
    integer(int32) :: ierr
    
    expr(:, 1) = [1.0_real64, 2.0_real64, ieee_value(0.0_real64, ieee_quiet_nan), 3.0_real64]  ! mean = (1+2+3)/3 = 2.0
    expr(:, 2) = [ieee_value(0.0_real64, ieee_quiet_nan), 5.0_real64, 7.0_real64, 9.0_real64]  ! mean = (5+7+9)/3 = 7.0
    expr(:, 3) = [10.0, 20.0, 30.0, 40.0]  ! mean = 25.0
    
    call compute_gene_means(n_genes, n_reps, expr, means, ierr)
    
    call assert_equal_int(ierr, ERR_OK, "test_compute_gene_means_with_nan: should succeed")
    call assert_equal_real(means(1), 2.0_real64, TOL, "test_compute_gene_means_with_nan: gene 1 mean")
    call assert_equal_real(means(2), 7.0_real64, TOL, "test_compute_gene_means_with_nan: gene 2 mean")
    call assert_equal_real(means(3), 25.0_real64, TOL, "test_compute_gene_means_with_nan: gene 3 mean")
  end subroutine test_compute_gene_means_with_nan

  ! Test case 3: compute_gene_means with all NaN values for a gene.
  subroutine test_compute_gene_means_all_nan()
    integer, parameter :: n_genes = 2, n_reps = 3
    real(real64) :: expr(n_reps, n_genes), means(n_genes)
    integer(int32) :: ierr
    
    expr(:, 1) = [1.0, 2.0, 3.0]  ! Normal gene
    expr(:, 2) = ieee_value(0.0_real64, ieee_quiet_nan)  ! All NaN gene
    
    call compute_gene_means(n_genes, n_reps, expr, means, ierr)
    
    call assert_equal_int(ierr, ERR_OK, "test_compute_gene_means_all_nan: should succeed")
    call assert_equal_real(means(1), 2.0_real64, TOL, "test_compute_gene_means_all_nan: gene 1 mean")
    call assert_true(ieee_is_nan(means(2)), "test_compute_gene_means_all_nan: gene 2 should be NaN")
  end subroutine test_compute_gene_means_all_nan

  ! Test case 4: compute_gene_means with invalid input.
  subroutine test_compute_gene_means_invalid_input()
    integer, parameter :: n_genes = 0, n_reps = 3, n_genes_neg = -1
    real(real64) :: expr(3, 1), means(1)
    integer(int32) :: ierr
    
    ! Test with zero genes
    call compute_gene_means(n_genes, n_reps, expr, means, ierr)
    call assert_not_equal_int(ierr, ERR_OK, "test_compute_gene_means_invalid_input: zero genes should fail")
    
    ! Test with negative genes
    call compute_gene_means(n_genes_neg, n_reps, expr, means, ierr)
    call assert_not_equal_int(ierr, ERR_OK, "test_compute_gene_means_invalid_input: negative genes should fail")
    
    ! Test with zero replicates
    call compute_gene_means(n_genes, 0, expr, means, ierr)
    call assert_not_equal_int(ierr, ERR_OK, "test_compute_gene_means_invalid_input: zero replicates should fail")
  end subroutine test_compute_gene_means_invalid_input

  ! --------------------------------------------------------------------------
  ! Test Cases for compute_residuals
  ! --------------------------------------------------------------------------

  ! Test case 5: Basic compute_residuals functionality.
  subroutine test_compute_residuals_basic()
    integer, parameter :: n_genes = 4, n_reps = 3
    real(real64) :: expr(n_reps, n_genes), means(n_genes), resid(n_reps, n_genes)
    real(real64) :: expected_resid(n_reps, n_genes)
    integer(int32) :: ierr
    
    expr = reshape([1.0, 2.0, 3.0,    &   ! Gene 1
                    4.0, 5.0, 6.0,    &   ! Gene 2
                    10.0, 20.0, 30.0, &   ! Gene 3
                    0.0, 0.0, 0.0],   &   ! Gene 4
                   [n_reps, n_genes])
    
    means = [2.0, 5.0, 20.0, 0.0]
    expected_resid = reshape([-1.0, 0.0, 1.0,     &   ! Gene 1 residuals
                              -1.0, 0.0, 1.0,     &   ! Gene 2 residuals
                              -10.0, 0.0, 10.0,   &   ! Gene 3 residuals
                              0.0, 0.0, 0.0],     &   ! Gene 4 residuals
                             [n_reps, n_genes])
    
    call compute_residuals(n_genes, n_reps, expr, means, resid, ierr)
    
    call assert_equal_int(ierr, ERR_OK, "test_compute_residuals_basic: should succeed")
    call assert_allclose_array_real(reshape(resid, [n_reps*n_genes]), &
                                    reshape(expected_resid, [n_reps*n_genes]), &
                                    n_reps*n_genes, 0.0_real64, TOL, &
                                    "test_compute_residuals_basic: residuals")
  end subroutine test_compute_residuals_basic

  ! Test case 6: compute_residuals with NaN values.
  subroutine test_compute_residuals_with_nan()
    integer, parameter :: n_genes = 2, n_reps = 4
    real(real64) :: expr(n_reps, n_genes), means(n_genes), resid(n_reps, n_genes)
    integer(int32) :: ierr
    
    expr(:, 1) = [1.0_real64, 2.0_real64, ieee_value(0.0_real64, ieee_quiet_nan), 3.0_real64]
    expr(:, 2) = [ieee_value(0.0_real64, ieee_quiet_nan), 5.0_real64, 7.0_real64, 9.0_real64]
    means = [2.0_real64, 7.0_real64]
    
    call compute_residuals(n_genes, n_reps, expr, means, resid, ierr)
    
    call assert_equal_int(ierr, ERR_OK, "test_compute_residuals_with_nan: should succeed")
    ! Check specific values
    call assert_equal_real(resid(1, 1), -1.0_real64, TOL, "test_compute_residuals_with_nan: resid(1,1)")
    call assert_equal_real(resid(2, 1), 0.0_real64, TOL, "test_compute_residuals_with_nan: resid(2,1)")
    call assert_true(ieee_is_nan(resid(3, 1)), "test_compute_residuals_with_nan: resid(3,1) should be NaN")
    call assert_equal_real(resid(4, 1), 1.0_real64, TOL, "test_compute_residuals_with_nan: resid(4,1)")
    
    call assert_true(ieee_is_nan(resid(1, 2)), "test_compute_residuals_with_nan: resid(1,2) should be NaN")
    call assert_equal_real(resid(2, 2), -2.0_real64, TOL, "test_compute_residuals_with_nan: resid(2,2)")
    call assert_equal_real(resid(3, 2), 0.0_real64, TOL, "test_compute_residuals_with_nan: resid(3,2)")
    call assert_equal_real(resid(4, 2), 2.0_real64, TOL, "test_compute_residuals_with_nan: resid(4,2)")
  end subroutine test_compute_residuals_with_nan

  ! Test case 7: compute_residuals with all NaN values.
  subroutine test_compute_residuals_all_nan()
    integer, parameter :: n_genes = 2, n_reps = 3
    real(real64) :: expr(n_reps, n_genes), means(n_genes), resid(n_reps, n_genes)
    integer(int32) :: ierr
    
    expr(:, 1) = [1.0, 2.0, 3.0]
    expr(:, 2) = ieee_value(0.0_real64, ieee_quiet_nan)  ! All NaN
    means = [2.0, 0.0]  ! Second mean is irrelevant
    
    call compute_residuals(n_genes, n_reps, expr, means, resid, ierr)
    
    call assert_equal_int(ierr, ERR_OK, "test_compute_residuals_all_nan: should succeed")
    ! All residuals for gene 2 should be NaN
    call assert_true(all(ieee_is_nan(resid(:, 2))), "test_compute_residuals_all_nan: all residuals for NaN gene should be NaN")
  end subroutine test_compute_residuals_all_nan

  ! Test case 8: compute_residuals with invalid input.
  subroutine test_compute_residuals_invalid_input()
    integer, parameter :: n_genes = 0, n_reps = 3
    real(real64) :: expr(3, 1), means(1), resid(3, 1)
    integer(int32) :: ierr
    
    ! Test with zero genes
    call compute_residuals(n_genes, n_reps, expr, means, resid, ierr)
    call assert_not_equal_int(ierr, ERR_OK, "test_compute_residuals_invalid_input: zero genes should fail")

  end subroutine test_compute_residuals_invalid_input

  ! --------------------------------------------------------------------------
  ! Test Cases for pool_means_alloc
  ! --------------------------------------------------------------------------

  ! Test case 9: Basic pool_means_alloc functionality.
  subroutine test_pool_means_alloc_basic()
    integer, parameter :: n_genes_S1 = 5, n_genes_S2 = 5, n_points = 3
    real(real64) :: mean_S1(n_genes_S1), mean_S2(n_genes_S2), x_star(n_points)
    integer(int32) :: N_pool, ierr
    
    mean_S1 = [1.0, 3.0, 5.0, 7.0, 9.0]
    mean_S2 = [2.0, 4.0, 6.0, 8.0, 10.0]
    
    call pool_means_alloc(n_genes_S1, mean_S1, n_genes_S2, mean_S2, n_points, N_pool, x_star, ierr)
    
    call assert_equal_int(ierr, ERR_OK, "test_pool_means_alloc_basic: should succeed")
    call assert_equal_int(N_pool, 10, "test_pool_means_alloc_basic: N_pool should be 10")
    
    ! Check that x_star contains quantiles from pooled data
    ! Pooled data: [1,2,3,4,5,6,7,8,9,10]
    ! For n_points=3, quantiles at positions: 10/4=2.5, 20/4=5.0, 30/4=7.5
    ! Floored: 2, 5, 7 -> values: 2, 5, 7 -> interpolation to 3.25, 5.5 and 7.75
    call assert_equal_real(x_star(1), 3.25_real64, TOL, "test_pool_means_alloc_basic: first quantile")
    call assert_equal_real(x_star(2), 5.5_real64, TOL, "test_pool_means_alloc_basic: second quantile")
    call assert_equal_real(x_star(3), 7.75_real64, TOL, "test_pool_means_alloc_basic: third quantile")
  end subroutine test_pool_means_alloc_basic

  ! Test case 10: pool_means_alloc with NaN values.
  subroutine test_pool_means_alloc_with_nan()
    integer, parameter :: n_genes_S1 = 4, n_genes_S2 = 4, n_points = 2
    real(real64) :: mean_S1(n_genes_S1), mean_S2(n_genes_S2), x_star(n_points)
    integer(int32) :: N_pool, ierr
    
    mean_S1 = [1.0_real64, ieee_value(0.0_real64, ieee_quiet_nan), 3.0_real64, 5.0_real64]
    mean_S2 = [2.0_real64, 4.0_real64, ieee_value(0.0_real64, ieee_quiet_nan), 6.0_real64]
    
    call pool_means_alloc(n_genes_S1, mean_S1, n_genes_S2, mean_S2, n_points, N_pool, x_star, ierr)
    
    call assert_equal_int(ierr, ERR_OK, "test_pool_means_alloc_with_nan: should succeed")
    call assert_equal_int(N_pool, 6, "test_pool_means_alloc_with_nan: N_pool should exclude NaN values")
    
    ! Pooled data (excluding NaN): [1,2,3,4,5,6]
    ! Values: 2.666, 4.3333 -> interpolation
    call assert_equal_real(x_star(1), 2.0_real64 + 2.0_real64/3.0_real64, TOL, "test_pool_means_alloc_with_nan: first quantile")
    call assert_equal_real(x_star(2), 4.0_real64 + 1.0_real64/3.0_real64, TOL, "test_pool_means_alloc_with_nan: second quantile")
  end subroutine test_pool_means_alloc_with_nan

  ! Test case 11: pool_means_alloc with single study.
  subroutine test_pool_means_alloc_single_study()
    integer, parameter :: n_genes_S1 = 5, n_genes_S2 = 1, n_points = 3
    real(real64) :: mean_S1(n_genes_S1), mean_S2(1), x_star(n_points)
    integer(int32) :: N_pool, ierr
    
    mean_S1 = [1.0, 2.0, 3.0, 4.0, 5.0]
    mean_S2 = [0.0]  ! Dummy
    
    call pool_means_alloc(n_genes_S1, mean_S1, n_genes_S2, mean_S2, n_points, N_pool, x_star, ierr)
    
    call assert_equal_int(ierr, ERR_OK, "test_pool_means_alloc_single_study: should succeed")
    ! Does not work with a single study
  end subroutine test_pool_means_alloc_single_study

  ! Test case 12: pool_means_alloc with invalid input.
  subroutine test_pool_means_alloc_invalid_input()
    integer, parameter :: n_genes_S1 = 0, n_genes_S2 = 5, n_points = 3
    real(real64) :: mean_S1(1), mean_S2(5), x_star(n_points)
    integer(int32) :: N_pool, ierr
    
    mean_S2 = [1.0, 2.0, 3.0, 4.0, 5.0]
    
    ! Test with zero genes in S1
    call pool_means_alloc(n_genes_S1, mean_S1, n_genes_S2, mean_S2, n_points, N_pool, x_star, ierr)
    call assert_not_equal_int(ierr, ERR_OK, "test_pool_means_alloc_invalid_input: zero genes in S1 should fail")
    
    ! Test with zero points
    call pool_means_alloc(5, mean_S2, n_genes_S2, mean_S2, 0, N_pool, x_star, ierr)
    call assert_not_equal_int(ierr, ERR_OK, "test_pool_means_alloc_invalid_input: zero points should fail")
  end subroutine test_pool_means_alloc_invalid_input

  subroutine test_construct_neighborhoods_basic()
      integer(int32), parameter :: n_points    = 2
      integer(int32), parameter :: n_genes_S   = 5
      integer(int32), parameter :: n_reps_S    = 3
      integer(int32), parameter :: n_neighbors = 2

      integer(int32) :: ierr
      real(real64) :: x_star(n_points)
      real(real64) :: mean_S(n_genes_S)
      real(real64) :: resid_S(n_reps_S, n_genes_S)
      real(real64) :: tmp_distances(n_genes_S)
      integer(int32) :: tmp_distances_perm(n_genes_S)
      real(real64) :: neighborhood_residuals(n_reps_S, n_neighbors, n_points)
      integer(int32) :: neighborhood_indices(n_neighbors, n_points)

      ! -----------------------------
      ! Inputs
      ! -----------------------------
      x_star = [ 2.0_real64, 10.0_real64 ]
      mean_S = [ 1.0, 2.5, 9.0, 10.5, 20.0 ]

      resid_S = reshape([ &
          1.0,  2.0,  3.0,  4.0,  5.0, &   ! rep 1
          10.0,20.0,30.0,40.0,50.0, &   ! rep 2
          -1.0,-2.0,-3.0,-4.0,-5.0  &    ! rep 3
      ], shape(resid_S))

      ! -----------------------------
      ! Call routine
      ! -----------------------------
      call construct_neighborhoods( &
          n_points, x_star, n_genes_S, mean_S, n_reps_S, resid_S, &
          tmp_distances, tmp_distances_perm, &
          neighborhood_residuals, neighborhood_indices, &
          n_neighbors, ierr )

      call assert_equal_int(ierr, ERR_OK, "ierr must be ERR_OK")

      ! -----------------------------
      ! Expected neighbors
      !
      ! For x_star(1)=2.0:
      !   distances = [1.0, 0.5, 7.0, 8.5, 18.0]
      !   sorted → gene 2, gene 1
      !   tie-breaking: ascending gene index  <-- IMPORTANT
      !
      ! For x_star(2)=10.0:
      !   distances = [9.0, 7.5, 1.0, 0.5, 10.0]
      !   sorted → gene 4, gene 3
      ! -----------------------------

      call assert_equal_array_int( neighborhood_indices(:,1), [2,1], n_neighbors, &
          "test_construct_neighborhoods_basic: Incorrect neighborhood indices for point 1" )

      call assert_equal_array_int( neighborhood_indices(:,2), [4,3], n_neighbors, &
          "test_construct_neighborhoods_basic: Incorrect neighborhood indices for point 2" )

      ! -----------------------------
      ! Expected residuals
      ! -----------------------------
      call assert_equal_array_real( neighborhood_residuals(:,1,1), resid_S(:, 2), n_reps_S, TOL, &
          "test_construct_neighborhoods_basic: Incorrect residuals(:,1,1)" )

      call assert_equal_array_real( neighborhood_residuals(:,2,1), resid_S(:, 1), n_reps_S, TOL, &
          "test_construct_neighborhoods_basic: Incorrect residuals(:,2,1)" )

      call assert_equal_array_real( neighborhood_residuals(:,1,2), resid_S(:, 4), n_reps_S, TOL, &
          "test_construct_neighborhoods_basic: Incorrect residuals(:,1,2)" )

      call assert_equal_array_real( neighborhood_residuals(:,2,2), resid_S(:, 3), n_reps_S, TOL, &
          "test_construct_neighborhoods_basic: Incorrect residuals(:,2,2)" )

      call construct_neighborhoods( &
          0_int32, x_star, n_genes_S, mean_S, n_reps_S, resid_S, &
          tmp_distances, tmp_distances_perm, &
          neighborhood_residuals, neighborhood_indices, &
          n_neighbors, ierr )

      call assert_equal_int(ierr, ERR_EMPTY_INPUT, "ierr must be ERR_EMPTY_INPUT")

      call construct_neighborhoods( &
          n_points, x_star, 0_int32, mean_S, n_reps_S, resid_S, &
          tmp_distances, tmp_distances_perm, &
          neighborhood_residuals, neighborhood_indices, &
          n_neighbors, ierr )

      call assert_equal_int(ierr, ERR_EMPTY_INPUT, "ierr must be ERR_EMPTY_INPUT")

      call construct_neighborhoods( &
          n_points, x_star, n_genes_S, mean_S, 0_int32, resid_S, &
          tmp_distances, tmp_distances_perm, &
          neighborhood_residuals, neighborhood_indices, &
          n_neighbors, ierr )

      call assert_equal_int(ierr, ERR_EMPTY_INPUT, "ierr must be ERR_EMPTY_INPUT")

      call construct_neighborhoods( &
          n_points, x_star, n_genes_S, mean_S, n_reps_S, resid_S, &
          tmp_distances, tmp_distances_perm, &
          neighborhood_residuals, neighborhood_indices, &
          0_int32, ierr )

      call assert_equal_int(ierr, ERR_EMPTY_INPUT, "ierr must be ERR_EMPTY_INPUT")
  end subroutine test_construct_neighborhoods_basic

  subroutine test_construct_neighborhoods_nan_means()
      integer(int32), parameter :: n_points    = 1
      integer(int32), parameter :: n_genes_S   = 4
      integer(int32), parameter :: n_reps_S    = 2
      integer(int32), parameter :: n_neighbors = 2

      integer(int32) :: ierr
      real(real64) :: x_star(n_points)
      real(real64) :: mean_S(n_genes_S)
      real(real64) :: resid_S(n_reps_S, n_genes_S)
      real(real64) :: tmp_distances(n_genes_S)
      integer(int32) :: tmp_distances_perm(n_genes_S)
      real(real64) :: neighborhood_residuals(n_reps_S, n_neighbors, n_points)
      integer(int32) :: neighborhood_indices(n_neighbors, n_points)

      x_star = [ 5.0_real64 ]
      mean_S = [ 4.0_real64, ieee_value(1.0_real64, ieee_quiet_nan), 6.0_real64, ieee_value(1.0_real64, ieee_quiet_nan) ]

      resid_S = reshape([ &
          1.0, 2.0, 3.0, 4.0, &
          10.0,20.0,30.0,40.0 &
      ], shape(resid_S))

      call construct_neighborhoods( &
          n_points, x_star, n_genes_S, mean_S, n_reps_S, resid_S, &
          tmp_distances, tmp_distances_perm, &
          neighborhood_residuals, neighborhood_indices, &
          n_neighbors, ierr )

      call assert_equal_int(ierr, ERR_OK, "ierr must be ERR_OK")

      ! Only genes 1 and 3 are valid (non-NaN)
      call assert_equal_array_int( neighborhood_indices(:,1), [1,3], n_neighbors, &
          "test_construct_neighborhoods_nan_means: NaN mean handling incorrect" )
  end subroutine test_construct_neighborhoods_nan_means
end module mod_test_jensen_shannon