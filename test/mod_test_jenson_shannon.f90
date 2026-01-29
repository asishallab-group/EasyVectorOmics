! filepath: test/mod_test_jenson_shannon.f90
!> Unit test suite for Jensen-Shannon divergence and related statistical routines.
module mod_test_jenson_shannon
  use asserts
  use tox_jenson_shannon_test
  use tox_errors, only: ERR_INVALID_INPUT, ERR_EMPTY_INPUT
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

contains

  !> Get array of all available tests.
  function get_all_tests() result(all_tests)
    type(test_case) :: all_tests(18)
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
    all_tests(14) = test_case("test_construct_neighborhoods_with_nan", test_construct_neighborhoods_with_nan)
    all_tests(15) = test_case("test_construct_neighborhoods_kx_limits", test_construct_neighborhoods_kx_limits)
    all_tests(16) = test_case("test_construct_neighborhoods_small_dataset", test_construct_neighborhoods_small_dataset)
    all_tests(17) = test_case("test_construct_neighborhoods_large_kx", test_construct_neighborhoods_large_kx)
    all_tests(18) = test_case("test_construct_neighborhoods_edge_cases", test_construct_neighborhoods_edge_cases)
  end function get_all_tests

  !> Run all Jensen-Shannon tests.
  subroutine run_all_tests_jenson_shannon()
    type(test_case) :: all_tests(18)
    integer(int32) :: i
    all_tests = get_all_tests()
    do i = 1, size(all_tests)
      call all_tests(i)%test_proc()
      print *, trim(all_tests(i)%name), " passed."
    end do
    print *, "All Jensen-Shannon tests passed successfully."
  end subroutine run_all_tests_jenson_shannon

  !> Run specific tests by name.
  subroutine run_named_tests_jenson_shannon(test_names)
    character(len=*), intent(in) :: test_names(:)
    type(test_case) :: all_tests(18)
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
  end subroutine run_named_tests_jenson_shannon

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
    
    call assert_equal_int(ierr, 0, "test_compute_gene_means_basic: should succeed")
    call assert_allclose_array_real(means, expected_means, n_genes, 0.0_real64, &
                                    1e-9_real64, "test_compute_gene_means_basic: means")
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
    
    call assert_equal_int(ierr, 0, "test_compute_gene_means_with_nan: should succeed")
    call assert_equal_real(means(1), 2.0_real64, 1e-9_real64, "test_compute_gene_means_with_nan: gene 1 mean")
    call assert_equal_real(means(2), 7.0_real64, 1e-9_real64, "test_compute_gene_means_with_nan: gene 2 mean")
    call assert_equal_real(means(3), 25.0_real64, 1e-9_real64, "test_compute_gene_means_with_nan: gene 3 mean")
  end subroutine test_compute_gene_means_with_nan

  ! Test case 3: compute_gene_means with all NaN values for a gene.
  subroutine test_compute_gene_means_all_nan()
    integer, parameter :: n_genes = 2, n_reps = 3
    real(real64) :: expr(n_reps, n_genes), means(n_genes)
    integer(int32) :: ierr
    
    expr(:, 1) = [1.0, 2.0, 3.0]  ! Normal gene
    expr(:, 2) = ieee_value(0.0_real64, ieee_quiet_nan)  ! All NaN gene
    
    call compute_gene_means(n_genes, n_reps, expr, means, ierr)
    
    call assert_equal_int(ierr, 0, "test_compute_gene_means_all_nan: should succeed")
    call assert_equal_real(means(1), 2.0_real64, 1e-9_real64, "test_compute_gene_means_all_nan: gene 1 mean")
    call assert_true(ieee_is_nan(means(2)), "test_compute_gene_means_all_nan: gene 2 should be NaN")
  end subroutine test_compute_gene_means_all_nan

  ! Test case 4: compute_gene_means with invalid input.
  subroutine test_compute_gene_means_invalid_input()
    integer, parameter :: n_genes = 0, n_reps = 3, n_genes_neg = -1
    real(real64) :: expr(3, 1), means(1)
    integer(int32) :: ierr
    
    ! Test with zero genes
    call compute_gene_means(n_genes, n_reps, expr, means, ierr)
    call assert_true(ierr /= 0, "test_compute_gene_means_invalid_input: zero genes should fail")
    
    ! Test with negative genes
    call compute_gene_means(n_genes_neg, n_reps, expr, means, ierr)
    call assert_true(ierr /= 0, "test_compute_gene_means_invalid_input: negative genes should fail")
    
    ! Test with zero replicates
    call compute_gene_means(n_genes, 0, expr, means, ierr)
    call assert_true(ierr /= 0, "test_compute_gene_means_invalid_input: zero replicates should fail")
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
    
    call assert_equal_int(ierr, 0, "test_compute_residuals_basic: should succeed")
    call assert_allclose_array_real(reshape(resid, [n_reps*n_genes]), &
                                    reshape(expected_resid, [n_reps*n_genes]), &
                                    n_reps*n_genes, 0.0_real64, 1e-9_real64, &
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
    
    call assert_equal_int(ierr, 0, "test_compute_residuals_with_nan: should succeed")
    ! Check specific values
    call assert_equal_real(resid(1, 1), -1.0_real64, 1e-9_real64, "test_compute_residuals_with_nan: resid(1,1)")
    call assert_equal_real(resid(2, 1), 0.0_real64, 1e-9_real64, "test_compute_residuals_with_nan: resid(2,1)")
    call assert_true(ieee_is_nan(resid(3, 1)), "test_compute_residuals_with_nan: resid(3,1) should be NaN")
    call assert_equal_real(resid(4, 1), 1.0_real64, 1e-9_real64, "test_compute_residuals_with_nan: resid(4,1)")
    
    call assert_true(ieee_is_nan(resid(1, 2)), "test_compute_residuals_with_nan: resid(1,2) should be NaN")
    call assert_equal_real(resid(2, 2), -2.0_real64, 1e-9_real64, "test_compute_residuals_with_nan: resid(2,2)")
    call assert_equal_real(resid(3, 2), 0.0_real64, 1e-9_real64, "test_compute_residuals_with_nan: resid(3,2)")
    call assert_equal_real(resid(4, 2), 2.0_real64, 1e-9_real64, "test_compute_residuals_with_nan: resid(4,2)")
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
    
    call assert_equal_int(ierr, 0, "test_compute_residuals_all_nan: should succeed")
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
    call assert_true(ierr /= 0, "test_compute_residuals_invalid_input: zero genes should fail")

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
    
    call assert_equal_int(ierr, 0, "test_pool_means_alloc_basic: should succeed")
    call assert_equal_int(N_pool, 10, "test_pool_means_alloc_basic: N_pool should be 10")
    
    ! Check that x_star contains quantiles from pooled data
    ! Pooled data: [1,2,3,4,5,6,7,8,9,10]
    ! For n_points=3, quantiles at positions: 10/4=2.5, 20/4=5.0, 30/4=7.5
    ! Floored: 2, 5, 7 -> values: 2, 5, 7 -> interpolation to 3.25, 5.5 and 7.75
    call assert_equal_real(x_star(1), 3.25_real64, 1e-9_real64, "test_pool_means_alloc_basic: first quantile")
    call assert_equal_real(x_star(2), 5.5_real64, 1e-9_real64, "test_pool_means_alloc_basic: second quantile")
    call assert_equal_real(x_star(3), 7.75_real64, 1e-9_real64, "test_pool_means_alloc_basic: third quantile")
  end subroutine test_pool_means_alloc_basic

  ! Test case 10: pool_means_alloc with NaN values.
  subroutine test_pool_means_alloc_with_nan()
    integer, parameter :: n_genes_S1 = 4, n_genes_S2 = 4, n_points = 2
    real(real64) :: mean_S1(n_genes_S1), mean_S2(n_genes_S2), x_star(n_points)
    integer(int32) :: N_pool, ierr
    
    mean_S1 = [1.0_real64, ieee_value(0.0_real64, ieee_quiet_nan), 3.0_real64, 5.0_real64]
    mean_S2 = [2.0_real64, 4.0_real64, ieee_value(0.0_real64, ieee_quiet_nan), 6.0_real64]
    
    call pool_means_alloc(n_genes_S1, mean_S1, n_genes_S2, mean_S2, n_points, N_pool, x_star, ierr)
    
    call assert_equal_int(ierr, 0, "test_pool_means_alloc_with_nan: should succeed")
    call assert_equal_int(N_pool, 6, "test_pool_means_alloc_with_nan: N_pool should exclude NaN values")
    
    ! Pooled data (excluding NaN): [1,2,3,4,5,6]
    ! Values: 2.666, 4.3333 -> interpolation
    call assert_equal_real(x_star(1), 2.0_real64 + 2.0_real64/3.0_real64, 1e-9_real64, "test_pool_means_alloc_with_nan: first quantile")
    call assert_equal_real(x_star(2), 4.0_real64 + 1.0_real64/3.0_real64, 1e-9_real64, "test_pool_means_alloc_with_nan: second quantile")
  end subroutine test_pool_means_alloc_with_nan

  ! Test case 11: pool_means_alloc with single study.
  subroutine test_pool_means_alloc_single_study()
    integer, parameter :: n_genes_S1 = 5, n_genes_S2 = 0, n_points = 3
    real(real64) :: mean_S1(n_genes_S1), mean_S2(1), x_star(n_points)
    integer(int32) :: N_pool, ierr
    
    mean_S1 = [1.0, 2.0, 3.0, 4.0, 5.0]
    mean_S2 = [0.0]  ! Dummy
    
    call pool_means_alloc(n_genes_S1, mean_S1, n_genes_S2, mean_S2, n_points, N_pool, x_star, ierr)
    
    call assert_equal_int(ierr, ERR_INVALID_INPUT, "test_pool_means_alloc_single_study: should succeed")
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
    call assert_true(ierr /= 0, "test_pool_means_alloc_invalid_input: zero genes in S1 should fail")
    
    ! Test with zero points
    call pool_means_alloc(5, mean_S2, n_genes_S2, mean_S2, 0, N_pool, x_star, ierr)
    call assert_true(ierr /= 0, "test_pool_means_alloc_invalid_input: zero points should fail")
  end subroutine test_pool_means_alloc_invalid_input

  ! --------------------------------------------------------------------------
  ! Test Cases for construct_neighborhoods
  ! --------------------------------------------------------------------------

  ! Test case 13: Basic construct_neighborhoods functionality.
  subroutine test_construct_neighborhoods_basic()
    integer, parameter :: n_points = 2, n_genes_S = 10, n_reps_S = 3
    integer(int32) :: N_pool, k_x
    real(real64) :: x_star(n_points), mean_S(n_genes_S), resid_S(n_reps_S, n_genes_S)
    real(real64) :: distances(n_genes_S)
    integer(int32) :: work_indices(n_genes_S)
    real(real64) :: neighborhood_residuals(n_points, n_reps_S * 1000)
    integer(int32) :: neighborhood_indices(n_points, 1000)
    integer(int32) :: i, ierr, neighborhood_size
    
    x_star = [5.0, 15.0]
    mean_S = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
    
    do i = 1, n_genes_S
      resid_S(1, i) = real(i, real64)
      resid_S(2, i) = -real(i, real64)
      resid_S(3, i) = 0.0_real64
    end do
    
    N_pool = n_genes_S * 2

    neighborhood_size = 10
    
    call construct_neighborhoods(n_points, x_star, n_genes_S, mean_S, n_reps_S, resid_S, &
                                 distances, work_indices, N_pool, k_x, &
                                 neighborhood_residuals, neighborhood_indices, ierr, neighborhood_size)
    
    call assert_true(k_x == 10, "test_construct_neighborhoods_basic: k_x within limits")
    call assert_true(any(neighborhood_indices(1, 1:k_x) == 5), "test_construct_neighborhoods_basic: gene 5 should be nearest")
    call assert_true(any(neighborhood_indices(2, 1:k_x) == 10), "test_construct_neighborhoods_basic: gene 10 should be nearest")
    call assert_true(.not. all(ieee_is_nan(neighborhood_residuals)), "test_construct_neighborhoods_basic: residuals should exist")
  end subroutine test_construct_neighborhoods_basic

  ! Test case 14: construct_neighborhoods with NaN mean values.
  subroutine test_construct_neighborhoods_with_nan()
    integer, parameter :: n_points = 2, n_genes_S = 8, n_reps_S = 2
    integer(int32) :: N_pool, k_x
    real(real64) :: x_star(n_points), mean_S(n_genes_S), resid_S(n_reps_S, n_genes_S)
    real(real64) :: distances(n_genes_S)
    integer(int32) :: work_indices(n_genes_S)
    real(real64) :: neighborhood_residuals(n_points, n_reps_S * 1000)
    integer(int32) :: neighborhood_indices(n_points, 1000)
    integer(int32) :: i, ierr
    
    x_star = [3.0_real64, 7.0_real64]
    mean_S = [1.0_real64, 2.0_real64, ieee_value(0.0_real64, ieee_quiet_nan), 4.0_real64, &
              5.0_real64, ieee_value(0.0_real64, ieee_quiet_nan), 7.0_real64, 8.0_real64]
    
    do i = 1, n_genes_S
      if (.not. ieee_is_nan(mean_S(i))) then
        resid_S(:, i) = [real(i, real64), -real(i, real64)]
      else
        resid_S(:, i) = ieee_value(0.0_real64, ieee_quiet_nan)
      end if
    end do
    
    N_pool = n_genes_S * 2
    
    call construct_neighborhoods(n_points, x_star, n_genes_S, mean_S, n_reps_S, resid_S, &
                                 distances, work_indices, N_pool, k_x, &
                                 neighborhood_residuals, neighborhood_indices, ierr)
    
    call assert_true(.not. any(neighborhood_indices == 3), "test_construct_neighborhoods_with_nan: gene 3 (NaN) should not be in neighborhood")
    call assert_true(.not. any(neighborhood_indices == 6), "test_construct_neighborhoods_with_nan: gene 6 (NaN) should not be in neighborhood")
  end subroutine test_construct_neighborhoods_with_nan

  ! Test case 15: Test k_x calculation limits.
  subroutine test_construct_neighborhoods_kx_limits()
    integer, parameter :: n_points = 10, n_genes_S = 1000, n_reps_S = 3
    integer(int32) :: N_pool, k_x
    real(real64) :: x_star(n_points), mean_S(n_genes_S), resid_S(n_reps_S, n_genes_S)
    real(real64) :: distances(n_genes_S)
    integer(int32) :: work_indices(n_genes_S)
    real(real64) :: neighborhood_residuals(n_points, n_reps_S * 1000)
    integer(int32) :: neighborhood_indices(n_points, 1000)
    integer(int32) :: i, ierr
    
    do i = 1, n_genes_S
      mean_S(i) = real(i, real64) / 10.0
      resid_S(:, i) = [1.0, -1.0, 0.0]
    end do
    x_star = [(real(i, real64), i = 1, n_points)]
    
    ! Test lower bound: N_pool/(2*n_points) = 1000/(2*10) = 50, but min is 100
    N_pool = 1000
    call construct_neighborhoods(n_points, x_star, n_genes_S, mean_S, n_reps_S, resid_S, &
                                 distances, work_indices, N_pool, k_x, &
                                 neighborhood_residuals, neighborhood_indices, ierr)
    call assert_equal_int(k_x, 100, "test_construct_neighborhoods_kx_limits: lower bound k_x=100")
    
    ! Test upper bound: N_pool/(2*n_points) = 1000000/(2*10) = 50000, but max is 1000
    N_pool = 1000000
    call construct_neighborhoods(n_points, x_star, n_genes_S, mean_S, n_reps_S, resid_S, &
                                 distances, work_indices, N_pool, k_x, &
                                 neighborhood_residuals, neighborhood_indices, ierr)
    call assert_equal_int(k_x, 1000, "test_construct_neighborhoods_kx_limits: upper bound k_x=1000")
  end subroutine test_construct_neighborhoods_kx_limits

  ! Test case 16: construct_neighborhoods with small dataset.
  subroutine test_construct_neighborhoods_small_dataset()
    integer, parameter :: n_points = 3, n_genes_S = 50, n_reps_S = 2
    integer(int32) :: N_pool, k_x
    real(real64) :: x_star(n_points), mean_S(n_genes_S), resid_S(n_reps_S, n_genes_S)
    real(real64) :: distances(n_genes_S)
    integer(int32) :: work_indices(n_genes_S)
    real(real64) :: neighborhood_residuals(n_points, n_reps_S * 1000)
    integer(int32) :: neighborhood_indices(n_points, 1000)
    integer(int32) :: i, ierr
    
    do i = 1, n_genes_S
      mean_S(i) = real(i, real64)
      resid_S(:, i) = [real(i, real64), -real(i, real64)]
    end do
    x_star = [10.0, 25.0, 40.0]
    
    N_pool = n_genes_S * 2
    
    call construct_neighborhoods(n_points, x_star, n_genes_S, mean_S, n_reps_S, resid_S, &
                                 distances, work_indices, N_pool, k_x, &
                                 neighborhood_residuals, neighborhood_indices, ierr)
    
    call assert_true(k_x >= 100, "test_construct_neighborhoods_small_dataset: k_x >= 100")
    call assert_true(.not. all(neighborhood_indices == -1), "test_construct_neighborhoods_small_dataset: indices should not be all -1")
  end subroutine test_construct_neighborhoods_small_dataset

  ! Test case 17: construct_neighborhoods with k_x larger than number of genes.
  subroutine test_construct_neighborhoods_large_kx()
    integer, parameter :: n_points = 1, n_genes_S = 50, n_reps_S = 2
    integer(int32) :: N_pool, k_x
    real(real64) :: x_star(n_points), mean_S(n_genes_S), resid_S(n_reps_S, n_genes_S)
    real(real64) :: distances(n_genes_S)
    integer(int32) :: work_indices(n_genes_S)
    real(real64) :: neighborhood_residuals(n_points, n_reps_S * 1000)
    integer(int32) :: neighborhood_indices(n_points, 1000)
    integer(int32) :: i, ierr
    
    do i = 1, n_genes_S
      mean_S(i) = real(i, real64)
      resid_S(:, i) = [1.0, -1.0]
    end do
    x_star = [25.0]
    
    N_pool = 200  ! N_pool/(2*n_points) = 200/2 = 100, which is > n_genes_S
    
    call construct_neighborhoods(n_points, x_star, n_genes_S, mean_S, n_reps_S, resid_S, &
                                 distances, work_indices, N_pool, k_x, &
                                 neighborhood_residuals, neighborhood_indices, ierr)
    
    call assert_equal_int(k_x, 100, "test_construct_neighborhoods_large_kx: k_x should be 100")
    call assert_true(all(neighborhood_indices(1, 51:1000) == -1), "test_construct_neighborhoods_large_kx: indices beyond n_genes_S should be -1")
  end subroutine test_construct_neighborhoods_large_kx

  ! Test case 18: Edge cases for construct_neighborhoods.
  subroutine test_construct_neighborhoods_edge_cases()
    integer, parameter :: n_points = 1, n_genes_S = 1, n_reps_S = 1
    integer(int32) :: N_pool, k_x
    real(real64) :: x_star(n_points), mean_S(n_genes_S), resid_S(n_reps_S, n_genes_S)
    real(real64) :: distances(n_genes_S)
    integer(int32) :: work_indices(n_genes_S)
    real(real64) :: neighborhood_residuals(n_points, n_reps_S * 1000)
    integer(int32) :: neighborhood_indices(n_points, 1000)
    integer(int32) :: ierr
    
    x_star = [5.0]
    mean_S = [5.0]  ! Exact match
    resid_S = reshape([2.5], [n_reps_S, n_genes_S])  ! Proper 2D array
    
    N_pool = 100
    
    call construct_neighborhoods(n_points, x_star, n_genes_S, mean_S, n_reps_S, resid_S, &
                                 distances, work_indices, N_pool, k_x, &
                                 neighborhood_residuals, neighborhood_indices, ierr)
    
    call assert_equal_int(neighborhood_indices(1, 1), 1, "test_construct_neighborhoods_edge_cases: single gene should be selected")
    call assert_equal_real(neighborhood_residuals(1, 1), 2.5_real64, 1e-9_real64, "test_construct_neighborhoods_edge_cases: residual should be extracted")
    call assert_true(all(neighborhood_indices(1, 2:1000) == -1), "test_construct_neighborhoods_edge_cases: remaining indices should be -1")
  end subroutine test_construct_neighborhoods_edge_cases

end module mod_test_jenson_shannon