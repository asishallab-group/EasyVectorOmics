! filepath: test/mod_test_normalization_pipeline.f90
!> Unit test suite for normalization_pipeline routine.
module mod_test_normalization_pipeline
  use asserts
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use tox_normalization
  use test_suite, only: test_case
  implicit none
  public


contains

  !> Get array of all available tests.
  function get_all_tests_normalization_pipeline() result(all_tests)
    type(test_case),allocatable :: all_tests(:)
    allocate(all_tests(4))
    all_tests(1) = test_case("test_pipeline_basic", test_pipeline_basic)
    all_tests(2) = test_case("test_pipeline_edge_cases", test_pipeline_edge_cases)
    all_tests(3) = test_case("test_pipeline_vs_manual", test_pipeline_vs_manual)
    all_tests(4) = test_case("test_pipeline_empty_matrix", test_pipeline_empty_matrix)
  end function get_all_tests_normalization_pipeline
 
  

  !> Basic test: pipeline with small matrix, no fold change
  subroutine test_pipeline_basic()
    integer(int32), parameter :: n_genes = 10, n_tissues = 6, n_grps = 2, max_stack = 20
    real(real64), dimension(n_genes * n_tissues) :: input_matrix, buf_stddev, buf_quant
    real(real64), dimension(n_genes * n_grps) :: buf_avg, buf_log
    real(real64), dimension(n_genes) :: temp_col, rank_means, loess_x, loess_y, yhat_global
    integer(int32), dimension(n_genes) :: perm, indices_used
    integer(int32), dimension(max_stack) :: stack_left, stack_right
    integer(int32), dimension(n_grps) :: group_s, group_c
    real(real64) :: span
    integer(int32) :: degree, ierr, i, j

    do j = 1, n_tissues
      do i = 1, n_genes
        input_matrix((j - 1) * n_genes + i) = dble(i) + dble(j) * 0.25d0
      end do
    end do
    group_s = [1,4]
    group_c = [3,3]
    span = 0.75d0
    degree = 2

    call normalization_pipeline(n_genes, n_tissues, input_matrix, buf_stddev, buf_quant, buf_avg, buf_log, &
                                temp_col, rank_means, perm, stack_left, stack_right, max_stack, &
                                group_s, group_c, n_grps, loess_x, loess_y, indices_used, yhat_global, &
                                span, degree, use_quantile=1, ierr=ierr)
    call assert_equal_int(ierr, 0, "normalization_pipeline returned error")
    call assert_no_nan_real(buf_log, n_genes * n_grps, "test_pipeline_basic: NaN in output")
    call assert_equal_int(size(buf_log), n_genes * n_grps, "test_pipeline_basic: output size incorrect")
    call assert_true(all(buf_log >= 0.0d0), "test_pipeline_basic: output should be non-negative")
  end subroutine test_pipeline_basic

  !> Edge case test: zero input, check output shape and values
  subroutine test_pipeline_edge_cases()
    integer(int32), parameter :: n_genes = 10, n_tissues = 6, n_grps = 2, max_stack = 20
    real(real64), dimension(n_genes * n_tissues) :: input_matrix, buf_stddev, buf_quant
    real(real64), dimension(n_genes * n_grps) :: buf_avg, buf_log
    real(real64), dimension(n_genes) :: temp_col, rank_means, loess_x, loess_y, yhat_global
    integer(int32), dimension(n_genes) :: perm, indices_used
    integer(int32), dimension(max_stack) :: stack_left, stack_right
    integer(int32), dimension(n_grps) :: group_s, group_c
    real(real64) :: span
    integer(int32) :: degree, ierr

    input_matrix = 0.0d0
    group_s = [1,4]
    group_c = [3,3]
    span = 0.75d0
    degree = 2

    call normalization_pipeline(n_genes, n_tissues, input_matrix, buf_stddev, buf_quant, buf_avg, buf_log, &
                                temp_col, rank_means, perm, stack_left, stack_right, max_stack, &
                                group_s, group_c, n_grps, loess_x, loess_y, indices_used, yhat_global, &
                                span, degree, use_quantile=1, ierr=ierr)
    call assert_true(ierr /= 0, "test_pipeline_edge_cases: expected error for zero-variance input")
  end subroutine test_pipeline_edge_cases

  !> Test: pipeline output matches manual stepwise normalization (no fold change)
  subroutine test_pipeline_vs_manual()
    integer(int32), parameter :: n_genes = 10, n_tissues = 6, n_grps = 2, max_stack = 20
    real(real64), dimension(n_genes * n_tissues) :: input_matrix, buf_stddev, buf_quant
    real(real64), dimension(n_genes * n_grps) :: buf_avg, buf_log, manual_out
    real(real64), dimension(n_genes * n_grps) :: buf_avg_no_quant, buf_log_no_quant, manual_out_no_quant
    real(real64), dimension(n_genes) :: temp_col, rank_means, loess_x, loess_y, yhat_global
    integer(int32), dimension(n_genes) :: perm, indices_used
    integer(int32), dimension(max_stack) :: stack_left, stack_right
    integer(int32), dimension(n_grps) :: group_s, group_c
    real(real64), dimension(n_genes, n_tissues) :: input_2d, stddev_2d
    real(real64) :: span
    integer(int32) :: degree, ierr, i, j

    ! Build varied data with non-zero variance per gene
    do j = 1, n_tissues
      do i = 1, n_genes
        input_matrix((j - 1) * n_genes + i) = dble(i) * 2.0d0 + dble(j) * 0.5d0
      end do
    end do
    group_s = [1,4]
    group_c = [3,3]
    span = 0.75d0
    degree = 2

    ! Manual stepwise normalization
    do j = 1, n_tissues
      do i = 1, n_genes
        input_2d(i, j) = input_matrix((j - 1) * n_genes + i)
      end do
    end do

    call normalize_by_std_dev(n_genes, n_tissues, input_2d, stddev_2d, &
                              loess_x, loess_y, indices_used, yhat_global, &
                              span, degree, ierr)
    call assert_equal_int(ierr, 0, "normalize_by_std_dev returned error")

    do j = 1, n_tissues
      do i = 1, n_genes
        buf_stddev((j - 1) * n_genes + i) = stddev_2d(i, j)
      end do
    end do

    call quantile_normalization(n_genes, n_tissues, buf_stddev, buf_quant, temp_col, rank_means, perm, stack_left, stack_right, max_stack, ierr)
    call assert_equal_int(ierr, 0, "quantile_normalization returned error")
    call calc_tiss_avg(n_genes, n_grps, group_s, group_c, buf_quant, buf_avg, ierr)
    call assert_equal_int(ierr, 0, "calc_tiss_avg returned error")
    call log2_transformation(n_genes, n_grps, buf_avg, manual_out, ierr)
    call assert_equal_int(ierr, 0, "log2_transformation returned error")

    ! Manual stepwise normalization without quantile
    call calc_tiss_avg(n_genes, n_grps, group_s, group_c, buf_stddev, buf_avg_no_quant, ierr)
    call assert_equal_int(ierr, 0, "calc_tiss_avg (no quantile) returned error")
    call log2_transformation(n_genes, n_grps, buf_avg_no_quant, manual_out_no_quant, ierr)
    call assert_equal_int(ierr, 0, "log2_transformation (no quantile) returned error")

    ! Pipeline normalization
    call normalization_pipeline(n_genes, n_tissues, input_matrix, buf_stddev, buf_quant, buf_avg, buf_log, &
                                temp_col, rank_means, perm, stack_left, stack_right, max_stack, &
                                group_s, group_c, n_grps, loess_x, loess_y, indices_used, yhat_global, &
                                span, degree, use_quantile=1, ierr=ierr)
    call assert_equal_int(ierr, 0, "normalization_pipeline returned error")
    call assert_true(all(abs(buf_log - manual_out) < 1.0d-12), "test_pipeline_vs_manual: pipeline and manual outputs differ")

    call normalization_pipeline(n_genes, n_tissues, input_matrix, buf_stddev, buf_quant, buf_avg_no_quant, buf_log_no_quant, &
                  temp_col, rank_means, perm, stack_left, stack_right, max_stack, &
                  group_s, group_c, n_grps, loess_x, loess_y, indices_used, yhat_global, &
                  span, degree, use_quantile=0, ierr=ierr)
    call assert_equal_int(ierr, 0, "normalization_pipeline (no quantile) returned error")
    call assert_true(all(abs(buf_log_no_quant - manual_out_no_quant) < 1.0d-12), "test_pipeline_vs_manual: pipeline and manual outputs differ (no quantile)")
  end subroutine test_pipeline_vs_manual

  !> Test with empty input matrix.
  subroutine test_pipeline_empty_matrix()
    integer(int32) :: n_genes, n_tissues, n_grps, max_stack, degree, ierr
    real(real64) :: span
    real(real64), dimension(0) :: input_matrix, buf_stddev, buf_quant, buf_avg, buf_log
    real(real64), dimension(0) :: temp_col, rank_means, loess_x, loess_y, yhat_global
    integer(int32), dimension(0) :: perm, stack_left, stack_right, group_s, group_c, indices_used
    n_genes = 0; n_tissues = 0; n_grps = 0; max_stack = 0
    span = 0.75d0
    degree = 2
    call normalization_pipeline(n_genes, n_tissues, input_matrix, buf_stddev, buf_quant, buf_avg, buf_log, &
                                temp_col, rank_means, perm, stack_left, stack_right, max_stack, &
                                group_s, group_c, n_grps, loess_x, loess_y, indices_used, yhat_global, &
                                span, degree, ierr=ierr)
    call assert_equal_int(ierr, 202, "normalization_pipeline should return error for empty input")
    ! No further assertion needed: just check no crash
  end subroutine test_pipeline_empty_matrix

end module mod_test_normalization_pipeline
