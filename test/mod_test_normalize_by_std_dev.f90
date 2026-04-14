! filepath: test/mod_test_normalize_by_std_dev.f90
!> Unit test suite for normalize_by_std_dev routine.
module mod_test_normalize_by_std_dev
  use asserts
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use tox_normalization
  use test_suite
  implicit none
  public

contains

  !> Get array of all available tests.
  function get_all_tests_normalize_by_std_dev() result(all_tests)
    type(test_case), allocatable :: all_tests(:)
    allocate(all_tests(2))
    
    all_tests(1) = test_case("test_loess_normalization_outlier_correction", test_loess_normalization_outlier_correction)
    all_tests(2) = test_case("test_loess_zero_variance_handling", test_loess_zero_variance_handling)

  end function get_all_tests_normalize_by_std_dev

  !> Main test: Verifies that an SD outlier is corrected by the global curve.
  subroutine test_loess_normalization_outlier_correction()
    integer(int32), parameter :: ng = 20, nt = 10
    real(real64) :: mat(ng, nt), res(ng, nt)
    real(real64) :: lx(ng), ly(ng), yhat(ng)
    integer(int32) :: idx(ng), ierr, i
    real(real64) :: span = 0.75_real64
    integer(int32) :: deg = 1

    ! 1. Create data where SD = Mean (Perfect relationship)
    do i = 1, ng
       ! Gene i has mean i, and values fluctuate to give SD approximately i
       mat(i, :) = real(i, real64) 
       mat(i, 1) = mat(i, 1) + real(i, real64) * 0.5_real64
       mat(i, 2) = mat(i, 2) - real(i, real64) * 0.5_real64
    end do

    ! 2. Introduce an OUTLIER in gene 10
    ! Assign it a huge variance that breaks the linear trend
    mat(10, 1) = 1000.0_real64 

    call normalize_by_std_dev(ng, nt, mat, res, lx, ly, idx, yhat, span, deg, ierr)

    call assert_equal_int(ierr, 0, "LOESS normalization failed")

    ! 3. Verification:
    ! In RMS, gene 10 would be almost 0 due to its huge denominator.
    ! In LOESS, its fitted_sd (yhat) should be close to 10 (the trend), not 300+.
    call assert_true(yhat(10) < 50.0_real64, "LOESS failed to suppress SD outlier")
    call assert_no_nan_real(res, ng*nt, "NaNs in LOESS result")
  end subroutine test_loess_normalization_outlier_correction

    !> Verifies that genes with SD=0 do not break the routine.
    subroutine test_loess_zero_variance_handling()
        integer(int32), parameter :: ng = 10, nt = 5
        real(real64) :: mat(ng, nt), res(ng, nt)
        real(real64) :: lx(ng), ly(ng), yhat(ng)
        integer(int32) :: idx(ng), ierr, i
        
        ! 1. Inicializar todo constante (SD = 0)
        mat = 1.0_real64 

        ! 2. Darle a los primeros 7 genes algo de variabilidad y medias distintas
        !    Usamos el índice 'i' para que las medias (X) sean [1.1, 1.2, 1.3...]
        do i = 1, 7
        mat(i, :) = 1.0_real64 + real(i, real64) * 0.1_real64
        mat(i, 1) = mat(i, 1) + 0.5_real64  ! Genera una SD > 0
        end do

        ! 3. El gen 10 sigue siendo constante (mat(10, :) = 1.0)
        
        call normalize_by_std_dev(ng, nt, mat, res, lx, ly, idx, yhat, 0.75d0, 1, ierr)
        
        ! Verificamos que no haya error
        call assert_equal_int(ierr, 0, "LOESS failed even with valid points")
        
        ! El gen 10 (inválido) debe haber quedado intacto (1.0)
        call assert_equal_real(res(10,1), 1.0_real64, 1d-12, "Zero variance gene altered")
    end subroutine test_loess_zero_variance_handling

end module mod_test_normalize_by_std_dev