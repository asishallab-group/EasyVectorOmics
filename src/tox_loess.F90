module tox_loess
  use safeguard
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use tox_errors, only: set_ok, set_err, is_err, validate_dimension_size, validate_in_range_real, validate_all_in_range_int, ERR_INVALID_INPUT
  
  ! ---- LOESS netlib externals ----
  ! ============================================================
  ! LOESS Subroutines from Netlib
  ! ============================================================
  ! These subroutines are part of the LOESS implementation provided by the Netlib library.
  ! They are used to perform various operations such as decomposition, fitting, evaluation,
  ! and robust weight computation for LOESS models. The subroutines are interfaced here
  ! to integrate the Netlib routines into the Tensor Omics project.
  !
  ! Subroutine Descriptions:
  ! - lowesd: Initializes workspace arrays and performs LOESS decomposition.
  ! - lowesb: Fits the LOESS model and computes the diagonal elements of the hat matrix.
  ! - lowese: Evaluates the LOESS model and computes the smoothed response variable array.
  ! - lowesw: Computes robust weights for LOESS using residuals.
  ! ============================================================
  interface
    ! ============================================================
    ! Subroutine: lowesd
    ! ============================================================
    !> Perform LOESS decomposition.
    !| This subroutine computes the decomposition of the LOESS model.
    !| It initializes the integer and real workspace arrays based on the input parameters.
    subroutine lowesd(mode, iv, liv, lv, wv, dim, n, span, degree, nvmax, setlf)
      use, intrinsic :: iso_fortran_env, only: real64, int32
      integer(int32), intent(in) :: mode
      !| Mode of operation
      integer(int32), intent(in) :: liv
      !| Length of the integer workspace array
      integer(int32), intent(in) :: lv
      !| Length of the real workspace array
      integer(int32), intent(in) :: dim
      !| Dimensionality of the data
      integer(int32), intent(in) :: n
      !| Number of data points
      integer(int32), intent(in) :: degree
      !| Degree of the LOESS polynomial
      integer(int32), intent(in) :: nvmax
      !| Maximum neighborhood size
      real(real64), intent(in) :: span
      !| Smoothing parameter for LOESS
      logical, intent(in) :: setlf
      !| Save matrix factorization flag
      integer(int32), intent(inout) :: iv(liv)
      !| Integer workspace array
      real(real64), intent(inout) :: wv(lv)
      !| Real workspace array
    end subroutine lowesd

    ! ============================================================
    ! Subroutine: lowesb
    ! ============================================================
    !> Perform LOESS fitting.
    !| This subroutine fits the LOESS model to the data.
    !| It uses the input predictor and response variables to compute the diagonal elements of the hat matrix.
    subroutine lowesb(x, y, w, diagl, infl, iv, liv, lv, wv)
      import :: real64, int32
      real(real64), intent(in) :: x(*)
      !| Predictor variable array
      real(real64), intent(in) :: y(*)
      !| Response variable array
      real(real64), intent(in) :: w(*)
      !| Weight array for data points
      real(real64), intent(out) :: diagl(*)
      !| Diagonal elements of the hat matrix
      logical, intent(in) :: infl
      !| Influence calculation flag
      integer(int32), intent(in) :: iv(*)
      !| Integer workspace array
      integer(int32), intent(in) :: liv
      !| Length of the integer workspace array
      integer(int32), intent(in) :: lv
      !| Length of the real workspace array
      real(real64), intent(inout) :: wv(*)
      !| Real workspace array
    end subroutine lowesb

    ! ============================================================
    ! Subroutine: lowese
    ! ============================================================
    !> Perform LOESS evaluation.
    !| This subroutine evaluates the LOESS model.
    !| It computes the smoothed response variable array based on the input predictor variables.
    subroutine lowese(iv, liv, lv, wv, n, z, s)
      import :: real64, int32
      integer(int32), intent(in) :: liv
      !| Length of the integer workspace array
      integer(int32), intent(in) :: lv
      !| Length of the real workspace array
      integer(int32), intent(in) :: n
      !| Number of data points
      integer(int32), intent(inout) :: iv(*)
      !| Integer workspace array
      real(real64), intent(inout) :: wv(*)
      !| Real workspace array
      real(real64), intent(in) :: z(*)
      !| Predictor variable array
      real(real64), intent(out) :: s(*)
      !| Smoothed response variable array
    end subroutine lowese

    ! ============================================================
    ! Subroutine: lowesw
    ! ============================================================
    !> Compute robust weights for LOESS.
    !| This subroutine computes the robust weights for the LOESS model.
    !| It uses the residuals to update the weights for robust fitting.
    subroutine lowesw(res, n, rw, pi)
      use, intrinsic :: iso_fortran_env, only: real64, int32
      real(real64), intent(in) :: res(*)
      !| Residuals array
      integer(int32), intent(in) :: n
      !| Number of data points
      real(real64), intent(out) :: rw(*)
      !| Robust weights array
      integer(int32), intent(out) :: pi(*)
      !| Permutation indices array
    end subroutine lowesw
  end interface

contains

  ! ============================================================
  ! Recommend workspace sizes based on Netlib exact formulas
  ! ============================================================
  !> Recommend workspace sizes based on Netlib exact formulas.
  !| Computes the required sizes for integer and real workspace arrays.
  !| These sizes depend on the dimensionality of the data and the maximum neighborhood size.
  subroutine tox_loess_required_workspace(d, nvmax, liv, lv, setlf)
    !! Recommend sizes for integer pool iv(:) and real pool v(:).
    !! nvmax is the maximum neighborhood size (usually n_train).
    !! setlf indicates if matrix factorizations are saved.
    integer(int32), intent(in) :: d
    !| Dimensionality of the data
    integer(int32), intent(in) :: nvmax
    !| Maximum neighborhood size
    logical, intent(in) :: setlf
    !| Save matrix factorization flag
    integer(int32), intent(out) :: liv
    !| Length of the integer workspace array
    integer(int32), intent(out) :: lv
    !| Length of the real workspace array

    ! Adjusted formulas to ensure sufficient workspace
    ! liv = 50 + (2**d + 6) * nvmax + 100
    ! if (setlf) liv = liv + nvmax

    ! lv = 50 + (3 * d + 4) * nvmax + nvmax + 100
    ! if (setlf) lv = lv + (d + 1) * nvmax
    liv  = max(50000, 500 * nvmax)
    lv   = max(500000, 5000 * nvmax)
  end subroutine tox_loess_required_workspace

  ! ============================================================
  ! Plain LOESS fitting
  ! ============================================================
  !> Perform plain LOESS fitting.
  !| Fits a LOESS model to the data using the specified smoothing parameter.
  !| Outputs the smoothed response variable array.
  subroutine loess_fit_plain(n, x, y, w, z, span, degree, nvmax, infl, setlf, iv, liv, wv, lv, diagl, yhat, ierr)
    integer(int32), intent(in) :: n
    !| Total number of data points
    integer(int32), intent(in) :: degree
    !| Degree of the LOESS polynomial
    integer(int32), intent(in) :: nvmax
    !| Maximum neighborhood size
    integer(int32), intent(in) :: liv
    !| Length of the integer workspace array
    integer(int32), intent(in) :: lv
    !| Length of the real workspace array

    real(real64), intent(in) :: x(n)
    !| Predictor variable array
    real(real64), intent(in) :: y(n)
    !| Response variable array
    real(real64), intent(in) :: w(n)
    !| Weight array for data points
    real(real64), intent(in) :: z(n,1)
    !| Additional predictor variable array
    real(real64), intent(in) :: span
    !| Smoothing parameter for LOESS

    logical, intent(in) :: infl
    !| Influence calculation flag
    logical, intent(in) :: setlf
    !| Save matrix factorization flag

    integer(int32), intent(inout) :: iv(liv)
    !| Integer workspace array
    real(real64), intent(inout) :: wv(lv)
    !| Real workspace array
    real(real64), intent(inout) :: diagl(n)
    !| Diagonal elements of the hat matrix

    real(real64), intent(out) :: yhat(n)
    !| Smoothed response variable array
    integer(int32), intent(out) :: ierr
    !! Error code

    ! Initialize error code
    call set_ok(ierr)

    ! Validate inputs
    call validate_dimension_size(n, ierr)
    call validate_in_range_real(span, ierr, min=0.0_real64)
    if (is_err(ierr)) return

    ! Validate indices
    call validate_all_in_range_int(iv, liv, ierr, min=1, max=n)
    if (is_err(ierr)) return

    ! Validate workspace sizes
    if (liv < 50000 .or. lv < 500000) then
      call set_err(ierr, ERR_INVALID_INPUT)
      return
    end if

    call lowesd(106, iv, liv, lv, wv, 1, n, span, degree, nvmax, setlf)
    call lowesb(x, y, w, diagl, infl, iv, liv, lv, wv)
    call lowese(iv, liv, lv, wv, n, z, yhat)
  end subroutine loess_fit_plain

  ! ============================================================
  ! Robust LOESS fitting (optimized)
  ! ============================================================
  !> Perform robust LOESS fitting with bisquare reweighting.
  !| Fits a LOESS model to the data using robust iterations to handle outliers.
  !| Outputs the smoothed response variable array.
  subroutine loess_fit_robust(n, x, y, w, z, span, degree, nvmax, infl, setlf, n_iters, iv, liv, wv, lv, diagl, rw, ww, res, pi, yhat, ierr)
    integer(int32), intent(in) :: n
    !| Total number of data points
    integer(int32), intent(in) :: degree
    !| Degree of the LOESS polynomial
    integer(int32), intent(in) :: nvmax
    !| Maximum neighborhood size
    integer(int32), intent(in) :: n_iters
    !| Number of robust iterations
    integer(int32), intent(in) :: liv
    !| Length of the integer workspace array
    integer(int32), intent(in) :: lv
    !| Length of the real workspace array

    real(real64), intent(in) :: x(n)
    !| Predictor variable array
    real(real64), intent(in) :: y(n)
    !| Response variable array
    real(real64), intent(in) :: w(n)
    !| Weight array for data points
    real(real64), intent(in) :: z(n,1)
    !| Additional predictor variable array
    real(real64), intent(in) :: span
    !| Smoothing parameter for LOESS

    logical, intent(in) :: infl
    !| Influence calculation flag
    logical, intent(in) :: setlf
    !| Save matrix factorization flag

    integer(int32), intent(inout) :: iv(liv)
    !| Integer workspace array
    real(real64), intent(inout) :: wv(lv)
    !| Real workspace array
    real(real64), intent(inout) :: diagl(n)
    !| Diagonal elements of the hat matrix
    real(real64), intent(inout) :: rw(n)
    !| Robust weights array
    real(real64), intent(inout) :: ww(n)
    !| Working weights array
    real(real64), intent(inout) :: res(n)
    !| Residuals array
    integer(int32), intent(inout) :: pi(n)
    !| Permutation indices array

    real(real64), intent(out) :: yhat(n)
    !| Smoothed response variable array
    integer(int32), intent(out) :: ierr
    !! Error code

    integer :: it, i, dim_val

    ! Initialize error code
    call set_ok(ierr)

    ! Validate inputs
    call validate_dimension_size(n, ierr)
    call validate_in_range_real(span, ierr, min=0.0_real64)
    if (is_err(ierr)) return

    ! Validate indices
    call validate_all_in_range_int(iv, liv, ierr, min=1, max=n)
    if (is_err(ierr)) return

    ! Validate workspace sizes
    if (liv < 50000 .or. lv < 500000) then
      call set_err(ierr, ERR_INVALID_INPUT)
      return
    end if

    rw = 1.0_real64
    dim_val = 1

    do it = 1, n_iters
      do i = 1, n
        ww(i) = w(i) * rw(i)
      end do

      call lowesd(106, iv, liv, lv, wv, dim_val, n, span, degree, nvmax, setlf)
      call lowesb(x, y, ww, diagl, infl, iv, liv, lv, wv)
      call lowese(iv, liv, lv, wv, n, z, yhat)

      do i = 1, n
        res(i) = y(i) - yhat(i)
      end do

      call lowesw(res, n, rw, pi)
    end do
  end subroutine loess_fit_robust

end module tox_loess
