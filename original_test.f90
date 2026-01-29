program loess_from_tsv_netlib
  use f42_utils, only: sort_real

  implicit none
  integer, parameter :: dp = kind(1.0d0)

  character(len=*), parameter :: in_tsv  = "material/kallisto_sex_data.tsv"
  character(len=*), parameter :: out_csv = "gene_mean_std_loess.csv"

  ! -----------------------------
  ! CONFIG: spans + robust toggle
  ! -----------------------------
  integer, parameter :: n_spans = 5
  real(dp), parameter :: spans(n_spans) = (/ 0.20_dp, 0.40_dp, 0.60_dp, 0.80_dp, 0.90_dp /)

  logical, parameter :: ROBUST = .false.          ! <--- CAMBIA AQUÍ
  integer, parameter :: N_ROBUST_ITERS = 3        ! típico 2..4

  ! LOESS params
  integer, parameter :: d = 1
  integer :: degree, nvmax
  logical :: infl, setlf

  ! Data
  integer :: n_samples, n_genes_total, n_valid
  character(len=256), allocatable :: gene_id(:)
  real(dp), allocatable :: x(:), y(:)        ! valid only
  real(dp), allocatable :: w(:), z(:,:)      ! loess inputs
  real(dp), allocatable :: yhat(:,:)         ! (n_valid, n_spans)
  integer, allocatable :: order(:), stack_left(:), stack_right(:)

  ! LOESS workspace
  integer :: liv, lv, i
  integer, allocatable :: iv(:)
  real(dp), allocatable :: wv(:)
  real(dp), allocatable :: diagl(:)

  ! Robust working arrays
  real(dp), allocatable :: rw(:), ww(:), res(:)
  integer, allocatable  :: pi(:)

  integer :: sidx

  ! ---- LOESS netlib externals ----
  external :: lowesd, lowesb, lowese, lowesw

  ! -----------------------------
  ! Basic settings
  ! -----------------------------
  degree = 2
  infl   = .false.     ! infl = trace(L)/influence extras; NO es robust
  setlf  = .false.

  call scan_tsv_dims(in_tsv, n_samples, n_genes_total)
  print *, "TSV:", trim(in_tsv)
  print *, "  n_samples =", n_samples
  print *, "  n_genes_total =", n_genes_total

  allocate(gene_id(n_genes_total))
  call count_valid_genes(in_tsv, n_samples, n_genes_total, gene_id, n_valid)
  print *, "  n_valid (std>0) =", n_valid
  if (n_valid < 5) stop "Too few valid points for LOESS."

  allocate(x(n_valid), y(n_valid))
  call fill_xy(in_tsv, n_samples, n_genes_total, gene_id, n_valid, x, y)

  ! Sort by x for nicer output + more stable curves
  allocate(order(n_valid))
  allocate(stack_left(n_valid), stack_right(n_valid))
  order = [(i, i=1,n_valid)]
  call sort_real(x, order, stack_left, stack_right)
  x = x(order)
  y = y(order)

  allocate(w(n_valid), z(n_valid,1), yhat(n_valid, n_spans))
  w = 1.0_dp
  z(:,1) = x

  ! conservative workspace sizes
  liv  = max(50000, 200*n_valid)
  lv   = max(500000, 500*n_valid)
  allocate(iv(liv), wv(lv), diagl(n_valid))

  nvmax = max(200, 4*n_valid)

  ! robust buffers
  allocate(rw(n_valid), ww(n_valid), res(n_valid), pi(n_valid))
  rw = 1.0_dp
  ww = 1.0_dp
  res = 0.0_dp
  pi = [(sidx, sidx=1,n_valid)]   ! 1..n

  print *, "Running LOESS netlib..."
  print *, "  ROBUST =", ROBUST, "  (iters=", N_ROBUST_ITERS, ")"
  do sidx = 1, n_spans
    print '(A,F5.2)', "  span = ", spans(sidx)

    if (.not.ROBUST) then
      call loess_fit_predict_plain(n_valid, x, y, w, z, spans(sidx), degree, nvmax, infl, setlf, &
                                   iv, liv, wv, lv, diagl, yhat(:,sidx))
    else
      call loess_fit_predict_robust(n_valid, x, y, w, z, spans(sidx), degree, nvmax, infl, setlf, &
                                    N_ROBUST_ITERS, iv, liv, wv, lv, diagl, &
                                    rw, ww, res, pi, yhat(:,sidx))
    end if
  end do

  call write_csv(out_csv, x, y, yhat, spans, n_valid, n_spans)
  print *, "Wrote:", trim(out_csv)

contains

  ! ============================================================
  ! Plain LOESS
  ! ============================================================
  subroutine loess_fit_predict_plain(n, x, y, w, z, span, degree, nvmax, infl, setlf, &
                                     iv, liv, wv, lv, diagl, yhat)
    integer, intent(in) :: n, degree, nvmax, liv, lv
    real(dp), intent(in) :: x(n), y(n), w(n), z(n,1), span
    logical, intent(in) :: infl, setlf
    integer, intent(inout) :: iv(liv)
    real(dp), intent(inout) :: wv(lv)
    real(dp), intent(inout) :: diagl(n)
    real(dp), intent(out) :: yhat(n)

    call lowesd(106, iv, liv, lv, wv, 1, n, span, degree, nvmax, setlf)
    print *, "lowesd called with: liv=", liv, ", lv=", lv, ", n=", n, ", span=", span, ", degree=", degree, ", nvmax=", nvmax
    call lowesb(x, y, w, diagl, infl, iv, liv, lv, wv)
    print *, "lowesb called with: size(x)=", size(x), ", size(y)=", size(y), ", size(w)=", size(w), ", size(diagl)=", size(diagl), ", liv=", liv, ", lv=", lv
    call lowese(iv, liv, lv, wv, n, z, yhat)
    print *, "lowese called with: liv=", liv, ", lv=", lv, ", n=", n, ", size(z)=", size(z), ", size(yhat)=", size(yhat)
  end subroutine

  ! ============================================================
  ! Robust LOESS (bisquare reweighting via lowesw)
  ! ============================================================
  subroutine loess_fit_predict_robust(n, x, y, w, z, span, degree, nvmax, infl, setlf, &
                                      n_iters, iv, liv, wv, lv, diagl, rw, ww, res, pi, yhat)
    integer, intent(in) :: n, degree, nvmax, n_iters, liv, lv
    real(dp), intent(in) :: x(n), y(n), w(n), z(n,1), span
    logical, intent(in) :: infl, setlf
    integer, intent(inout) :: iv(liv)
    real(dp), intent(inout) :: wv(lv)
    real(dp), intent(inout) :: diagl(n)
    real(dp), intent(inout) :: rw(n), ww(n), res(n)
    integer, intent(inout) :: pi(n)
    real(dp), intent(out) :: yhat(n)

    integer :: it, i

    rw = 1.0_dp

    do it = 1, n_iters
      do i=1,n
        ww(i) = w(i) * rw(i)
      end do

      ! IMPORTANT: re-init each iteration (fresh state)
      call lowesd(106, iv, liv, lv, wv, 1, n, span, degree, nvmax, setlf)
      call lowesb(x, y, ww, diagl, infl, iv, liv, lv, wv)
      call lowese(iv, liv, lv, wv, n, z, yhat)

      do i=1,n
        res(i) = y(i) - yhat(i)
      end do

      call lowesw(res, n, rw, pi)
    end do
  end subroutine

  ! ============================================================
  ! TSV parsing + helpers (tus rutinas casi iguales)
  ! ============================================================
  subroutine scan_tsv_dims(fname, n_samples, n_genes)
    character(len=*), intent(in) :: fname
    integer, intent(out) :: n_samples, n_genes
    integer :: u, ios
    character(len=100000) :: line
    integer :: ncols

    open(newunit=u, file=fname, status='old', action='read', iostat=ios)
    if (ios /= 0) stop "Cannot open TSV"

    read(u, '(A)', iostat=ios) line
    if (ios /= 0) stop "Cannot read header"

    ncols = 1 + count_tabs(line)
    n_samples = max(0, ncols - 1)

    n_genes = 0
    do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) > 0) n_genes = n_genes + 1
    end do
    close(u)
  end subroutine

  integer function count_tabs(line)
    character(len=*), intent(in) :: line
    integer :: k
    count_tabs = 0
    do k=1,len_trim(line)
      if (line(k:k) == char(9)) count_tabs = count_tabs + 1
    end do
  end function

  subroutine count_valid_genes(fname, n_samples, n_genes_total, gene_id, n_valid)
    character(len=*), intent(in) :: fname
    integer, intent(in) :: n_samples, n_genes_total
    character(len=256), intent(out) :: gene_id(n_genes_total)
    integer, intent(out) :: n_valid

    integer :: u, ios, i
    character(len=100000) :: line
    character(len=256) :: gid
    real(dp) :: mean, std

    n_valid = 0
    open(newunit=u, file=fname, status='old', action='read', iostat=ios)
    if (ios /= 0) stop "Cannot open TSV"
    read(u, '(A)') line ! header

    do i=1,n_genes_total
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      call compute_mean_std_from_line(line, n_samples, mean, std, gid)
      gene_id(i) = gid
      if (std > 0.0_dp .and. mean /= 0.0_dp) n_valid = n_valid + 1
    end do
    close(u)
  end subroutine

  subroutine fill_xy(fname, n_samples, n_genes_total, gene_id, n_valid, x, y)
    character(len=*), intent(in) :: fname
    integer, intent(in) :: n_samples, n_genes_total, n_valid
    character(len=256), intent(in) :: gene_id(n_genes_total)
    real(dp), intent(out) :: x(n_valid), y(n_valid)

    integer :: u, ios, i, j
    character(len=100000) :: line
    character(len=256) :: gid
    real(dp) :: mean, std

    j = 0
    open(newunit=u, file=fname, status='old', action='read', iostat=ios)
    if (ios /= 0) stop "Cannot open TSV"
    read(u, '(A)') line ! header

    do i=1,n_genes_total
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit

      call compute_mean_std_from_line(line, n_samples, mean, std, gid)
      if (std > 0.0_dp .and. mean /= 0.0_dp) then
        j = j + 1
        x(j) = mean
        y(j) = std
      end if
      if (j == n_valid) exit
    end do
    close(u)
  end subroutine

  subroutine compute_mean_std_from_line(line_in, n_samples, mean, std, gid)
    character(len=*), intent(in) :: line_in
    integer, intent(in) :: n_samples
    real(dp), intent(out) :: mean, std
    character(len=256), intent(out) :: gid

    character(len=100000) :: tmp
    integer :: pos, k, ios, nobs
    character(len=256) :: tok
    real(dp) :: x, m, s, delta

    tmp = line_in
    call tabs_to_spaces(tmp)

    pos = 1
    call next_token(tmp, pos, gid)

    m = 0.0_dp
    s = 0.0_dp
    nobs = 0

    do k=1,n_samples
      call next_token(tmp, pos, tok)
      if (len_trim(tok) == 0) exit
      read(tok, *, iostat=ios) x
      if (ios /= 0) cycle
      nobs = nobs + 1
      delta = x - m
      m = m + delta / real(nobs, dp)
      s = s + delta * (x - m)
    end do

    if (nobs >= 2) then
      mean = m
      std  = sqrt(s / real(nobs - 1, dp))
    else
      mean = 0.0_dp
      std  = 0.0_dp
    end if
  end subroutine

  subroutine tabs_to_spaces(str)
    character(len=*), intent(inout) :: str
    integer :: i, n
    n = len_trim(str)
    do i=1,n
      if (str(i:i) == char(9)) str(i:i) = ' '
    end do
  end subroutine

  subroutine next_token(str, pos, tok)
    character(len=*), intent(in) :: str
    integer, intent(inout) :: pos
    character(len=256), intent(out) :: tok
    integer :: n, i, j

    n = len_trim(str)
    tok = ""

    do while (pos <= n .and. str(pos:pos) == ' ')
      pos = pos + 1
    end do
    if (pos > n) return

    i = pos
    do while (pos <= n .and. str(pos:pos) /= ' ')
      pos = pos + 1
    end do
    j = pos - 1
    tok = str(i:min(j, i+255))
  end subroutine

  subroutine argsort_real(a, idx, n)
    real(dp), intent(in) :: a(n)
    integer, intent(out) :: idx(n)
    integer, intent(in) :: n
    integer :: i, j, t
    do i=1,n
      idx(i) = i
    end do
    do i=2,n
      t = idx(i)
      j = i - 1
      do while (j >= 1 .and. a(idx(j)) > a(t))
        idx(j+1) = idx(j)
        j = j - 1
      end do
      idx(j+1) = t
    end do
  end subroutine

  subroutine apply_order_inplace(x, y, idx, n)
    real(dp), intent(inout) :: x(n), y(n)
    integer, intent(in) :: idx(n), n
    real(dp), allocatable :: xt(:), yt(:)
    integer :: i
    allocate(xt(n), yt(n))
    do i=1,n
      xt(i) = x(idx(i))
      yt(i) = y(idx(i))
    end do
    x = xt
    y = yt
    deallocate(xt, yt)
  end subroutine

  subroutine write_csv(fname, x, y, yhat, spans, n, ns)
    character(len=*), intent(in) :: fname
    real(dp), intent(in) :: x(n), y(n), yhat(n,ns), spans(ns)
    integer, intent(in) :: n, ns
    integer :: u, i, s
    character(len=64) :: col

    open(newunit=u, file=fname, status='replace', action='write')
    write(u,'(A)', advance='no') 'x_mean,y_std'
    do s=1,ns
      write(col,'(A,F0.3)') ',y_smooth_span_', spans(s)
      call sanitize(col)
      write(u,'(A)', advance='no') trim(col)
    end do
    write(u,'(A)') ""

    do i=1,n
      write(u,'(ES24.15, A, ES24.15)', advance='no') x(i), ',', y(i)
      do s=1,ns
        write(u,'(A, ES24.15)', advance='no') ',', yhat(i,s)
      end do
      write(u,'(A)') ""
    end do
    close(u)
  end subroutine

  subroutine sanitize(s)
    character(len=*), intent(inout) :: s
    integer :: i
    do i=1,len_trim(s)
      if (s(i:i) == '.') s(i:i) = 'p'
    end do
  end subroutine

end program loess_from_tsv_netlib
