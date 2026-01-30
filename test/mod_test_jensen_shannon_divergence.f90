!> Unit test suite for tox_jensen_shannon_divergence routine.
module mod_test_jensen_shannon_divergence
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use tox_jensen_shannon_divergence
    use tox_errors
    use f42_utils, only: above, below
    implicit none

    ! Abstract interface for all test procedures
    abstract interface
        subroutine test_interface()
        end subroutine test_interface
    end interface

    ! Type to hold test name and procedure pointer
    type :: test_case
        character(len=128) :: name
        procedure(test_interface), pointer, nopass :: test_proc => null()
    end type test_case

    real(real64), parameter :: TOL = 1d-12

contains

    !> Get array of all available tests.
    function get_all_tests() result(all_tests)
        type(test_case) :: all_tests(4)
        all_tests(1) = test_case("test_determine_shared_residual_range", test_determine_shared_residual_range)
        all_tests(2) = test_case("test_build_residual_histograms", test_build_residual_histograms)
        all_tests(3) = test_case("test_compute_divergence_per_reference_point", test_compute_divergence_per_reference_point)
        all_tests(4) = test_case("test_compute_weighted_global_divergence", test_compute_weighted_global_divergence)
    end function get_all_tests

    !> Run all tox_jensen_shannon_divergence tests.
    subroutine run_all_tests_tox_jensen_shannon_divergence
        type(test_case), allocatable :: all_tests(:)
        integer(int32) :: i

        all_tests = get_all_tests()

        do i = 1, size(all_tests)
            call all_tests(i)%test_proc()
            print "(' ',A,' passed.')", trim(all_tests(i)%name)
        end do
        print *, "All tox_jensen_shannon_divergence tests passed successfully."
    end subroutine run_all_tests_tox_jensen_shannon_divergence

    !> Run specific tox_jensen_shannon_divergence tests by name.
    subroutine run_named_tests_tox_jensen_shannon_divergence(test_names)
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
                    print "(' ',A,' passed.')", trim(test_names(i))
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                print *, "Unknown test: ", trim(test_names(i))
            end if
        end do
    end subroutine run_named_tests_tox_jensen_shannon_divergence

    subroutine test_determine_shared_residual_range
        integer(int32), parameter :: n_residuals = 4, n_neighbors = 3
        real(real64), dimension(n_residuals, n_neighbors) :: S1, S2
        real(real64) :: R
        integer(int32) :: ierr
        real(real64) :: q

        ! ============================================================
        ! Test 1 — Basic correctness with simple values
        ! ============================================================
        !
        ! S1 abs values: [1,2,3,4,  5,6,7,8,  9,10,11,12]
        ! S2 abs values: [2,4,6,8,  1,3,5,7,  9, 0, 1, 2]
        !
        ! Combined abs pool (24 values):
        !   [1,2,3,4,5,6,7,8,9,10,11,12,
        !    2,4,6,8,1,3,5,7,9,0,1,2]
        !
        ! Sorted:
        !   [0,1,1,1,2,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,11,12]
        !
        ! 95% quantile → 0.95 * (24 - 1)) + 1 = 22.85
        ! sorted(22) = 10
        ! sorted(23) = 11
        !
        ! Expected R = 10 + 0.85 * (11-10) = 10.85
        !
        S1 = reshape([ &
            1,2,3,4,  5,6,-7,8,  9,10,11,12 ], [n_residuals,n_neighbors])
        S2 = reshape([ &
            2,-4,6,8,  1,3,5,7,  9,0,1,2 ], [n_residuals,n_neighbors])

        call determine_shared_residual_range_alloc(S1, S2, n_residuals, n_neighbors, R, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_determine_shared_residual_range: Test 1: ierr should be OK")
        call assert_equal_real(R, 10.85_real64, TOL, "test_determine_shared_residual_range: Test 1: R should be 10.85")

        ! ============================================================
        ! Test 2 — Custom quantile (50%)
        ! ============================================================
        !
        ! Median of sorted array above = 0.5 * (sorted(12) + sorted(13)) = 5.5
        !
        q = 50.0_real64
        call determine_shared_residual_range_alloc(S1, S2, n_residuals, n_neighbors, R, ierr, q)
        call assert_equal_int(ierr, ERR_OK, "test_determine_shared_residual_range: Test 2: ierr should be OK")
        call assert_equal_real(R, 5.0_real64, TOL, "test_determine_shared_residual_range: Test 2: R should be 5.0")

        ! ============================================================
        ! Test 3 — Quantile < 0 → error
        ! ============================================================
        q = below(0.0_real64)
        call determine_shared_residual_range_alloc(S1, S2, n_residuals, n_neighbors, R, ierr, q)
        call assert_equal_int(ierr, ERR_INVALID_INPUT, "test_determine_shared_residual_range: Test 3: ierr should be INVALID_INPUT")

        ! ============================================================
        ! Test 4 — Quantile > 100 → error
        ! ============================================================
        q = above(100.0_real64)
        call determine_shared_residual_range_alloc(S1, S2, n_residuals, n_neighbors, R, ierr, q)
        call assert_equal_int(ierr, ERR_INVALID_INPUT, "test_determine_shared_residual_range: Test 4: ierr should be INVALID_INPUT")

        ! ============================================================
        ! Test 5 — NaNs must be ignored
        ! ============================================================
        !
        ! Replace some values with NaN; remaining values should determine R.
        !
        S1 = reshape([ &
            -1, 2, 3, 4, &
            5, 6, -7, 8, &
            9, 10, -11, 12 ], [n_residuals,n_neighbors])
        S2 = S1
        S1(1,1) = ieee_value(1.0_real64, ieee_quiet_nan)
        S2(4,3) = ieee_value(1.0_real64, ieee_quiet_nan)

        ! Pool now excludes two NaNs → 22 values
        ! sorted = [1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12]
        ! 95% quantile → 0.95*(22-1)+1=20.95
        ! sorted(20) = 11
        ! sorted(21) = 11
        ! -> R = 11
        !
        call determine_shared_residual_range_alloc(S1, S2, n_residuals, n_neighbors, R, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_determine_shared_residual_range: Test 5: ierr should be OK")
        call assert_equal_real(R, 11.0_real64, TOL, "test_determine_shared_residual_range: Test 5: R should ignore NaNs")

        ! ============================================================
        ! Test 6 — All zeros
        ! ============================================================
        S1 = 0.0_real64
        S2 = 0.0_real64
        call determine_shared_residual_range_alloc(S1, S2, n_residuals, n_neighbors, R, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_determine_shared_residual_range: Test 6: ierr should be OK")
        call assert_equal_real(R, 0.0_real64, TOL, "test_determine_shared_residual_range: Test 6: R should be zero")

        ! ============================================================
        ! Test 7 — Single residual (n_residuals=1, n_neighbors=1)
        ! ============================================================
        S1 = 3.0_real64
        S2 = -4.0_real64
        ! sorted = [3, 4]
        ! rank = 0.95 * (2-1) + 1 = 1.95
        ! R = 3 + (4-3)*0.95 = 3.95
        call determine_shared_residual_range_alloc(S1, S2, 1_int32, 1_int32, R, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_determine_shared_residual_range: Test 6: ierr should be OK")
        call assert_equal_real(R, 3.95_real64, TOL, "test_determine_shared_residual_range: Test 7: R should be 3.95")

    end subroutine test_determine_shared_residual_range

    subroutine test_build_residual_histograms
        integer(int32), parameter :: n_residuals = 6, n_neighbors = 3
        integer(int32), parameter :: n_bins = 4
        real(real64), dimension(n_residuals, n_neighbors) :: E
        real(real64), dimension(n_neighbors, n_bins) :: pmf, expected_pmf
        integer(int32), dimension(n_neighbors, n_bins) :: counts, expected_counts
        integer(int32), dimension(n_neighbors) :: included
        real(real64) :: R
        integer(int32) :: ierr

        ! ============================================================
        ! Test 1 — Simple symmetric case, no NaNs
        ! ============================================================
        !
        ! R = 2, M = 4 → bin width w = 1
        ! Bins: [-2,-1), [-1,0), [0,1), [1,2]
        !
        ! Residuals for each neighbor j:
        ! j=1: [-2, -0.5, 0.2, 1.7, 0.9, -1.2]
        ! j=2: [0, 0, 0, 0, 0, 0]
        ! j=3: [2.5, -3.0, 1.2, 0.4, -0.1, 0.0]  (clamping applies)
        !
        R = 2.0_real64

        E(:,1) = [-2.0, -0.5, 0.2, 1.7, 0.9, -1.2]
        E(:,2) = 0.0_real64
        E(:,3) = [2.5, -3.0, 1.2, 0.4, -0.1, 0.0]   ! will clamp to [-2,2]

        call build_residual_histograms(E, n_residuals, n_neighbors, R, n_bins, &
                                       counts, pmf, included, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_build_residual_histograms: Test 1: ierr should be OK")

        ! j = 1
        ! Values fall into bins:
        ! [-2,-1): -2, -1.2 → 2
        ! [-1,0): -0.5 → 1
        ! [0,1): 0.2, 0.9 → 2
        ! [1,2]: 1.7 → 1
        ! 
        ! j = 2 — all zeros → all in bin [0,1)
        ! 
        ! j = 3 — clamping:
        ! 2.5 → 2
        ! -3 → -2
        ! bins:
        ! [-2,-1): -2 → 1
        ! [-1,0): -0.1 → 1
        ! [0,1): 0.4, 0.0 → 2
        ! [1,2]: 1.2, 2 → 2
        expected_counts = reshape([&
            2, 0, 1,&
            1, 0, 1,&
            2, 6, 2,&
            1, 0, 2&
        ], [n_neighbors, n_bins])
        expected_pmf = reshape([&
            0.3333333333333333_real64,  0.0_real64, 0.16666666666666666_real64,&
            0.16666666666666666_real64, 0.0_real64, 0.16666666666666666_real64,&
            0.3333333333333333_real64,  1.0_real64, 0.3333333333333333_real64,&
            0.16666666666666666_real64, 0.0_real64, 0.3333333333333333_real64&
        ], [n_neighbors, n_bins])

        call assert_equal_array_int(counts, expected_counts, size(counts, kind=int32), "test_build_residual_histograms: Test 1: counts don't match")
        call assert_equal_array_real(pmf, expected_pmf, size(pmf, kind=int32), TOL, "test_build_residual_histograms: Test 1: pmf don't match")
        call assert_equal_int(included(1), 6, "test_build_residual_histograms: Test 1: included row 1")
        call assert_equal_int(included(2), 6, "test_build_residual_histograms: Test 1: included row 2")
        call assert_equal_int(included(3), 6, "test_build_residual_histograms: Test 1: included row 3")

        ! ============================================================
        ! Test 2 — NaNs must be ignored
        ! ============================================================
        E = 0.0_real64
        E(1,1) = ieee_value(1.0_real64, ieee_quiet_nan)
        E(3,1) = ieee_value(1.0_real64, ieee_quiet_nan)
        E(2,2) = ieee_value(1.0_real64, ieee_quiet_nan)
        E(6,3) = ieee_value(1.0_real64, ieee_quiet_nan)

        call build_residual_histograms(E, n_residuals, n_neighbors, R, n_bins, &
                                       counts, pmf, included, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_build_residual_histograms: Test 2: ierr should be OK")

        expected_counts = reshape([&
            0, 0, 0,&
            0, 0, 0,&
            4, 5, 5,&
            0, 0, 0&
        ], [n_neighbors, n_bins])
        expected_pmf = reshape([&
            0.0_real64,  0.0_real64, 0.0_real64,&
            0.0_real64,  0.0_real64, 0.0_real64,&
            1.0_real64,  1.0_real64, 1.0_real64,&
            0.0_real64,  0.0_real64, 0.0_real64&
        ], [n_neighbors, n_bins])

        call assert_equal_array_int(counts, expected_counts, size(counts, kind=int32), "test_build_residual_histograms: Test 2: counts don't match")
        call assert_equal_array_real(pmf, expected_pmf, size(pmf, kind=int32), TOL, "test_build_residual_histograms: Test 2: pmf don't match")
        call assert_equal_int(included(1), 4, "test_build_residual_histograms: Test 2: included row 1")
        call assert_equal_int(included(2), 5, "test_build_residual_histograms: Test 2: included row 2")
        call assert_equal_int(included(3), 5, "test_build_residual_histograms: Test 2: included row 3")

        ! ============================================================
        ! Test 3 — All NaN → pmf = 0, counts = 0, included = 0
        ! ============================================================
        E = ieee_value(1.0_real64, ieee_quiet_nan)

        call build_residual_histograms(E, n_residuals, n_neighbors, R, n_bins, &
                                       counts, pmf, included, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_build_residual_histograms: Test 3: ierr should be OK")

        expected_counts = 0.0_real64
        expected_pmf = 0.0_real64
        call assert_equal_array_int(counts, expected_counts, size(counts, kind=int32), "test_build_residual_histograms: Test 3: counts don't match")
        call assert_equal_array_real(pmf, expected_pmf, size(pmf, kind=int32), TOL, "test_build_residual_histograms: Test 3: pmf don't match")
        call assert_equal_int(included(1), 0, "test_build_residual_histograms: Test 2: included row 1")
        call assert_equal_int(included(2), 0, "test_build_residual_histograms: Test 2: included row 2")
        call assert_equal_int(included(3), 0, "test_build_residual_histograms: Test 2: included row 3")

        ! ============================================================
        ! Test 4 — Residuals exactly on boundaries
        ! ============================================================
        !
        ! R = 2, bins = 4, width = 1
        ! Values: -2, -1, 0, 1, 2
        !
        E = reshape([ -2.0, -1.0, 0.0, 1.0, 2.0, 0.0, &
                      -2.0, -1.0, 0.0, 1.0, 2.0, 0.0, &
                      -2.0, -1.0, 0.0, 1.0, 2.0, 0.0 ], [n_residuals,n_neighbors])

        call build_residual_histograms(E, n_residuals, n_neighbors, R, n_bins, &
                                       counts, pmf, included, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_build_residual_histograms: Test 4: ierr should be OK")

        ! Expected binning:
        ! -2 → bin 1
        ! -1 → bin 2
        !  0 → bin 3
        !  1 → bin 4
        !  2 → bin 4 (right boundary included)
        !
        expected_counts = reshape([&
            1, 1, 1,&
            1, 1, 1,&
            2, 2, 2,&
            2, 2, 2&
        ], [n_neighbors, n_bins])
        expected_pmf = reshape([&
            0.16666666666666666_real64, 0.16666666666666666_real64, 0.16666666666666666_real64,&
            0.16666666666666666_real64, 0.16666666666666666_real64, 0.16666666666666666_real64,&
            0.3333333333333333_real64, 0.3333333333333333_real64, 0.3333333333333333_real64,&
            0.3333333333333333_real64, 0.3333333333333333_real64, 0.3333333333333333_real64&
        ], [n_neighbors, n_bins])

        call assert_equal_array_int(counts, expected_counts, size(counts, kind=int32), "test_build_residual_histograms: Test 3: counts don't match")
        call assert_equal_array_real(pmf, expected_pmf, size(pmf, kind=int32), TOL, "test_build_residual_histograms: Test 3: pmf don't match")
        call assert_equal_int(included(1), 6, "test_build_residual_histograms: Test 3: included row 1")
        call assert_equal_int(included(2), 6, "test_build_residual_histograms: Test 3: included row 2")
        call assert_equal_int(included(3), 6, "test_build_residual_histograms: Test 3: included row 3")

    end subroutine test_build_residual_histograms

    subroutine test_compute_divergence_per_reference_point
        integer(int32), parameter :: n_neighbors = 3, n_bins = 4
        real(real64), dimension(n_neighbors, n_bins) :: p, q
        real(real64), dimension(n_neighbors) :: jsd, expected_jsd
        integer(int32) :: ierr
        real(real64) :: tol

        ! ============================================================
        ! Test 1 — Identical PMFs → JSD = 0
        ! ============================================================
        p = reshape([0.1, 0.2, 0.3, 0.4, &
                     0.25,0.25,0.25,0.25, &
                     1.0, 0.0, 0.0, 0.0], [n_neighbors,n_bins])
        q = p

        call compute_divergence_per_reference_point(p, q, n_neighbors, n_bins, jsd, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_compute_divergence_per_reference_point: Test 1: ierr OK")

        expected_jsd = 0.0_real64
        call assert_equal_array_real(jsd, expected_jsd, size(jsd, kind=int32), TOL, "test_compute_divergence_per_reference_point: Test 1: identical PMFs → JSD=0")

        ! ============================================================
        ! Test 2 — Completely disjoint PMFs → JSD = log(2)
        ! ============================================================
        !
        ! p = [1,0,0,0]
        ! q = [0,1,0,0]
        !
        ! mix = [0.5,0.5,0,0]
        !
        ! KL(p||mix) = 1 * log(1/0.5) = log(2)
        ! KL(q||mix) = log(2)
        !
        ! JSD = 0.5*(log2 + log2) = log(2)
        !
        p = 0.0_real64
        q = 0.0_real64
        p(1,1) = 1.0
        q(1,2) = 1.0

        call compute_divergence_per_reference_point(p, q, n_neighbors, n_bins, jsd, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_compute_divergence_per_reference_point: Test 2: ierr OK")

        expected_jsd = 0.0_real64
        expected_jsd(1) = log(2.0_real64)
        call assert_equal_array_real(jsd, expected_jsd, size(jsd, kind=int32), TOL, "test_compute_divergence_per_reference_point: Test 2: disjoint PMFs → JSD=log(2)")
        call assert_equal_real(jsd(2), 0.0_real64, TOL, "test_compute_divergence_per_reference_point: Test 2: rows 2,3 are zero PMFs → JSD=0")
        call assert_equal_real(jsd(3), 0.0_real64, TOL, "test_compute_divergence_per_reference_point: Test 2: rows 2,3 are zero PMFs → JSD=0")

        ! ============================================================
        ! Test 3 — Partially overlapping PMFs (analytic check)
        ! ============================================================
        !
        ! p = [0.5, 0.5, 0, 0]
        ! q = [0.0, 1.0, 0, 0]
        !
        ! mix = [0.25, 0.75, 0, 0]
        !
        ! KL(p||mix) = 0.5*log(0.5/0.25) + 0.5*log(0.5/0.75)
        !            = 0.5*log(2) + 0.5*log(2/3)
        !
        ! KL(q||mix) = 1.0*log(1.0/0.75)
        !
        ! JSD = 0.5*(KL_p + KL_q)
        !
        p = 0.0_real64
        q = 0.0_real64
        p(1,:) = [0.5, 0.5, 0.0, 0.0]
        q(1,:) = [0.0, 1.0, 0.0, 0.0]

        call compute_divergence_per_reference_point(p, q, n_neighbors, n_bins, jsd, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_compute_divergence_per_reference_point: Test 3: ierr OK")

        ! Compute expected value analytically
        expected_jsd = 0.0_real64
        expected_jsd(1) = 0.5_real64 * ( &
            0.5_real64*log(2.0_real64) + 0.5_real64*log(2.0_real64/3.0_real64) &
          + log(1.0_real64/0.75_real64) )
        call assert_equal_array_real(jsd, expected_jsd, size(jsd, kind=int32), TOL, "test_compute_divergence_per_reference_point: Test 3: analytic partial-overlap JSD")

        ! ============================================================
        ! Test 4 — Zero-probability bins handled correctly
        ! ============================================================
        !
        ! p = [1,0,0,0]
        ! q = [0,0,1,0]
        !
        ! mix = [0.5,0,0.5,0]
        !
        ! KL(p||mix) = log(1/0.5) = log(2)
        ! KL(q||mix) = log(1/0.5) = log(2)
        !
        ! JSD = log(2)
        !
        p = 0.0_real64
        q = 0.0_real64
        p(1,1) = 1.0
        q(1,3) = 1.0

        call compute_divergence_per_reference_point(p, q, n_neighbors, n_bins, jsd, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_compute_divergence_per_reference_point: Test 4: ierr OK")

        expected_jsd = 0.0_real64
        expected_jsd(1) = log(2.0_real64)
        call assert_equal_array_real(jsd, expected_jsd, size(jsd, kind=int32), TOL, "test_compute_divergence_per_reference_point: Test 4: zero-probability bins handled correctly")

        ! ============================================================
        ! Test 5 — Multiple neighbors, mixed patterns
        ! ============================================================
        p = reshape([ &
            0.2,0.3,0.5,&
            0.0,1.0,0.0,&
            0.0,0.0,0.25,&
            0.25,0.25,0.25&
        ], [n_neighbors,n_bins])

        q = reshape([ &
            0.2,0.3,0.5,&
            0.0,0.0,1.0,&
            0.0,0.0,0.25,&
            0.25,0.25,0.25&
        ], [n_neighbors,n_bins])

        call compute_divergence_per_reference_point(p, q, n_neighbors, n_bins, jsd, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_compute_divergence_per_reference_point: Test 5: ierr OK")

        expected_jsd = 0.0_real64
        expected_jsd(2) = 0.5 * (1.0_real64 * log(1.0_real64 / 0.5_real64))
        expected_jsd(3) = 0.5 * (1.0_real64 * log(1.0_real64 / 0.5_real64))
        call assert_equal_array_real(jsd, expected_jsd, size(jsd, kind=int32), TOL, "test_compute_divergence_per_reference_point: Test 5: mixed patterns")
    end subroutine test_compute_divergence_per_reference_point

    subroutine test_compute_weighted_global_divergence
        integer(int32), parameter :: n_neighbors = 4
        real(real64), dimension(n_neighbors) :: jsd
        integer(int32), dimension(n_neighbors) :: n1, n2
        real(real64), dimension(n_neighbors) :: w, expected_weights
        real(real64) :: global_jsd
        integer(int32) :: ierr
        real(real64) :: expected, tol

        ! ============================================================
        ! Test 1 — Simple case: equal sample counts → uniform weights
        ! ============================================================
        !
        ! jsd = [0.1, 0.2, 0.3, 0.4]
        ! n1 = [5,5,5,5]
        ! n2 = [5,5,5,5]
        !
        ! n_j = [10,10,10,10]
        ! T = 40
        ! w = [0.25,0.25,0.25,0.25]
        !
        ! global_jsd = 0.25*(0.1+0.2+0.3+0.4) = 0.25
        !
        jsd = [0.1_real64, 0.2_real64, 0.3_real64, 0.4_real64]
        n1  = [5.0_real64,5.0_real64,5.0_real64,5.0_real64]
        n2  = [5.0_real64,5.0_real64,5.0_real64,5.0_real64]

        call compute_weighted_global_divergence(jsd, n_neighbors, n1, n2, &
                                                global_jsd, w, ierr)

        expected_weights = 0.25_real64
        call assert_equal_int(ierr, ERR_OK, "test_compute_weighted_global_divergence: Test 1: ierr OK")
        call assert_equal_array_real(w, expected_weights, size(w, kind=int32), TOL, "test_compute_weighted_global_divergence: Test 1: uniform weights")
        call assert_equal_real(global_jsd, 0.25_real64, TOL, "test_compute_weighted_global_divergence: Test 1: global JSD")

        ! ============================================================
        ! Test 2 — Unequal sample counts → weighted average
        ! ============================================================
        !
        ! jsd = [1.0, 2.0, 3.0, 4.0]
        ! n1 = [10, 20, 30, 40]
        ! n2 = [ 0, 10, 10, 10]
        !
        ! n_j = [10, 30, 40, 50]
        ! T = 130
        ! w = [10/130, 30/130, 40/130, 50/130]
        !
        ! global_jsd = sum_j w(j)*jsd(j)
        !
        jsd = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
        n1  = [10.0_real64,20.0_real64,30.0_real64,40.0_real64]
        n2  = [ 0.0_real64,10.0_real64,10.0_real64,10.0_real64]

        call compute_weighted_global_divergence(jsd, n_neighbors, n1, n2, &
                                                global_jsd, w, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_compute_weighted_global_divergence: Test 2: ierr OK")
        
        expected_weights = [10.0_real64, 30.0_real64, 40.0_real64, 50.0_real64] / 130.0_real64
        call assert_equal_array_real(w, expected_weights, size(w, kind=int32), TOL, "test_compute_weighted_global_divergence: Test 2: wweights")

        expected = (1.0_real64*(10.0_real64/130.0_real64) + &
                    2.0_real64*(30.0_real64/130.0_real64) + &
                    3.0_real64*(40.0_real64/130.0_real64) + &
                    4.0_real64*(50.0_real64/130.0_real64))

        call assert_equal_real(global_jsd, expected, TOL, "test_compute_weighted_global_divergence: Test 2: weighted global JSD")

        ! ============================================================
        ! Test 3 — Some neighborhoods have zero samples → weight = 0
        ! ============================================================
        !
        ! jsd = [0.5, 1.0, 2.0, 4.0]
        ! n1 = [0, 10, 0, 5]
        ! n2 = [0,  0, 0, 5]
        !
        ! n_j = [0,10,0,10]
        ! T = 20
        ! w = [0, 10/20, 0, 10/20] = [0,0.5,0,0.5]
        !
        ! global_jsd = 0.5*1.0 + 0.5*4.0 = 2.5
        !
        jsd = [0.5_real64, 1.0_real64, 2.0_real64, 4.0_real64]
        n1  = [0.0_real64,10.0_real64,0.0_real64,5.0_real64]
        n2  = [0.0_real64, 0.0_real64,0.0_real64,5.0_real64]

        call compute_weighted_global_divergence(jsd, n_neighbors, n1, n2, &
                                                global_jsd, w, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_compute_weighted_global_divergence: Test 3: ierr OK")

        expected_weights = [0.0_real64,0.5_real64,0.0_real64,0.5_real64]
        call assert_equal_array_real(w, expected_weights, size(w, kind=int32), TOL, "test_compute_weighted_global_divergence: Test 3: weights with zero-sample neighborhoods")
        call assert_equal_real(global_jsd, 2.5_real64, TOL, "test_compute_weighted_global_divergence: Test 3: weighted global JSD")

        ! ============================================================
        ! Test 4 — All neighborhoods have zero samples → weights=0, global JSD=0
        ! ============================================================
        !
        ! jsd = [1,2,3,4]
        ! n1 = [0,0,0,0]
        ! n2 = [0,0,0,0]
        !
        ! n_j = [0,0,0,0]
        ! T = 0 → weights = 0
        ! global_jsd = 0
        !
        jsd = [1.0_real64,2.0_real64,3.0_real64,4.0_real64]
        n1  = 0
        n2  = 0

        call compute_weighted_global_divergence(jsd, n_neighbors, n1, n2, &
                                                global_jsd, w, ierr)

        expected_weights = 0.0_real64
        call assert_equal_int(ierr, ERR_OK, "test_compute_weighted_global_divergence: Test 4: ierr OK")
        call assert_no_nan_real(w, size(w, kind=int32), "test_compute_weighted_global_divergence: Test 4: all weights not NaN")
        call assert_equal_array_real(w, expected_weights, size(w, kind=int32), TOL, "test_compute_weighted_global_divergence: Test 4: all weights zero")
        call assert_equal_real(global_jsd, 0.0_real64, TOL, "test_compute_weighted_global_divergence: Test 4: global JSD zero when no samples")

        ! ============================================================
        ! Test 5 — Mixed jsd values, mixed sample counts
        ! ============================================================
        !
        ! jsd = [0.0, 0.5, 1.0, 2.0]
        ! n1 = [5, 0, 10, 5]
        ! n2 = [5, 5,  0, 5]
        !
        ! n_j = [10,5,10,10]
        ! T = 35
        ! w = [10/35, 5/35, 10/35, 10/35]
        !
        ! global_jsd = sum_j w(j)*jsd(j)
        !
        jsd = [0.0_real64, 0.5_real64, 1.0_real64, 2.0_real64]
        n1  = [5.0_real64,0.0_real64,10.0_real64,5.0_real64]
        n2  = [5.0_real64,5.0_real64, 0.0_real64,5.0_real64]

        call compute_weighted_global_divergence(jsd, n_neighbors, n1, n2, &
                                                global_jsd, w, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_compute_weighted_global_divergence: Test 5: ierr OK")

        expected_weights = [10.0_real64, 5.0_real64, 10.0_real64, 10.0_real64] / 35.0_real64

        call assert_equal_array_real(w, expected_weights, size(w, kind=int32), TOL, "test_compute_weighted_global_divergence: Test 5: all weights zero")

        expected = (0.0_real64*(10.0_real64/35.0_real64) + &
                    0.5_real64*( 5.0_real64/35.0_real64) + &
                    1.0_real64*(10.0_real64/35.0_real64) + &
                    2.0_real64*(10.0_real64/35.0_real64))

        call assert_equal_real(global_jsd, expected, TOL, "test_compute_weighted_global_divergence: Test 5: weighted global JSD")

    end subroutine test_compute_weighted_global_divergence

end module mod_test_jensen_shannon_divergence
