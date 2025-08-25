
!> @brief Main program to run the entire test suite.
PROGRAM run_tests
    USE, INTRINSIC :: iso_fortran_env, ONLY: INT32
    ! Import all test suites
    USE mod_test_gene_centroids 

    IMPLICIT NONE

    ABSTRACT INTERFACE
        SUBROUTINE run_all_tests_interface()
        END SUBROUTINE
        SUBROUTINE run_named_tests_interface(test_names)
            CHARACTER(LEN=*), INTENT(IN) :: test_names(:)
        END SUBROUTINE
    END INTERFACE

    TYPE :: test_suite
        CHARACTER(LEN=64) :: name
        PROCEDURE(run_all_tests_interface), POINTER, NOPASS :: run_all => NULL()
        PROCEDURE(run_named_tests_interface), POINTER, NOPASS :: run_named => NULL()
    END TYPE test_suite

    TYPE(test_suite), ALLOCATABLE :: available_suites(:)
    INTEGER :: i, n_args
    CHARACTER(LEN=256) :: arg
    LOGICAL :: suite_found
    CHARACTER(LEN=64), ALLOCATABLE :: test_names_to_run(:)

    CALL initialize_suites()
    n_args = command_argument_count()

    IF (n_args == 0) THEN
        WRITE(*, '(A)') "Running all tests from all suites..."
        WRITE(*,*)
        DO i = 1, SIZE(available_suites)
            CALL available_suites(i)%run_all()
        END DO
        WRITE(*, '(A)') "All tests completed."
        STOP 0
    END IF

    IF (n_args >= 1) THEN
        CALL get_command_argument(1, arg)
        suite_found = .FALSE.
        DO i = 1, SIZE(available_suites)
            IF (TRIM(arg) == TRIM(available_suites(i)%name)) THEN
                suite_found = .TRUE.
                IF (n_args == 1) THEN
                    CALL available_suites(i)%run_all()
                ELSE IF (n_args == 2) THEN
                    CALL get_command_argument(2, arg)
                    CALL parse_test_names(arg, test_names_to_run)
                    CALL available_suites(i)%run_named(test_names_to_run)
                END IF
                EXIT
            END IF
        END DO
        IF (.NOT. suite_found) THEN
            WRITE(*,*) "ERROR: Test suite '", TRIM(arg), "' not found."
            STOP 1
        END IF
    END IF

CONTAINS

    SUBROUTINE add_suite(name, run_all_proc, run_named_proc)
        CHARACTER(LEN=*), INTENT(IN) :: name
        PROCEDURE(run_all_tests_interface) :: run_all_proc
        PROCEDURE(run_named_tests_interface) :: run_named_proc
        TYPE(test_suite), ALLOCATABLE :: new_suites(:)
        INTEGER :: old_size
        old_size = 0
        IF (ALLOCATED(available_suites)) old_size = SIZE(available_suites)
        ALLOCATE(new_suites(old_size + 1))
        IF (old_size > 0) new_suites(1:old_size) = available_suites
        new_suites(old_size + 1)%name = name
        new_suites(old_size + 1)%run_all => run_all_proc
        new_suites(old_size + 1)%run_named => run_named_proc
        CALL move_alloc(new_suites, available_suites)
    END SUBROUTINE add_suite

    SUBROUTINE initialize_suites()
        ! The USE statement at the top of the program makes the procedures available,
        ! so the EXTERNAL statement is not needed and has been removed.
        CALL add_suite("gene_centroids", run_all_tests_gene_centroids, run_named_tests_gene_centroids) 
    END SUBROUTINE initialize_suites

    SUBROUTINE parse_test_names(names_str, names_arr)
        CHARACTER(LEN=*), INTENT(IN) :: names_str
        CHARACTER(LEN=64), ALLOCATABLE, INTENT(OUT) :: names_arr(:)
        CHARACTER(LEN=LEN(names_str)) :: temp_str
        INTEGER :: count, i, start
        temp_str = TRIM(names_str)
        count = 1
        DO i = 1, LEN_TRIM(temp_str); IF (temp_str(i:i) == ',') count = count + 1; END DO
        ALLOCATE(names_arr(count))
        start = 1; count = 1
        DO i = 1, LEN_TRIM(temp_str)
            IF (temp_str(i:i) == ',') THEN
                names_arr(count) = temp_str(start:i-1)
                start = i + 1; count = count + 1
            END IF
        END DO
        names_arr(count) = temp_str(start:)
    END SUBROUTINE parse_test_names

END PROGRAM run_tests
