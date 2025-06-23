program test_normalize_by_std_dev
  use, intrinsic :: iso_fortran_env, only: real64
  use testdrive
  use tox_normalization_mod
  implicit none

  call run_tests([ &
    new_unittest("handles large numbers", test_normalizes_large_numbers) &
  ])

contains

  pure function itoa(i) result(str)
    integer, intent(in) :: i
    character(len=12) :: str
    write(str,'(I0)') i
  end function itoa

  subroutine test_normalizes_large_numbers(this)
    use testdrive, only: assert_true

    type(unittest_type) :: this
    real(real64) :: mat(2,2), result(2,2), expected(2,2), std_dev(2)
    integer :: i, j
    character(len=32) :: msg

    mat = reshape([1e6_real64, 2e6_real64, 1e6_real64, 2e6_real64], [2,2])
    call normalize_by_std_dev_core(2, 2, mat, result)
    do i = 1, 2
      std_dev(i) = sqrt((mat(i,1)**2 + mat(i,2)**2)/2.0_real64)
      do j = 1, 2
        expected(i,j) = mat(i,j) / std_dev(i)
        write(msg, '(A,I0,A,I0,A)') "Large number at (", i, ",", j, ")"
        call assert_true(abs(result(i,j) - expected(i,j)) < 1e-12, msg)

      end do
    end do
  end subroutine

end program test_normalize_by_std_dev
