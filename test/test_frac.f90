program frac_test
   use precision, only: wp
   use stdlib_error, only: check
   use frac, only: frac_diff_coeff, frac_diff
   implicit none

   real(wp) :: coeff(3)
   real(wp) :: data(5)
   real(wp) :: diff(5)

   data = [0.0_wp, 1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp]

   call frac_diff_coeff(coeff, 0.0_wp)
   call frac_diff(coeff, data, diff)
   call check(maxval(abs(diff - [0.0_wp, 1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp])) < 1e-5_wp)

   call frac_diff_coeff(coeff, 1.0_wp)
   call frac_diff(coeff, data, diff)
   call check(maxval(abs(diff - [0.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp])) < 1e-5_wp)

   call frac_diff_coeff(coeff, -1.0_wp)
   call frac_diff(coeff, data, diff)
   call check(maxval(abs(diff - [0.0_wp, 1.0_wp, 3.0_wp, 6.0_wp, 9.0_wp])) < 1e-5_wp)
end program frac_test
