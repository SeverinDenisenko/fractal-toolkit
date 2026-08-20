program autoreg_test
   use precision, only: wp
   use autoreg, only: ar_coeff
   use stdlib_error, only: check
   implicit none

   real(wp) :: coeff(3)
   real(wp) :: data(7)

   data = [0, 1, 2, 3, 4, 5, 6]
   call ar_coeff(coeff, data)
   call check(maxval(abs(coeff - [1.33333333_wp, 0.33333333_wp, -0.66666667_wp])) < 1e-5_wp)

   data = [1, -1, 1, -1, 1, -1, 1]
   call ar_coeff(coeff, data)
   call check(maxval(abs(coeff - [-0.33333333_wp,  0.33333333_wp, -0.33333333_wp])) < 1e-5_wp)
end program autoreg_test
