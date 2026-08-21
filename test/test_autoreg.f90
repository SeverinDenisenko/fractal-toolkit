program autoreg_test
   use precision, only: wp
   use autoreg, only: yw_ar_coeff, burg_ar_coeff, ar_predict
   use stdlib_error, only: check
   implicit none

   real(wp) :: coeff(3)
   real(wp) :: data(16)

   data(:) = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
   call yw_ar_coeff(coeff, data)
   call check(abs(ar_predict(coeff, data, size(data) + 1) - 16.0_wp) < 1e-5_wp)
   call check(maxval(abs(coeff - [1.33333333_wp, 0.33333333_wp, -0.66666667_wp])) < 1e-5_wp)

   data(:) = [1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1]
   call yw_ar_coeff(coeff, data)
   call check(abs(ar_predict(coeff, data, size(data) + 1) - 1.0_wp) < 1e-5_wp)
   call check(maxval(abs(coeff - [-0.33333333_wp,  0.33333333_wp, -0.33333333_wp])) < 1e-5_wp)

   data(:) = [1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1]
   call burg_ar_coeff(coeff, data)
   call check(abs(ar_predict(coeff, data, size(data) + 1) - 1.0_wp) < 1e-5_wp)
   call check(maxval(abs(coeff - [-1.0_wp,  0.0_wp, 0.0_wp])) < 1e-5_wp)
end program autoreg_test
