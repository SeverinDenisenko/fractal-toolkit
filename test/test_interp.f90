program interp_test
   use precision, only: wp
   use interp, only: cubic_resample
   use stdlib_error, only: check
   implicit none

   real(wp) :: x(4)
   real(wp) :: y(4)
   real(wp) :: x1(4)
   real(wp) :: y1(4)

   x = [0, 1, 2, 3]
   y = [0, 1, 2, 3]
   x1 = [0, 1, 2, 3]
   y1 = [0, 0, 0, 0]
   call cubic_resample(x, y, x1, y1)
   call check(sum(abs(y1 - [0, 1, 2, 3])) < 1e-5)

   x = [0, 1, 2, 3]
   y = [0, 1, 2, 3]
   x1 = [0.0, 0.25, 0.5, 0.75]
   y1 = [0, 0, 0, 0]
   call cubic_resample(x, y, x1, y1)
   call check(sum(abs(y1 - [0.0, 0.25, 0.5, 0.75])) < 1e-5)

   x = [0, 1, 2, 3]
   y = [0, 1, 8, 27]
   x1 = [0.25, 0.5, 0.75, 1.0]
   y1 = [0, 0, 0, 0]
   call cubic_resample(x, y, x1, y1)
   call check(sum(abs(y1 - [0.015625, 0.125, 0.421875, 1.0])) < 1e-5)

   x = [10, 11, 20, 100]
   y = [100, 121, 400, 10000]
   x1 = [14, 16, 18, 20]
   y1 = [0, 0, 0, 0]
   call cubic_resample(x, y, x1, y1)
   call check(sum(abs(y1 - [196, 256, 324, 400])) < 1e-5)
end program interp_test
