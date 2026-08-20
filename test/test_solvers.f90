program solvers_test
   use precision, only: wp
   use stdlib_error, only: check
   use solvers, only: linregress, powerregress
   implicit none

   real(wp) :: x(3)
   real(wp) :: y(3)
   real(wp) :: x1(10)
   real(wp) :: y1(10)
   real(wp) :: a, k, c
   real(wp) :: sigma2, a_err, k_err, c_err

   x1(:) = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
   y1(:) = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
   call linregress(x1, y1, k, a, sigma2, k_err, a_err)
   call check(abs(k - 1.0_wp) < 1e-5_wp)
   call check(abs(a - 0.0_wp) < 1e-5_wp)
   call check(abs(k_err - 0.0_wp) < 1e-5_wp)
   call check(abs(a_err - 0.0_wp) < 1e-5_wp)

   x(:) = [0, 1, 2]
   y(:) = [1, 2, 3]
   call linregress(x, y, k, a)
   call check(abs(k - 1.0_wp) < 1e-5_wp)
   call check(abs(a - 1.0_wp) < 1e-5_wp)

   x(:) = [0, 1, 2]
   y(:) = [1, 3, 5]
   call linregress(x, y, k, a)
   call check(abs(k - 2.0_wp) < 1e-5_wp)
   call check(abs(a - 1.0_wp) < 1e-5_wp)

   x(:) = [1, 2, 3]
   y(:) = [1, 4, 9]
   call powerregress(x, y, a, c, sigma2, a_err, c_err)
   call check(abs(sigma2 - 0.0_wp) < 1e-5_wp)
   call check(abs(c - 1.0_wp) < 1e-5_wp)
   call check(abs(a - 2.0_wp) < 1e-5_wp)
   call check(abs(a_err - 0.0_wp) < 1e-5_wp)
   call check(abs(c_err - 0.0_wp) < 1e-5_wp)
end program solvers_test
