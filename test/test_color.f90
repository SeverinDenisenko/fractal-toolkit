program color_test
   use precision, only: wp
   use hurst, only: estimate_hurst_berg
   use generators, only: generate_color
   use stdlib_error, only: check
   implicit none

   real(wp), allocatable :: series(:)
   integer :: n, m
   real(wp) :: H, a, H_err, a_err, sigma2

   n = 10000
   m = 100

   allocate(series(n))

   call generate_color(series, 3.0_wp)
   call estimate_hurst_berg(series, m, H, a, H_err, a_err, sigma2)
   call check(abs(H - 1.0_wp) < 1e-2)
   call check(abs(a - 3.0_wp) < 1e-2)
   call check(abs(H_err - 0.0_wp) < 1e-2)
   call check(abs(a_err - 0.0_wp) < 1e-2)
   call check(abs(sigma2 - 0.0_wp) < 1e-1)

   call generate_color(series, 2.0_wp)
   call estimate_hurst_berg(series, m, H, a, H_err, a_err, sigma2)
   call check(abs(H - 0.5_wp) < 1e-2)
   call check(abs(a - 2.0_wp) < 1e-2)
   call check(abs(H_err - 0.0_wp) < 1e-2)
   call check(abs(a_err - 0.0_wp) < 1e-2)
   call check(abs(sigma2 - 0.0_wp) < 1e-1)

   call generate_color(series, 0.0_wp)
   call estimate_hurst_berg(series, m, H, a, H_err, a_err, sigma2)
   call check(abs(H - 0.5_wp) < 1e-2)
   call check(abs(a - 0.0_wp) < 1e-2)
   call check(abs(H_err - 0.0_wp) < 1e-2)
   call check(abs(a_err - 0.0_wp) < 1e-2)
   call check(abs(sigma2 - 0.0_wp) < 1e-1)
end program color_test
