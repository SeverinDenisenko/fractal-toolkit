program color_test
   use precision, only: wp
   use hurst, only: estimate_hurst_berg, estimate_hurst_yw, estimate_hurst_rs, estimate_hurst_lssd
   use generators, only: generate_color, generate_fgn
   use stdlib_error, only: check
   implicit none

   real(wp), allocatable :: series(:)
   integer :: n, m
   real(wp) :: H, a, H_err, a_err, sigma2

   n = 2**13
   m = 100

   allocate(series(n))

   call generate_color(series, 0.0_wp)
   call estimate_hurst_berg(series, m, H, a, H_err, a_err, sigma2)
   call check(abs(H - 0.5_wp) < 1e-1)
   call check(abs(a - 0.0_wp) < 1e-1)
   call check(abs(H_err - 0.0_wp) < 1e-2)
   call check(abs(a_err - 0.0_wp) < 1e-2)
   call check(abs(sigma2 - 0.0_wp) < 1e-1)

   call generate_color(series, 1.0_wp)
   call estimate_hurst_yw(series, m, H, a, H_err, a_err, sigma2)
   call check(abs(H - 1.0_wp) < 1e-1)
   call check(abs(a - 1.0_wp) < 1e-1)
   call check(abs(H_err - 0.0_wp) < 1e-2)
   call check(abs(a_err - 0.0_wp) < 1e-2)
   call check(abs(sigma2 - 0.0_wp) < 1e-1)

   call generate_color(series, 2.0_wp)
   call estimate_hurst_berg(series, m, H, a, H_err, a_err, sigma2)
   call check(abs(H - 0.5_wp) < 1e-1)
   call check(abs(a - 2.0_wp) < 1e-1)
   call check(abs(H_err - 0.0_wp) < 1e-2)
   call check(abs(a_err - 0.0_wp) < 1e-2)
   call check(abs(sigma2 - 0.0_wp) < 1e-1)

   call generate_color(series, 3.0_wp)
   call estimate_hurst_yw(series, m, H, a, H_err, a_err, sigma2)
   call check(abs(H - 1.0_wp) < 1e-1)
   call check(abs(a - 3.0_wp) < 1e-1)
   call check(abs(H_err - 0.0_wp) < 1e-2)
   call check(abs(a_err - 0.0_wp) < 1e-2)
   call check(abs(sigma2 - 0.0_wp) < 1e-1)

   call generate_color(series, 0.0_wp)
   call estimate_hurst_rs(series, H, H_err, sigma2)
   call check(abs(H - 0.5_wp) < 1e-1)
   call check(abs(H_err - 0.0_wp) < 1e-1)
   call check(abs(sigma2 - 0.0_wp) < 1e-1)

   call generate_color(series, 1.0_wp)
   call estimate_hurst_rs(series, H, H_err, sigma2)
   call check(abs(H - 1.0_wp) < 1e-1)
   call check(abs(H_err - 0.0_wp) < 1e-1)
   call check(abs(sigma2 - 0.0_wp) < 1e-1)

   call generate_color(series, 0.0_wp)
   call estimate_hurst_lssd(series, 1, 100, H, H_err)
   call check(abs(H - 0.5_wp) < 1e-1)
   call check(abs(H_err - 0.0_wp) < 1e-1)

   call generate_color(series, 1.0_wp)
   call estimate_hurst_lssd(series, 3, 1000, H, H_err)
   call check(abs(H - 1.0_wp) < 1e-2)
   call check(abs(H_err - 0.0_wp) < 1e-2)

   call generate_fgn(series, 0.5_wp)
   call estimate_hurst_lssd(series, 3, 100, H, H_err)
   call check(abs(H - 0.5_wp) < 1e-2)
   call check(abs(H_err - 0.0_wp) < 1e-2)

   call generate_fgn(series, 0.001_wp)
   call estimate_hurst_lssd(series, 3, 50, H, H_err)
   call check(abs(H - 0.001_wp) < 1e-1)
   call check(abs(H_err - 0.0_wp) < 1e-2)

   deallocate(series)
end program color_test
