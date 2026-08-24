program stat_test
   use precision, only: wp
   use stat, only: fd_bins, percentile, fit_normal
   use stdlib_error, only: check
   use generators, only: generate_gauss
   implicit none

   real(wp) :: pct(2)
   real(wp) :: data(1000)
   real(wp) :: mu, sigma, se_mu, se_sigma, ks, cvm, sk, ku

   call check(fd_bins([1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp, 7.0_wp, 8.0_wp, 9.0_wp, 10.0_wp]) == 3)
   call check(fd_bins([1.0_wp, 2.0_wp]) == 2)

   call percentile([1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp, 7.0_wp, 8.0_wp, 9.0_wp], [75.0_wp, 25.0_wp], pct)
   call check(abs(7.0_wp - pct(1)) < 1e-5)
   call check(abs(3.0_wp - pct(2)) < 1e-5)

   call generate_gauss(data, 1.0_wp, 2.0_wp)
   call fit_normal(data, mu, sigma, se_mu, se_sigma, ks, cvm, sk, ku)
   call check(abs(mu - 1.0_wp) < 1e-1)
   call check(abs(sigma - 2.0_wp) < 1e-1)
   call check(abs(se_mu - 0.0_wp) < 1e-1)
   call check(abs(se_sigma - 0.0_wp) < 1e-1)
   call check(abs(ks - 0.0_wp) < 1e-1)
   call check(abs(cvm - 0.0_wp) < 1e-1)
   call check(abs(sk - 0.0_wp) < 1e-1)
   call check(abs(ku - 0.0_wp) < 1e-0)

end program stat_test
