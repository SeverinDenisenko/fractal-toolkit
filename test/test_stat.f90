program stat_test
   use precision, only: wp
   use stdlib_error, only: check
   use stat, only: friedman_diaconis_bins, percentile
   implicit none

   real(wp) :: pct(2)

   call check(friedman_diaconis_bins([1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp, 7.0_wp, 8.0_wp, 9.0_wp, 10.0_wp]) == 3)
   call check(friedman_diaconis_bins([1.0_wp, 2.0_wp]) == 2)

   call percentile([1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp, 7.0_wp, 8.0_wp, 9.0_wp], [75.0_wp, 25.0_wp], pct)
   call check(abs(7.0_wp - pct(1)) < 1e-5)
   call check(abs(3.0_wp - pct(2)) < 1e-5)

end program stat_test
