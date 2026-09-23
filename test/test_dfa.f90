program dfa_test
   use precision, only: wp
   use hurst, only: dfa_chart, dfa_chart_size, estimate_hurst_dfa
   use generators, only: generate_gauss, generate_fgn
   use stdlib_error, only: check
   implicit none

   real(wp), allocatable :: series(:), F(:), N(:)
   real(wp) :: H, H_err, sigma2, Hmean
   integer :: nlen, n_sizes

   call check(dfa_chart_size(8) == 2)
   call check(dfa_chart_size(16) == 3)
   call check(dfa_chart_size(1024) == 9)

   nlen = 16
   n_sizes = dfa_chart_size(nlen)
   allocate(series(nlen), F(n_sizes), N(n_sizes))

   call generate_gauss(series, 0.0_wp, 1.0_wp, 42)
   call dfa_chart(series, F, N)
   call check(size(N) == n_sizes)
   call check(maxval(abs(N - [4.0_wp, 8.0_wp, 16.0_wp])) < 1e-5_wp)
   call check(all(F > 0.0_wp))
   call check(all(F < huge(1.0_wp)))

   deallocate(series, F, N)

   nlen = 2048
   allocate(series(nlen))

   call generate_gauss(series, 0.0_wp, 1.0_wp, 42)
   call estimate_hurst_dfa(series, H, H_err, sigma2)
   call check(abs(H - 0.5_wp) < 0.1_wp)

   call check_hurst_average(series, 0.3_wp, 0.05_wp, Hmean)
   call check_hurst_average(series, 0.5_wp, 0.05_wp, Hmean)
   call check_hurst_average(series, 0.7_wp, 0.05_wp, Hmean)

   deallocate(series)
contains
   subroutine check_hurst_average(series, Htrue, tol, Hmean)
      real(wp), intent(inout) :: series(:)
      real(wp), intent(in) :: Htrue, tol
      real(wp), intent(out) :: Hmean

      real(wp) :: H, H_err, sigma2
      integer :: i

      Hmean = 0.0_wp
      do i = 1, 5
         call generate_fgn(series, Htrue, seed_in=i)
         call estimate_hurst_dfa(series, H, H_err, sigma2)
         Hmean = Hmean + H
      end do
      Hmean = Hmean / 5.0_wp
      call check(abs(Hmean - Htrue) < tol)
   end subroutine check_hurst_average
end program dfa_test