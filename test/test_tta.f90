program tta_test
   use precision, only: wp
   use hurst, only: tta_chart, tta_chart_size, estimate_hurst_tta
   use generators, only: generate_gauss, generate_fgn
   use stdlib_error, only: check
   implicit none

   real(wp), allocatable :: series(:), S(:), T(:)
   real(wp) :: H, H_err, sigma2, Hmean
   integer :: nlen, n_sizes
   integer :: i

   call check(tta_chart_size(21) == 10)
   call check(tta_chart_size(20) == 9)
   call check(tta_chart_size(1024) == 10)

   nlen = 2048
   n_sizes = tta_chart_size(nlen)
   allocate(series(nlen), S(n_sizes), T(n_sizes))

   call generate_gauss(series, 0.0_wp, 1.0_wp, 42)
   call tta_chart(series, S, T)
   call check(size(T) == n_sizes)
   call check(maxval(abs(T - [(real(i, wp), i = 1, n_sizes)])) < 1e-5_wp)
   call check(all(S > 0.0_wp))
   call check(all(S < huge(1.0_wp)))

   call estimate_hurst_tta(series, H, H_err, sigma2)
   call check(abs(H - 0.5_wp) < 0.1_wp)

   call check_hurst_average(series, 0.3_wp, 0.05_wp, Hmean)
   call check_hurst_average(series, 0.5_wp, 0.05_wp, Hmean)
   call check_hurst_average(series, 0.7_wp, 0.05_wp, Hmean)

   deallocate(series, S, T)
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
         call estimate_hurst_tta(series, H, H_err, sigma2)
         Hmean = Hmean + H
      end do
      Hmean = Hmean / 5.0_wp
      call check(abs(Hmean - Htrue) < tol)
   end subroutine check_hurst_average
end program tta_test