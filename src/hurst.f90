module hurst
   use precision, only: wp
   use solvers, only: powerregress
   use spectra, only: berg_psd, yw_psd, psd_size
   use stdlib_error, only: check
   implicit none

contains
   ! Compute Hurst expoenent from PSD slope `a` where PSD~f^a
   real(wp) function slope_to_hurst(a) result(H)
      real(wp), intent(in) :: a

      if (a < 1.0_wp) then
         H = (a + 1.0_wp) / 2.0_wp
      else
         H = (a - 1.0_wp) / 2.0_wp
      end if
   end function slope_to_hurst

   subroutine estimate_hurst_psd(f, P, H, a, H_err, a_err, sigma2)
      real(wp), intent(in) :: f(:), P(:)
      real(wp), intent(out) :: H, a, H_err, a_err, sigma2

      real(wp) :: c, c_err

      call powerregress(f(2:), P(2:), a, c, sigma2, a_err, c_err)

      a = -a
      H = slope_to_hurst(a)
      H_err = a_err / 2.0_wp
   end subroutine estimate_hurst_psd

   subroutine estimate_hurst_berg(series, m, H, a, H_err, a_err, sigma2)
      real(wp), intent(in) :: series(:)
      integer, intent(in) :: m
      real(wp), intent(out) :: H, a, H_err, a_err, sigma2

      real(wp), allocatable :: P(:), f(:)
      integer :: n

      n = size(series)

      allocate(f(psd_size(n)), P(psd_size(n)))

      call berg_psd(f, P, series, 1.0_wp, m)
      call estimate_hurst_psd(f, P, H, a, H_err, a_err, sigma2)

      deallocate(f, P)
   end subroutine estimate_hurst_berg

   subroutine estimate_hurst_yw(series, m, H, a, H_err, a_err, sigma2)
      real(wp), intent(in) :: series(:)
      integer, intent(in) :: m
      real(wp), intent(out) :: H, a, H_err, a_err, sigma2

      real(wp), allocatable :: P(:), f(:)
      integer :: n

      n = size(series)

      allocate(f(psd_size(n)), P(psd_size(n)))

      call yw_psd(f, P, series, 1.0_wp, m)
      call estimate_hurst_psd(f, P, H, a, H_err, a_err, sigma2)

      deallocate(f, P)
   end subroutine estimate_hurst_yw

   ! TODO
   ! subroutine rs_analysis(series, H, a, H_err, a_err, sigma2)
   !    real(wp), intent(in) :: series(:)
   !    real(wp), intent(out) :: H, a, H_err, a_err, sigma2
   ! end subroutine rs_analysis
end module hurst
