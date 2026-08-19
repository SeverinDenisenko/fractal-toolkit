module hurst
   use precision, only: wp
   use stdlib_error, only: check
   use solvers, only: powerregress
   use spectra, only: berg_psd, berg_psd_size
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

   subroutine estimate_hurst_berg(series, m, H, a, H_err, a_err, sigma2)
      real(wp), intent(in) :: series(:)
      integer, intent(in) :: m
      real(wp), intent(out) :: H, H_err
      real(wp), intent(out) :: a, a_err
      real(wp), intent(out) :: sigma2

      real(wp), allocatable :: P(:)
      real(wp), allocatable :: f(:)
      real(wp) :: c, c_err
      integer :: n

      n = size(series)

      allocate(f(berg_psd_size(n)), P(berg_psd_size(n)))

      call berg_psd(f, P, series, 1.0_wp, m)
      call powerregress(f(2:), P(2:), a, c, sigma2, a_err, c_err)

      a = -a
      H = slope_to_hurst(a)
      H_err = a_err / 2.0_wp
   end subroutine estimate_hurst_berg
end module hurst
