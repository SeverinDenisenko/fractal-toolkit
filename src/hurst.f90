module hurst
   use precision, only: wp
   use solvers, only: powerregress
   use spectra, only: berg_psd, yw_psd, psd_size
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

   subroutine estimate_hurst_psd(f, P, H, a, H_err, a_err, sigma2, ierr)
      real(wp), intent(in) :: f(:), P(:)
      real(wp), intent(out) :: H, a, H_err, a_err, sigma2
      integer, intent(out), optional :: ierr

      real(wp) :: c, c_err

      H = 0.0_wp
      a = 0.0_wp
      H_err = 0.0_wp
      a_err = 0.0_wp
      sigma2 = 0.0_wp

      call powerregress(f(2:), P(2:), a, c, sigma2, a_err, c_err, ierr=ierr)
      if (present(ierr) .and. ierr /= 0) return

      a = -a
      H = slope_to_hurst(a)
      H_err = a_err / 2.0_wp
   end subroutine estimate_hurst_psd

   subroutine estimate_hurst_berg(series, m, H, a, H_err, a_err, sigma2, ierr)
      real(wp), intent(in) :: series(:)
      integer, intent(in) :: m
      real(wp), intent(out) :: H, a, H_err, a_err, sigma2
      integer, intent(out), optional :: ierr

      real(wp), allocatable :: P(:), f(:)
      integer :: n

      n = size(series)

      allocate(f(psd_size(n)), P(psd_size(n)))

      call berg_psd(f, P, series, 1.0_wp, m, ierr=ierr)
      if (present(ierr) .and. ierr /= 0) return

      call estimate_hurst_psd(f, P, H, a, H_err, a_err, sigma2, ierr=ierr)

      deallocate(f, P)
   end subroutine estimate_hurst_berg

   subroutine estimate_hurst_yw(series, m, H, a, H_err, a_err, sigma2, ierr)
      real(wp), intent(in) :: series(:)
      integer, intent(in) :: m
      real(wp), intent(out) :: H, a, H_err, a_err, sigma2
      integer, intent(out), optional :: ierr

      real(wp), allocatable :: P(:), f(:)
      integer :: n

      n = size(series)

      allocate(f(psd_size(n)), P(psd_size(n)))

      call yw_psd(f, P, series, 1.0_wp, m, ierr=ierr)
      if (present(ierr) .and. ierr /= 0) return

      call estimate_hurst_psd(f, P, H, a, H_err, a_err, sigma2, ierr=ierr)

      deallocate(f, P)
   end subroutine estimate_hurst_yw

   ! TODO
   ! subroutine rs_analysis(series, H, a, H_err, a_err, sigma2)
   !    real(wp), intent(in) :: series(:)
   !    real(wp), intent(out) :: H, a, H_err, a_err, sigma2
   ! end subroutine rs_analysis
end module hurst
