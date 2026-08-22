module hurst
   use precision, only: wp
   use solvers, only: powerregress
   use spectra, only: berg_psd, yw_psd, psd_size
   use stat, only: mean, variance
   use checks, only: check
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

   integer function rs_chart_size(n) result(m)
      integer :: n
      m = n / 2 - 1
   end function rs_chart_size

   subroutine rs_chart(series, RS, N, ierr)
      real(wp), intent(in) :: series(:)
      real(wp), intent(out) :: RS(:), N(:)
      integer, intent(out), optional :: ierr

      integer :: i, j, a, n_curr, n_min, n_max, begin, end
      real(wp), allocatable :: cumdiv(:), stddiv(:)
      real(wp) :: R

      n_min = 2
      n_max = size(series) / 2
      N = [(i, i = n_min, n_max)]

      if(check(rs_chart_size(size(series)) == size(N), msg="rs_chart: size mismatch", ierr=ierr)) return
      if(check(rs_chart_size(size(series)) == size(RS), msg="rs_chart: size mismatch", ierr=ierr)) return

      do j = 1,size(N)
         n_curr = n_min + j - 1
         a = size(series) / n_curr
         allocate(cumdiv(a), stddiv(a))

         do i = 1, a
            begin = 1 + (i - 1) * n_curr
            end = begin + n_curr - 1

            if (end .gt. size(series)) exit

            cumdiv(i) = sum(series(begin:end) - mean(series(begin:end)))
            stddiv(i) = sqrt(variance(series(begin:end)))
         end do

         R = maxval(cumdiv) - minval(cumdiv)
         RS(j) = sum(R / stddiv) / a

         deallocate(cumdiv, stddiv)
      end do

   end subroutine rs_chart

   subroutine rs_analysis(series, H, H_err, sigma2, ierr)
      real(wp), intent(in) :: series(:)
      real(wp), intent(out) :: H, H_err, sigma2
      integer, intent(out), optional :: ierr

      integer :: m
      real(wp), allocatable :: RS(:), N(:)
      real(wp) :: c, c_err
      m = rs_chart_size(size(series))
      allocate(RS(m), N(m))

      call rs_chart(series, RS, N, ierr)
      if (present(ierr) .and. ierr /= 0) return

      call powerregress(RS, N, H, c, sigma2, H_err, c_err, ierr=ierr)
      if (present(ierr) .and. ierr /= 0) return

      deallocate(RS, N)
   end subroutine rs_analysis
end module hurst
