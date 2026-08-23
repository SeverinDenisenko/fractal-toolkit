module hurst
   use precision, only: wp
   use solvers, only: powerregress
   use spectra, only: berg_psd, yw_psd, psd_size
   use stat, only: mean, variance
   use checks, only: check
   use math, only: log2
   implicit none

   private
   public :: slope_to_hurst, estimate_hurst_psd, estimate_hurst_berg, estimate_hurst_yw
   public :: rs_chart_size, rs_chart, rs_analysis

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
      m = int(log2(real(n, wp)))
   end function rs_chart_size

   subroutine rs_process_chunk(data, n, n_real, range_val, std_dev, mean_val)
      real(wp), intent(in) :: data(:)
      integer, intent(in) :: n
      real(wp), intent(in) :: n_real
      real(wp), intent(out) :: range_val, std_dev, mean_val

      real(wp) :: sum_val, sum_sq_val, cum_sum, cum_max, cum_min
      real(wp) :: variance_val
      integer :: k

      sum_val = sum(data)
      mean_val = sum_val / n_real

      cum_sum = 0.0_wp
      cum_max = -huge(1.0_wp)
      cum_min = huge(1.0_wp)
      sum_val = 0.0_wp
      sum_sq_val = 0.0_wp

      do k = 1, n
         cum_sum = cum_sum + (data(k) - mean_val)
         cum_max = max(cum_max, cum_sum)
         cum_min = min(cum_min, cum_sum)

         sum_val = sum_val + (data(k) - mean_val)
         sum_sq_val = sum_sq_val + (data(k) - mean_val)**2
      end do

      range_val = cum_max - cum_min

      if (n > 1) then
         variance_val = (sum_sq_val - sum_val**2 / n_real) / (n_real - 1.0_wp)
         std_dev = sqrt(max(variance_val, 0.0_wp))
      else
         std_dev = 0.0_wp
      endif
   end subroutine rs_process_chunk

   subroutine rs_chart(series, RS, N, ierr)
      real(wp), intent(in) :: series(:)
      real(wp), intent(out) :: RS(:), N(:)
      integer, intent(out), optional :: ierr

      integer :: i, j, n_chunks, n_curr
      integer :: start_idx, end_idx
      real(wp), allocatable :: cumdiv(:), stddiv(:)
      real(wp) :: R, mean_val
      real(wp) :: n_curr_real

      N = [(2 ** i, i = 1, rs_chart_size(size(series)))]

      if(check(size(N) == rs_chart_size(size(series)), msg="rs_chart: size mismatch in N", ierr=ierr)) return
      if(check(size(RS) == rs_chart_size(size(series)), msg="rs_chart: size mismatch in RS", ierr=ierr)) return

      do j = 1, size(N)
         n_curr = int(N(j))
         n_curr_real = real(n_curr, wp)
         n_chunks = size(series) / n_curr

         allocate(cumdiv(n_chunks), stddiv(n_chunks))

         do i = 1, n_chunks
            start_idx = (i - 1) * n_curr + 1
            end_idx = start_idx + n_curr - 1

            if (end_idx > size(series)) exit

            call rs_process_chunk(series(start_idx:end_idx), n_curr, n_curr_real, R, stddiv(i), mean_val)

            if (stddiv(i) > 0.0_wp) then
               cumdiv(i) = R / stddiv(i)
            else
               cumdiv(i) = 0.0_wp
            endif
         end do

         RS(j) = sum(cumdiv) / real(n_chunks, wp)

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

      call powerregress(N, RS, H, c, sigma2, H_err, c_err, ierr=ierr)
      if (present(ierr) .and. ierr /= 0) return

      deallocate(RS, N)
   end subroutine rs_analysis
end module hurst
