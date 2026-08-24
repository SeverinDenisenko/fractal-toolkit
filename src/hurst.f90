module hurst
   use precision, only: wp
   use solvers, only: powerregress
   use spectra, only: berg_psd, yw_psd, psd_size
   use stat, only: mean, variance
   use checks, only: check
   use math, only: log2
   use stdlib_optval, only: optval
   use stdlib_math, only: swap
   implicit none

   private
   public :: slope_to_hurst, estimate_hurst_psd, estimate_hurst_berg, estimate_hurst_yw
   public :: rs_chart_size, rs_chart, estimate_hurst_rs, estimate_hurst_lssd

contains
   ! Compute Hurst expoenent from PSD slope `a` where PSD~1/f^a
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

   subroutine estimate_hurst_rs(series, H, H_err, sigma2, ierr)
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
   end subroutine estimate_hurst_rs

   real(wp) function cm_lssd(m, n, H) result(c)
      integer, intent(in) :: m, n
      real(wp), intent(in) :: H

      real(wp) :: u
      u = real(n, wp) / real(m, wp)
      c = sqrt((u - u ** (2.0_wp * H - 1.0_wp)) / (u - 0.5_wp))
   end function cm_lssd

   real(wp) function dm_h_lssd(m, n, H) result(d)
      integer, intent(in) :: m, n
      real(wp), intent(in) :: H

      real(wp) :: u
      u = real(n, wp) / real(m, wp)
      d = log(real(m, kind=wp)) + log(u) / (1.0_wp - u ** (2.0_wp - 2.0_wp * H))
   end function dm_h_lssd

   ! Contractive Mapping for the LSSD method
   real(wp) function ctm_lssd(H, n, p, q, T, S) result(g)
      real(wp), intent(in) :: H
      integer, intent(in) :: n, p, q
      integer, intent(in) :: T(:)
      real(wp), intent(in) :: S(:)

      integer :: m, mmax, idx
      real(wp) :: a11, a12, a21, a22, b1, b2, sm, cm, dm, u

      mmax = size(T)
      a11 = 0
      a12 = 0
      a21 = 0
      a22 = 0
      b1 = 0
      b2 = 0

      do idx = 1,mmax
         m = T(idx)
         sm = S(idx)
         cm = cm_lssd(m, n, H)
         dm = dm_h_lssd(m, n, H)
         u = real(m, wp) ** p
         a11 = a11 + 1.0_wp / u
         a12 = a12 + log(real(m, wp)) / u
         a21 = a21 + dm / u
         a22 = a22 + dm * log(real(m, wp)) / u
         b1 = b1 + (log(sm) - log(cm)) / u
         b2 = b2 + dm * (log(sm) - log(cm)) / u
      end do

      g = (a11 * (b2 - H ** q) - a21 * b1) / (a11 * a22 - a21 * a12)
   end function ctm_lssd

   ! Least Squares via Standard Deviation (Koutsoyiannis estimator)
   subroutine estimate_hurst_lssd(X, p, q, H, H_err, eps, maxiter, ierr)
      real(wp), intent(in) :: X(:)
      integer, intent(in) :: p, q
      real(wp), intent(out) :: H, H_err
      real(wp), intent(in), optional :: eps
      integer, intent(in), optional :: maxiter
      integer, intent(inout), optional :: ierr

      integer :: n, m, idx, k, i, mmax, begin, end, curriter
      real(wp), allocatable :: S(:), Z(:)
      integer, allocatable :: T(:)
      real(wp) :: H0

      H = 0.0_wp
      H_err = 0.0_wp

      n = size(X)
      if(check(n >= 10, msg="estimate_hurst_lssd: series too short", ierr=ierr)) return
      if(check(p >= 0 .and. p <= 10, msg="estimate_hurst_lssd: invalid weight p", ierr=ierr)) return
      if(check(q >= 1 .and. q <= 10000, msg="estimate_hurst_lssd: invalid penalty q", ierr=ierr)) return

      mmax = n / 10
      allocate(S(mmax), T(mmax))

      T = [(m, m = 1, mmax)]
      do idx = 1, size(T)
         m = T(idx)
         k = n / m
         allocate(Z(k))
         do i = 1, k
            begin = (i - 1) * m + 1
            end = (i - 1) * m + m
            Z(i) = sum(X(begin:end))
         end do
         S(idx) = sqrt(variance(Z))
         deallocate(Z)
      end do

      if(check(all(S > 0.0_wp), msg="estimate_hurst_lssd: non-positive std deviation", ierr=ierr)) return

      curriter = 0
      H = 0.5_wp
      do
         H0 = H
         H = ctm_lssd(H0, n, p, q, T, S)
         if(check(H > 0.0_wp .and. H < 1.0_wp, msg="estimate_hurst_lssd: H out of (0,1)", ierr=ierr)) return
         if (abs(H0 - H) < optval(eps, 0.01_wp)) then
            exit
         end if
         curriter = curriter + 1
         if(check(curriter < optval(maxiter, 1000), msg="estimate_hurst_lssd: max iteration", ierr=ierr)) return
      end do
      H_err = abs(H - H0)

      deallocate(T, S)
   end subroutine estimate_hurst_lssd
end module hurst
