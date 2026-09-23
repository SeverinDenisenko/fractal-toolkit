module hurst
   use precision, only: wp
   use solvers, only: powerregress
   use spectra, only: berg_psd, complex_yw_psd, psd_size
   use stat, only: mean, variance
   use checks, only: check
   use math, only: log2
   use stdlib_optval, only: optval
   use stdlib_math, only: swap
   implicit none

   private
   public :: slope_to_hurst, estimate_hurst_psd, estimate_hurst_berg, estimate_hurst_yw
   public :: rs_chart_size, rs_chart, estimate_hurst_rs, estimate_hurst_lssd
   public :: complex_estimate_hurst_psd, complex_estimate_hurst_yw
   public :: complex_rs_chart, complex_estimate_hurst_rs
   public :: dfa_chart_size, dfa_chart, estimate_hurst_dfa
   public :: tta_chart_size, tta_chart, estimate_hurst_tta

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

   subroutine complex_estimate_hurst_psd(f, P, H, a, H_err, a_err, sigma2, ierr)
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
   end subroutine complex_estimate_hurst_psd

   subroutine estimate_hurst_psd(f, P, H, a, H_err, a_err, sigma2, ierr)
      real(wp), intent(in) :: f(:), P(:)
      real(wp), intent(out) :: H, a, H_err, a_err, sigma2
      integer, intent(out), optional :: ierr

      call complex_estimate_hurst_psd(f, P, H, a, H_err, a_err, sigma2, ierr=ierr)
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

   subroutine complex_estimate_hurst_yw(series, m, H, a, H_err, a_err, sigma2, ierr)
      complex(wp), intent(in) :: series(:)
      integer, intent(in) :: m
      real(wp), intent(out) :: H, a, H_err, a_err, sigma2
      integer, intent(out), optional :: ierr

      real(wp), allocatable :: P(:), f(:)
      integer :: n

      n = size(series)

      allocate(f(psd_size(n)), P(psd_size(n)))

      call complex_yw_psd(f, P, series, 1.0_wp, m, ierr=ierr)
      if (present(ierr) .and. ierr /= 0) return

      call complex_estimate_hurst_psd(f, P, H, a, H_err, a_err, sigma2, ierr=ierr)

      deallocate(f, P)
   end subroutine complex_estimate_hurst_yw

   subroutine estimate_hurst_yw(series, m, H, a, H_err, a_err, sigma2, ierr)
      real(wp), intent(in) :: series(:)
      integer, intent(in) :: m
      real(wp), intent(out) :: H, a, H_err, a_err, sigma2
      integer, intent(out), optional :: ierr

      complex(wp), allocatable :: cseries(:)

      allocate(cseries(size(series)))
      cseries = cmplx(series, 0.0_wp, kind=wp)
      call complex_estimate_hurst_yw(cseries, m, H, a, H_err, a_err, sigma2, ierr=ierr)
      deallocate(cseries)
   end subroutine estimate_hurst_yw

   integer function rs_chart_size(n) result(m)
      integer :: n
      m = int(log2(real(n, wp)))
   end function rs_chart_size

   subroutine complex_rs_process_chunk(data, n, n_real, range_val, std_dev, mean_val)
      complex(wp), intent(in) :: data(:)
      integer, intent(in) :: n
      real(wp), intent(in) :: n_real
      complex(wp), intent(out) :: mean_val
      real(wp), intent(out) :: range_val, std_dev

      complex(wp), allocatable :: cum(:)
      real(wp) :: d2, r2, sum_sq_val, variance_val, dx, dy
      complex(wp) :: sum_dev
      integer :: k, j

      mean_val = sum(data) / n_real

      ! Cumulative deviations C_k = sum_{m<=k} (data(m) - mean_val)
      allocate(cum(n))
      cum(1) = data(1) - mean_val
      do k = 2, n
         cum(k) = cum(k - 1) + (data(k) - mean_val)
      end do

      ! Range as the diameter of the cumulative walk in the complex plane.
      ! For real data this reduces exactly to max(C_k) - min(C_k).
      r2 = 0.0_wp
      do j = 1, n - 1
         do k = j + 1, n
            dx = real(cum(k) - cum(j), wp)
            dy = aimag(cum(k) - cum(j))
            d2 = dx * dx + dy * dy
            if (d2 > r2) r2 = d2
         end do
      end do
      range_val = sqrt(r2)

      sum_dev = (0.0_wp, 0.0_wp)
      sum_sq_val = 0.0_wp
      do k = 1, n
         sum_dev = sum_dev + (data(k) - mean_val)
         sum_sq_val = sum_sq_val + abs(data(k) - mean_val)**2
      end do

      if (n > 1) then
         variance_val = (sum_sq_val - abs(sum_dev)**2 / n_real) / (n_real - 1.0_wp)
         std_dev = sqrt(max(variance_val, 0.0_wp))
      else
         std_dev = 0.0_wp
      endif

      deallocate(cum)
   end subroutine complex_rs_process_chunk

   subroutine complex_rs_chart(series, RS, N, ierr)
      complex(wp), intent(in) :: series(:)
      real(wp), intent(out) :: RS(:), N(:)
      integer, intent(out), optional :: ierr

      integer :: i, j, n_chunks, n_curr
      integer :: start_idx, end_idx
      real(wp), allocatable :: cumdiv(:), stddiv(:)
      complex(wp) :: mean_val
      real(wp) :: R
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

            call complex_rs_process_chunk(series(start_idx:end_idx), n_curr, n_curr_real, R, stddiv(i), mean_val)

            if (stddiv(i) > 0.0_wp) then
               cumdiv(i) = R / stddiv(i)
            else
               cumdiv(i) = 0.0_wp
            endif
         end do

         RS(j) = sum(cumdiv) / real(n_chunks, wp)

         deallocate(cumdiv, stddiv)
      end do
   end subroutine complex_rs_chart

   subroutine rs_chart(series, RS, N, ierr)
      real(wp), intent(in) :: series(:)
      real(wp), intent(out) :: RS(:), N(:)
      integer, intent(out), optional :: ierr

      complex(wp), allocatable :: cseries(:)

      allocate(cseries(size(series)))
      cseries = cmplx(series, 0.0_wp, kind=wp)
      call complex_rs_chart(cseries, RS, N, ierr=ierr)
      deallocate(cseries)
   end subroutine rs_chart

   subroutine complex_estimate_hurst_rs(series, H, H_err, sigma2, ierr)
      complex(wp), intent(in) :: series(:)
      real(wp), intent(out) :: H, H_err, sigma2
      integer, intent(out), optional :: ierr

      integer :: m
      real(wp), allocatable :: RS(:), N(:)
      real(wp) :: c, c_err
      m = rs_chart_size(size(series))
      allocate(RS(m), N(m))

      call complex_rs_chart(series, RS, N, ierr)
      if (present(ierr) .and. ierr /= 0) return

      call powerregress(N, RS, H, c, sigma2, H_err, c_err, ierr=ierr)
      if (present(ierr) .and. ierr /= 0) return

      deallocate(RS, N)
   end subroutine complex_estimate_hurst_rs

   subroutine estimate_hurst_rs(series, H, H_err, sigma2, ierr)
      real(wp), intent(in) :: series(:)
      real(wp), intent(out) :: H, H_err, sigma2
      integer, intent(out), optional :: ierr

      complex(wp), allocatable :: cseries(:)

      allocate(cseries(size(series)))
      cseries = cmplx(series, 0.0_wp, kind=wp)
      call complex_estimate_hurst_rs(cseries, H, H_err, sigma2, ierr=ierr)
      deallocate(cseries)
   end subroutine estimate_hurst_rs

   integer function dfa_chart_size(n) result(m)
      integer, intent(in) :: n

      m = int(log2(real(n, wp))) - 1
   end function dfa_chart_size

   ! Detrended Fluctuation Analysis (Peng et al.): S(m) ~ c * m^H
   ! Box sizes N = 4, 8, ..., ~2^floor(log2(n))
   subroutine dfa_chart(series, F, N, ierr)
      real(wp), intent(in) :: series(:)
      real(wp), intent(out) :: F(:), N(:)
      integer, intent(out), optional :: ierr

      integer :: i, itau, j, k, m, npts, n_sizes
      integer :: start_idx, end_idx
      real(wp), allocatable :: profile(:)
      real(wp) :: sm1, sqsum, delta, xbar
      real(wp) :: sy, sxy, sy2, b, a, ss_res, s_tau

      npts = size(series)
      n_sizes = dfa_chart_size(npts)

      if (check(n_sizes > 0, msg="dfa_chart: series too short", ierr=ierr)) return
      if (check(size(N) == n_sizes, msg="dfa_chart: size mismatch in N", ierr=ierr)) return
      if (check(size(F) == n_sizes, msg="dfa_chart: size mismatch in F", ierr=ierr)) return

      N = [(2.0_wp ** i, i = 2, n_sizes + 1)]

      allocate(profile(npts))
      xbar = mean(series)
      profile(1) = series(1) - xbar
      do i = 2, npts
         profile(i) = profile(i - 1) + series(i) - xbar
      end do

      do j = 1, n_sizes
         m = int(N(j))
         k = npts / m
         sm1 = real(m, wp) * real(m + 1, wp) / 2.0_wp
         sqsum = real(m, wp) * real(m + 1, wp) * real(2 * m + 1, wp) / 6.0_wp
         delta = real(m, wp) * sqsum - sm1 * sm1

         s_tau = 0.0_wp
         do itau = 1, k
            start_idx = (itau - 1) * m + 1
            end_idx = start_idx + m - 1

            sy = sum(profile(start_idx:end_idx))
            sy2 = sum(profile(start_idx:end_idx)**2)
            sxy = 0.0_wp
            do i = 1, m
               sxy = sxy + real(i, wp) * profile(start_idx + i - 1)
            end do

            b = (real(m, wp) * sxy - sm1 * sy) / delta
            a = (sy - b * sm1) / real(m, wp)
            ss_res = sy2 - a * sy - b * sxy
            s_tau = s_tau + sqrt(max(ss_res, 0.0_wp) / real(m - 1, wp))
         end do

         F(j) = s_tau / real(k, wp)
      end do

      deallocate(profile)
   end subroutine dfa_chart

   subroutine estimate_hurst_dfa(series, H, H_err, sigma2, ierr)
      real(wp), intent(in) :: series(:)
      real(wp), intent(out) :: H, H_err, sigma2
      integer, intent(out), optional :: ierr

      integer :: m
      real(wp), allocatable :: F(:), N(:)
      real(wp) :: c, c_err

      m = dfa_chart_size(size(series))
      allocate(F(m), N(m))

      call dfa_chart(series, F, N, ierr)
      if (present(ierr) .and. ierr /= 0) return

      call powerregress(N, F, H, c, sigma2, H_err, c_err, ierr=ierr)
      if (present(ierr) .and. ierr /= 0) return

      deallocate(F, N)
   end subroutine estimate_hurst_dfa

   ! Triangles Total Areas (TTA) method (Lotfalinezhad & Maleki, 2020):
   ! profile Y_i = sum_{t<=i} (X_t - xbar), total triangle area over lag tau
   ! S(tau) = (tau/2) * sum_i |Y_{i+2tau} - 2Y_{i+tau} + Y_i| ~ c * tau^H
   integer function tta_chart_size(n) result(m)
      integer, intent(in) :: n

      m = min(10, (n - 1) / 2)
   end function tta_chart_size

   subroutine tta_chart(series, S, T, ierr)
      real(wp), intent(in) :: series(:)
      real(wp), intent(out) :: S(:), T(:)
      integer, intent(out), optional :: ierr

      integer :: i, tau, npts, n_lags, k
      real(wp), allocatable :: Y(:)
      real(wp) :: xbar, stau

      npts = size(series)
      n_lags = tta_chart_size(npts)

      if (check(n_lags > 0, msg="tta_chart: series too short", ierr=ierr)) return
      if (check(size(T) == n_lags, msg="tta_chart: size mismatch in T", ierr=ierr)) return
      if (check(size(S) == n_lags, msg="tta_chart: size mismatch in S", ierr=ierr)) return

      T = [(real(tau, wp), tau = 1, n_lags)]

      allocate(Y(npts))
      xbar = mean(series)
      Y(1) = series(1) - xbar
      do i = 2, npts
         Y(i) = Y(i - 1) + series(i) - xbar
      end do

      do tau = 1, n_lags
         k = (npts - 1) / (2 * tau)
         stau = 0.0_wp
         do i = 1, k
            stau = stau + abs(Y(2 * (i - 1) * tau + 1 + 2 * tau) - 2.0_wp * Y(2 * (i - 1) * tau + 1 + tau) + Y(2 * (i - 1) * tau + 1))
         end do
         S(tau) = real(tau, wp) * stau / 2.0_wp
      end do

      deallocate(Y)
   end subroutine tta_chart

   subroutine estimate_hurst_tta(series, H, H_err, sigma2, ierr)
      real(wp), intent(in) :: series(:)
      real(wp), intent(out) :: H, H_err, sigma2
      integer, intent(out), optional :: ierr

      integer :: m
      real(wp), allocatable :: S(:), T(:)
      real(wp) :: c, c_err

      m = tta_chart_size(size(series))
      allocate(S(m), T(m))

      call tta_chart(series, S, T, ierr)
      if (present(ierr) .and. ierr /= 0) return

      call powerregress(T, S, H, c, sigma2, H_err, c_err, ierr=ierr)
      if (present(ierr) .and. ierr /= 0) return

      deallocate(S, T)
   end subroutine estimate_hurst_tta

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
