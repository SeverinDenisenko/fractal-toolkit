module spectra
   use precision, only: wp
   use int_math, only: up2power
   use stdlib_linalg, only: solve_lstsq
   use stdlib_error, only: check
   use stdlib_stats, only: var
   use stdlib_constants, only: PI_dp
   implicit none

contains
   ! Calculate AR coefficients `phi` for time series `S` by solving linear system of equatins
   subroutine berg_ar_coeff(phi, S)
      real(wp), intent(out) :: phi(:)
      real(wp), intent(in) :: S(:)

      integer :: n, p, i, j
      real(wp), allocatable :: X(:,:)

      p = size(phi)
      n = size(S) - p

      allocate(X(n, p))

      do i = 1, n
         do j = 1, p
            X(i, j) = S(i + p - j)
         end do
      end do

      call check(size(S(p+1:)) == n)
      call solve_lstsq(X, S(p+1:), x=phi)

      deallocate(X)
   end subroutine berg_ar_coeff

   ! Compute prediction of timeseries `S` for AR model with coefficients `phi` at index `i` from previous values (not including `i`)
   real(wp) function ar_predict(phi, S, i) result(x)
      real(wp), intent(in) :: phi(:)
      real(wp), intent(in) :: S(:)
      integer, intent(in) :: i

      x = dot_product(phi, S(i - 1:i - size(phi):-1))
   end function ar_predict

   ! Compute optimal frequencies amount for time series with length `n`
   integer function berg_psd_size(n) result(m)
      integer, intent(in) :: n
      m = n / 2 + 1
   end function berg_psd_size

   ! Calculate normalized (0 to 1) frequencies from 0 up to 2 * Nyquist frequency for time step `dt`
   subroutine freq_full(f, dt)
      real(wp), intent(out) :: f(:)
      real(wp), intent(in) :: dt

      integer :: i, n
      real(wp) :: step

      n = size(f)
      step = 1.0_wp / dt
      f = [(step * i / n, i = 0, n - 1)]
   end subroutine freq_full

   ! Calculate normalized (0 to 1) frequencies from 0 up to Nyquist frequency (including) for time step `dt`
   subroutine freq_nyquist(f, dt)
      real(wp), intent(out) :: f(:)
      real(wp), intent(in) :: dt

      integer :: i, n
      real(wp) :: nyquist

      n = size(f)
      nyquist = 1.0_wp / dt / 2.0_wp
      f = [(nyquist * i / (n - 1), i = 0, n - 1)]
   end subroutine freq_nyquist

   ! Calculate frequency response of an AR filter `phi` on frequencies `f` (0 to 1)
   subroutine ar_freq_response(phi, f, H)
      real(wp), intent(in) :: phi(:)
      real(wp), intent(in) :: f(:)
      complex(wp), intent(out) :: H(:)

      integer :: p, n
      integer :: i, k
      real(wp) :: omega
      complex(wp) :: denominator

      p = size(phi)
      n = size(f)

      call check(n == size(H), msg="ar_freq_response: size mismatch")

      do i = 1, n
         omega = f(i) * 2.0_wp * PI_dp

         denominator = cmplx(1.0_wp, 0.0_wp, kind=wp)
         do k = 1, p
            denominator = denominator - phi(k) * exp(-cmplx(0.0_wp, 1.0_wp, kind=wp) * k * omega)
         end do

         H(i) = cmplx(1.0_wp, 0.0_wp, kind=wp) / denominator
      end do
   end subroutine ar_freq_response

   ! Calculate power spectrum density of a series `S` with even time step `dt` using Berg method of order `m`
   subroutine berg_psd(f, P, S, dt, m)
      real(wp), intent(out) :: P(:)
      real(wp), intent(out) :: f(:)
      real(wp), intent(in) :: S(:)
      real(wp), intent(in) :: dt
      integer, intent(in) :: m

      integer :: i
      real(wp), allocatable :: phi(:)
      real(wp), allocatable :: predict(:)
      complex(wp), allocatable :: H(:)
      real(wp) :: sigma2

      call check(size(S) > m, msg="berg_psd: size(S) must be larger then m")
      call check(size(f) == berg_psd_size(size(S)), msg="berg_psd: size missmatch")
      call check(size(P) == berg_psd_size(size(S)), msg="berg_psd: size missmatch")

      allocate(predict(size(S) - m))
      allocate(phi(m))
      allocate(H(size(P)))

      ! Calculate AR
      call berg_ar_coeff(phi, S)

      ! Calculate the prediction and prediction error
      do i=1,size(S)-m
         predict(i) = ar_predict(phi, S, i + m)
      end do
      sigma2 = var(S(m+1:) - predict(:))

      ! Calculate the frequency response
      call freq_nyquist(f, dt)
      call ar_freq_response(phi, f, H)

      ! Calculate the PSD
      P = sigma2 * abs(H) ** 2
      P = P / (f(2) - f(1))
      P = P / sum(P) * var(S, corrected=.false.)
   end subroutine berg_psd

end module spectra
