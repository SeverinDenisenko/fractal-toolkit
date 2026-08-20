module spectra
   use precision, only: wp
   use integers, only: up2power
   use autoreg, only: ar_coeff, ar_predict, ar_freq_response
   use stdlib_error, only: check
   use stdlib_stats, only: var
   implicit none

contains
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
      call ar_coeff(phi, S)

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

      deallocate(predict, phi, H)
   end subroutine berg_psd
end module spectra
