module spectra
   use precision, only: wp
   use integers, only: up2power
   use autoreg, only: yw_ar_coeff, burg_ar_coeff, ar_predict, ar_freq_response
   use stat, only: variance
   use checks, only: check
   implicit none

contains
   ! Compute optimal frequencies amount for time series with length `n`
   integer function psd_size(n) result(m)
      integer, intent(in) :: n
      m = n / 2 + 1
   end function psd_size

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

   subroutine ar_psd(f, P, S, dt, phi, ierr)
      real(wp), intent(out) :: P(:)
      real(wp), intent(out) :: f(:)
      real(wp), intent(in) :: S(:)
      real(wp), intent(in) :: phi(:)
      real(wp), intent(in) :: dt
      integer, intent(out), optional :: ierr

      complex(wp), allocatable :: H(:)

      if (check(size(S) > size(phi), msg="ar_psd: size(S) must be larger then size(phi)", ierr=ierr)) return
      if (check(size(f) == psd_size(size(S)), msg="ar_psd: size missmatch", ierr=ierr)) return
      if (check(size(P) == psd_size(size(S)), msg="ar_psd: size missmatch", ierr=ierr)) return

      allocate(H(size(P)))

      ! Calculate the frequency response
      call freq_nyquist(f, dt)
      call ar_freq_response(phi, f, H, ierr=ierr)
      if (present(ierr) .and. ierr /= 0) return

      ! Calculate the PSD
      P = abs(H) ** 2
      P = P / (f(2) - f(1))
      P = P / sum(P) * variance(S, i=0.0_wp)

      deallocate(H)
   end subroutine ar_psd

   ! Calculate power spectrum density of a series `S` with even time step `dt` using Berg method of order `m`
   subroutine berg_psd(f, P, S, dt, m, ierr)
      real(wp), intent(out) :: f(:), P(:)
      real(wp), intent(in) :: S(:), dt
      integer, intent(in) :: m
      integer, intent(out), optional :: ierr

      real(wp), allocatable :: phi(:)

      allocate(phi(m))

      call burg_ar_coeff(phi, S)
   
      call ar_psd(f, P, S, dt, phi, ierr=ierr)
      if (present(ierr) .and. ierr /= 0) return

      deallocate(phi)
   end subroutine berg_psd

   ! Calculate power spectrum density of a series `S` with even time step `dt` using Yule-Walker method of order `m`
   subroutine yw_psd(f, P, S, dt, m, ierr)
      real(wp), intent(out) :: f(:), P(:)
      real(wp), intent(in) :: S(:), dt
      integer, intent(in) :: m
      integer, intent(out), optional :: ierr

      real(wp), allocatable :: phi(:)

      allocate(phi(m))

      call yw_ar_coeff(phi, S, ierr=ierr)
      if (present(ierr) .and. ierr /= 0) return

      call ar_psd(f, P, S, dt, phi, ierr=ierr)
      if (present(ierr) .and. ierr /= 0) return

      deallocate(phi)
   end subroutine yw_psd
end module spectra
