module spectra
   use precision, only: wp
   use autoreg, only: complex_yw_ar_coeff, complex_burg_ar_coeff, complex_ar_freq_response
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

   ! Variance of complex series `S` with correction `i` added to sample size
   ! (same convention as `stat::variance`)
   real(wp) function complex_variance(S, i) result(sigma2)
      complex(wp), intent(in) :: S(:)
      real(wp), intent(in) :: i

      complex(wp) :: mu
      real(wp) :: d

      d = real(size(S), wp) + i
      mu = sum(S) / size(S)
      sigma2 = sum(abs(S - mu) ** 2) / d
   end function complex_variance

   subroutine complex_ar_psd(f, P, S, dt, phi, ierr)
      real(wp), intent(out) :: P(:)
      real(wp), intent(out) :: f(:)
      complex(wp), intent(in) :: S(:)
      complex(wp), intent(in) :: phi(:)
      real(wp), intent(in) :: dt
      integer, intent(out), optional :: ierr

      complex(wp), allocatable :: H(:)

      if (check(size(S) > size(phi), msg="complex_ar_psd: size(S) must be larger then size(phi)", ierr=ierr)) return
      if (check(size(f) == psd_size(size(S)), msg="complex_ar_psd: size missmatch", ierr=ierr)) return
      if (check(size(P) == psd_size(size(S)), msg="complex_ar_psd: size missmatch", ierr=ierr)) return

      allocate(H(size(P)))

      ! Calculate the frequency response
      call freq_nyquist(f, dt)
      call complex_ar_freq_response(phi, f, H, ierr=ierr)
      if (present(ierr) .and. ierr /= 0) return

      ! Calculate the PSD
      P = abs(H) ** 2
      P = P / (f(2) - f(1))
      P = P / sum(P) * complex_variance(S, 0.0_wp)

      deallocate(H)
   end subroutine complex_ar_psd

   subroutine ar_psd(f, P, S, dt, phi, ierr)
      real(wp), intent(out) :: P(:), f(:)
      real(wp), intent(in) :: S(:), phi(:), dt
      integer, intent(out), optional :: ierr

      complex(wp), allocatable :: cS(:), cphi(:)

      allocate(cS(size(S)), cphi(size(phi)))
      cS = cmplx(S, 0.0_wp, kind=wp)
      cphi = cmplx(phi, 0.0_wp, kind=wp)
      call complex_ar_psd(f, P, cS, dt, cphi, ierr=ierr)
      deallocate(cS, cphi)
   end subroutine ar_psd

   ! Calculate power spectrum density of a complex series `S` with even time step `dt` using Berg method of order `m`
   subroutine complex_berg_psd(f, P, S, dt, m, ierr)
      real(wp), intent(out) :: f(:), P(:)
      complex(wp), intent(in) :: S(:)
      real(wp), intent(in) :: dt
      integer, intent(in) :: m
      integer, intent(out), optional :: ierr

      complex(wp), allocatable :: phi(:)

      allocate(phi(m))

      call complex_burg_ar_coeff(phi, S)

      call complex_ar_psd(f, P, S, dt, phi, ierr=ierr)
      if (present(ierr) .and. ierr /= 0) return

      deallocate(phi)
   end subroutine complex_berg_psd

   ! Calculate power spectrum density of a series `S` with even time step `dt` using Berg method of order `m`.
   ! Real data is treated as complex with zero imaginary part.
   subroutine berg_psd(f, P, S, dt, m, ierr)
      real(wp), intent(out) :: f(:), P(:)
      real(wp), intent(in) :: S(:), dt
      integer, intent(in) :: m
      integer, intent(out), optional :: ierr

      complex(wp), allocatable :: cS(:)

      allocate(cS(size(S)))
      cS = cmplx(S, 0.0_wp, kind=wp)
      call complex_berg_psd(f, P, cS, dt, m, ierr=ierr)
      deallocate(cS)
   end subroutine berg_psd

   ! Calculate power spectrum density of a complex series `S` with even time step `dt` using Yule-Walker method of order `m`
   subroutine complex_yw_psd(f, P, S, dt, m, ierr)
      real(wp), intent(out) :: f(:), P(:)
      complex(wp), intent(in) :: S(:)
      real(wp), intent(in) :: dt
      integer, intent(in) :: m
      integer, intent(out), optional :: ierr

      complex(wp), allocatable :: phi(:)

      allocate(phi(m))

      call complex_yw_ar_coeff(phi, S, ierr=ierr)
      if (present(ierr) .and. ierr /= 0) return

      call complex_ar_psd(f, P, S, dt, phi, ierr=ierr)
      if (present(ierr) .and. ierr /= 0) return

      deallocate(phi)
   end subroutine complex_yw_psd

   ! Calculate power spectrum density of a series `S` with even time step `dt` using Yule-Walker method of order `m`.
   ! Real data is treated as complex with zero imaginary part.
   subroutine yw_psd(f, P, S, dt, m, ierr)
      real(wp), intent(out) :: f(:), P(:)
      real(wp), intent(in) :: S(:), dt
      integer, intent(in) :: m
      integer, intent(out), optional :: ierr

      complex(wp), allocatable :: cS(:)

      allocate(cS(size(S)))
      cS = cmplx(S, 0.0_wp, kind=wp)
      call complex_yw_psd(f, P, cS, dt, m, ierr=ierr)
      deallocate(cS)
   end subroutine yw_psd
end module spectra
