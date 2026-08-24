module generators
   use precision, only: wp
   use stat, only: mean
   use frac, only: frac_diff_simple
   use fourier, only: rfft1d, irfft1d, fft1d, ifft1d
   use integers, only: up2power
   use checks, only: check
   use stdlib_random, only: random_seed
   use stdlib_stats_distribution_normal, only: rvs_normal
   use stdlib_stats_distribution_uniform, only: rvs_uniform
   use stdlib_optval, only: optval
   implicit none

   integer :: default_seed = 42
   real(wp), parameter :: default_sigma = 1.0_wp
   real(wp), parameter :: default_mu = 0.0_wp
contains
   subroutine generate_white(array, mu, sigma, seed)
      real(wp), intent(out) :: array(:)
      real(wp), optional, intent(in) :: mu
      real(wp), optional, intent(in) :: sigma
      integer, optional, intent(in) :: seed

      integer :: n

      n = size(array)
      call random_seed(optval(seed, default_seed), default_seed)

      array(:) = rvs_uniform(optval(mu, default_mu) - optval(sigma, default_sigma) / 2.0_wp, optval(mu, default_mu) + optval(sigma, default_sigma) / 2.0_wp, n)
   end subroutine generate_white

   subroutine generate_gauss(array, mu, sigma, seed)
      real(wp), intent(out) :: array(:)
      real(wp), optional, intent(in) :: mu
      real(wp), optional, intent(in) :: sigma
      integer, optional, intent(in) :: seed

      integer :: n

      n = size(array)

      call random_seed(optval(seed, default_seed), default_seed)

      array(:) = rvs_normal(optval(mu, default_mu), optval(sigma, default_sigma), n)
   end subroutine generate_gauss

   ! Produces series by integrating fgn fractionaly
   subroutine generate_gauss_integrate(series, intorder, seed_in, ierr)
      real(wp), intent(out) :: series(:)
      real(wp), intent(in) :: intorder
      integer, optional, intent(in) :: seed_in
      integer, intent(out), optional :: ierr

      integer :: n
      real(wp), allocatable :: fgn(:)

      n = size(series)
      allocate(fgn(n))

      call generate_gauss(fgn, 0.0_wp, 1.0_wp, seed_in)
      call frac_diff_simple(-intorder, fgn, series, ierr=ierr)

      deallocate(fgn)
   end subroutine generate_gauss_integrate

   ! Produces series by integrating white noise fractionally
   subroutine generate_white_integrate(series, intorder, seed_in, ierr)
      real(wp), intent(out) :: series(:)
      real(wp), intent(in) :: intorder
      integer, optional, intent(in) :: seed_in
      integer, intent(out), optional :: ierr

      integer :: n
      real(wp), allocatable :: noise(:)

      n = size(series)
      allocate(noise(n))

      call generate_white(noise, 0.0_wp, 1.0_wp, seed_in)
      call frac_diff_simple(-intorder, noise, series, ierr=ierr)

      deallocate(noise)
   end subroutine generate_white_integrate

   ! Produces series with PSD~1/f^a
   subroutine generate_color(series, a, seed_in, ierr)
      real(wp), intent(out) :: series(:)
      real(wp), intent(in) :: a
      integer, optional, intent(in) :: seed_in
      integer, intent(out), optional :: ierr

      real(wp), allocatable :: white(:), f(:)
      complex(wp), allocatable :: white_fft(:)
      real(wp), allocatable :: S(:)
      integer :: n, i

      n = size(series)
      n = up2power(n, ierr=ierr)
      if (present(ierr) .and. ierr /= 0) return

      allocate(white(n), white_fft(n/2), S(n/2), f(n/2))

      call generate_white(white, default_mu, default_sigma, seed_in)

      call rfft1d(white, white_fft, ierr=ierr)
      if (present(ierr) .and. ierr /= 0) return

      S(:) = [0.0_wp, 1.0_wp / [((real(i, wp) / n) ** a, i = 1, n/2 - 1)] / n]
      S(:) = sqrt(S)
      S(:) = S / sqrt(mean(S**2))
      white_fft(:) = white_fft * S
      
      call irfft1d(white, white_fft, ierr=ierr)
      if (present(ierr) .and. ierr /= 0) return

      series = white(1:size(series))

      deallocate(white, white_fft, S, f)
   end subroutine generate_color

   subroutine generate_fgn(series, H, seed_in, ierr)
      real(wp), intent(out) :: series(:)
      real(wp), intent(in) :: H
      integer, optional, intent(in) :: seed_in
      integer, intent(out), optional :: ierr

      real(wp), allocatable :: m(:), n(:), rho(:)
      complex(wp), allocatable :: rhofft(:), V(:), W(:)
      integer :: k, j, l

      if(check(H > 0.0_wp, msg="generate_fgn: invalid H", ierr=ierr)) return
      if(check(H < 1.0_wp, msg="generate_fgn: invalid H", ierr=ierr)) return

      l = up2power(size(series))
      allocate(m(l), n(l), rho(l))
      allocate(rhofft(l * 2), V(l * 2), W(l * 2))

      ! The autocorrelation function of the FGN sequence
      do k = 0,size(rho)-1
         rho(k+1) = 0.5_wp * (abs(k - 1) ** (2.0_wp * H) - 2.0_wp * k ** (2.0_wp * H) + (k + 1) ** (2.0_wp * H))
      end do

      rhofft = [rho(1:size(rho)), 0.0_wp, rho(size(rho):2:-1)]
      call fft1d(rhofft, ierr=ierr)
      if (present(ierr) .and. ierr /= 0) return

      ! Eigenvalues of the correlation sequence
      V = sqrt(rhofft)

      call generate_gauss(m, default_mu, default_sigma, seed_in)
      call generate_gauss(n, default_mu, default_sigma, default_seed)

      W(:) = 0.0_wp
      W(1) = V(1) / sqrt(2.0_wp * l) * m(1)
      do j = 2,l
         W(j) = V(j) / sqrt(4.0_wp * l) * (m(j) + cmplx(0.0_wp, 1.0_wp, kind=wp) * n(j))
         W(j + l) = V(j + l) / sqrt(4.0_wp * l) * (m(l - j + 2) - cmplx(0.0_wp, 1.0_wp, kind=wp) * n(l - j + 2))
      end do
      W(l + 1) = V(l + 1) / sqrt(2.0_wp * l) * n(1)

      call ifft1d(W, ierr=ierr)
      if (present(ierr) .and. ierr /= 0) return

      series = l ** (-H) * real(W(1:size(series)))

      deallocate(m, n, rho, rhofft, V, W)
   end subroutine generate_fgn
end module generators
