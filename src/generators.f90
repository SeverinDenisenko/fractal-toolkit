module generators
   use precision, only: wp
   use stat, only: mean
   use frac, only: frac_diff_simple
   use fourier, only: rfft1d, irfft1d
   use integers, only: up2power
   use stdlib_random, only: random_seed
   use stdlib_stats_distribution_normal, only: rvs_normal
   use stdlib_stats_distribution_uniform, only: rvs_uniform
   use stdlib_optval, only: optval
   implicit none

   integer, parameter :: default_seed = 42
   real(wp), parameter :: default_sigma = 1.0_wp
   real(wp), parameter :: default_mu = 0.0_wp
contains
   subroutine generate_white(array, mu, sigma, seed)
      real(wp), intent(out) :: array(:)
      real(wp), optional, intent(in) :: mu
      real(wp), optional, intent(in) :: sigma
      integer, optional, intent(in) :: seed

      integer :: seed_get
      integer :: n

      n = size(array)
      call random_seed(optval(seed, default_seed), seed_get)

      array(:) = rvs_uniform(optval(mu, default_mu) - optval(sigma, default_sigma) / 2.0_wp, optval(mu, default_mu) + optval(sigma, default_sigma) / 2.0_wp, n)
   end subroutine generate_white

   subroutine generate_fgn(array, mu, sigma, seed)
      real(wp), intent(out) :: array(:)
      real(wp), optional, intent(in) :: mu
      real(wp), optional, intent(in) :: sigma
      integer, optional, intent(in) :: seed

      integer :: seed_get
      integer :: n

      n = size(array)

      call random_seed(optval(seed, default_seed), seed_get)

      array(:) = rvs_normal(optval(mu, default_mu), optval(sigma, default_sigma), n)
   end subroutine generate_fgn

   ! Produces series by integrating fgn fractionaly
   subroutine generate_fgn_integrate(series, intorder, seed_in)
      real(wp), intent(out) :: series(:)
      real(wp), intent(in) :: intorder
      integer, optional, intent(in) :: seed_in

      integer :: n
      real(wp), allocatable :: fgn(:)

      n = size(series)
      allocate(fgn(n))

      call generate_fgn(fgn, 0.0_wp, 1.0_wp, seed_in)
      call frac_diff_simple(-intorder, fgn, series)

      deallocate(fgn)
   end subroutine generate_fgn_integrate

   ! Produces series by integrating white noise fractionally
   subroutine generate_white_integrate(series, intorder, seed_in)
      real(wp), intent(out) :: series(:)
      real(wp), intent(in) :: intorder
      integer, optional, intent(in) :: seed_in

      integer :: n
      real(wp), allocatable :: noise(:)

      n = size(series)
      allocate(noise(n))

      call generate_white(noise, 0.0_wp, 1.0_wp, seed_in)
      call frac_diff_simple(-intorder, noise, series)

      deallocate(noise)
   end subroutine generate_white_integrate

   ! Produces series with PSD~1/f^a
   subroutine generate_color(series, a, seed_in)
      real(wp), intent(out) :: series(:)
      real(wp), intent(in) :: a
      integer, optional, intent(in) :: seed_in

      real(wp), allocatable :: white(:), f(:)
      complex(wp), allocatable :: white_fft(:)
      real(wp), allocatable :: S(:)
      integer :: n, i

      n = size(series)
      n = up2power(n)

      allocate(white(n), white_fft(n/2), S(n/2), f(n/2))

      call generate_white(white, default_mu, default_sigma, seed_in)

      call rfft1d(white, white_fft)
      S(:) = [0.0_wp, 1.0_wp / [((real(i, wp) / n) ** a, i = 1, n/2 - 1)] / n]
      S(:) = sqrt(S)
      S(:) = S / sqrt(mean(S**2))
      white_fft(:) = white_fft * S
      call irfft1d(white, white_fft)
      series = white(1:size(series))

      deallocate(white, white_fft, S, f)
   end subroutine generate_color
end module generators
