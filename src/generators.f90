module generators
   use precision, only: wp
   use stdlib_random, only: random_seed
   use stdlib_stats_distribution_normal, only: rvs_normal
   use stdlib_stats_distribution_uniform, only: rvs_uniform
   use stdlib_optval, only: optval
   use frac, only: frac_diff_simple
   implicit none

   integer :: default_seed = 42
   real(wp) :: default_sigma = 1.0_wp
   real(wp) :: default_mu = 0.0_wp
contains
   subroutine generte_white(array, mu, sigma, seed)
      real(wp), intent(out) :: array(:)
      real(wp), optional, intent(in) :: mu
      real(wp), optional, intent(in) :: sigma
      integer, optional, intent(in) :: seed

      integer :: seed_get
      integer :: n

      n = size(array)
      call random_seed(optval(seed, default_seed), seed_get)

      array(:) = rvs_uniform(optval(mu, default_mu) - optval(sigma, default_sigma) / 2.0_wp, optval(mu, default_mu) + optval(sigma, default_sigma) / 2.0_wp, n)
   end subroutine generte_white

   subroutine generte_fgn(array, mu, sigma, seed)
      real(wp), intent(out) :: array(:)
      real(wp), optional, intent(in) :: mu
      real(wp), optional, intent(in) :: sigma
      integer, optional, intent(in) :: seed

      integer :: seed_get
      integer :: n

      n = size(array)

      call random_seed(optval(seed, default_seed), seed_get)

      array(:) = rvs_normal(optval(mu, default_mu), optval(sigma, default_sigma), n)
   end subroutine generte_fgn

   ! Produces series by integrating fgn fractionaly
   subroutine generte_fgn_integrate(series, intorder, seed_in)
      real(wp), intent(out) :: series(:)
      real(wp), intent(in) :: intorder
      integer, optional, intent(in) :: seed_in

      integer :: n
      real(wp), allocatable :: fgn(:)

      n = size(series)
      allocate(fgn(n))

      call generte_fgn(fgn, 0.0_wp, 1.0_wp, seed_in)
      call frac_diff_simple(-intorder, fgn, series)

      deallocate(fgn)
   end subroutine generte_fgn_integrate

   ! Produces series by integrating white noise fractionally
   subroutine generte_white_integrate(series, intorder, seed_in)
      real(wp), intent(out) :: series(:)
      real(wp), intent(in) :: intorder
      integer, optional, intent(in) :: seed_in

      integer :: n
      real(wp), allocatable :: noise(:)

      n = size(series)
      allocate(noise(n))

      call generte_white(noise, 0.0_wp, 1.0_wp, seed_in)
      call frac_diff_simple(-intorder, noise, series)

      deallocate(noise)
   end subroutine generte_white_integrate
end module generators
