program main
   use stdlib_error, only: check
   use gen_series, only: gen_series_fgn_integrate, gen_series_white, gen_series_fgn, gen_series_white_integrate
   use spectra, only: berg_psd, berg_psd_size
   use fit, only: powerregress
   use precision, only: wp
   implicit none

   integer :: n, m, t
   real(wp), allocatable :: data(:)
   real(wp), allocatable :: P(:)
   real(wp), allocatable :: f(:)
   real(wp) :: a, c, sigma2, H

   n = 10000
   m = 200
   t = 10
   allocate(data(n), P(berg_psd_size(n)), f(berg_psd_size(n)))

   call gen_series_white_integrate(data, 0.500_wp)
   call berg_psd(f, P, data, 1.0_wp, m)

   call powerregress(f(t:), P(t:), a, c, sigma2)
   H = (-a + 1) / 2

   write (*,*) -a, H, sigma2
end program main
