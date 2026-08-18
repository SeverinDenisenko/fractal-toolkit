program main_berg
   use stdlib_error, only: check
   use generators, only: generte_white_integrate, generte_fgn, generte_fgn_integrate, generte_white
   use spectra, only: berg_psd, berg_psd_size
   use solvers, only: powerregress
   use io, only: write_table_file
   use hurst, only: slope_to_hurst
   use precision, only: wp
   implicit none

   integer :: n, m, t, n1
   real(wp), allocatable :: data(:)
   real(wp), allocatable :: P(:)
   real(wp), allocatable :: f(:)
   real(wp) :: a, c, sigma2, H

   n = 10000
   m = n / 20
   n1 = berg_psd_size(n)
   t = n / 100
   allocate(data(n), P(n1), f(n1))

   call generte_fgn_integrate(data, 0.3_wp)
   call berg_psd(f, P, data, 1.0_wp, m)

   call powerregress(f(t:n1-t), P(t:n1-t), a, c, sigma2)

   H = slope_to_hurst(-a)
   write (*,*) H, sigma2

   call write_table_file('output.dat', f, P)
end program main_berg
