program main
   use stdlib_error, only: check
   use gen_series, only: gen_series_white_integrate, gen_series_fgn, gen_series_fgn_integrate, gen_series_white
   use spectra, only: berg_psd, berg_psd_size
   use fit, only: powerregress
   use precision, only: wp
   implicit none

   integer :: n, m, t, i, n1
   real(wp), allocatable :: data(:)
   real(wp), allocatable :: P(:)
   real(wp), allocatable :: f(:)
   real(wp) :: a, c, sigma2, H

   n = 10000
   m = n / 20
   n1 = berg_psd_size(n)
   t = n / 100
   allocate(data(n), P(n1), f(n1))

   call gen_series_fgn_integrate(data, 0.5_wp)
   call berg_psd(f, P, data, 1.0_wp, m)

   call powerregress(f(t:n1-t), P(t:n1-t), a, c, sigma2)

   a = -a
   if (a < 2) then
      H = (a + 1) / 2
   else
      H = (a - 1) / 2
   end if

   write (*,*) H, sigma2

   open(unit=10, file='output.dat', status='replace')
   do i = t, size(f)-t
      write(10, '(ES20.12, 2X, ES20.12)') f(i), P(i)
   end do
   close(10)
end program main
