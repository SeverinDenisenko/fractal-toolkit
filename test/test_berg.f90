program berg_test
   use precision, only: wp
   use spectra, only: berg_psd, psd_size
   use stdlib_error, only: check
   implicit none

   real(wp), allocatable :: P(:)
   real(wp), allocatable :: f(:)
   integer :: n
   real(wp) :: data2(8)

   data2(:) = [1, -1, 1, -1, 1, -1, 1, -1]
   n = psd_size(size(data2))
   call check(n == 5)
   allocate(P(n))
   allocate(f(n))
   call berg_psd(f, P, data2, 1.0_wp, 4)
   call check(maxval(abs(f - [0.000_wp, 0.125_wp, 0.250_wp, 0.375_wp, 0.500_wp])) < 1e-5_wp)
   call check(maxval(abs(P - [0.000_wp, 0.000_wp, 0.000_wp, 0.000_wp, 1.000_wp])) < 1e-5_wp)

   deallocate(P, f)
end program berg_test
