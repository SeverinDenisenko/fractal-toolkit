program spectra_test
   use precision, only: wp
   use constants, only: pi
   use spectra, only: psd_size, berg_psd, yw_psd, complex_berg_psd, complex_yw_psd
   use stdlib_error, only: check
   implicit none

   integer, parameter :: n = 16
   real(wp) :: data2(8)
   real(wp), allocatable :: P(:), Pr(:), f(:)
   complex(wp) :: cdata(n)
   integer :: i, ns, nloc

   data2(:) = [1, -1, 1, -1, 1, -1, 1, -1]
   ns = psd_size(size(data2))
   allocate(P(ns), Pr(ns), f(ns))

   call complex_berg_psd(f, P, cmplx(data2, 0.0_wp, kind=wp), 1.0_wp, 4)
   call check(maxval(abs(f - [0.000_wp, 0.125_wp, 0.250_wp, 0.375_wp, 0.500_wp])) < 1e-5_wp)
   call check(maxval(abs(P - [0.000_wp, 0.000_wp, 0.000_wp, 0.000_wp, 1.000_wp])) < 1e-5_wp)

   call berg_psd(f, Pr, data2, 1.0_wp, 4)
   call check(maxval(abs(P - Pr)) < 1e-12_wp)

   call complex_yw_psd(f, P, cmplx(data2, 0.0_wp, kind=wp), 1.0_wp, 4)
   call yw_psd(f, Pr, data2, 1.0_wp, 4)
   call check(maxval(abs(P - Pr)) < 1e-12_wp)

   deallocate(P, Pr, f)

   do i = 0, n - 1
      cdata(i + 1) = exp(cmplx(0.0_wp, 2.0_wp * pi * 0.2_wp * i, kind=wp))
   end do

   ns = psd_size(n)
   allocate(P(ns), f(ns))

   call complex_yw_psd(f, P, cdata, 1.0_wp, 1)
   nloc = maxloc(P, 1)
   call check(nloc > 1 .and. nloc < ns)
   call check(maxval(P) > 0.5_wp)

   call complex_berg_psd(f, P, cdata, 1.0_wp, 1)
   nloc = maxloc(P, 1)
   call check(nloc > 1 .and. nloc < ns)
   call check(maxval(P) > 0.5_wp)

   deallocate(P, f)

   print '(a)', "spectra_test passed"
end program spectra_test