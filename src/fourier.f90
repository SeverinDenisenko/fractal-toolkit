module fourier
   use precision, only: wp
   use integers, only: is2power, up2power
   use constants, only: pi
   use complex, only: zroots
   use check_mod, only: check
   use stdlib_math, only: swap
   use stdlib_optval, only: optval
   implicit none

contains
   subroutine fft2d(data, ierr)
      complex(wp), intent(inout) :: data(:,:)
      integer, intent(out), optional :: ierr

      call sfft2d(data, 1, ierr=ierr)
   end subroutine fft2d

   subroutine ifft2d(data, ierr)
      complex(wp), intent(inout) :: data(:,:)
      integer, intent(out), optional :: ierr

      call sfft2d(data, -1, ierr=ierr)
   end subroutine ifft2d

   subroutine fft1d(data, ierr)
      complex(wp), intent(inout) :: data(:)
      integer, intent(out), optional :: ierr

      call sfft1d(data, 1, ierr=ierr)
   end subroutine fft1d

   subroutine ifft1d(data, ierr)
      complex(wp), intent(inout) :: data(:)
      integer, intent(out), optional :: ierr

      call sfft1d(data, -1, ierr=ierr)
   end subroutine ifft1d

   subroutine rfft1d(data, zdata, ierr)
      real(wp), intent(inout) :: data(:)
      complex(wp), intent(inout) :: zdata(:)
      integer, intent(out), optional :: ierr

      call srfft1d(data, 1, zdata, ierr=ierr)
   end subroutine rfft1d

   subroutine irfft1d(data, zdata, ierr)
      real(wp), intent(inout) :: data(:)
      complex(wp), intent(inout) :: zdata(:)
      integer, intent(out), optional :: ierr

      call srfft1d(data, -1, zdata, ierr=ierr)
   end subroutine irfft1d

   subroutine sfft1d(data, isign, ierr)
      complex(wp), intent(inout) :: data(:)
      integer, intent(in) :: isign
      integer, intent(out), optional :: ierr
      complex(wp), allocatable :: tmp(:,:)
      integer :: n
      n = size(data)
      if (check(is2power(n), 'sfft1d: n must be a power of 2', ierr=ierr)) return

      allocate(tmp(1, size(data)))
      tmp(1, 1:size(data)) = data

      call sfft2d(tmp, isign, ierr=ierr)
      if (present(ierr) .and. ierr /= 0) return

      data = tmp(1, 1:size(data))
      deallocate(tmp)
   end subroutine sfft1d

   subroutine sfft2d(data, isign, ierr)
      complex(wp), intent(inout) :: data(:,:)
      integer, intent(in) :: isign
      integer, intent(out), optional :: ierr

      integer :: n, i, istep, j, m, mmax, n2
      real(wp) :: theta
      complex(wp) :: temp(size(data, 1))
      complex(wp) :: w, ws, wwp

      n = size(data,2)

      if (check(is2power(n), 'fft: n must be a power of 2', ierr=ierr)) return

      n2 = n / 2
      j = n2
      do i = 1, n - 2
         if (j > i) then
            call swap(data(:,j+1), data(:,i+1))
         end if
         m = n2
         do
            if (m < 2 .or. j < m) exit
            j = j - m
            m = m / 2
         end do
         j = j + m
      end do
      mmax = 1
      do
         if (n <= mmax) exit
         istep = 2 * mmax
         theta = pi / (isign * mmax)
         wwp = cmplx(-2.0_wp * sin(0.5_wp * theta)**2, sin(theta), kind=wp)
         w = cmplx(1.0_wp, 0.0_wp, kind=wp)
         do m = 1, mmax
            ws = w
            do i = m, n, istep
               j = i + mmax
               temp(:) = ws * data(:,j)
               data(:,j) = data(:,i) - temp
               data(:,i) = data(:,i) + temp
            end do
            w = w * wwp + w
         end do
         mmax = istep
      end do
   end subroutine sfft2d

   subroutine srfft1d(data, isign, zdata, ierr)
      real(wp), intent(inout) :: data(:)
      integer, intent(in) :: isign
      complex(wp), optional, target :: zdata(:)
      integer, intent(out), optional :: ierr
   
      integer :: n, nh, nq
      complex(wp) :: w(size(data)/4)
      complex(wp), dimension(size(data)/4-1) :: h1, h2
      complex(wp), pointer :: cdata(:)
      complex(wp) :: z
      real(wp) :: c1 = 0.5_wp, c2
    
      n = size(data)
      if (check(is2power(n), 'srfft: n must be a power of 2', ierr=ierr)) return
      nh = n / 2
      nq = n / 4
    
      if (present(zdata)) then
         if (check(n / 2 == size(zdata), msg='srfft', ierr=ierr)) return
         cdata => zdata
         if (isign == 1) then
            cdata = cmplx(data(1:n-1:2), data(2:n:2), kind=wp)
         end if
      else
         allocate(cdata(n / 2))
         cdata = cmplx(data(1:n-1:2), data(2:n:2), kind=wp)
      end if
   
      if (isign == 1) then
         c2 = -0.5_wp
         call sfft1d(cdata, +1)
      else
         c2 = 0.5_wp
      end if

      w(:) = zroots(sign(n, isign), n/4)
      w(:) = cmplx(-aimag(w), real(w), kind=wp)
      h1(:) = c1 * (cdata(2:nq) + conjg(cdata(nh:nq+2:-1)))
      h2(:) = c2 * (cdata(2:nq) - conjg(cdata(nh:nq+2:-1)))
      cdata(2:nq) = h1 + w(2:nq) * h2
      cdata(nh:nq+2:-1) = conjg(h1 - w(2:nq) * h2)
      z = cdata(1)

      if (isign == 1) then
         cdata(1) = cmplx(real(z) + aimag(z), real(z) - aimag(z), kind=wp)
      else
         cdata(1) = cmplx(c1 * (real(z) + aimag(z)), c1 * (real(z) - aimag(z)), kind=wp)
         call sfft1d(cdata, -1)
      end if

      if (present(zdata)) then
         if (isign /= 1) then
            data(1:n-1:2) = real(cdata)
            data(2:n:2) = aimag(cdata)
         end if
      else
         data(1:n-1:2) = real(cdata)
         data(2:n:2) = aimag(cdata)
         deallocate(cdata)
      end if
   end subroutine srfft1d
end module fourier
