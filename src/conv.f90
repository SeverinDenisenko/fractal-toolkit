module conv
   use stdlib_error, only: check
   use precision, only: wp
   use intergers, only: up2power
   use fourier, only: fft1d, ifft1d
   implicit none

contains
   subroutine conv1d_full(y, x, b)
      real(wp), intent(in) :: x(:)
      real(wp), intent(in) :: b(:)
      real(wp), intent(out) :: y(:)

      integer :: nx, nb, ny, n
      complex(wp), allocatable :: ax(:), ab(:)

      nx = size(x)
      nb = size(b)
      ny = nx + nb - 1

      call check(ny == size(y), "conv1d_full: size mismatch")

      n = up2power(ny)

      allocate(ax(n), ab(n))

      ax = (0.0_wp, 0.0_wp)
      ab = (0.0_wp, 0.0_wp)

      ax(1:nx) = cmplx(x, 0.0_wp, kind=wp)
      ab(1:nb) = cmplx(b, 0.0_wp, kind=wp)

      call fft1d(ax)
      call fft1d(ab)

      ax = ax * ab

      call ifft1d(ax)

      y(1:ny) = real(ax(1:ny), kind=wp) / real(n, wp)

      deallocate(ax, ab)
   end subroutine conv1d_full

   subroutine conv1d_same(y, x, b)
      real(wp), intent(in) :: x(:)
      real(wp), intent(in) :: b(:)
      real(wp), intent(out) :: y(:)

      integer :: nx, nb, ny, n
      complex(wp), allocatable :: ax(:), ab(:)

      nx = size(x)
      nb = size(b)
      ny = nx

      call check(ny == size(y), "conv1d_same: size mismatch")

      n = up2power(nx + nb - 1)

      allocate(ax(n), ab(n))

      ax = (0.0_wp, 0.0_wp)
      ab = (0.0_wp, 0.0_wp)

      ax(1:nx) = cmplx(x, 0.0_wp, kind=wp)
      ab(1:nb) = cmplx(b, 0.0_wp, kind=wp)

      call fft1d(ax)
      call fft1d(ab)

      ax = ax * ab

      call ifft1d(ax)

      y(1:ny) = real(ax(1:ny), kind=wp) / real(n, wp)

      deallocate(ax, ab)
   end subroutine conv1d_same
end module conv
