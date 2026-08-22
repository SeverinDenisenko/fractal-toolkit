module conv
   use precision, only: wp
   use integers, only: up2power
   use fourier, only: fft1d, ifft1d
   use check_mod, only: check
   implicit none

 contains
   subroutine conv1d_full(y, x, b, ierr)
      real(wp), intent(in) :: x(:)
      real(wp), intent(in) :: b(:)
      real(wp), intent(out) :: y(:)
      integer, intent(out), optional :: ierr

      integer :: nx, nb, ny, n
      complex(wp), allocatable :: ax(:), ab(:)

      nx = size(x)
      nb = size(b)
      ny = nx + nb - 1

      if (check(ny == size(y), "conv1d_full: size mismatch", ierr=ierr)) return

      n = up2power(ny, ierr=ierr)
      if (present(ierr) .and. ierr /= 0) return

      allocate(ax(n), ab(n))

      ax(:) = (0.0_wp, 0.0_wp)
      ab(:) = (0.0_wp, 0.0_wp)

      ax(1:nx) = cmplx(x, 0.0_wp, kind=wp)
      ab(1:nb) = cmplx(b, 0.0_wp, kind=wp)

      call fft1d(ax, ierr=ierr)
      if (present(ierr) .and. ierr /= 0) return

      call fft1d(ab, ierr=ierr)
      if (present(ierr) .and. ierr /= 0) return

      ax(:) = ax * ab

      call ifft1d(ax, ierr=ierr)
      if (present(ierr) .and. ierr /= 0) return

      y(1:ny) = real(ax(1:ny), kind=wp) / real(n, wp)

      deallocate(ax, ab)
   end subroutine conv1d_full

   subroutine conv1d_same(y, x, b, ierr)
      real(wp), intent(in) :: x(:)
      real(wp), intent(in) :: b(:)
      real(wp), intent(out) :: y(:)
      integer, intent(out), optional :: ierr

      integer :: nx, nb, ny, n
      complex(wp), allocatable :: ax(:), ab(:)

      nx = size(x)
      nb = size(b)
      ny = nx

      if (check(ny == size(y), "conv1d_same: size mismatch", ierr=ierr)) return

      n = up2power(nx + nb - 1, ierr=ierr)
      if (present(ierr) .and. ierr /= 0) return

      allocate(ax(n), ab(n))

      ax(:) = (0.0_wp, 0.0_wp)
      ab(:) = (0.0_wp, 0.0_wp)

      ax(1:nx) = cmplx(x, 0.0_wp, kind=wp)
      ab(1:nb) = cmplx(b, 0.0_wp, kind=wp)

      call fft1d(ax, ierr=ierr)
      if (present(ierr) .and. ierr /= 0) return
      call fft1d(ab, ierr=ierr)
      if (present(ierr) .and. ierr /= 0) return

      ax(:) = ax * ab

      call ifft1d(ax, ierr=ierr)
      if (present(ierr) .and. ierr /= 0) return

      y(1:ny) = real(ax(1:ny), kind=wp) / real(n, wp)

      deallocate(ax, ab)
   end subroutine conv1d_same
end module conv
