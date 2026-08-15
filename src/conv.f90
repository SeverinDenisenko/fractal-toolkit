module conv
   use stdlib_error, only: check
   use precision, only: wp
   implicit none

contains
   subroutine conv1d_full(y, x, b)
      real(wp), intent(in) :: x(:)
      real(wp), intent(in) :: b(:)
      real(wp), intent(out) :: y(:)

      integer :: ib, ix, iy, nx, nb, ny

      y(:) = 0

      nx = size(x)
      nb = size(b)
      ny = nx + nb - 1

      call check(ny == size(y), "conv1d_full: size missmatch")

      do ix=1, nx
         do ib=1, nb
            iy = ix + ib - 1
            y(iy) = y(iy) + x(ix) * b(ib)
         end do
      end do
   end subroutine conv1d_full

   subroutine conv1d_same(y, x, b)
      real(wp), intent(in) :: x(:)
      real(wp), intent(in) :: b(:)
      real(wp), intent(out) :: y(:)

      integer :: ib, ix, iy, nx, nb, ny

      y(:) = 0

      nx = size(x)
      nb = size(b)
      ny = nx

      call check(ny == size(y), "conv1d_same: size missmatch")

      do ix=1, nx
         do ib=1, nb
            iy = ix + ib - 1
            if (iy > ny) then
               exit
            end if
            y(iy) = y(iy) + x(ix) * b(ib)
         end do
      end do
   end subroutine conv1d_same
end module conv
