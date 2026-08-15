module frac
   use stdlib_error, only: check
   use precision, only: wp
   use conv, only: conv1d_same
   implicit none

contains
   ! Calculate Maclaurin Series coefficients for fractional difference
   ! Will fill `coeff` array with coefficients
   subroutine frac_diff_coeff(coeff, order)
      real(wp), intent(in) :: order
      real(wp), intent(out) :: coeff(:)

      integer :: len
      integer :: i

      len = size(coeff)

      coeff(1) = 1.0_wp
      do i = 1, len - 1
         coeff(i + 1) = -coeff(i) / i * (order - i + 1)
      end do
   end subroutine frac_diff_coeff

   ! Calculate fractional difference
   subroutine frac_diff(coeff, series_in, series_out)
      real(wp), intent(in) :: coeff(:)
      real(wp), intent(in) :: series_in(:)
      real(wp), intent(out) :: series_out(:)

      call check(size(series_in) == size(series_out), msg = "frac_diff: series must be equal length")

      call conv1d_same(series_out, series_in, coeff)
   end subroutine frac_diff

   ! Calculate fractional difference
   subroutine frac_diff_simple(order, series_in, series_out)
      real(wp), intent(in) :: order
      real(wp), intent(in) :: series_in(:)
      real(wp), intent(out) :: series_out(:)

      integer :: len
      real(wp), allocatable :: coeff(:)

      len = size(series_in)
      allocate(coeff(len))

      call frac_diff_coeff(coeff, order)
      call frac_diff(coeff, series_in, series_out)

      deallocate(coeff)
   end subroutine frac_diff_simple
end module frac
