module stat
   use precision, only: wp
   use stdlib_error, only: check
   use stdlib_sorting, only: sort
   implicit none

contains
   subroutine percentile(series, q, pct)
      real(wp), intent(in) :: series(:)
      real(wp), intent(in) :: q(:)
      real(wp), intent(out) :: pct(:)

      real(wp), allocatable :: sorted(:)
      integer :: n, i
      real(8) :: rank, fraction
      integer :: lower_idx, upper_idx

      n = size(series)
      allocate(sorted(n))
      sorted(:) = series
      call sort(sorted)

      do i=1,size(q)
         rank = (n - 1) * q(i) / 100.0_wp
         lower_idx = floor(rank) + 1
         upper_idx = min(lower_idx + 1, n)
         fraction = rank - floor(rank)

         if (lower_idx < n) then
            pct(i) = series(lower_idx) + fraction * (series(upper_idx) - series(lower_idx))
         else
            pct(i) = series(n)
         end if
      end do

      call check(size(q) == size(pct), msg='percentile: size mismatch')

      deallocate(sorted)
   end subroutine percentile

   ! Friedman-Diaconis method for estimating optimal bins count for historgam
   integer function fd_bins(series) result(bins)
      real(wp), intent(in) :: series(:)
      real(wp) :: pct(2), iqr, range, h
      
      call percentile(series, [25.0_wp, 75.0_wp], pct)

      iqr = pct(2) - pct(1)
      h = 2.0_wp * iqr / (size(series) ** (1.0_wp / 3.0_wp))
      range = maxval(series) - minval(series)

      bins = ceiling(range / h)
   end function fd_bins
end module stat
