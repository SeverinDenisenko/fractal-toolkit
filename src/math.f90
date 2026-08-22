module math
   use precision, only: wp
   implicit none

contains
   real(wp) function log2(a) result(x)
      real(wp), intent(in) :: a

      x = logn(a, 2.0_wp)
   end function log2

   real(wp) function logn(a, n) result(x)
      real(wp), intent(in) :: a, n

      x = log(a) / log(n)
   end function logn
end module math
