module complex
   use precision, only: wp
   use constants, only: twopi
   implicit none

contains
   function zroots(n,nn)
      integer, intent(in) :: n,nn
      complex(wp) :: zroots(nn)

      integer :: k
      real(wp) :: theta

      zroots(1) = 1.0
      theta = twopi / n

      k = 1
      do
         if (k >= nn) then
            exit
         end if
         zroots(k+1) = cmplx(cos(k * theta), sin(k * theta), wp)
         zroots(k+2:min(2*k,nn)) = zroots(k + 1) * zroots(2:min(k, nn-k))
         k = 2 * k
      end do
   end function zroots
end module complex
