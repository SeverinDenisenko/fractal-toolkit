module constants
   use precision, only: wp
   use stdlib_constants, only: PI_dp
   implicit none

   real(wp), parameter :: pi = PI_dp
   real(wp), parameter :: twopi = PI_dp * 2.0_wp
end module constants
