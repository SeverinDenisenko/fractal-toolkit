module version
   implicit none

   integer, parameter :: max_ver_len = 16
contains
   subroutine ver(v)
      character(len=*) :: v
      v = '1.0'
   end subroutine ver
end module version
