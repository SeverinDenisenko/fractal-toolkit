program version_test
   use version, only: ver
   use stdlib_error, only: check
   implicit none

   character(len=16) :: v

   call ver(v)
   call check(v == '1.0')
end program version_test
