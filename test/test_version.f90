program version_test
   use version, only: max_ver_len, ver
   use stdlib_error, only: check
   implicit none

   character(len=max_ver_len) :: v

   call ver(v)
   call check(v == '1.0')
end program version_test
