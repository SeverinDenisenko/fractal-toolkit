module integers
   use stdlib_error, only: check
   implicit none

contains
   integer function up2power(n) result(x)
      integer, intent(in) :: n

      call check(n > 0, msg="up2power: n must be positive")
      x = n - 1
      x = ior(x, ishft(x, -1))
      x = ior(x, ishft(x, -2))
      x = ior(x, ishft(x, -4))
      x = ior(x, ishft(x, -8))
      x = ior(x, ishft(x, -16))
      x = x + 1
   end function up2power

   logical function is2power(n) result(res)
      integer, intent(in) :: n

      res = (n /= 0) .and. (iand(n, n - 1) == 0)
   end function is2power
end module integers
