module integers
   use checks, only: check
   implicit none

 contains
   integer function up2power(n, ierr) result(x)
      integer, intent(in) :: n
      integer, intent(out), optional :: ierr

      x = 0
      if (check(n > 0, msg="up2power: n must be positive", ierr=ierr)) return
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
