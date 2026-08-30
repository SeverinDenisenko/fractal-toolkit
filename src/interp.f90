module interp
   use precision, only: wp
   use checks, only: check
   implicit none

   private
   public :: cubic_resample, cubic_solve, cubic_interp

contains
   subroutine cubic_resample(X, Y, X1, Y1, ierr)
      real(wp), intent(in) :: X(:), Y(:), X1(:)
      real(wp), intent(out) :: Y1(:)
      integer, intent(out), optional :: ierr

      real(wp), allocatable :: b(:), c(:), d(:)
      integer :: n

      n = size(X)
      if(check(n > 3, msg="cubic_resample: n must be bigger then 3", ierr=ierr)) return
      if(check(n == size(Y), msg="cubic_resample: size mismatch", ierr=ierr)) return
      if(check(size(X1) == size(Y1), msg="cubic_resample: size mismatch", ierr=ierr)) return
      allocate(b(n), c(n), d(n))

      call cubic_solve(X, Y, b, c, d, ierr)

      call cubic_interp(X, Y, b, c, d, X1, Y1, ierr)

      deallocate(b, c, d)
   end subroutine cubic_resample

   subroutine cubic_solve(X, Y, b, c, d, ierr)
      real(wp), intent(in) :: X(:), Y(:)
      real(wp), intent(out) :: b(:), c(:), d(:)
      integer, intent(out), optional :: ierr

      integer :: n, i
      real(wp), allocatable :: main_diag(:), sub_diag(:), super_diag(:), rhs(:), h(:)

      n = size(X)
      if(check(n > 3, msg="cubic_resample: n must be bigger then 3", ierr=ierr)) return
      if(check(n == size(Y), msg="cubic_resample: size mismatch", ierr=ierr)) return
      if(check(n == size(b), msg="cubic_resample: size mismatch", ierr=ierr)) return
      if(check(n == size(c), msg="cubic_resample: size mismatch", ierr=ierr)) return
      if(check(n == size(d), msg="cubic_resample: size mismatch", ierr=ierr)) return

      allocate(main_diag(n), sub_diag(n-1), super_diag(n-1), rhs(n), h(n-1))

      ! Compute and store the spacings
      do i = 1, n-1
         h(i) = X(i+1) - X(i)
      end do

      call construct_tridiagonal(X, Y, h, main_diag, sub_diag, super_diag, rhs)
      call solve_tridiagonal(main_diag, sub_diag, super_diag, rhs, c)
      call compute_coefficients(X, Y, h, c, b, d)

      deallocate(main_diag, sub_diag, super_diag, rhs, h)
   end subroutine cubic_solve

   subroutine construct_tridiagonal(X, Y, h, main_diag, sub_diag, super_diag, rhs)
      real(wp), intent(in) :: X(:), Y(:), h(:)
      real(wp), intent(out) :: main_diag(:), sub_diag(:), super_diag(:), rhs(:)

      integer :: n, i

      n = size(X)

      sub_diag = h
      super_diag = h

      ! Build the right-hand side and main diagonal
      rhs(1) = 0.0_wp
      rhs(2) = (Y(2) - Y(1)) / h(1)

      do i = 2, n-1
         main_diag(i) = 2.0_wp * (h(i-1) + h(i))
         rhs(i+1) = (Y(i+1) - Y(i)) / h(i)
         rhs(i) = rhs(i+1) - rhs(i)
      end do

      ! Boundary conditions
      main_diag(1) = -h(1)
      main_diag(n) = -h(n-1)
      rhs(n) = 0.0_wp

      ! Natural spline boundary conditions
      rhs(1) = rhs(3) / (X(4) - X(2)) - rhs(2) / (X(3) - X(1))
      rhs(n) = rhs(n-1) / (X(n) - X(n-2)) - rhs(n-2) / (X(n-1) - X(n-3))
      rhs(1) = rhs(1) * h(1)**2 / (X(4) - X(1))
      rhs(n) = -rhs(n) * h(n-1)**2 / (X(n) - X(n-3))
   end subroutine construct_tridiagonal

   subroutine solve_tridiagonal(main_diag, sub_diag, super_diag, rhs, solution)
      real(wp), intent(inout) :: main_diag(:)
      real(wp), intent(in) :: sub_diag(:), super_diag(:)
      real(wp), intent(inout) :: rhs(:)
      real(wp), intent(out) :: solution(:)

      integer :: n, i, j
      real(wp) :: factor

      n = size(main_diag)

      solution = rhs

      ! Forward elimination (Thomas algorithm)
      do i = 2, n
         factor = sub_diag(i-1) / main_diag(i-1)
         main_diag(i) = main_diag(i) - factor * super_diag(i-1)
         solution(i) = solution(i) - factor * solution(i-1)
      end do

      ! Back substitution
      solution(n) = solution(n) / main_diag(n)
      do j = 1, n-1
         i = n - j
         solution(i) = (solution(i) - super_diag(i) * solution(i+1)) / main_diag(i)
      end do
   end subroutine solve_tridiagonal

   subroutine compute_coefficients(X, Y, h, c, b, d)
      real(wp), intent(in) :: X(:), Y(:), h(:)
      real(wp), intent(inout) :: c(:)
      real(wp), intent(out) :: b(:), d(:)

      integer :: n, i

      n = size(X)

      b(n) = (Y(n) - Y(n-1)) / h(n-1) + h(n-1) * (c(n-1) + 2.0_wp * c(n))
      do i = 1, n-1
         b(i) = (Y(i+1) - Y(i)) / h(i) - h(i) * (c(i+1) + 2.0_wp * c(i))
         d(i) = (c(i+1) - c(i)) / h(i)
         c(i) = 3.0_wp * c(i)
      end do

      c(n) = 3.0_wp * c(n)
      d(n) = d(n-1)
   end subroutine compute_coefficients

   subroutine cubic_interp(X, Y, b, c, d, X1, Y1, ierr)
      real(wp), intent(in) :: X(:), Y(:), b(:), c(:), d(:), X1(:)
      real(wp), intent(out) :: Y1(:)
      integer, intent(out), optional :: ierr

      integer :: j, k, n, n1
      real(wp), allocatable :: uu(:)
      real(wp) :: dx

      n = size(X)
      n1 = size(X1)
      if(check(n == size(Y), msg="cubic_interp: size mismatch", ierr=ierr)) return
      if(check(n1 == size(Y1), msg="cubic_interp: size mismatch", ierr=ierr)) return
      if(check(n == size(Y), msg="cubic_interp: size mismatch", ierr=ierr)) return
      if(check(n == size(b), msg="cubic_interp: size mismatch", ierr=ierr)) return
      if(check(n == size(c), msg="cubic_interp: size mismatch", ierr=ierr)) return
      if(check(n == size(d), msg="cubic_interp: size mismatch", ierr=ierr)) return
      allocate(uu(n))

      do j = 1,n1
         uu = x - X1(j)
         where(uu .gt. 0)
            uu = 0
         else where
            uu = 1
         end where
         k = int(sum(uu))

         if (k .eq. 0) then
            Y1(j) = y(1)
         elseif (k .eq. n) then
            Y1(j) = y(n)
         else
            dx = X1(j) - x(k)
            Y1(j) = y(k) + dx * (b(k) + dx * (c(k) + dx * d(k)))
         end if
      end do
   end subroutine cubic_interp
end module interp
