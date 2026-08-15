module fit
   use precision, only: wp
   use stdlib_linalg, only: solve_lstsq
   use stdlib_error, only: check
   use stdlib_stats, only: var
   implicit none

contains
   ! Solves y = x * k + a
   subroutine linregress(x, y, k, a, sigma2)
      real(wp), intent(in) :: x(:)
      real(wp), intent(in) :: y(:)
      real(wp), intent(out) :: k
      real(wp), intent(out) :: a
      real(wp), optional, intent(out) :: sigma2

      real(wp), allocatable :: M(:,:)
      real(wp) :: solution(2)
      integer :: n

      call check(size(x) == size(y), msg="linregress: size missmatch")
      n = size(x)
      allocate(M(n, 2))

      M(:, 1) = x(:)
      M(:, 2) = 1.0_wp

      call solve_lstsq(M, y, x=solution)

      k = solution(1)
      a = solution(2)

      if(present(sigma2)) then
         sigma2 = sum(((k * x + a) - y)**2) / (n - 2)
      end if
   end subroutine linregress

   ! Solves y = c * x^a
   subroutine powerregress(x, y, a, c, sigma2)
      real(wp), intent(in) :: x(:)
      real(wp), intent(in) :: y(:)
      real(wp), intent(out) :: a
      real(wp), intent(out) :: c
      real(wp), optional, intent(out) :: sigma2

      real(wp), allocatable :: M(:,:)
      real(wp) :: solution(2)
      real(wp), allocatable :: log_x(:), log_y(:)
      integer :: n

      call check(size(x) == size(y), msg="powerregress: size mismatch")
      n = size(x)

      call check(all(x > 0.0_wp), msg="powerregress: x must be positive")
      call check(all(y > 0.0_wp), msg="powerregress: y must be positive")

      allocate(M(n, 2))
      allocate(log_x(n), log_y(n))

      log_x = log(x)
      log_y = log(y)

      M(:, 1) = log_x(:)
      M(:, 2) = 1.0_wp

      call solve_lstsq(M, log_y, x=solution)

      a = solution(1)
      c = exp(solution(2))

      if(present(sigma2)) then
         sigma2 = sum((a * log_x + solution(2) - log_y)**2) / (n - 2)
      end if
   end subroutine powerregress
end module fit
