module solvers
   use precision, only: wp
   use stdlib_linalg, only: solve_lstsq, solve, inv
   use check_mod, only: check
   implicit none

 contains
   ! Compute variation for residuals with `m` freedom degrees
   real(wp) function residuals_var(x, y, m, ierr) result(sigma2)
      real(wp), intent(in) :: x(:), y(:)
      integer, intent(in) :: m
      integer, intent(out), optional :: ierr

      integer :: n

      sigma2 = 0.0_wp
      if (check(size(x) == size(y), msg="residuals_var: size mismatch", ierr=ierr)) return
      n = size(x)

      sigma2 = sum((x - y)**2)/(n - m)
   end function residuals_var

   ! Solves y = x * k + a
   subroutine linregress(x, y, k, a, sigma2, k_err, a_err, ierr)
      real(wp), intent(in) :: x(:), y(:)
      real(wp), intent(out) :: k, a
      real(wp), optional, intent(out) :: sigma2, k_err, a_err
      integer, intent(out), optional :: ierr

      real(wp), allocatable :: M(:,:)
      real(wp) :: solution(2), cov_matrix(2,2)
      integer :: n

      k = 0.0_wp
      a = 0.0_wp
      if (check(size(x) == size(y), msg="linregress: size mismatch", ierr=ierr)) return

      n = size(x)
      allocate(M(n, 2))

      M(:, 1) = x(:)
      M(:, 2) = 1.0_wp

      call solve_lstsq(M, y, x=solution)

      k = solution(1)
      a = solution(2)

      if (present(sigma2)) then
         sigma2 = residuals_var(k * x + a, y, 2)

         cov_matrix(:,:) = sigma2 * inv(matmul(transpose(M), M))

         if (present(k_err)) then
            k_err = sqrt(cov_matrix(1,1))
         end if

          if (present(a_err)) then
             a_err = sqrt(cov_matrix(2,2))
          end if
      end if

      deallocate(M)
   end subroutine linregress

   ! Solves y = c * x^a
   subroutine powerregress(x, y, a, c, sigma2, a_err, c_err, ierr)
      real(wp), intent(in) :: x(:)
      real(wp), intent(in) :: y(:)
      real(wp), intent(out) :: a, c
      real(wp), optional, intent(out) :: sigma2, a_err, c_err
      integer, intent(out), optional :: ierr

      real(wp), allocatable :: M(:,:)
      real(wp) :: solution(2), cov_matrix(2,2)
      integer :: n

      a = 0.0_wp
      c = 0.0_wp
      if (check(size(x) == size(y), msg="powerregress: size mismatch", ierr=ierr)) return
      n = size(x)

      if (check(all(x > 0.0_wp), msg="powerregress: x must be positive", ierr=ierr)) return
      if (check(all(y > 0.0_wp), msg="powerregress: y must be positive", ierr=ierr)) return

      allocate(M(n, 2))

      M(:, 1) = log(x)
      M(:, 2) = 1.0_wp

      call solve_lstsq(M, log(y), x=solution)

      a = solution(1)
      c = exp(solution(2))

      if (present(sigma2)) then
         sigma2 = residuals_var(a * log(x) + log(c), log(y), 2)

         cov_matrix(:,:) = sigma2 * inv(matmul(transpose(M), M))

         if (present(a_err)) then
            a_err = sqrt(cov_matrix(1,1))
         end if

          if (present(c_err)) then
             c_err = sqrt(cov_matrix(2,2))
          end if
      end if

      deallocate(M)
   end subroutine powerregress
end module solvers
