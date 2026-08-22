module autoreg
   use precision, only: wp
   use constants, only: pi
   use stdlib_linalg, only: solve_lstsq
   use check_mod, only: check
   implicit none

contains
   ! Calculate AR coefficients `phi` for time series `S` by solving Yule-Walker equations
   subroutine yw_ar_coeff(phi, S, ierr)
      real(wp), intent(out) :: phi(:)
      real(wp), intent(in) :: S(:)
      integer, intent(out), optional :: ierr

      integer :: n, p, i, j
      real(wp), allocatable :: X(:,:)

      p = size(phi)
      n = size(S) - p

      if (check(size(S(p+1:)) == n, msg="yw_ar_coeff: size mismatch", ierr=ierr)) return

      allocate(X(n, p))

      do i = 1, n
         do j = 1, p
            X(i, j) = S(i + p - j)
         end do
      end do

      call solve_lstsq(X, S(p+1:), x=phi)

      deallocate(X)
   end subroutine yw_ar_coeff

   ! Calculate AR coefficients `phi` for time series `S` by using Burg method
   ! AR coefficients are not negated and one is not present in the front
   subroutine burg_ar_coeff(phi, S)
      real(wp), intent(out) :: phi(:)
      real(wp), intent(in) :: S(:)

      integer :: n, p, k, i
      real(wp), allocatable :: f(:), b(:) ! Forward and backward prediction errors
      real(wp) :: ref_coef, f_old, num, den
      real(wp), allocatable :: a(:)

      p = size(phi)
      n = size(S)

      allocate(f(n), b(n))
      f = S
      b = S

      allocate(a(0:p))
      a(1:p) = 0.0_wp
      a(0) = 1.0_wp

      do k = 1, p
         ! Reflection coefficient
         num = sum(f(k+1:n)*b(k:n-1))
         den = sum(f(k+1:n)**2 + b(k:n-1)**2)
         if (den < tiny(den)) then
            ref_coef = 0.0_wp
         else
            ref_coef = -2.0_wp * num / den
         end if

         ! Levinson recursion
         do i = 1, k-1
            a(i) = a(i) + ref_coef * a(k-i)
         end do
         a(k) = ref_coef

         ! Update prediction errors
         if (k < p) then
            do i = n, k+1, -1
               f_old = f(i)
               f(i)  = f_old + ref_coef * b(i-1)
               b(i)  = b(i-1) + ref_coef * f_old
            end do
         end if
      end do

      phi = -a(1:p)

      deallocate(f, b, a)
   end subroutine burg_ar_coeff

   ! Compute prediction of timeseries `S` for AR model with coefficients `phi` at index `i` from previous values (not including `i`)
   real(wp) function ar_predict(phi, S, i) result(x)
      real(wp), intent(in) :: phi(:)
      real(wp), intent(in) :: S(:)
      integer, intent(in) :: i

      x = dot_product(phi, S(i - 1:i - size(phi):-1))
   end function ar_predict

   ! Calculate frequency response `H` of an AR filter `phi` on frequencies `f` (0 to 1)
   subroutine ar_freq_response(phi, f, H, ierr)
      real(wp), intent(in) :: phi(:)
      real(wp), intent(in) :: f(:)
      complex(wp), intent(out) :: H(:)
      integer, intent(out), optional :: ierr

      integer :: p, n
      integer :: i, k
      real(wp) :: omega
      complex(wp) :: D

      p = size(phi)
      n = size(f)

      if (check(n == size(H), msg="ar_freq_response: size mismatch", ierr=ierr)) return

      do i = 1, n
         omega = f(i) * 2.0_wp * pi

         D = cmplx(1.0_wp, 0.0_wp, kind=wp)
         do k = 1, p
            D = D - phi(k) * exp(-cmplx(0.0_wp, 1.0_wp, kind=wp) * k * omega)
         end do

         H(i) = cmplx(1.0_wp, 0.0_wp, kind=wp) / D
      end do
   end subroutine ar_freq_response
end module autoreg
