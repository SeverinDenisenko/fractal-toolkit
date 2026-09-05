module autoreg
   use precision, only: wp
   use constants, only: pi
   use stdlib_linalg, only: solve_lstsq
   use checks, only: check
   implicit none

contains
   ! Calculate AR coefficients `phi` for complex time series `S` by solving Yule-Walker equations.
   ! The complex least-squares problem is recast into an equivalent real one
   !   [Xr  -Xi] [pr]   [br]
   !   [Xi   Xr] [pi] = [bi],
   ! which yields the same minimum-norm solution and reduces exactly to the
   ! real Yule-Walker solve for real input.
   subroutine complex_yw_ar_coeff(phi, S, ierr)
      complex(wp), intent(out) :: phi(:)
      complex(wp), intent(in) :: S(:)
      integer, intent(out), optional :: ierr

      integer :: n, p, i, j
      complex(wp), allocatable :: X(:,:)
      real(wp), allocatable :: A(:,:), b(:), z(:)

      p = size(phi)
      n = size(S) - p

      if (check(size(S(p+1:)) == n, msg="yw_ar_coeff: size mismatch", ierr=ierr)) return
      if (check(n > 0, msg="yw_ar_coeff: series too short for the given order", ierr=ierr)) return

      allocate(X(n, p))
      do i = 1, n
         do j = 1, p
            X(i, j) = S(i + p - j)
         end do
      end do

      allocate(A(2*n, 2*p), b(2*n), z(2*p))
      A(1:n, 1:p) = real(X, wp)
      A(1:n, p+1:2*p) = -aimag(X)
      A(n+1:2*n, 1:p) = aimag(X)
      A(n+1:2*n, p+1:2*p) = real(X, wp)

      b(1:n) = real(S(p+1:), wp)
      b(n+1:2*n) = aimag(S(p+1:))

      call solve_lstsq(A, b, x=z)

      phi = cmplx(z(1:p), z(p+1:2*p), kind=wp)

      deallocate(X, A, b, z)
   end subroutine complex_yw_ar_coeff

   ! Calculate AR coefficients `phi` for real time series `S` by solving Yule-Walker equations.
   ! Real data is treated as complex with zero imaginary part.
   subroutine yw_ar_coeff(phi, S, ierr)
      real(wp), intent(out) :: phi(:)
      real(wp), intent(in) :: S(:)
      integer, intent(out), optional :: ierr

      complex(wp), allocatable :: cphi(:), cS(:)

      allocate(cphi(size(phi)), cS(size(S)))

      cS = cmplx(S, 0.0_wp, kind=wp)

      call complex_yw_ar_coeff(cphi, cS, ierr=ierr)
      if (present(ierr) .and. ierr /= 0) return

      phi = real(cphi, wp)

      deallocate(cphi, cS)
   end subroutine yw_ar_coeff

   ! Calculate AR coefficients `phi` for complex time series `S` by using Burg method
   ! AR coefficients are not negated and one is not present in the front
   subroutine complex_burg_ar_coeff(phi, S)
      complex(wp), intent(out) :: phi(:)
      complex(wp), intent(in) :: S(:)

      integer :: n, p, k, i
      complex(wp), allocatable :: f(:), b(:) ! Forward and backward prediction errors
      complex(wp) :: ref_coef, f_old, num
      real(wp) :: den
      complex(wp), allocatable :: a(:)

      p = size(phi)
      n = size(S)

      allocate(f(n), b(n))
      f = S
      b = S

      allocate(a(0:p))
      a(1:p) = (0.0_wp, 0.0_wp)
      a(0) = (1.0_wp, 0.0_wp)

      do k = 1, p
         ! Reflection coefficient
         num = sum(f(k+1:n) * conjg(b(k:n-1)))
         den = sum(abs(f(k+1:n))**2 + abs(b(k:n-1))**2)
         if (den < tiny(den)) then
            ref_coef = (0.0_wp, 0.0_wp)
         else
            ref_coef = -2.0_wp * num / den
         end if

         ! Levinson recursion
         do i = 1, k-1
            a(i) = a(i) + ref_coef * conjg(a(k-i))
         end do
         a(k) = ref_coef

         ! Update prediction errors
         if (k < p) then
            do i = n, k+1, -1
               f_old = f(i)
               f(i)  = f_old + ref_coef * b(i-1)
               b(i)  = b(i-1) + conjg(ref_coef) * f_old
            end do
         end if
      end do

      phi = -a(1:p)

      deallocate(f, b, a)
   end subroutine complex_burg_ar_coeff

   ! Calculate AR coefficients `phi` for real time series `S` by using Burg method.
   ! Real data is treated as complex with zero imaginary part.
   ! AR coefficients are not negated and one is not present in the front
   subroutine burg_ar_coeff(phi, S)
      real(wp), intent(out) :: phi(:)
      real(wp), intent(in) :: S(:)

      complex(wp), allocatable :: cphi(:), cS(:)

      allocate(cphi(size(phi)), cS(size(S)))
      cS = cmplx(S, 0.0_wp, kind=wp)
      call complex_burg_ar_coeff(cphi, cS)
      phi = real(cphi, wp)
      deallocate(cphi, cS)
   end subroutine burg_ar_coeff

   ! Compute prediction of complex timeseries `S` for AR model with coefficients `phi` at index `i`
   complex(wp) function complex_ar_predict(phi, S, i) result(x)
      complex(wp), intent(in) :: phi(:)
      complex(wp), intent(in) :: S(:)
      integer, intent(in) :: i

      x = sum(phi * S(i - 1:i - size(phi):-1))
   end function complex_ar_predict

   ! Compute prediction of timeseries `S` for AR model with coefficients `phi` at index `i` from previous values (not including `i`).
   ! Real data is treated as complex with zero imaginary part.
   real(wp) function ar_predict(phi, S, i) result(x)
      real(wp), intent(in) :: phi(:)
      real(wp), intent(in) :: S(:)
      integer, intent(in) :: i

      complex(wp), allocatable :: cphi(:), cS(:)

      allocate(cphi(size(phi)), cS(size(S)))
      cphi = cmplx(phi, 0.0_wp, kind=wp)
      cS = cmplx(S, 0.0_wp, kind=wp)
      x = real(complex_ar_predict(cphi, cS, i), wp)
      deallocate(cphi, cS)
   end function ar_predict

   ! Calculate frequency response `H` of an AR filter `phi` on frequencies `f` (0 to 1)
   ! Real data is treated as complex with zero imaginary part.
   subroutine ar_freq_response(phi, f, H, ierr)
      real(wp), intent(in) :: phi(:)
      real(wp), intent(in) :: f(:)
      complex(wp), intent(out) :: H(:)
      integer, intent(out), optional :: ierr

      complex(wp), allocatable :: cphi(:)

      allocate(cphi(size(phi)))
      cphi = cmplx(phi, 0.0_wp, kind=wp)
      call complex_ar_freq_response(cphi, f, H, ierr=ierr)
      deallocate(cphi)
   end subroutine ar_freq_response

   ! Calculate frequency response `H` of a complex AR filter `phi` on frequencies `f` (0 to 1)
   subroutine complex_ar_freq_response(phi, f, H, ierr)
      complex(wp), intent(in) :: phi(:)
      real(wp), intent(in) :: f(:)
      complex(wp), intent(out) :: H(:)
      integer, intent(out), optional :: ierr

      integer :: p, n
      integer :: i, k
      real(wp) :: omega
      complex(wp) :: D

      p = size(phi)
      n = size(f)

      if (check(n == size(H), msg="complex_ar_freq_response: size mismatch", ierr=ierr)) return

      do i = 1, n
         omega = f(i) * 2.0_wp * pi

         D = cmplx(1.0_wp, 0.0_wp, kind=wp)
         do k = 1, p
            D = D - phi(k) * exp(-cmplx(0.0_wp, 1.0_wp, kind=wp) * k * omega)
         end do

         H(i) = cmplx(1.0_wp, 0.0_wp, kind=wp) / D
      end do
   end subroutine complex_ar_freq_response
end module autoreg
