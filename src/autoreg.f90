module autoreg
   use precision, only: wp
   use constants, only: pi
   use stdlib_linalg, only: solve_lstsq
   use stdlib_error, only: check
   implicit none

contains
   ! Calculate AR coefficients `phi` for time series `S` by solving linear system of equatins with Hankel matrix
   subroutine ar_coeff(phi, S)
      real(wp), intent(out) :: phi(:)
      real(wp), intent(in) :: S(:)

      integer :: n, p, i, j
      real(wp), allocatable :: X(:,:)

      p = size(phi)
      n = size(S) - p

      allocate(X(n, p))

      do i = 1, n
         do j = 1, p
            X(i, j) = S(i + p - j)
         end do
      end do

      call check(size(S(p+1:)) == n, msg="berg_ar_coeff: size mismatch")
      call solve_lstsq(X, S(p+1:), x=phi)

      deallocate(X)
   end subroutine ar_coeff

   ! Compute prediction of timeseries `S` for AR model with coefficients `phi` at index `i` from previous values (not including `i`)
   real(wp) function ar_predict(phi, S, i) result(x)
      real(wp), intent(in) :: phi(:)
      real(wp), intent(in) :: S(:)
      integer, intent(in) :: i

      x = dot_product(phi, S(i - 1:i - size(phi):-1))
   end function ar_predict

   ! Calculate frequency response `H` of an AR filter `phi` on frequencies `f` (0 to 1)
   subroutine ar_freq_response(phi, f, H)
      real(wp), intent(in) :: phi(:)
      real(wp), intent(in) :: f(:)
      complex(wp), intent(out) :: H(:)

      integer :: p, n
      integer :: i, k
      real(wp) :: omega
      complex(wp) :: D

      p = size(phi)
      n = size(f)

      call check(n == size(H), msg="ar_freq_response: size mismatch")

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
