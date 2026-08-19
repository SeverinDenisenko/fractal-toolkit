module fourier
   use precision, only: wp
   use intergers, only: is2power
   use stdlib_error, only: check
   use stdlib_math, only: swap
   use stdlib_constants, only: PI_dp
   implicit none

contains
   subroutine fft(data)
      complex(wp), dimension(:,:), intent(inout) :: data
      call sfft(data, 1)
   end subroutine fft

   subroutine ifft(data)
      complex(wp), dimension(:,:), intent(inout) :: data
      call sfft(data, -1)
   end subroutine ifft

   subroutine sfft(data,isign)
      implicit none
      complex(wp), dimension(:,:), intent(inout) :: data
      integer, intent(in) :: isign

      integer :: n, i, istep, j, m, mmax, n2
      real(wp) :: theta
      complex(wp), dimension(size(data,1)) :: temp
      complex(wp) :: w, wwp
      complex(wp) :: ws
    
      n = size(data,2)

      call check(is2power(n), 'fft: n must be a power of 2')
   
      n2 = n / 2
      j = n2
      do i=1, n-2
         if (j > i) then
            call swap(data(:,j+1), data(:,i+1))
         end if
         m = n2
         do
            if (m < 2 .or. j < m) exit
            j = j - m
            m = m / 2
         end do
         j = j + m
      end do
      mmax = 1
      do
         if (n <= mmax) exit
         istep = 2 * mmax
         theta = PI_dp / (isign * mmax)
         wwp = cmplx(-2.0_wp * sin(0.5_wp * theta)**2, sin(theta), kind=wp)
         w = cmplx(1.0_wp, 0.0_wp, kind=wp)
         do m = 1,mmax
            ws=w
            do i = m,n,istep
               j = i + mmax
               temp = ws * data(:,j)
               data(:,j) = data(:,i) - temp
               data(:,i) = data(:,i) + temp
            end do
            w = w * wwp + w
         end do
         mmax = istep
      end do
   end subroutine sfft
end module fourier
