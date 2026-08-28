module stat
   use precision, only: wp
   use checks, only: check
   use stdlib_sorting, only: sort
   use stdlib_optval, only: optval
   implicit none

contains
   real(wp) function mean(S) result(mu)
      real(wp), intent(in) :: S(:)

      mu = sum(S) / size(S)
   end function mean

   ! i = -1.0: eliminates bias
   ! i = +0.0: the variance of the sample
   ! i = +1.0: minimizes mean squared error for the normal distribution
   ! i = +1.5: mostly eliminates bias in unbiased estimation of standard deviation for the normal distribution
   real(wp) function variance(S, i) result(sigma2)
      real(wp), intent(in) :: S(:)
      real(wp), optional, intent(in) :: i
      real(wp) :: d

      d = real(size(S), wp) + optval(i, -1.0_wp)
      sigma2 = sum((S - mean(S))**2) / d
   end function variance

   subroutine percentile(S, q, pct, ierr)
      real(wp), intent(in) :: S(:)
      real(wp), intent(in) :: q(:)
      real(wp), intent(out) :: pct(:)
      integer, intent(out), optional :: ierr

      real(wp), allocatable :: sorted(:)
      integer :: n, i
      real(8) :: rank, fraction
      integer :: lower_idx, upper_idx

      if (check(size(q) == size(pct), msg='percentile: size mismatch', ierr=ierr)) return

      n = size(S)
      allocate(sorted(n))
      sorted(:) = S
      call sort(sorted)

      do i=1,size(q)
         rank = (n - 1) * q(i) / 100.0_wp
         lower_idx = floor(rank) + 1
         upper_idx = min(lower_idx + 1, n)
         fraction = rank - floor(rank)

         if (lower_idx < n) then
            pct(i) = S(lower_idx) + fraction * (S(upper_idx) - S(lower_idx))
         else
            pct(i) = S(n)
         end if
      end do

      deallocate(sorted)
   end subroutine percentile

   ! Friedman-Diaconis method for estimating optimal bins count for historgam
   integer function fd_bins(S) result(bins)
      real(wp), intent(in) :: S(:)
      real(wp) :: pct(2), iqr, range, h

      call percentile(S, [25.0_wp, 75.0_wp], pct)

      iqr = pct(2) - pct(1)
      h = 2.0_wp * iqr / (size(S) ** (1.0_wp / 3.0_wp))
      range = maxval(S) - minval(S)

      bins = ceiling(range / h)
   end function fd_bins

   ! Measures if the tail is longer on the right (positive) or left (negative)
   real(wp) function skewness(S)
      real(wp), intent(in) :: S(:)
      real(wp) :: mu, n

      n = real(size(S), wp)
      mu = mean(S)

      skewness = (sum((S - mu) ** 3) / n) / (sum((S - mu) ** 2) / n) ** (3.0_wp / 2.0_wp)
   end function skewness

   ! Measures if the distribution has heavier tails (positive) or is flatter (negative) than a normal distribution
   real(wp) function kurtosis(S)
      real(wp), intent(in) :: S(:)
      real(wp) :: mu, n

      n = real(size(S), wp)
      mu = mean(S)

      kurtosis = (sum((S - mu) ** 4) / n) / (sum((S - mu) ** 2) / n) ** 2 - 3.0_wp
   end function kurtosis

   elemental real(wp) function normal_cdf(x, mu, sigma) result(cdf)
      real(wp), intent(in) :: x, mu, sigma

      cdf = 0.5 * (1 + erf((x - mu) / (sigma * sqrt(2.0_wp))))
   end function normal_cdf

   ! Kolmogorov-Smirnov Statistic
   ! Accepts sorted in ascending order series
   real(wp) function ks_statistic(S) result(D)
      real(wp), intent(in) :: S(:)
      real(wp) :: sigma, mu
      integer :: i, n

      n = size(S)
      sigma = sqrt(variance(S))
      mu = mean(S)

      D = maxval(abs([(real(i - 1, wp) / n, i = 1, n)] - normal_cdf(S, mu, sigma)))
   end function ks_statistic

   ! Cramer-von Mises Criterion
   ! Accepts sorted in ascending order series
   real(wp) function cvm_criterion(S) result(W2)
      real(wp), intent(in) :: S(:)
      real(wp) :: sigma, mu
      integer :: i, n

      n = size(S)
      sigma = sqrt(variance(S, i=0.0_wp))
      mu = mean(S)

      W2 = 1.0_wp / (12.0_wp * n) + sum(([((2.0_wp * i - 1.0_wp) / (2.0_wp * n), i=1, n)] - normal_cdf(S, mu, sigma)) ** 2)
   end function cvm_criterion

   subroutine fit_normal(S, mu, sigma, se_mu, se_sigma, ks, cvm, sk, ku)
      real(wp), intent(in) :: S(:)
      real(wp), intent(out) :: mu, sigma
      real(wp), optional :: se_mu, se_sigma
      real(wp), optional :: ks, cvm
      real(wp), optional :: sk, ku

      integer :: n
      real(wp), allocatable :: Ss(:)

      n = size(S)

      mu = mean(S)
      sigma = sqrt(variance(S))

      if(present(se_mu))then
         se_mu = sigma / sqrt(real(n, wp))
      end if
      if(present(se_sigma))then
         se_sigma = sigma / sqrt(2.0_wp * (n - 1_wp))
      end if

      if(present(ks) .or. present(cvm)) then
         allocate(Ss(n))
         Ss(:) = S
         call sort(Ss)

         if(present(ks)) then
            ks = ks_statistic(Ss)
         end if
         if(present(cvm)) then
            cvm = cvm_criterion(Ss)
         end if

         deallocate(Ss)
      end if

      if(present(sk)) then
         sk = skewness(S)
      end if
      if(present(ks)) then
         ku = kurtosis(S)
      end if
   end subroutine fit_normal
end module stat
