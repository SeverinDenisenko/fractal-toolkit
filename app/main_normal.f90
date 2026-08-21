program main_normal
   use io, only: read_single_column
   use version, only: max_ver_len, ver
   use precision, only: wp
   use stat, only: fit_normal
   implicit none

   character(len=128) :: arg
   character(len=128) :: input
   character(len=max_ver_len) :: v
   real(wp), allocatable :: series(:)
   real(wp) :: mu, sigma, se_mu, se_sigma, ks, cvm, sk, ku
   integer :: i

   input = 'input.dat'

   i = 1
   do while (i <= command_argument_count())
      call get_command_argument(i, arg)

      select case (arg)
       case ('-v', '--version')
         call ver(v)
         print '(2a)', 'version ', v
         stop
       case ('-h', '--help')
         call print_help()
         stop
       case ('-i', '--input')
         if (i >= command_argument_count()) then
            print '(a)', 'Error: missing value for -i'
            stop 1
         end if
         i = i + 1
         call get_command_argument(i, arg)
         read(arg, *) input
       case default
         print '(a,a,/)', 'Unrecognized command-line option: ', arg
         call print_help()
         stop 1
      end select
      i = i + 1
   end do

   call read_single_column(input, series)

   call fit_normal(series, mu, sigma, se_mu, se_sigma, ks, cvm, sk, ku)

   print '("                             Mean = ", ES10.3)', mu
   print '("               Standard Deviation = ", ES10.3)', sigma
   print '("              Mean Standard Error = ", ES10.3)', se_mu
   print '("Standard Deviation Standard Error = ", ES10.3)', se_sigma
   print '("     Kolmogorov-Smirnov Statistic = ", ES10.3)', ks
   print '("       Cramer-von Mises Criterion = ", ES10.3)', cvm
   print '("                         Skewness = ", ES10.3)', sk
   print '("                  Excess Kurtosis = ", ES10.3)', ku
contains
   subroutine print_help()
      print '(a)', 'usage: normal [options]'
      print '(a)', ''
      print '(a)', 'cmdline options:'
      print '(a)', ''
      print '(a)', '  -v, --version     print version information and exit'
      print '(a)', '  -h, --help        print usage information and exit'
      print '(a)', '  -i, --input       select input file'
   end subroutine print_help
end program main_normal
